"""
Cross-sensor validation: does Landsat 8 see the same meltwater as Sentinel-2?

Independent satellite, independent optics, independent processing chain, same
ground on the same day. Agreement is evidence the detector measures the
surface rather than an artefact of one sensor. Disagreement localises the
fault.

The same detection algorithm is applied to both sensors - NDWI threshold, a
brightness floor, and the NDSI cloud halo - so that any difference is
attributable to the sensor rather than to the method. Thresholds are
expressed in reflectance rather than in raw counts, which makes them
sensor-neutral.

Four things that would silently corrupt this comparison, all found the hard
way while building it:

1. Landsat Collection-2 *Level-2 surface reflectance is invalid over snow and
   ice*. Measured over this study area, 96.9% of its pixels fall outside the
   product's own documented valid range, with green reflectance reading 1.22
   where reflectance cannot exceed 1.0. The LaSRC atmospheric correction is
   not designed for bright cryospheric surfaces. Using it would have produced
   an official-looking validation built on unphysical numbers. This module
   therefore uses Level-1 and computes top-of-atmosphere reflectance directly.

2. Sentinel-2 L2A is bottom-of-atmosphere while Landsat Level-1 is
   top-of-atmosphere. These are different physical quantities. NDWI is a
   normalised ratio and is far more robust to this than raw reflectance is,
   but it is not immune, so this is a stated caveat rather than a solved
   problem - see validate_processing_level.py, which measures the size of the
   effect directly.

3. Landsat is 30 m against Sentinel-2's 10 m and cannot resolve the narrow
   supraglacial channels that carry much of the signal, so it is expected to
   under-detect. That bias is measured, not treated as error.

4. Landsat pixels are warped onto the Sentinel-2 grid using georeferencing
   read from the imagery itself. A first attempt hardcoded the tile origin and
   was wrong by two pixels, which would have shifted every mask by 20 m and
   looked like genuine sensor disagreement rather than a bug.

Data routing, after three dead ends: the AWS copy is Requester Pays, USGS
LandsatLook returns an authentication page, and Planetary Computer's Level-1
collection holds only Landsat 1-5 MSS. Google's public `gcp-public-data-landsat`
bucket serves Landsat 8 Collection-1 Level-1 anonymously and is used here. It
covers 2013 to end-2021, which overlaps four of this project's melt seasons.
Scene discovery still uses STAC; only the pixels come from Google.

Run:  python src/validate_landsat.py
"""

import csv
import glob
import re
from datetime import date

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import rasterio
import requests
from pystac_client import Client
from rasterio.enums import Resampling
from rasterio.vrt import WarpedVRT

import melt

PC_STAC = "https://planetarycomputer.microsoft.com/api/stac/v1"
GCS_BUCKET = "gcp-public-data-landsat"
GCS_LIST = f"https://storage.googleapis.com/storage/v1/b/{GCS_BUCKET}/o"
GCS_DATA = f"https://storage.googleapis.com/{GCS_BUCKET}"

MAX_DAY_GAP = 1
MAX_LS_SCENE_CLOUD = 40.0  # coarse prefilter; the AOI is screened properly later
MIN_SHARED_FRAC = 0.50  # both sensors must see at least half the study area

# Below this, an area ratio is noise over noise and is reported as n/a.
RATIO_FLOOR_KM2 = 0.1

# Reflectance-unit equivalents of the pipeline's DN thresholds, so the same
# rule can be applied to a sensor with different scaling.
S2_SCALE = 1.0 / 10000.0
BRIGHTNESS_FLOOR_REFL = melt.BRIGHTNESS_FLOOR * S2_SCALE

# Landsat 8 Collection-1 BQA bits.
BQA_FILL = 0
BQA_CLOUD = 4
BQA_SHADOW_CONF = 7  # 2-bit field, 3 == high confidence
BQA_CIRRUS_CONF = 11  # 2-bit field, 3 == high confidence

OUT_CSV = melt.ROOT / "output" / "validation_landsat.csv"
OUT_PNG = melt.ROOT / "output" / "validation_landsat.png"


# --- Landsat access ----------------------------------------------------------

def gcs_scene_dir(path, row, acq_date):
    """Locate the Collection-1 scene folder for a path/row/acquisition date."""
    prefix = f"LC08/01/{path:03d}/{row:03d}/"
    params = {"prefix": prefix, "delimiter": "/", "maxResults": 1000}
    token = None
    stamp = acq_date.strftime("%Y%m%d")
    while True:
        if token:
            params["pageToken"] = token
        j = requests.get(GCS_LIST, params=params, timeout=60).json()
        for p in j.get("prefixes", []):
            name = p.rstrip("/").rsplit("/", 1)[-1]
            parts = name.split("_")
            if len(parts) > 3 and parts[3] == stamp:
                return p.rstrip("/"), name
        token = j.get("nextPageToken")
        if not token:
            return None, None


def read_mtl(scene_dir, scene_name):
    """Radiometric rescaling coefficients and sun elevation from the MTL."""
    txt = requests.get(f"{GCS_DATA}/{scene_dir}/{scene_name}_MTL.txt", timeout=60).text
    out = {}
    for line in txt.splitlines():
        if "=" not in line:
            continue
        k, v = line.split("=", 1)
        out[k.strip()] = v.strip().strip('"')
    return out


def landsat_toa(scene_dir, scene_name, band, mtl, ref_item):
    """One band as top-of-atmosphere reflectance on the Sentinel-2 grid.

    TOA reflectance = (MULT * DN + ADD) / sin(sun elevation). The sun-angle
    term is common to all bands and cancels inside NDWI, but it is applied
    anyway so the brightness floor means what it says.
    """
    n = int(re.sub(r"\D", "", band))
    mult = float(mtl[f"REFLECTANCE_MULT_BAND_{n}"])
    add = float(mtl[f"REFLECTANCE_ADD_BAND_{n}"])
    sun = np.radians(float(mtl["SUN_ELEVATION"]))

    dn = _warp_to_aoi(f"{GCS_DATA}/{scene_dir}/{scene_name}_{band}.TIF",
                      ref_item, Resampling.bilinear).astype("f4")
    refl = (mult * dn + add) / np.sin(sun)
    return np.where(dn > 0, refl, np.nan)


def _warp_to_aoi(href, ref_item, resampling):
    crs, _ = melt.tile_georeference(ref_item)
    dst = melt.aoi_transform(ref_item)
    with rasterio.open(href) as src:
        with WarpedVRT(src, crs=crs, transform=dst, width=melt.WIN_SIZE,
                       height=melt.WIN_SIZE, resampling=resampling) as vrt:
            return vrt.read(1)


def landsat_scene(scene_dir, scene_name, ref_item):
    """Green, NIR, SWIR as TOA reflectance plus a validity mask."""
    mtl = read_mtl(scene_dir, scene_name)
    green = landsat_toa(scene_dir, scene_name, "B3", mtl, ref_item)
    nir = landsat_toa(scene_dir, scene_name, "B5", mtl, ref_item)
    swir = landsat_toa(scene_dir, scene_name, "B6", mtl, ref_item)

    bqa = _warp_to_aoi(f"{GCS_DATA}/{scene_dir}/{scene_name}_BQA.TIF",
                       ref_item, Resampling.nearest)
    bad = ((bqa >> BQA_FILL) & 1).astype(bool)
    bad |= ((bqa >> BQA_CLOUD) & 1).astype(bool)
    bad |= ((bqa >> BQA_SHADOW_CONF) & 3) == 3
    bad |= ((bqa >> BQA_CIRRUS_CONF) & 3) == 3

    valid = np.isfinite(green) & np.isfinite(nir) & ~bad
    return green, nir, swir, valid


# --- Shared detection --------------------------------------------------------

def detect_reflectance(green, nir, swir, valid):
    """The pipeline's detector, expressed in reflectance so it is sensor-neutral."""
    with np.errstate(invalid="ignore", divide="ignore"):
        ndwi = np.where(valid, (green - nir) / np.maximum(green + nir, 1e-6), np.nan)
        ndsi = np.where(valid, (green - swir) / np.maximum(green + swir, 1e-6), np.nan)

    ponds = (np.nan_to_num(ndwi, nan=-9) > melt.NDWI_THRESHOLD)
    ponds &= valid & (np.nan_to_num(green, nan=0) > BRIGHTNESS_FLOOR_REFL)

    # Same cloud-halo logic as the pipeline, at native resolution here.
    from scipy import ndimage

    low = valid & (np.nan_to_num(ndsi, nan=1.0) < melt.NDSI_ICE_MIN)
    px = int(round(melt.CLOUD_HALO_KM * 1000 / melt.PIXEL_M))
    # Dilate on a decimated grid; a 100-pixel structuring element is ruinous.
    d = 10
    low_small = low[::d, ::d]
    halo_small = ndimage.binary_dilation(low_small, iterations=max(px // d, 1))
    halo = np.kron(halo_small, np.ones((d, d), bool))[:low.shape[0], :low.shape[1]]

    return ponds & ~halo, valid & ~halo


def s2_scene(item):
    bands = melt.load_scene(item, bands=("green", "nir"))
    green = bands["green"] * S2_SCALE
    nir = bands["nir"] * S2_SCALE
    swir = melt.read_coarse(item, "swir16")
    off = melt.boa_offset(item)
    swir = np.where(swir > 0, (swir + off) * S2_SCALE, np.nan)
    scale = melt.WIN_SIZE // swir.shape[0]
    swir = np.kron(swir, np.ones((scale, scale), "f4"))

    valid = (bands["green"] > 0) & (bands["nir"] > 0)
    valid &= ~melt.reject_mask(item)
    return green, nir, swir, valid


# --- Pairing and comparison --------------------------------------------------

def find_pairs():
    rows = []
    for p in sorted(glob.glob(str(melt.ROOT / "output" / "season_*.csv"))):
        with open(p, newline="") as f:
            rows.extend(list(csv.DictReader(f)))
    rows.sort(key=lambda r: r["date"])
    if not rows:
        return []

    ref_item = melt.get_item(rows[0]["item_id"])
    crs, _ = melt.tile_georeference(ref_item)
    tr = melt.aoi_transform(ref_item)
    from rasterio.warp import transform_bounds

    x0, y0 = tr * (0, 0)
    x1, y1 = tr * (melt.WIN_SIZE, melt.WIN_SIZE)
    w, s, e, n = transform_bounds(crs, "EPSG:4326", min(x0, x1), min(y0, y1),
                                  max(x0, x1), max(y0, y1))

    cat = Client.open(PC_STAC)
    ls = list(cat.search(collections=["landsat-c2-l2"], bbox=[w, s, e, n],
                         datetime="2013-01-01/2021-12-31",
                         query={"eo:cloud_cover": {"lt": MAX_LS_SCENE_CLOUD}},
                         limit=100).items())
    print(f"[stac] {len(ls)} Landsat scenes over the AOI in the "
          f"Collection-1 era, under {MAX_LS_SCENE_CLOUD:.0f}% scene cloud")

    pairs, seen = [], set()
    for r in rows:
        d = date.fromisoformat(r["date"])
        best = None
        for i in ls:
            gap = abs((i.datetime.date() - d).days)
            if gap > MAX_DAY_GAP:
                continue
            key = (gap, i.properties.get("eo:cloud_cover") or 100)
            if best is None or key < best[0]:
                best = (key, i)
        if best and (r["date"], best[1].id) not in seen:
            seen.add((r["date"], best[1].id))
            pairs.append((r, best[1]))
    return pairs


RES_RATIO = 3  # 10 m Sentinel-2 -> 30 m Landsat


def _block_mean(a, k=RES_RATIO):
    h = a.shape[0] // k * k
    with np.errstate(invalid="ignore"):
        return np.nanmean(a[:h, :h].reshape(h // k, k, h // k, k), axis=(1, 3))


def _block_any(a, k=RES_RATIO):
    h = a.shape[0] // k * k
    return a[:h, :h].reshape(h // k, k, h // k, k).any(axis=(1, 3))


def compare(row, ls_item):
    s2_item = melt.get_item(row["item_id"])
    path = int(ls_item.properties["landsat:wrs_path"])
    wrow = int(ls_item.properties["landsat:wrs_row"])
    scene_dir, scene_name = gcs_scene_dir(path, wrow, ls_item.datetime.date())
    if scene_dir is None:
        return None, "no Collection-1 scene on GCS"

    lg, ln, lsw, lv = landsat_scene(scene_dir, scene_name, s2_item)
    ls_ponds, lv = detect_reflectance(lg, ln, lsw, lv)

    sg, sn, ssw, sv = s2_scene(s2_item)
    s2_ponds, sv = detect_reflectance(sg, sn, ssw, sv)

    both = sv & lv
    if both.mean() < MIN_SHARED_FRAC:
        return None, f"only {both.mean()*100:.0f}% shared clear ground"

    a, b = s2_ponds & both, ls_ponds & both
    inter, union = int((a & b).sum()), int((a | b).sum())

    # Matched-resolution comparison. Sentinel-2 reflectance is averaged to
    # Landsat's 30 m grid *before* detection, which is the only like-for-like
    # test: at 10 m the two sensors are not measuring the same quantity,
    # because a 10 m channel is a sub-pixel mixture to Landsat and its NDWI is
    # diluted below threshold no matter how good either instrument is.
    sg3, sn3 = _block_mean(np.where(sv, sg, np.nan)), _block_mean(np.where(sv, sn, np.nan))
    v3 = np.isfinite(sg3) & np.isfinite(sn3)
    with np.errstate(invalid="ignore"):
        ndwi3 = (sg3 - sn3) / np.maximum(sg3 + sn3, 1e-6)
    a3 = (np.nan_to_num(ndwi3, nan=-9) > melt.NDWI_THRESHOLD) & v3
    a3 &= np.nan_to_num(sg3, nan=0) > BRIGHTNESS_FLOOR_REFL
    both3 = _block_any(both) & v3
    a3, b3 = a3 & both3, _block_any(ls_ponds) & both3
    inter3, union3 = int((a3 & b3).sum()), int((a3 | b3).sum())

    # Positional tolerance: does each detection have a counterpart within one
    # Landsat pixel? Separates "wrong place" from "slightly offset", which
    # matters because inter-sensor co-registration is itself ~10-15 m.
    from scipy import ndimage

    near_b = ndimage.binary_dilation(b, iterations=RES_RATIO)
    near_a = ndimage.binary_dilation(a, iterations=RES_RATIO)
    matched = int((a & near_b).sum()) + int((b & near_a).sum())
    tol = matched / max(int(a.sum()) + int(b.sum()), 1)

    px30 = (melt.PIXEL_M * RES_RATIO) ** 2 / 1e6
    return {
        "date_s2": row["date"],
        "date_ls": ls_item.datetime.date().isoformat(),
        "s2_id": row["item_id"],
        "ls_id": scene_name,
        "day_gap": abs((ls_item.datetime.date() - date.fromisoformat(row["date"])).days),
        "shared_frac": round(float(both.mean()), 4),
        "s2_km2": round(int(a.sum()) * melt.PIXEL_KM2, 4),
        "ls_km2": round(int(b.sum()) * melt.PIXEL_KM2, 4),
        "iou_native": round(inter / union, 4) if union else float("nan"),
        "s2_km2_at30": round(int(a3.sum()) * px30, 4),
        "ls_km2_at30": round(int(b3.sum()) * px30, 4),
        "iou_matched": round(inter3 / union3, 4) if union3 else float("nan"),
        "agree_within_1px": round(tol, 4),
        "ls_recall_of_s2": round(inter / int(a.sum()), 4) if a.sum() else float("nan"),
    }, None


def main():
    pairs = find_pairs()
    print(f"[pairs] {len(pairs)} same-day Sentinel-2 / Landsat pairs "
          f"within {MAX_DAY_GAP} day(s)\n")

    results = []
    for row, ls_item in pairs:
        try:
            res, why = compare(row, ls_item)
        except Exception as e:
            print(f"  {row['date']}  ERROR {type(e).__name__}: {str(e)[:70]}")
            continue
        if res is None:
            print(f"  {row['date']}  skip ({why})")
            continue
        results.append(res)
        # Both sensors near zero is agreement, not a ratio; printing one is
        # meaningless and invites reading noise as a 36-million-fold error.
        if min(res["s2_km2_at30"], res["ls_km2_at30"]) < RATIO_FLOOR_KM2:
            r30 = "  n/a"
        else:
            r30 = f"{res['ls_km2_at30'] / res['s2_km2_at30']:6.2f}"
        print(f"  {res['date_s2']} vs {res['date_ls']} (gap {res['day_gap']}d, "
              f"shared {res['shared_frac']*100:3.0f}%)  "
              f"S2 {res['s2_km2']:6.2f} LS {res['ls_km2']:6.2f} km2 | "
              f"IoU {res['iou_native']:.2f} native, {res['iou_matched']:.2f} matched | "
              f"within-1px {res['agree_within_1px']:.2f} | LS/S2@30m {r30}")

    if not results:
        print("\nNo comparable pairs.")
        return

    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(results[0].keys()))
        w.writeheader()
        w.writerows(results)

    s2 = np.array([r["s2_km2_at30"] for r in results])
    ls = np.array([r["ls_km2_at30"] for r in results])
    iou_n = np.array([r["iou_native"] for r in results])
    iou_m = np.array([r["iou_matched"] for r in results])
    tol = np.array([r["agree_within_1px"] for r in results])
    ratio = ls / np.maximum(s2, 1e-9)

    # Same-day, well-covered pairs are the only ones that test the detector
    # rather than the weather; a one-day gap lets ponds freeze, drain or fill.
    best = [i for i, r in enumerate(results)
            if r["day_gap"] == 0 and r["shared_frac"] >= 0.8]

    print(f"\n  pairs compared              {len(results)}")
    print(f"  median IoU, native res      {np.median(iou_n):.3f}")
    print(f"  median IoU, matched 30 m    {np.median(iou_m):.3f}")
    print(f"  median agreement within 1px {np.median(tol):.3f}")
    print(f"  median LS/S2 area @30 m     {np.median(ratio):.2f}")
    if len(s2) > 2:
        print(f"  Pearson r (areas)           {np.corrcoef(s2, ls)[0, 1]:.3f}")
        rk = lambda v: np.argsort(np.argsort(v))
        print(f"  Spearman (rank)             {np.corrcoef(rk(s2), rk(ls))[0, 1]:.3f}")
    if best:
        print(f"\n  same-day pairs with >=80% shared ground: {len(best)}")
        print(f"    median IoU matched        {np.median(iou_m[best]):.3f}")
        print(f"    median within-1px         {np.median(tol[best]):.3f}")
        print(f"    median LS/S2 area         {np.median(ratio[best]):.2f}")

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(14, 5.8))
    lim = float(max(s2.max(), ls.max())) * 1.15 or 1.0
    ax.plot([0, lim], [0, lim], "--", color="#888", lw=1, label="1:1")
    for i, r in enumerate(results):
        same = r["day_gap"] == 0 and r["shared_frac"] >= 0.8
        ax.scatter(s2[i], ls[i], s=80 if same else 45,
                   color="#0077b6" if same else "#b8c4cc",
                   edgecolor="white", zorder=3)
        ax.annotate(r["date_s2"][:7], (s2[i], ls[i]), textcoords="offset points",
                    xytext=(6, -9), fontsize=7.5, color="#555")
    ax.scatter([], [], s=80, color="#0077b6", label="same day, >=80% shared")
    ax.scatter([], [], s=45, color="#b8c4cc", label="1-day gap or partial")
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.set_xlabel("Sentinel-2 meltwater (km$^2$, resampled to 30 m)")
    ax.set_ylabel("Landsat 8 meltwater (km$^2$, 30 m)")
    ax.set_title("Same ground, same day, independent satellite")
    ax.grid(alpha=0.3)
    ax.legend(loc="upper left", fontsize=8.5)

    x = np.arange(len(results))
    ax2.bar(x - 0.2, iou_n, width=0.4, color="#b8c4cc", label="IoU at native res")
    ax2.bar(x + 0.2, iou_m, width=0.4, color="#0077b6", label="IoU at matched 30 m")
    ax2.plot(x, tol, "o-", color="#c1121f", ms=4, lw=1.2,
             label="agreement within 1 Landsat pixel")
    ax2.set_xticks(x)
    ax2.set_xticklabels([r["date_s2"] for r in results], rotation=45,
                        ha="right", fontsize=8)
    ax2.set_ylabel("agreement")
    ax2.set_ylim(0, 1)
    ax2.set_title("Resolution explains most of the disagreement")
    ax2.grid(alpha=0.3, axis="y")
    ax2.legend(fontsize=8)

    fig.suptitle("Cross-sensor validation - George VI Ice Shelf", fontsize=13)
    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=130, bbox_inches="tight")
    print(f"\n[write] {OUT_CSV.relative_to(melt.ROOT)}")
    print(f"[write] {OUT_PNG.relative_to(melt.ROOT)}")


if __name__ == "__main__":
    main()

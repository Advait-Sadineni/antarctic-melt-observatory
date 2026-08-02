"""M4 Task 4b: depth maps + Monte Carlo volumes for GVI pilot seasons.

    python scripts/depth_product.py 2019-20

Recipe (spec v1.1): peak-season production scene -> dual-band depths
(g from output/depth/calibration.json; red clipped to its <=3 m validity,
green secondary) -> validity = PONDED cells of the fused melt-state record
(state_<season>.npz, 120 m) upsampled to 30 m -> volume with Monte Carlo
uncertainty. Writes output/depth/product_<season>.json + depth GeoTIFF.
Gates evaluated here: centre-deeper-than-rim gradient; volume recorded for
the Corr et al. (2022) x2 comparison.
"""
import json
import sys
from pathlib import Path

import numpy as np
import rasterio
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
import depth
import icesat2
import melt
import shelf

RED_ZMAX = 3.0          # red-band validity (Pope-lineage; deeper is green's job)
WATER_TILES = ("19DEA", "19CEV")
R_INF = 0.05


def state_validity_mask(season, aoi_shape, tile_item):
    """PONDED-ever cells from the fused state record, nearest-upsampled from
    the 120 m SAR grid onto this scene's 30 m AOI pixels."""
    import sar
    z = np.load(sar.OUT / f"state_{season}_t2.npz")
    ponded = z["days_ponded"] > 0
    grid = shelf.build_fixed_grid()
    tr4, gw4, gh4, _ = sar.sar_grid(grid)
    from rasterio.enums import Resampling
    from rasterio.warp import reproject
    crs, _ = melt.tile_georeference(tile_item)
    dst = np.zeros(aoi_shape, "u1")
    reproject(source=ponded.astype("u1"), destination=dst,
              src_transform=tr4, src_crs=shelf.GRID_CRS,
              dst_transform=melt.aoi_transform(tile_item), dst_crs=crs,
              resampling=Resampling.nearest)
    return dst.astype(bool)


def centre_rim_gate(zmap, ponds, min_pond_px=100):
    """Physics sanity: pond centres must be deeper than pond rims."""
    labels, n = ndimage.label(ponds & np.isfinite(zmap))
    deeper = total = 0
    for i, sl in enumerate(ndimage.find_objects(labels), start=1):
        if sl is None:
            continue
        m = labels[sl] == i
        if m.sum() < min_pond_px:
            continue
        er = ndimage.binary_erosion(m, iterations=3)
        rim = m & ~ndimage.binary_erosion(m)
        if not er.any() or not rim.any():
            continue
        zc = np.nanmean(zmap[sl][er])
        zr = np.nanmean(zmap[sl][rim])
        total += 1
        deeper += bool(zc > zr)
    return deeper, total


def run(season):
    calib = json.loads((icesat2.OUT / "calibration.json").read_text())
    bands_fit = calib.get("bands_tight") or calib["bands"]  # decision 0003:
    g_red = bands_fit["red"]["g"]                           # tight-window fit
    g_green = bands_fit["green"]["g"]

    y1 = int(season.split("-")[0]) + 1
    out_rows = []
    for tile in WATER_TILES:
        it = shelf.production_scene(tile, f"{y1}-01-01", f"{y1}-02-28")
        if it is None:
            continue
        bands = melt.load_scene(it, bands=("green", "nir", "red"))
        reject = melt.reject_mask(it) | melt.cloud_mask(it)["mask"]
        _, ponds, _ = melt.detect(bands["green"], bands["nir"], reject,
                                  red=bands["red"])
        _, ad_red = depth.rim_albedo(bands["red"], ponds)
        _, ad_green = depth.rim_albedo(bands["green"], ponds)
        zr = depth.depth_single(bands["red"], np.nanmean(ad_red[ponds]), R_INF, g_red)
        zr[zr > RED_ZMAX] = np.nan
        zg = depth.depth_single(bands["green"], np.nanmean(ad_green[ponds]), R_INF, g_green)
        z, uncertain = depth.dual_band_merge(zr, zg)
        z_final = np.where(np.isfinite(z), z,
                           np.where(np.isfinite(zr), zr, zg))  # merged, else best single
        valid = state_validity_mask(season, ponds.shape, it) & ponds
        z_final = np.where(valid, z_final, np.nan).astype("f4")

        deeper, total = centre_rim_gate(z_final, valid)
        px_area = 100.0   # 10 m pixels

        def z_fn(g_shift, ad_shift, rinf, _b=bands, _ad=ad_red, _v=valid):
            zz = depth.depth_single(_b["red"],
                                    np.nanmean(_ad[ponds]) + ad_shift,
                                    rinf, g_red * (1 + g_shift))
            zz[zz > RED_ZMAX] = np.nan
            return np.where(_v, zz, np.nan)

        mc = depth.mc_volume(z_fn, px_area, n=500,
                             g_sigma=0.2,   # calibration posterior (decision 0003)
                             ad_sigma=0.02, rinf_range=(0.03, 0.07))
        crs, _ = melt.tile_georeference(it)
        prof = dict(driver="GTiff", height=z_final.shape[0], width=z_final.shape[1],
                    count=1, dtype="float32", crs=crs,
                    transform=melt.aoi_transform(it),
                    tiled=True, blockxsize=512, blockysize=512,
                    compress="deflate", predictor=2, nodata=np.nan)
        tifp = icesat2.OUT / f"depth_{season}_{tile}.tif"
        with rasterio.open(tifp, "w", **prof) as dst:
            dst.write(z_final, 1)
        row = {"tile": tile, "scene": it.datetime.date().isoformat(),
               "ponded_valid_km2": round(float(np.isfinite(z_final).sum()) * px_area / 1e6, 2),
               "depth_mean_m": round(float(np.nanmean(z_final)), 2)
                   if np.isfinite(z_final).any() else None,
               "depth_p90_m": round(float(np.nanpercentile(z_final, 90)), 2)
                   if np.isfinite(z_final).any() else None,
               "volume_km3": mc, "uncertain_frac": round(float(
                   (uncertain & valid).sum() / max(valid.sum(), 1)), 3),
               "gate_centre_deeper": {"deeper": deeper, "ponds": total}}
        out_rows.append(row)
        print(f"  [{season} {tile}] {row}", flush=True)

    outp = icesat2.OUT / f"product_{season}.json"
    outp.write_text(json.dumps({"season": season, "g_red": g_red,
                                "g_green": g_green, "tiles": out_rows}, indent=1))
    print(f"[{season}] product written", flush=True)


if __name__ == "__main__":
    run(sys.argv[1] if len(sys.argv) > 1 else "2020-21")

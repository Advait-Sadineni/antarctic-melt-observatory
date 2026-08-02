"""M4 Task 4a: calibrate the attenuation coefficient g against ICESat-2.

For every QC-passing crossover depth: reload the exact optical scene the pond
came from, take the pond-masked reflectance around the photon segment (red and
green), the pond's rim albedo, and form X = ln(A_d-R_inf) - ln(R_w-R_inf).
Then g is a closed-form least-squares fit (z = X/g) with a 3-sigma trim.

    python scripts/depth_calibrate.py        -> output/depth/calibration.json
"""
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from rasterio.transform import rowcol
from rasterio.warp import transform as warp_transform
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
import depth
import icesat2
import melt
import shelf

R_INF = 0.05
RMSE_GATE_M = 0.5
N_GATE = 30
SHELF_CTX = {  # tile -> set_shelf args (None = default George VI context)
    "41CPV": ("Amery", ["41CPV"], "reference/boundaries/Amery.geojson"),
}


def rows_by_scene():
    groups = defaultdict(list)
    for f in sorted((icesat2.OUT).glob("crossovers_*.json")):
        for r in json.loads(f.read_text()):
            if not r.get("qc_pass"):
                continue
            tile = r["pond"].split("-")[0]
            groups[(tile, r["scene_date"])].append(r)
    return groups


def scene_for(tile, day):
    ctx = SHELF_CTX.get(tile)
    if ctx is not None:
        name, tiles, bnd = ctx
        shelf.set_shelf(name, tiles, str(Path(bnd).resolve()))
    src = melt.get_source()
    items = src.search(query=src.tile_query(tile),
                       datetime=f"{day}/{day}")
    if not items:
        return None
    it = items[0]
    row, col, size = shelf._tile_window(it)
    if size < 256:
        return None
    melt.set_aoi(tile, row, col, size)
    return it


def sample(band_img, ad_img, ponds, rr, cc, halfwin=2):
    """Pond-masked window means of reflectance and rim albedo at (rr, cc)."""
    r0, r1 = max(0, rr - halfwin), rr + halfwin + 1
    c0, c1 = max(0, cc - halfwin), cc + halfwin + 1
    m = ponds[r0:r1, c0:c1]
    if not m.any():
        return None, None
    return (float(np.nanmean(band_img[r0:r1, c0:c1][m])),
            float(np.nanmean(ad_img[r0:r1, c0:c1][m])))


def main():
    groups = rows_by_scene()
    print(f"[calib] {sum(len(v) for v in groups.values())} QC depths in "
          f"{len(groups)} scene groups", flush=True)
    pairs = []          # (site, X_red, X_green, z)
    for (tile, day), rows in sorted(groups.items()):
        try:
            it = scene_for(tile, day)
        except Exception as ex:
            print(f"  [skip {tile} {day}] {type(ex).__name__}", flush=True)
            continue
        if it is None:
            print(f"  [no scene] {tile} {day}", flush=True)
            continue
        bands = melt.load_scene(it, bands=("green", "nir", "red"))
        reject = melt.reject_mask(it) | melt.cloud_mask(it)["mask"]
        _, ponds, _ = melt.detect(bands["green"], bands["nir"], reject,
                                  red=bands["red"])
        crs, _ = melt.tile_georeference(it)
        tr = melt.aoi_transform(it)
        _, ad_red = depth.rim_albedo(bands["red"], ponds)
        _, ad_green = depth.rim_albedo(bands["green"], ponds)
        site = "amery" if tile in SHELF_CTX else "george_vi"
        n0 = len(pairs)
        for r in rows:
            xs, ys = warp_transform("EPSG:4326", crs, [r["lon"]], [r["lat"]])
            rr, cc = rowcol(tr, xs[0], ys[0])
            if not (0 <= rr < ponds.shape[0] and 0 <= cc < ponds.shape[1]):
                continue
            rw_r, adr = sample(bands["red"], ad_red, ponds, rr, cc)
            rw_g, adg = sample(bands["green"], ad_green, ponds, rr, cc)
            if rw_r is None:
                continue
            Xr = depth.calib_X(np.array([rw_r]), adr, R_INF)[0]
            Xg = depth.calib_X(np.array([rw_g]), adg, R_INF)[0]
            pairs.append((site, float(Xr), float(Xg), r["depth_m"],
                          r.get("gap_days")))
        print(f"  {tile} {day}: +{len(pairs)-n0} pairs (total {len(pairs)})",
              flush=True)

    (icesat2.OUT / "calibration_pairs.json").write_text(json.dumps(
        [{"site": s, "X_red": xr, "X_green": xg, "z": z_, "gap_days": gp}
         for s, xr, xg, z_, gp in pairs], indent=1))
    out = {"r_inf": R_INF, "pairs_total": len(pairs), "bands": {},
           "bands_tight": {}, "sites": {}}
    Xr = np.array([p[1] for p in pairs])
    Xg = np.array([p[2] for p in pairs])
    z = np.array([p[3] for p in pairs])
    gap = np.array([p[4] if p[4] is not None else 99 for p in pairs])
    tight = gap <= 3          # the spec's original crossover window
    for band, X in (("red", Xr), ("green", Xg)):
        # red is fit inside its physical validity (<= 3 m; Pope-lineage)
        dom = (z <= 3.0) if band == "red" else np.ones(len(z), bool)
        fit = depth.fit_g(X[dom], z[dom])
        fit["gate_rmse_lt_0.5"] = bool(fit["rmse_m"] < RMSE_GATE_M)
        fit["gate_n_ge_30"] = bool(fit["n"] >= N_GATE)
        out["bands"][band] = fit
        if tight.any():
            ft = depth.fit_g(X[dom & tight], z[dom & tight])
            ft["gate_rmse_lt_0.5"] = bool(ft["rmse_m"] < RMSE_GATE_M)
            ft["gate_n_ge_30"] = bool(ft["n"] >= N_GATE)
            out["bands_tight"][band] = ft
    for site in ("george_vi", "amery"):
        sel = np.array([p[0] == site for p in pairs])
        if sel.sum() >= 5:
            out["sites"][site] = depth.fit_g(Xr[sel & (z <= 3.0)], z[sel & (z <= 3.0)])
    (icesat2.OUT / "calibration.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1), flush=True)


if __name__ == "__main__":
    main()

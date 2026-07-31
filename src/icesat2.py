"""M4: ICESat-2 ATL03 crossover harvest for depth calibration (spec v1.1 §3).

The laser altimeter's photons return from BOTH the pond surface and the pond
bottom; the elevation histogram shows two peaks whose separation (x0.752
refraction correction, Parrish et al. 2019) is a direct depth measurement -
the ground truth our attenuation model calibrates against. Strong beams only.

Networked pieces (search/download) run through earthaccess with the netrc
credentials already configured; the peak-splitter is pure and offline-tested.
"""
import json
from datetime import date, timedelta
from pathlib import Path

import numpy as np

import melt

REFRACTION = 0.752
MIN_BOTTOM_PHOTONS = 50
MIN_SEPARATION_M = 0.3
BIN_M = 0.10
OUT = melt.ROOT / "output" / "depth"
STRONG_SUFFIX = {"forward": ("gt1l", "gt2l", "gt3l"),
                 "backward": ("gt1r", "gt2r", "gt3r")}


# --- pure (offline-tested) ----------------------------------------------------

def photon_peak_depth(elevations):
    """Depth from a photon elevation histogram: surface peak = strongest bin;
    bottom peak = strongest bin >= MIN_SEPARATION_M below it. None if the
    bottom return is too weak (< MIN_BOTTOM_PHOTONS in its peak +/- 1 bin)."""
    e = np.asarray(elevations, "f8")
    e = e[np.isfinite(e)]
    if e.size < MIN_BOTTOM_PHOTONS * 2:
        return None
    lo, hi = np.percentile(e, [1, 99])
    bins = np.arange(lo, hi + BIN_M, BIN_M)
    hist, edges = np.histogram(e, bins=bins)
    if hist.size < 5:
        return None
    i_surf = int(np.argmax(hist))
    surf = 0.5 * (edges[i_surf] + edges[i_surf + 1])
    below = edges[:-1] < (surf - MIN_SEPARATION_M)
    if not below.any():
        return None
    masked = np.where(below, hist, 0)
    i_bot = int(np.argmax(masked))
    support = hist[max(0, i_bot - 1): i_bot + 2].sum()
    if support < MIN_BOTTOM_PHOTONS:
        return None
    bottom = 0.5 * (edges[i_bot] + edges[i_bot + 1])
    return float((surf - bottom) * REFRACTION)


# --- networked ---------------------------------------------------------------

def search_atl03(bbox, day, pad_days=3):
    """ATL03 granules over bbox within +/-pad_days of `day`."""
    import earthaccess
    earthaccess.login(strategy="netrc")
    t0 = (day - timedelta(days=pad_days)).isoformat()
    t1 = (day + timedelta(days=pad_days)).isoformat()
    return earthaccess.search_data(short_name="ATL03",
                                   bounding_box=tuple(bbox),
                                   temporal=(t0, t1))


def granule_pond_depths(h5_path, pond_bboxes):
    """Depths from one downloaded ATL03 granule: for each strong beam, photons
    falling inside each pond lon/lat bbox -> peak-split depth. Returns
    [(pond_idx, lat, lon, depth_m), ...]."""
    import h5py
    results = []
    with h5py.File(h5_path, "r") as f:
        orient = int(f["orbit_info/sc_orient"][0])  # 1=forward, 0=backward
        beams = STRONG_SUFFIX["forward" if orient == 1 else "backward"]
        for b in beams:
            try:
                lat = f[f"{b}/heights/lat_ph"][:]
                lon = f[f"{b}/heights/lon_ph"][:]
                h = f[f"{b}/heights/h_ph"][:]
            except KeyError:
                continue
            for i, (w, s, e, n) in enumerate(pond_bboxes):
                sel = (lon >= w) & (lon <= e) & (lat >= s) & (lat <= n)
                if sel.sum() < MIN_BOTTOM_PHOTONS * 2:
                    continue
                d = photon_peak_depth(h[sel])
                if d is not None:
                    results.append((i, float(lat[sel].mean()),
                                    float(lon[sel].mean()), d))
    return results


def harvest(season, pond_geojson, out_name=None):
    """Crossover harvest for one season: for each (pond bbox, clear-scene date)
    pair in pond_geojson, search/download ATL03 and extract depths. Writes
    output/depth/crossovers_<season>.json."""
    import earthaccess
    OUT.mkdir(parents=True, exist_ok=True)
    gj = json.loads(Path(pond_geojson).read_text())
    rows = []
    for feat in gj["features"]:
        w, s, e, n = feat["bbox"]
        day = date.fromisoformat(feat["properties"]["scene_date"])
        grans = search_atl03([w, s, e, n], day)
        if not grans:
            continue
        files = earthaccess.download(grans, str(OUT / "atl03"))
        for fp in files:
            for idx, la, lo, d in granule_pond_depths(fp, [(w, s, e, n)]):
                rows.append({"pond": feat["properties"]["pond_id"],
                             "scene_date": feat["properties"]["scene_date"],
                             "lat": la, "lon": lo, "depth_m": d})
    out = OUT / (out_name or f"crossovers_{season}.json")
    out.write_text(json.dumps(rows, indent=1))
    print(f"[icesat2] {season}: {len(rows)} crossover depths -> {out.name}")
    return rows

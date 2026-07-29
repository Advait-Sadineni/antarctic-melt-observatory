"""Sentinel-1 wet-snow detection for the observatory (M3 spec, phase A).

Wet snow collapses radar backscatter: a pixel is WET when its gamma0 falls more
than THRESH_DB below that pixel's own frozen-winter baseline. Baselines are
per-pixel medians built PER RELATIVE ORBIT (incidence geometry differs by
track; mixing tracks smears the reference). Everything lands on the same fixed
EPSG:3031 grid the optical record uses, so the two sensors are pixel-comparable
by construction.

Networked entry points (winter_baseline, season_composite) are driven by
scripts/sar_gates.py; the math (to_db, _median_stack, wet_mask,
CompositeAccumulator) is pure and offline-tested.
"""
import json
from datetime import date

import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.vrt import WarpedVRT

import melt
import shelf

THRESH_DB = 3.0                       # spec: sweep 2-4 dB at gate time
OUT = melt.ROOT / "output" / "sar"


# --- pure math (offline-tested) ----------------------------------------------

def to_db(gamma0):
    """Backscatter power -> dB; non-positive (nodata) -> NaN."""
    g = np.asarray(gamma0, "f4")
    with np.errstate(divide="ignore", invalid="ignore"):
        db = 10.0 * np.log10(np.where(g > 0, g, np.nan))
    return db.astype("f4")


def _median_stack(arrays):
    """Per-pixel median across scenes, NaN-aware (nodata never votes)."""
    with np.errstate(all="ignore"):
        return np.nanmedian(np.stack(arrays), axis=0).astype("f4")


def wet_mask(db, baseline, thresh=THRESH_DB):
    """(wet, observed): wet where observed AND >= thresh dB below baseline.
    NaN in either scene or baseline means not-observed, never wet."""
    obs = ~np.isnan(db) & ~np.isnan(baseline)
    wet = np.zeros(db.shape, bool)
    wet[obs] = db[obs] < (baseline[obs] - thresh)
    return wet, obs


class CompositeAccumulator:
    """Per-pixel season bookkeeping: wet-observation count, observation count,
    first/last wet day-of-season (-1 = never)."""

    def __init__(self, shape):
        self.wet_days = np.zeros(shape, "u2")
        self.n_obs = np.zeros(shape, "u2")
        self.first_wet = np.full(shape, -1, "i2")
        self.last_wet = np.full(shape, -1, "i2")

    def add(self, wet, obs, day):
        self.wet_days += wet.astype("u2")
        self.n_obs += obs.astype("u2")
        newly = wet & (self.first_wet < 0)
        self.first_wet[newly] = day
        self.last_wet[wet] = day


# --- scene access (networked) -------------------------------------------------

def s1_source():
    from core.pc import Sentinel1Source
    return Sentinel1Source()


def list_scenes(source, bbox, start, end):
    """Scenes grouped by relative orbit, each list date-sorted."""
    items = source.search(collections=[source.COLLECTION], bbox=bbox,
                          datetime=f"{start}/{end}", limit=100)
    by_orbit = {}
    for it in items:
        by_orbit.setdefault(source.relative_orbit(it), []).append(it)
    for orb in by_orbit:
        by_orbit[orb].sort(key=lambda x: x.datetime)
    return by_orbit


def read_db_on_grid(item, grid_tr, gw, gh, source):
    """gamma0 dB warped onto the fixed grid; NaN outside the swath."""
    with rasterio.open(source.band_href(item, "hh")) as src:
        with WarpedVRT(src, crs=shelf.GRID_CRS, transform=grid_tr,
                       width=gw, height=gh,
                       resampling=Resampling.average, nodata=0) as vrt:
            g = vrt.read(1).astype("f4")
    return to_db(g)


def winter_baseline(year, orbit, grid, source, bbox, rebuild=False):
    """Per-pixel median dB over Jun-Aug of `year`, one orbit, disk-cached."""
    OUT.mkdir(parents=True, exist_ok=True)
    cache = OUT / f"baseline_{year}_{orbit}.npy"
    if cache.exists() and not rebuild:
        return np.load(cache)
    grid_tr, gw, gh, _ = grid
    scenes = list_scenes(source, bbox, f"{year}-06-01", f"{year}-08-31")
    stack = [read_db_on_grid(it, grid_tr, gw, gh, source)
             for it in scenes.get(orbit, [])]
    if len(stack) < 3:
        raise ValueError(f"only {len(stack)} winter scenes for orbit {orbit}")
    base = _median_stack(stack)
    np.save(cache, base)
    print(f"  [baseline] {year} orbit {orbit}: {len(stack)} scenes")
    return base


def season_composite(season, grid, source, bbox, thresh=THRESH_DB,
                     rebuild=False):
    """Nov-Mar wet/dry composite for one melt season on the fixed grid.

    Returns summary dict; saves npz (wet_days, n_obs, first_wet, last_wet) and
    a JSON summary. Baselines: preceding austral winter (Jun-Aug of the
    season's first year), per orbit."""
    OUT.mkdir(parents=True, exist_ok=True)
    tag = f"{season}_t{thresh:g}"
    npz = OUT / f"season_{tag}.npz"
    js = OUT / f"season_{tag}.json"
    if npz.exists() and not rebuild:
        return json.loads(js.read_text())

    y0 = int(season.split("-")[0])
    start, end = f"{y0}-11-01", f"{y0+1}-03-31"
    season_zero = date(y0, 11, 1).toordinal()
    grid_tr, gw, gh, shelfmask = grid
    acc = CompositeAccumulator((gh, gw))
    wet_max_km2, cell_km2 = 0.0, (shelf.GRID_RES ** 2) / 1e6

    by_orbit = list_scenes(source, bbox, start, end)
    used = skipped = 0
    for orbit, items in sorted(by_orbit.items()):
        try:
            base = winter_baseline(y0, orbit, grid, source, bbox)
        except ValueError as e:
            print(f"  [skip orbit {orbit}] {e}")
            skipped += len(items)
            continue
        for it in items:
            try:
                db = read_db_on_grid(it, grid_tr, gw, gh, source)
            except Exception as e:
                print(f"  [skip scene] {type(e).__name__}: {str(e)[:50]}")
                skipped += 1
                continue
            wet, obs = wet_mask(db, base, thresh)
            day = it.datetime.date().toordinal() - season_zero
            acc.add(wet, obs, day)
            used += 1
            wet_max_km2 = max(wet_max_km2,
                              float((wet & shelfmask).sum()) * cell_km2)

    on = shelfmask & (acc.n_obs > 0)
    summary = {
        "season": season, "thresh_db": thresh,
        "scenes_used": used, "scenes_skipped": skipped,
        "wet_extent_max_km2": round(wet_max_km2, 1),
        "shelf_observed_km2": round(float(on.sum()) * cell_km2, 1),
        "median_first_wet_day": int(np.median(acc.first_wet[on & (acc.first_wet >= 0)]))
            if np.any(on & (acc.first_wet >= 0)) else None,
        "median_last_wet_day": int(np.median(acc.last_wet[on & (acc.last_wet >= 0)]))
            if np.any(on & (acc.last_wet >= 0)) else None,
    }
    np.savez_compressed(npz, wet_days=acc.wet_days, n_obs=acc.n_obs,
                        first_wet=acc.first_wet, last_wet=acc.last_wet)
    js.write_text(json.dumps(summary, indent=1))
    print(f"  [composite] {season} t={thresh}: {summary}")
    return summary

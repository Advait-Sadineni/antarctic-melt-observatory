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
SAR_SCALE = 4                         # SAR grid = fixed grid at 4x (120 m).
# Wet/dry is a regional signal; full 30 m warps would cost ~1 GB per scene
# across ~200 scenes. Same origin and shelf mask, so pixels stay comparable;
# the fusion phase may revisit resolution with better plumbing.
MAX_BASELINE_SCENES = 12              # evenly-spaced cap; median is stable by ~10
OUT = melt.ROOT / "output" / "sar"


def sar_grid(grid):
    """Derive the 120 m SAR grid from the fixed 30 m grid (same origin)."""
    grid_tr, gw, gh, shelfmask = grid
    tr = grid_tr * rasterio.Affine.scale(SAR_SCALE)
    gw4, gh4 = gw // SAR_SCALE, gh // SAR_SCALE
    sm = shelfmask[:gh4 * SAR_SCALE, :gw4 * SAR_SCALE]
    sm = sm.reshape(gh4, SAR_SCALE, gw4, SAR_SCALE).any(axis=(1, 3))
    return tr, gw4, gh4, sm


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


def save_checkpoint(path, acc, done_ids, used, skipped, wet_max):
    """Composite progress -> disk, atomically (write tmp, replace). A composite
    is hours of streamed scenes; interruptions (reboots, gaming pauses) must
    cost minutes, not the whole run."""
    tmp = path.with_suffix(".tmp.npz")
    np.savez_compressed(tmp, wet_days=acc.wet_days, n_obs=acc.n_obs,
                        first_wet=acc.first_wet, last_wet=acc.last_wet,
                        done_ids=np.array(sorted(done_ids)),
                        counters=np.array([used, skipped], "i8"),
                        wet_max=np.array([wet_max], "f8"))
    tmp.replace(path)


def load_checkpoint(path, shape):
    """Resume state from save_checkpoint, or fresh state if absent/corrupt.
    Returns (acc, done_ids, used, skipped, wet_max)."""
    if path.exists():
        try:
            z = np.load(path, allow_pickle=False)
            acc = CompositeAccumulator(shape)
            acc.wet_days, acc.n_obs = z["wet_days"], z["n_obs"]
            acc.first_wet, acc.last_wet = z["first_wet"], z["last_wet"]
            return (acc, set(z["done_ids"].tolist()),
                    int(z["counters"][0]), int(z["counters"][1]),
                    float(z["wet_max"][0]))
        except Exception as e:
            print(f"  [ckpt unreadable, starting fresh] {type(e).__name__}")
    return CompositeAccumulator(shape), set(), 0, 0, 0.0


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


def read_db_on_grid(item, grid_tr, gw, gh, source, retries=2):
    """gamma0 dB warped onto the fixed grid; NaN outside the swath.

    Retries once on transient blob/HTTP failures (a single flaky tile read
    must never kill a multi-hour run); raises only after all attempts fail."""
    last = None
    for _ in range(retries):
        try:
            with rasterio.open(source.band_href(item, "hh")) as src:
                with WarpedVRT(src, crs=shelf.GRID_CRS, transform=grid_tr,
                               width=gw, height=gh,
                               resampling=Resampling.average, nodata=0) as vrt:
                    g = vrt.read(1).astype("f4")
            return to_db(g)
        except Exception as e:      # rasterio wraps CPLE errors variously
            last = e
    raise last


def winter_baseline(year, orbit, grid, source, bbox, rebuild=False):
    """Per-pixel median dB over Jun-Aug of `year`, one orbit, disk-cached."""
    OUT.mkdir(parents=True, exist_ok=True)
    cache = OUT / f"baseline_{year}_{orbit}.npy"
    if cache.exists() and not rebuild:
        return np.load(cache)
    grid_tr, gw, gh, _ = grid
    scenes = list_scenes(source, bbox, f"{year}-06-01", f"{year}-08-31")
    items = scenes.get(orbit, [])
    if len(items) > MAX_BASELINE_SCENES:   # evenly spaced subsample, RAM-bounded
        idx = np.linspace(0, len(items) - 1, MAX_BASELINE_SCENES).astype(int)
        items = [items[i] for i in idx]
    stack = []
    for it in items:                        # per-scene tolerance: skip, never die
        try:
            stack.append(read_db_on_grid(it, grid_tr, gw, gh, source))
        except Exception as e:
            print(f"  [baseline skip] {type(e).__name__}: {str(e)[:60]}")
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
    cell_km2 = (shelf.GRID_RES ** 2) / 1e6
    ckpt = OUT / f"ckpt_{tag}.npz"
    acc, done, used, skipped, wet_max_km2 = load_checkpoint(ckpt, (gh, gw))
    if done:
        print(f"  [resume] {tag}: {len(done)} scenes already composited")

    by_orbit = list_scenes(source, bbox, start, end)
    since_save = 0
    for orbit, items in sorted(by_orbit.items()):
        items = [it for it in items if it.id not in done]
        if not items:
            continue
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
            done.add(it.id)
            wet_max_km2 = max(wet_max_km2,
                              float((wet & shelfmask).sum()) * cell_km2)
            since_save += 1
            if since_save >= 15:
                save_checkpoint(ckpt, acc, done, used, skipped, wet_max_km2)
                since_save = 0

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
    ckpt.unlink(missing_ok=True)
    print(f"  [composite] {season} t={thresh}: {summary}")
    return summary

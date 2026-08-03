"""M3 Phase B: fuse radar wetness with optical ponds into melt states.

State encoding (uint8): 0=UNOBSERVED, 1=DRY, 2=WET, 3=PONDED.
PONDED requires radar-wet AND optical pond evidence - by construction the
classifier cannot emit PONDED where radar says dry; optical-pond-but-radar-dry
is a counted CONFLICT (the diagnostic of both systems), never a state.
"""
import numpy as np

UNOBSERVED, DRY, WET, PONDED = 0, 1, 2, 3


def classify_state(wet, obs, pond_evidence):
    """Vectorized state table; see module docstring."""
    st = np.zeros(wet.shape, "u1")
    st[obs] = DRY
    st[obs & wet] = WET
    st[obs & wet & pond_evidence] = PONDED
    return st


def conflict_mask(wet, obs, pond_evidence):
    """Optical says pond, radar (observed) says dry - count, don't classify."""
    return pond_evidence & obs & ~wet


def confidence(margin_db, gap_days):
    """Per-cell confidence in [0,1]: radar margin below baseline (saturating at
    6 dB) times optical-evidence freshness (e-folding 6 days)."""
    m = np.clip(np.asarray(margin_db, "f4") / 6.0, 0.0, 1.0)
    g = np.exp(-np.asarray(gap_days, "f4") / 6.0)
    return (m * g).astype("f4")


def nearest_pond_evidence(dates, masks, scene_date, window_days=6):
    """The pond mask nearest in time to scene_date, all-False if none within
    the window. dates: list[date]; masks: list[bool ndarray] (same shapes)."""
    return nearest_pond_with_gap(dates, masks, scene_date, window_days)[0]


def nearest_pond_with_gap(dates, masks, scene_date, window_days=6):
    """(mask, gap_days) of the temporally nearest pond evidence; (all-False,
    None) when nothing lies within the window - confidence needs the gap."""
    best, best_gap = None, window_days + 1
    for d, m in zip(dates, masks):
        gap = abs((d - scene_date).days)
        if gap <= window_days and gap < best_gap:
            best, best_gap = m, gap
    if best is None:
        return np.zeros(masks[0].shape, bool), None
    return best, best_gap


def split_conflict(conflict, day, first_wet, last_wet, margin=6):
    """Partition conflicts by the cell's own radar wet season: inside
    [first_wet - margin, last_wet + margin] they are true sensor
    disagreements (the error-like class the pilot gates); outside they are
    the lid / diurnal-refreeze signal (optics sees stored water, radar sees
    frozen surface - real physics, reported separately, never classified).
    Cells never radar-wet (first_wet < 0) count as off-season entirely."""
    ever = first_wet >= 0
    in_window = ever & (first_wet - margin <= day) & (day <= last_wet + margin)
    wet_part = conflict & in_window
    return wet_part, conflict & ~wet_part


class StateAccumulator:
    """Per-cell season bookkeeping over classified scenes: state-day counts,
    conflict days, confidence sums, ponded phenology (-1 = never)."""

    def __init__(self, shape):
        self.days_wet = np.zeros(shape, "u2")      # WET or PONDED
        self.days_ponded = np.zeros(shape, "u2")
        self.n_obs = np.zeros(shape, "u2")
        self.conflict_days = np.zeros(shape, "u2")
        self.conf_sum = np.zeros(shape, "f4")
        self.conf_n = np.zeros(shape, "u2")
        self.first_ponded = np.full(shape, -1, "i2")
        self.last_ponded = np.full(shape, -1, "i2")

    def add(self, state, conflict, conf, day):
        wet = state >= WET
        ponded = state == PONDED
        self.days_wet += wet.astype("u2")
        self.days_ponded += ponded.astype("u2")
        self.n_obs += (state != UNOBSERVED).astype("u2")
        self.conflict_days += conflict.astype("u2")
        self.conf_sum[wet] += conf[wet]
        self.conf_n += wet.astype("u2")
        newly = ponded & (self.first_ponded < 0)
        self.first_ponded[newly] = day
        self.last_ponded[ponded] = day


_CKPT_FIELDS = ("days_wet", "days_ponded", "n_obs", "conflict_days",
                "conf_sum", "conf_n", "first_ponded", "last_ponded")


def save_state_checkpoint(path, acc, done_ids, used, skipped,
                          wet_max=0.0, pond_max=0.0, conf_wet=0, conf_off=0):
    """Same interruption insurance the SAR composites carry (atomic replace)."""
    tmp = path.with_suffix(".tmp.npz")
    np.savez_compressed(tmp,
                        done_ids=np.array(sorted(done_ids)),
                        counters=np.array([used, skipped, conf_wet, conf_off], "i8"),
                        extents=np.array([wet_max, pond_max], "f8"),
                        **{f: getattr(acc, f) for f in _CKPT_FIELDS})
    tmp.replace(path)


def load_state_checkpoint(path, shape):
    """Resume state, or fresh if absent/corrupt:
    (acc, done, used, skipped, wet_max, pond_max, conf_wet, conf_off)."""
    if path.exists():
        try:
            z = np.load(path, allow_pickle=False)
            acc = StateAccumulator(shape)
            for f in _CKPT_FIELDS:
                setattr(acc, f, z[f])
            c = z["counters"]
            return (acc, set(z["done_ids"].tolist()),
                    int(c[0]), int(c[1]),
                    float(z["extents"][0]), float(z["extents"][1]),
                    int(c[2]) if len(c) > 2 else 0,
                    int(c[3]) if len(c) > 3 else 0)
        except Exception as e:
            print(f"  [state ckpt unreadable, fresh start] {type(e).__name__}")
    return StateAccumulator(shape), set(), 0, 0, 0.0, 0.0, 0, 0


# --- networked (driven by scripts/fusion_pilot.py) -----------------------------

MAX_OPTICAL_PER_TILE = 6
POND_FRACTION = 0.25          # optical water fraction that makes a 120 m cell a pond
EVIDENCE_WINDOW_DAYS = 6


def pond_series(season, sgrid, cache_root=None):
    """Optical pond evidence for a season on the SAR grid: [(date, mask), ...]
    from every GVI tile's clean pool over Nov-Mar, fraction >= POND_FRACTION
    per 120 m cell (the Gate-1 rule), same-date masks OR'd. Cached npz."""
    from pathlib import Path

    import melt
    import sar
    import shelf
    root = Path(cache_root) if cache_root else sar.OUT
    root.mkdir(parents=True, exist_ok=True)
    cache = root / f"ponds_{season}_core.npz"
    if cache.exists():
        z = np.load(cache, allow_pickle=False)
        from datetime import date as _date
        return [(_date.fromordinal(int(o)), m)
                for o, m in zip(z["ordinals"], z["masks"])]

    from rasterio.enums import Resampling
    from rasterio.warp import reproject
    y0 = int(season.split("-")[0])
    # melt-core window only: outside it optical "ponds" are dominated by
    # frozen lids (stored water under ice - real, but not PONDING evidence;
    # decision 0004). The March-lid class produced 24.7% false conflict.
    start, end = f"{y0}-12-01", f"{y0+1}-02-28"
    tr, gw4, gh4, sm = sgrid
    by_date = {}
    for tile in shelf.SHELF_TILES:
        pool = shelf.production_pool(tile, start, end,
                                     max_scenes=MAX_OPTICAL_PER_TILE)
        for it in pool:
            bands = melt.load_scene(it, bands=("green", "nir", "red"))
            reject = melt.reject_mask(it) | melt.cloud_mask(it)["mask"]
            _, ponds, _ = melt.detect(bands["green"], bands["nir"], reject,
                                      red=bands["red"])
            src_crs, _ = melt.tile_georeference(it)
            dst = np.zeros((gh4, gw4), "f4")
            reproject(source=ponds.astype("f4"), destination=dst,
                      src_transform=melt.aoi_transform(it), src_crs=src_crs,
                      dst_transform=tr, dst_crs=shelf.GRID_CRS,
                      resampling=Resampling.average)
            d = it.datetime.date()
            m = (dst >= POND_FRACTION) & sm
            by_date[d] = (by_date[d] | m) if d in by_date else m
            print(f"  [ponds] {tile} {d}: {int(m.sum())} cells", flush=True)
    series = sorted(by_date.items())
    if series:
        np.savez_compressed(cache,
                            ordinals=np.array([d.toordinal() for d, _ in series]),
                            masks=np.stack([m for _, m in series]))
    return series


def melt_state_season(season, grid, source, bbox, thresh, series=None):
    """The M3 product: walk the season's S1 scenes, classify each cell
    UNOBSERVED/DRY/WET/PONDED against nearest optical evidence, accumulate.
    Writes output/sar/state_<season>.npz + JSON summary. Checkpointed."""
    import json
    from datetime import date as _date

    import sar
    import shelf
    sar.OUT.mkdir(parents=True, exist_ok=True)
    tag = f"{season}_t{thresh:g}"
    npz = sar.OUT / f"state_{tag}.npz"
    js = sar.OUT / f"state_{tag}.json"
    if npz.exists() and js.exists():   # both, or rebuild: a kill can land
        return json.loads(js.read_text())   # between the two writes

    sgrid = sar.sar_grid(grid)
    tr, gw4, gh4, sm = sgrid
    series = series if series is not None else pond_series(season, sgrid)
    dates = [d for d, _ in series]
    masks = [m for _, m in series]
    print(f"  [state] {season}: {len(series)} optical evidence dates", flush=True)

    y0 = int(season.split("-")[0])
    start, end = f"{y0}-11-01", f"{y0+1}-03-31"
    season_zero = _date(y0, 11, 1).toordinal()
    cell_km2 = (shelf.GRID_RES * sar.SAR_SCALE) ** 2 / 1e6

    ckpt = sar.OUT / f"ckpt_state_{tag}.npz"
    (acc, done, used, skipped, wet_max, pond_max,
     conf_wet, conf_off) = load_state_checkpoint(ckpt, (gh4, gw4))
    if done:
        print(f"  [resume] state {tag}: {len(done)} scenes done", flush=True)

    # composite phenology for conflict phase-splitting (decision 0004)
    comp = np.load(sar.OUT / f"season_{tag}.npz")
    fw_c, lw_c = comp["first_wet"], comp["last_wet"]

    by_orbit = sar.list_scenes(source, bbox, start, end)
    since = 0
    for orbit, items in sorted(by_orbit.items()):
        items = [it for it in items if it.id not in done]
        if not items:
            continue
        try:
            base = sar.winter_baseline(y0, orbit, sgrid, source, bbox)
        except ValueError as e:
            print(f"  [skip orbit {orbit}] {e}")
            skipped += len(items)
            continue
        for it in items:
            try:
                db = sar.read_db_on_grid(it, tr, gw4, gh4, source)
            except Exception as e:
                print(f"  [skip scene] {type(e).__name__}: {str(e)[:50]}")
                skipped += 1
                continue
            wet, obs = sar.wet_mask(db, base, thresh)
            d = it.datetime.date()
            ev, gap = (nearest_pond_with_gap(dates, masks, d,
                                             EVIDENCE_WINDOW_DAYS)
                       if dates else (np.zeros((gh4, gw4), bool), None))
            st = classify_state(wet, obs, ev)
            with np.errstate(invalid="ignore"):
                margin = np.where(obs, base - db, 0.0)
            conf = np.clip(margin / 6.0, 0.0, 1.0).astype("f4")
            if gap is not None:
                pondcells = st == PONDED
                conf[pondcells] *= np.float32(np.exp(-gap / 6.0))
            cf = conflict_mask(wet, obs, ev)
            day_n = d.toordinal() - season_zero
            cw, co = split_conflict(cf, day_n, fw_c, lw_c)
            conf_wet += int(cw.sum())
            conf_off += int(co.sum())
            acc.add(st, cf, conf, day=day_n)
            used += 1
            done.add(it.id)
            wet_max = max(wet_max, float((st >= WET)[sm].sum()) * cell_km2)
            pond_max = max(pond_max, float((st == PONDED)[sm].sum()) * cell_km2)
            since += 1
            if since >= 15:
                save_state_checkpoint(ckpt, acc, done, used, skipped,
                                      wet_max, pond_max, conf_wet, conf_off)
                since = 0

    on = sm & (acc.n_obs > 0)
    ponded_ever = on & (acc.days_ponded > 0)
    ponded_days_total = int(acc.days_ponded[on].sum())
    wet_instances = conf_wet + ponded_days_total
    summary = {
        "season": season, "thresh_db": thresh,
        "scenes_used": used, "scenes_skipped": skipped,
        "optical_evidence_dates": len(series),
        "wet_extent_max_km2": round(wet_max, 1),
        "ponded_extent_max_km2": round(pond_max, 1),
        "ponded_ever_km2": round(float(ponded_ever.sum()) * cell_km2, 1),
        # gate class: disagreement while the cell was radar-wet-capable
        "conflict_rate": round(conf_wet / wet_instances, 4)
            if wet_instances else None,
        "conflict_days_wetseason": conf_wet,
        # lids / diurnal refreeze: optics sees stored water, radar frozen
        "conflict_days_offseason": conf_off,
        "mean_confidence_wet": round(float(
            (acc.conf_sum[on] / np.maximum(acc.conf_n[on], 1)).mean()), 3),
        "median_first_ponded_day": int(np.median(
            acc.first_ponded[on & (acc.first_ponded >= 0)]))
            if ponded_ever.any() else None,
    }
    np.savez_compressed(npz, **{f: getattr(acc, f) for f in _CKPT_FIELDS})
    js.write_text(json.dumps(summary, indent=1))
    ckpt.unlink(missing_ok=True)
    print(f"  [state] {season}: {summary}", flush=True)
    return summary

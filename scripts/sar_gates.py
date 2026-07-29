"""M3 validation gates (spec section 5) - networked, run in background.

    python scripts/sar_gates.py gate14     # baselines + Gate 1 sweep + Gate 4
    python scripts/sar_gates.py gate23     # 2019-20 & 2021-22 -> ranking + winter

Gate 1 (containment): >=90% of the blind-validated 19DEA pond pixels
(2021-01-24 scene) must be radar-WET on passes within +/-6 days.
Gate 4 (timing): median onset Nov-Dec, end Feb-Mar.
Verdicts accumulate in output/sar/gates.json.
"""
import json
import sys
from datetime import date, timedelta
from pathlib import Path

import numpy as np
from rasterio.enums import Resampling
from rasterio.warp import reproject

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
import melt
import sar
import shelf

GVI_BBOX = [-70.5, -72.6, -66.0, -70.0]
VALID_DATE = date(2021, 1, 24)          # blind-validated 19DEA scene date
GATES = sar.OUT / "gates.json"


def _save(key, value):
    sar.OUT.mkdir(parents=True, exist_ok=True)
    g = json.loads(GATES.read_text()) if GATES.exists() else {}
    g[key] = value
    GATES.write_text(json.dumps(g, indent=1))
    print(f"[gates] {key} = {value}")


def _validated_ponds_on_sargrid(sgrid):
    """The blind-validated 19DEA detection (2021-01-24) as a mask on the SAR
    grid: optical water fraction >= 0.25 of a 120 m cell."""
    it = melt.get_source().search(
        query=melt.get_source().tile_query("19DEA"),
        datetime="2021-01-24/2021-01-25")[0]
    row, col, size = shelf._tile_window(it)
    melt.set_aoi("19DEA", row, col, size)
    bands = melt.load_scene(it, bands=("green", "nir", "red"))
    reject = melt.reject_mask(it) | melt.cloud_mask(it)["mask"]
    _, ponds, _ = melt.detect(bands["green"], bands["nir"], reject, red=bands["red"])

    tr, gw4, gh4, sm = sgrid
    src_crs, _ = melt.tile_georeference(it)
    dst = np.zeros((gh4, gw4), "f4")
    reproject(source=ponds.astype("f4"), destination=dst,
              src_transform=melt.aoi_transform(it), src_crs=src_crs,
              dst_transform=tr, dst_crs=shelf.GRID_CRS,
              resampling=Resampling.average)
    return (dst >= 0.25) & sm, it


def gate14():
    grid = shelf.build_fixed_grid()
    sgrid = sar.sar_grid(grid)
    src = sar.s1_source()
    tr, gw4, gh4, sm = sgrid

    ponds, _ = _validated_ponds_on_sargrid(sgrid)
    print(f"[gate1] validated pond cells on SAR grid: {int(ponds.sum())}")

    # S1 passes within +/-6 days of the validated optical date
    t0, t1 = VALID_DATE - timedelta(days=6), VALID_DATE + timedelta(days=6)
    by_orbit = sar.list_scenes(src, GVI_BBOX, t0.isoformat(), t1.isoformat())
    print(f"[gate1] orbits in window: { {k: len(v) for k, v in by_orbit.items()} }")

    sweep = {}
    for thresh in (2.0, 3.0, 4.0):
        wet_any = np.zeros((gh4, gw4), bool)
        obs_any = np.zeros((gh4, gw4), bool)
        for orbit, items in by_orbit.items():
            try:
                base = sar.winter_baseline(2020, orbit, sgrid, src, GVI_BBOX)
            except ValueError as e:
                print(f"  [gate1 skip orbit {orbit}] {e}")
                continue
            for it in items:
                try:
                    db = sar.read_db_on_grid(it, tr, gw4, gh4, src)
                except Exception as e:
                    print(f"  [gate1 skip scene] {type(e).__name__}: {str(e)[:50]}")
                    continue
                wet, obs = sar.wet_mask(db, base, thresh)
                wet_any |= wet
                obs_any |= obs
        seen = ponds & obs_any
        frac = float((ponds & wet_any).sum() / max(seen.sum(), 1))
        sweep[str(thresh)] = {"containment": round(frac, 3),
                              "pond_cells_observed": int(seen.sum())}
        print(f"  [gate1] thresh {thresh} dB -> containment {frac:.3f}")

    best = max(sweep, key=lambda k: sweep[k]["containment"])
    _save("gate1_containment_sweep", sweep)
    _save("gate1_pass", sweep[best]["containment"] >= 0.90)
    _save("gate1_best_thresh_db", float(best))

    # Gate 4: full 2020-21 composite at the best threshold
    s = sar.season_composite("2020-21", sgrid, src, GVI_BBOX,
                             thresh=float(best))
    onset, endd = s["median_first_wet_day"], s["median_last_wet_day"]
    ok = (onset is not None and endd is not None
          and 0 <= onset <= 61 and 90 <= endd <= 151)   # Nov-Dec / Feb-Mar
    _save("gate4_timing", {"median_onset_day": onset, "median_end_day": endd,
                           "pass": bool(ok)})


def gate23():
    grid = shelf.build_fixed_grid()
    sgrid = sar.sar_grid(grid)
    src = sar.s1_source()
    g = json.loads(GATES.read_text())
    thresh = g.get("gate1_best_thresh_db", 3.0)

    ext = {}
    for season in ("2019-20", "2021-22"):
        ext[season] = sar.season_composite(season, sgrid, src, GVI_BBOX,
                                           thresh=thresh)["wet_extent_max_km2"]
    ext["2020-21"] = json.loads(
        (sar.OUT / f"season_2020-21_t{thresh:g}.json").read_text()
    )["wet_extent_max_km2"]
    ok = ext["2019-20"] > ext["2020-21"] > ext["2021-22"]
    _save("gate2_ranking", {"wet_extent_max_km2": ext, "pass": bool(ok)})

    # Gate 3: deep-winter sanity (Jun-Jul 2021), same baselines year 2021
    tr, gw4, gh4, sm = sgrid
    by_orbit = sar.list_scenes(src, GVI_BBOX, "2021-06-15", "2021-07-15")
    wet_frac = []
    for orbit, items in by_orbit.items():
        try:
            base = sar.winter_baseline(2021, orbit, sgrid, src, GVI_BBOX)
        except ValueError:
            continue
        for it in items[:3]:
            try:
                db = sar.read_db_on_grid(it, tr, gw4, gh4, src)
            except Exception:
                continue
            wet, obs = sar.wet_mask(db, base, thresh)
            on = sm & obs
            if on.sum():
                wet_frac.append(float((wet & sm).sum() / on.sum()))
    mid = float(np.median(wet_frac)) if wet_frac else None
    _save("gate3_winter", {"median_midwinter_wet_frac": mid,
                           "pass": bool(mid is not None and mid < 0.02)})


if __name__ == "__main__":
    {"gate14": gate14, "gate23": gate23}[sys.argv[1] if len(sys.argv) > 1 else "gate14"]()

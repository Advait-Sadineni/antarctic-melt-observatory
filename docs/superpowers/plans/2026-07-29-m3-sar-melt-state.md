# M3 SAR Melt-State Implementation Plan (Phase A: detection + gates)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Sentinel-1 wet-snow detection on the George VI grid (per-pixel winter
baseline, −3 dB rule), proven against the four validation gates from the M3 spec.

**Architecture:** A `Sentinel1Source` (M1 adapter subclass) feeds `src/sar.py`:
scenes grouped by relative orbit → per-pixel winter baselines (median dB, per
orbit, cached) → per-scene wet masks on the fixed EPSG:3031 GVI grid →
season composites (wet-extent, onset/end/duration). Gates run as scripts that
compare against the blind-validated optical record. **Fusion (spec §3) is
deliberately a SEPARATE follow-on plan**: its parameters (final threshold from
the 2–4 dB sweep, compositing bin width) are OUTPUTS of this phase's gates —
planning it before the gates would be fiction.

**Tech Stack:** Python 3.13, rasterio (WarpedVRT reprojection), numpy, pystac-client,
planetary-computer, pytest.

## Global Constraints

- 59 offline tests stay green after every task; new SAR math gets offline tests
  with synthetic arrays (no network in pytest).
- The fixed GVI grid (output/shelf/grid_fixed.npz) is the one grid; nothing
  bespoke. Region-agnostic: bbox/grid are parameters, never constants.
- Backscatter math in dB: `db = 10*log10(gamma0)`, guard `gamma0 <= 0` → NaN.
- Baselines are per (winter, relative_orbit); never mix orbits.
- Networked gate runs are background scripts writing logs + JSON verdicts under
  output/sar/; never inside pytest.
- George VI numbers stay frozen; SAR must never alter optical outputs.

## File Structure

- `src/core/sources.py` — add `Sentinel1Source` (10 lines; adapter seam)
- `src/sar.py` — all SAR science: scene listing, dB read-on-grid, baselines,
  wet masks, season composite
- `scripts/sar_gates.py` — gate runner (containment, ranking, winter, timing)
- `tests/test_sar.py` — offline synthetic tests
- `output/sar/` — baselines cache, composites, gate verdicts (gitignored)

---

### Task 1: Sentinel1Source adapter

**Files:** Modify `src/core/sources.py`; test `tests/test_sar.py`

**Interfaces:**
- Produces: `class Sentinel1Source(PlanetaryComputerSource)` importable from
  `core.sources` (re-export) with `COLLECTION="sentinel-1-rtc"`,
  `BAND_MAP={"hh": "hh"}`; `relative_orbit(item) -> int` reading
  `item.properties["sat:relative_orbit"]`.

- [ ] **Step 1: failing test**

```python
# tests/test_sar.py
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
import numpy as np
from core.sources import Sentinel1Source


class FakeS1Item:
    def __init__(self, orbit=65):
        self.id = "S1B_IW_GRDH_1SSH_20210131T073912_rtc"
        self.properties = {"sat:relative_orbit": orbit}
        class A:  # noqa
            href = "https://pc/hh.tif"
        self.assets = {"hh": A()}


def test_s1_source_band_and_orbit():
    src = Sentinel1Source.__new__(Sentinel1Source)
    it = FakeS1Item(orbit=65)
    assert src.band_href(it, "hh") == "https://pc/hh.tif"
    assert src.relative_orbit(it) == 65
    assert Sentinel1Source.COLLECTION == "sentinel-1-rtc"
```

- [ ] **Step 2:** `python -m pytest tests/test_sar.py -q` → FAIL (ImportError)
- [ ] **Step 3: implement** — in `src/core/pc.py`:

```python
class Sentinel1Source(PlanetaryComputerSource):
    """Sentinel-1 RTC (terrain-corrected gamma0). Not MGRS-tiled: search by
    bbox, group by relative orbit (radar geometry differs per track)."""
    COLLECTION = "sentinel-1-rtc"
    BAND_MAP = {"hh": "hh"}

    def relative_orbit(self, item):
        return int(item.properties["sat:relative_orbit"])
```

and in `src/core/sources.py` add at the bottom:

```python
def sentinel1_source():
    from .pc import Sentinel1Source
    return Sentinel1Source()
```

(the test imports `Sentinel1Source` from `core.sources`; add
`from .pc import Sentinel1Source  # noqa: F401` guarded in a function-level
import if circularity bites — pc imports sources, so re-export must be lazy:
implement test import target as `from core.pc import Sentinel1Source` instead
if needed and adjust the test import line accordingly.)

- [ ] **Step 4:** pytest tests/test_sar.py + full suite → PASS
- [ ] **Step 5:** commit `feat(m3): Sentinel1Source adapter`

### Task 2: dB read on the fixed grid

**Files:** Create `src/sar.py`; test `tests/test_sar.py`

**Interfaces:**
- Produces: `to_db(gamma0: ndarray) -> ndarray` (NaN where ≤0);
  `read_db_on_grid(item, grid_tr, gw, gh, source) -> ndarray` (float32 dB on
  grid, NaN outside swath) using `rasterio.vrt.WarpedVRT` against
  `source.band_href(item, "hh")`.

- [ ] **Step 1: failing test (to_db only — read is networked, not unit-tested)**

```python
def test_to_db_handles_zeros_and_scales():
    from sar import to_db
    g = np.array([[1.0, 0.1], [0.0, -1.0]], "f4")
    db = to_db(g)
    assert db[0, 0] == 0.0 and abs(db[0, 1] + 10.0) < 1e-5
    assert np.isnan(db[1, 0]) and np.isnan(db[1, 1])
```

- [ ] **Step 2:** FAIL (no module sar)
- [ ] **Step 3: implement `src/sar.py`** (module header + first functions):

```python
"""Sentinel-1 wet-snow detection for the observatory (M3 spec).

Wet snow collapses radar backscatter. Melt = gamma0 more than THRESH_DB below
that pixel's frozen-winter baseline, built per relative orbit."""
import numpy as np
import rasterio
from rasterio.vrt import WarpedVRT
from rasterio.enums import Resampling

import melt
import shelf

THRESH_DB = 3.0
OUT = melt.ROOT / "output" / "sar"


def to_db(gamma0):
    g = np.asarray(gamma0, "f4")
    with np.errstate(divide="ignore", invalid="ignore"):
        db = 10.0 * np.log10(np.where(g > 0, g, np.nan))
    return db.astype("f4")


def read_db_on_grid(item, grid_tr, gw, gh, source):
    """gamma0 dB warped onto the fixed grid; NaN outside the swath."""
    with rasterio.open(source.band_href(item, "hh")) as src:
        with WarpedVRT(src, crs=shelf.GRID_CRS, transform=grid_tr,
                       width=gw, height=gh,
                       resampling=Resampling.average, nodata=0) as vrt:
            g = vrt.read(1).astype("f4")
    return to_db(g)          # zeros (nodata) become NaN via to_db
```

- [ ] **Step 4:** pytest → PASS; commit `feat(m3): dB read on fixed grid`

### Task 3: winter baseline (per orbit, cached)

**Interfaces:**
- Produces: `winter_baseline(year, orbit, grid, source, bbox) -> ndarray`
  (per-pixel MEDIAN dB over Jun–Aug of `year`, that orbit only; cached to
  `output/sar/baseline_<year>_<orbit>.npy`); `list_scenes(source, bbox, start,
  end) -> dict[orbit, list[item]]`.

- [ ] **Step 1: failing test (median math + cache path, synthetic)**

```python
def test_baseline_median_ignores_nan(tmp_path, monkeypatch):
    import sar
    stack = [np.array([[1.0, np.nan]], "f4"),
             np.array([[3.0, 5.0]], "f4"),
             np.array([[2.0, np.nan]], "f4")]
    out = sar._median_stack(stack)
    assert out[0, 0] == 2.0 and out[0, 1] == 5.0
```

- [ ] **Step 2:** FAIL → **Step 3: implement**

```python
def _median_stack(arrays):
    return np.nanmedian(np.stack(arrays), axis=0).astype("f4")


def list_scenes(source, bbox, start, end):
    items = source.search(collections=[source.COLLECTION], bbox=bbox,
                          datetime=f"{start}/{end}", limit=100)
    by_orbit = {}
    for it in items:
        by_orbit.setdefault(source.relative_orbit(it), []).append(it)
    return by_orbit


def winter_baseline(year, orbit, grid, source, bbox, rebuild=False):
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
    return base
```

- [ ] **Step 4:** pytest full suite → PASS; commit `feat(m3): per-orbit winter baselines`

### Task 4: wet masks + season composite

**Interfaces:**
- Produces: `wet_mask(db, baseline, thresh=THRESH_DB) -> ndarray(bool)` (NaN-
  safe: NaN in either → False, tracked separately as observed mask);
  `season_composite(season, grid, source, bbox, thresh) -> dict` — iterates
  Nov 1 → Mar 31 scenes per orbit, accumulates per-pixel `wet_days`
  (count of wet observations), `n_obs`, `first_wet`/`last_wet` (day-of-season
  ints), saves `output/sar/season_<season>.npz` + JSON summary
  (wet-extent max km² on shelf mask, onset/end medians).

- [ ] **Step 1: failing tests (mask logic + composite math, synthetic)**

```python
def test_wet_mask_nan_safe():
    import sar
    db = np.array([[-20.0, -10.0], [np.nan, -14.0]], "f4")
    base = np.array([[-15.0, -12.0], [-13.0, np.nan]], "f4")
    wet, obs = sar.wet_mask(db, base)
    assert wet[0, 0] and not wet[0, 1]          # -20 < -18 yes; -10 > -15 no
    assert not wet[1, 0] and not wet[1, 1]      # NaN never wet
    assert obs[0, 0] and obs[0, 1] and not obs[1, 0] and not obs[1, 1]
```

- [ ] **Step 2:** FAIL → **Step 3: implement** (wet_mask + composite loop; the
  composite is a plain accumulation loop over `list_scenes` output, per orbit,
  using that orbit's baseline; day-of-season = (scene_date − Nov 1).days):

```python
def wet_mask(db, baseline, thresh=THRESH_DB):
    obs = ~np.isnan(db) & ~np.isnan(baseline)
    wet = np.zeros(db.shape, bool)
    wet[obs] = db[obs] < (baseline[obs] - thresh)
    return wet, obs
```

(season_composite: allocate int16 `first_wet` init −1, `last_wet` init −1,
uint16 `wet_days`, `n_obs`; for each scene in date order compute wet/obs,
`wet_days += wet`, `n_obs += obs`, set `first_wet` where `wet & (first_wet<0)`,
update `last_wet` where wet; after loop save npz + summary JSON with shelf-
masked wet-extent-max = max over scenes of wet&shelf area — track running max
during the loop.)

- [ ] **Step 4:** pytest → PASS; commit `feat(m3): wet masks + season composites`

### Task 5: gate runner + Gate 1/4 on 2020-21 (networked)

**Files:** Create `scripts/sar_gates.py`

- [ ] **Step 1: write the runner** — loads the fixed grid; builds 2020 winter
  baselines for every orbit seen; runs `season_composite("2020-21", ...)`;
  **Gate 1 (containment):** loads the blind-validated 19DEA pond pixels
  (regenerate detection for the validated scene via the optical pipeline,
  reproject to grid like tile_water_on_grid does), then reports the fraction
  of pond pixels with ≥1 wet observation within ±6 days of 2021-01-24;
  PASS ≥ 0.90. Runs the 2/3/4 dB sweep and reports each. **Gate 4 (timing):**
  median first_wet in Nov–Dec, median last_wet in Feb–Mar, wet fields
  spatially coherent (majority-filter change < 10%).
  Writes `output/sar/gates.json` verdicts.
- [ ] **Step 2: run in background**, review verdicts, iterate threshold choice
  from the sweep (this is the spec's declared tunable).
- [ ] **Step 3:** commit `feat(m3): gate runner; gate 1+4 verdicts`

### Task 6: Gates 2 + 3 (ranking + winter sanity)

- [ ] **Step 1:** composites for 2019-20 and 2021-22 (+ 2019 baselines);
  Gate 2: wet-extent-max ranks 2019-20 ≫ 2020-21 > 2021-22.
- [ ] **Step 2:** mid-winter composite (Jun–Jul 2021): wet fraction of shelf
  < 2% away from documented winter-storage areas; qualitative check of the
  Frontiers-2025 winter anomaly (visual figure to output/sar/).
- [ ] **Step 3:** all four verdicts into gates.json; commit
  `feat(m3): gates 2+3 — SAR detection validated on George VI`.

**Phase B (fusion + confidence + products) gets its own plan once gates.json
is green — its parameters are these gates' outputs.**

## Self-review notes
- Spec coverage: §2 detection → T1–T4; §5 gates → T5–T6; §3 fusion explicitly
  deferred to Phase B with reason; §4 scope enforced (GVI, 3 seasons).
- Types consistent: baselines/db float32 NaN-coded; masks bool; composite npz.
- No placeholders: every offline step has code; networked gate steps specify
  exact acceptance numbers.

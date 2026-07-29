# M1 Foundation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the validated detection pipeline backend-agnostic (Earth Search today, Microsoft Planetary Computer as the new cloud backend), pouring results into the observatory's day-one data model (COG + STAC), with George VI's frozen numbers as the regression gate.

**Architecture:** A tiny `SceneSource` adapter (search + band-href + tile-of) is introduced behind the existing code; `melt.py`/`shelf.py` keep their public functions so all 51 offline tests and the running scripts stay untouched. A `PlanetaryComputerSource` implements the same adapter with PC's asset names and signing. `run_season` gains a product writer that emits a water-fraction COG + STAC-style item JSON per (shelf, season) with reserved fields for M3/M4 (sensor confidence, depth, volume, uncertainty).

**Tech Stack:** Python 3.13, rasterio, pystac-client, planetary-computer (new dep), numpy, pytest.

## Global Constraints

- The 51 offline tests must pass after every task (`python -m pytest tests -q`).
- George VI numbers are frozen ground truth: 2020-21 = 142.48 km² (blind-validated), 2019-20 = 840.5, 2024-25 = 97.9, 2025-26 = 115.3. Any behavior change that moves them is a bug.
- Never modify running background processes' outputs (`output/shelf/history*.json` are written by them); new products go to `output/products/`.
- All new storage is open standards: COG for arrays, STAC-style JSON for metadata.
- Region-agnostic: nothing new may hard-code a shelf, tile list, or bbox (the `_tile_window` incident).
- Logical band names everywhere in core code: `red, green, blue, nir, swir16, scl`. Only adapters know provider asset keys.
- If `planetary-computer` is not installed: `pip install planetary-computer` (proactively, per user preference).

## File Structure

- `src/core/__init__.py` — package marker (empty)
- `src/core/sources.py` — `SceneSource` base + `EarthSearchSource` (current behavior, extracted)
- `src/core/pc.py` — `PlanetaryComputerSource`
- `src/core/products.py` — COG + STAC item writer (the data model)
- `src/melt.py` — gains `SOURCE` global + `set_source()`; asset/`Client.open` sites delegate to it
- `src/shelf.py`, `src/shelves.py` — their `Client.open` sites delegate to `melt.SOURCE.search`
- `scripts/regression_pc.py` — manual PC-vs-EarthSearch regression gate (networked, not pytest)
- `docs/decisions/0001-compute-backend.md` — decision log entry
- `tests/test_sources.py` — offline adapter tests (FakeItem pattern from `tests/test_melt.py`)

---

### Task 1: SceneSource + EarthSearchSource (extract current behavior)

**Files:**
- Create: `src/core/__init__.py`, `src/core/sources.py`
- Test: `tests/test_sources.py`

**Interfaces:**
- Produces (later tasks rely on these exact signatures):
  - `class SceneSource:` with methods
    `search(self, collections, query=None, datetime=None, bbox=None, ids=None, limit=100) -> list` (thin pystac wrapper returning items),
    `band_href(self, item, band: str) -> str`,
    `tile_of(self, item) -> str` (MGRS tile, e.g. `"19DEA"`).
  - `class EarthSearchSource(SceneSource):` with `STAC_API = "https://earth-search.aws.element84.com/v1"`; `BAND_MAP = {"red": "red", "green": "green", "blue": "blue", "nir": "nir", "swir16": "swir16", "scl": "scl"}`.

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_sources.py
"""Adapter tests - offline, FakeItem pattern from test_melt.py."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from core.sources import EarthSearchSource  # noqa: E402


class FakeAsset:
    def __init__(self, href):
        self.href = href


class FakeESItem:
    """Earth Search style: id S2B_19DEA_20210124_0_L2A, friendly asset keys."""
    def __init__(self):
        self.id = "S2B_19DEA_20210124_0_L2A"
        self.properties = {}
        self.assets = {k: FakeAsset(f"https://es/{k}.tif")
                       for k in ("red", "green", "blue", "nir", "swir16", "scl")}


def test_earthsearch_band_href_uses_friendly_names():
    src = EarthSearchSource.__new__(EarthSearchSource)  # no network in __init__
    it = FakeESItem()
    assert src.band_href(it, "green") == "https://es/green.tif"
    assert src.band_href(it, "swir16") == "https://es/swir16.tif"


def test_earthsearch_tile_of_parses_id():
    src = EarthSearchSource.__new__(EarthSearchSource)
    assert src.tile_of(FakeESItem()) == "19DEA"
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `python -m pytest tests/test_sources.py -q`
Expected: FAIL with `ModuleNotFoundError: No module named 'core'`

- [ ] **Step 3: Write the implementation**

```python
# src/core/__init__.py  (empty file)
```

```python
# src/core/sources.py
"""Scene-source adapters: the ONLY place that knows a provider's STAC API,
asset keys, or item-id conventions. Core code speaks logical band names
(red/green/blue/nir/swir16/scl) and calls search/band_href/tile_of."""


class SceneSource:
    STAC_API = None          # subclass sets
    COLLECTION = "sentinel-2-l2a"
    BAND_MAP = {}            # logical band -> provider asset key

    def __init__(self):
        from pystac_client import Client
        self.client = Client.open(self.STAC_API)

    def search(self, collections=None, query=None, datetime=None,
               bbox=None, ids=None, limit=100):
        kw = {"collections": collections or [self.COLLECTION], "limit": limit}
        if query is not None:
            kw["query"] = query
        if datetime is not None:
            kw["datetime"] = datetime
        if bbox is not None:
            kw["bbox"] = bbox
        if ids is not None:
            kw["ids"] = ids
        return list(self.client.search(**kw).items())

    def band_href(self, item, band):
        return item.assets[self.BAND_MAP[band]].href

    def tile_of(self, item):
        raise NotImplementedError

    def tile_query(self, tile):
        """Provider-specific STAC query fragment selecting one MGRS tile."""
        raise NotImplementedError


class EarthSearchSource(SceneSource):
    STAC_API = "https://earth-search.aws.element84.com/v1"
    BAND_MAP = {b: b for b in ("red", "green", "blue", "nir", "swir16", "scl")}

    def tile_of(self, item):
        # id like S2B_19DEA_20210124_0_L2A
        return item.id.split("_")[1]

    def tile_query(self, tile):
        return {"grid:code": {"eq": f"MGRS-{tile}"}}
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `python -m pytest tests/test_sources.py tests/test_melt.py -q`
Expected: all PASS (53)

- [ ] **Step 5: Commit**

```bash
git add src/core tests/test_sources.py
git commit -m "feat(m1): SceneSource adapter + EarthSearchSource (extracted current behavior)"
```

---

### Task 2: PlanetaryComputerSource

**Files:**
- Create: `src/core/pc.py`
- Test: `tests/test_sources.py` (append)

**Interfaces:**
- Consumes: `SceneSource` from Task 1.
- Produces: `class PlanetaryComputerSource(SceneSource)` with
  `STAC_API = "https://planetarycomputer.microsoft.com/api/stac/v1"`,
  `BAND_MAP = {"red": "B04", "green": "B03", "blue": "B02", "nir": "B08", "swir16": "B11", "scl": "SCL"}`,
  `tile_of()` reading `properties["s2:mgrs_tile"]`,
  `tile_query(tile)` returning `{"s2:mgrs_tile": {"eq": tile}}`,
  and hrefs signed via `planetary_computer.sign_inplace`.

- [ ] **Step 1: Ensure the dependency**

Run: `python -c "import planetary_computer; print(planetary_computer.__version__)"`
If it fails: `pip install planetary-computer` (then rerun to confirm).

- [ ] **Step 2: Write the failing tests**

```python
# append to tests/test_sources.py
from core.pc import PlanetaryComputerSource  # noqa: E402


class FakePCItem:
    """PC style: id S2B_MSIL2A_..._T19DEA_..., B-number asset keys."""
    def __init__(self):
        self.id = "S2B_MSIL2A_20210124T120000_R000_T19DEA_20210124T150000"
        self.properties = {"s2:mgrs_tile": "19DEA"}
        self.assets = {k: FakeAsset(f"https://pc/{k}.tif")
                       for k in ("B02", "B03", "B04", "B08", "B11", "SCL")}


def test_pc_band_href_maps_logical_to_b_numbers():
    src = PlanetaryComputerSource.__new__(PlanetaryComputerSource)
    it = FakePCItem()
    assert src.band_href(it, "green") == "https://pc/B03.tif"
    assert src.band_href(it, "scl") == "https://pc/SCL.tif"


def test_pc_tile_of_uses_property_not_id():
    src = PlanetaryComputerSource.__new__(PlanetaryComputerSource)
    assert src.tile_of(FakePCItem()) == "19DEA"


def test_pc_tile_query_shape():
    src = PlanetaryComputerSource.__new__(PlanetaryComputerSource)
    assert src.tile_query("19DEA") == {"s2:mgrs_tile": {"eq": "19DEA"}}
```

- [ ] **Step 3: Run tests to verify they fail**

Run: `python -m pytest tests/test_sources.py -q`
Expected: FAIL with `ModuleNotFoundError: No module named 'core.pc'`

- [ ] **Step 4: Write the implementation**

```python
# src/core/pc.py
"""Microsoft Planetary Computer adapter. Same Sentinel-2 archive, different
asset keys and item ids, and hrefs must be SAS-signed before rasterio reads
them - sign_inplace on the client handles that transparently."""
from .sources import SceneSource


class PlanetaryComputerSource(SceneSource):
    STAC_API = "https://planetarycomputer.microsoft.com/api/stac/v1"
    BAND_MAP = {"red": "B04", "green": "B03", "blue": "B02",
                "nir": "B08", "swir16": "B11", "scl": "SCL"}

    def __init__(self):
        import planetary_computer
        from pystac_client import Client
        self.client = Client.open(self.STAC_API,
                                  modifier=planetary_computer.sign_inplace)

    def tile_of(self, item):
        return item.properties["s2:mgrs_tile"]

    def tile_query(self, tile):
        return {"s2:mgrs_tile": {"eq": tile}}
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `python -m pytest tests/test_sources.py -q`
Expected: all PASS

- [ ] **Step 6: Commit**

```bash
git add src/core/pc.py tests/test_sources.py
git commit -m "feat(m1): PlanetaryComputerSource adapter (signed hrefs, B-number bands)"
```

---

### Task 3: Route melt/shelf/shelves through the active source

**Files:**
- Modify: `src/melt.py` (asset sites at lines ~193, ~274, ~344, ~355, ~464; `_tile_of` at ~180; `Client.open` at ~228, ~246)
- Modify: `src/shelf.py` (`Client.open` at ~161, ~452, ~488), `src/shelves.py` (~68)
- Test: `tests/test_sources.py` (append)

**Interfaces:**
- Consumes: `EarthSearchSource` (Task 1).
- Produces: `melt.SOURCE` (module global, default `EarthSearchSource()` built lazily), `melt.set_source(source)`, `melt.get_source()`.

- [ ] **Step 1: Write the failing test**

```python
# append to tests/test_sources.py
import melt  # noqa: E402


class RecordingSource:
    BAND_MAP = {b: b for b in ("red", "green", "blue", "nir", "swir16", "scl")}
    def __init__(self):
        self.calls = []
    def band_href(self, item, band):
        self.calls.append(("band_href", band))
        return item.assets[band].href
    def tile_of(self, item):
        return "19DEA"
    def search(self, **kw):
        self.calls.append(("search", kw)); return []
    def tile_query(self, tile):
        return {"grid:code": {"eq": f"MGRS-{tile}"}}


def test_set_source_swaps_and_tile_of_delegates():
    rec = RecordingSource()
    old = melt.get_source()
    try:
        melt.set_source(rec)
        assert melt.get_source() is rec
        assert melt._tile_of(FakeESItem()) == "19DEA"   # delegated, not id-parsed
    finally:
        melt.set_source(old)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest tests/test_sources.py -q`
Expected: FAIL with `AttributeError: module 'melt' has no attribute 'get_source'`

- [ ] **Step 3: Implement the routing**

In `src/melt.py` (near the STAC_API constant):

```python
_SOURCE = None

def get_source():
    """Active SceneSource. Defaults to Earth Search, constructed lazily so
    importing melt stays offline (tests never touch the network)."""
    global _SOURCE
    if _SOURCE is None:
        from core.sources import EarthSearchSource
        _SOURCE = EarthSearchSource()
    return _SOURCE

def set_source(source):
    global _SOURCE, _georef_cache
    _SOURCE = source
    _georef_cache = {}   # per-tile georefs come from hrefs; must not leak across sources
```

(Use the actual name of the existing georef cache dict — check `tile_georeference` at src/melt.py:186 and clear whatever cache it uses.)

Then replace, in `src/melt.py`:
- every `item.assets["<name>"].href` → `get_source().band_href(item, "<name>")` (5 sites; the variable-key sites `item.assets[asset].href` / `item.assets[name].href` become `get_source().band_href(item, asset)` / `(…, name)`)
- `_tile_of` body → `return get_source().tile_of(item)`
- both `Client.open(STAC_API)` sites → `get_source().search(...)` keeping identical arguments (the `ids=[item_id]` site becomes `get_source().search(ids=[item_id])[0]` — note `search` already returns a list).

In `src/shelf.py` and `src/shelves.py`, replace each `Client.open(melt.STAC_API).search(...)`/`cl.search(...)` with `melt.get_source().search(...)`, and replace every hard-coded `{"grid:code": {"eq": f"MGRS-{tile}"}}` query fragment with `melt.get_source().tile_query(tile)` so PC's different tile property works transparently.

- [ ] **Step 4: Run the full suite**

Run: `python -m pytest tests -q`
Expected: all PASS (55+). If a test constructs items whose id drives `_tile_of`, the delegation keeps Earth Search parsing — confirm no test regressions.

- [ ] **Step 5: Networked smoke (Earth Search unchanged)**

Run: `python -c "import sys; sys.path.insert(0,'src'); import melt, shelf; it = melt.get_source().search(query=melt.get_source().tile_query('19DEA'), datetime='2021-01-20/2021-01-28')[0]; print(it.id, melt._tile_of(it))"`
Expected: an item id + `19DEA`.

- [ ] **Step 6: Commit**

```bash
git add src/melt.py src/shelf.py src/shelves.py tests/test_sources.py
git commit -m "feat(m1): route all scene access through the active SceneSource"
```

---

### Task 4: Product writer — the day-one data model

**Files:**
- Create: `src/core/products.py`
- Modify: `src/shelf.py` `run_season` (after the total is computed, ~line 470)
- Test: `tests/test_sources.py` (append)

**Interfaces:**
- Consumes: nothing new.
- Produces: `write_product(shelf_name: str, season: str, water: np.ndarray, transform, crs: str, meta: dict, root: Path) -> Path` returning the item directory; layout `output/products/<shelf>/<season>/{water_fraction.tif, item.json}`.

- [ ] **Step 1: Write the failing test**

```python
# append to tests/test_sources.py
import json
import numpy as np
import rasterio
from rasterio.transform import from_origin
from core.products import write_product  # noqa: E402


def test_write_product_cog_and_item(tmp_path):
    water = np.zeros((40, 50), "f4"); water[10:14, 20:26] = 0.5
    tr = from_origin(-2_000_000, 1_000_000, 30, 30)
    d = write_product("george_vi", "2020-21", water, tr, "EPSG:3031",
                      {"area_km2": 142.48, "shelf_area_km2": 14390.6,
                       "obs_cloud": 1.4, "poorly_observed": False,
                       "sensor": "sentinel-2"}, root=tmp_path)
    with rasterio.open(d / "water_fraction.tif") as src:
        assert src.crs.to_string() == "EPSG:3031"
        assert src.transform == tr
        assert float(src.read(1)[12, 22]) == 0.5
    item = json.loads((d / "item.json").read_text())
    p = item["properties"]
    assert p["area_km2"] == 142.48 and p["sensor"] == "sentinel-2"
    # reserved fields exist and are null until M3/M4 fill them
    for k in ("volume_km3", "depth_mean_m", "uncertainty_km2", "sensor_confidence"):
        assert k in p and p[k] is None
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest tests/test_sources.py -q`
Expected: FAIL with `ModuleNotFoundError: No module named 'core.products'`

- [ ] **Step 3: Write the implementation**

```python
# src/core/products.py
"""The observatory's data model, v0: one product per (shelf, season).

A cloud-optimized GeoTIFF of per-cell water fraction on the fixed EPSG:3031
grid, plus a STAC-style item.json. Reserved fields (null for now) are the
full-vision schema from the design doc - M3 (sensor_confidence) and M4
(depth/volume/uncertainty) fill them in place; nothing ever migrates."""
import json
from pathlib import Path

import numpy as np
import rasterio

RESERVED = ("volume_km3", "depth_mean_m", "uncertainty_km2", "sensor_confidence")


def write_product(shelf_name, season, water, transform, crs, meta, root):
    d = Path(root) / shelf_name / season
    d.mkdir(parents=True, exist_ok=True)

    profile = dict(driver="GTiff", height=water.shape[0], width=water.shape[1],
                   count=1, dtype="float32", crs=crs, transform=transform,
                   tiled=True, blockxsize=512, blockysize=512,
                   compress="deflate", predictor=2)
    tif = d / "water_fraction.tif"
    with rasterio.open(tif, "w", **profile) as dst:
        dst.write(water.astype("float32"), 1)
        dst.build_overviews([2, 4, 8, 16], rasterio.enums.Resampling.average)

    props = {"shelf": shelf_name, "season": season, **meta}
    for k in RESERVED:
        props.setdefault(k, None)
    item = {"type": "Feature", "stac_version": "1.0.0",
            "id": f"{shelf_name}-{season}",
            "properties": props,
            "assets": {"water_fraction": {
                "href": "./water_fraction.tif",
                "type": "image/tiff; application=geotiff",
                "roles": ["data"]}}}
    (d / "item.json").write_text(json.dumps(item, indent=1))
    return d
```

- [ ] **Step 4: Run test to verify it passes**

Run: `python -m pytest tests/test_sources.py -q`
Expected: PASS

- [ ] **Step 5: Wire into run_season**

In `src/shelf.py`, at the end of `run_season` just before `return`:

```python
    try:
        from core.products import write_product
        write_product(CURRENT_SHELF, label, water, grid_tr, GRID_CRS,
                      {"area_km2": round(total, 2),
                       "shelf_area_km2": round(shelf_km2, 1),
                       "obs_cloud": obs_cloud, "poorly_observed": dom_poorly,
                       "sensor": "sentinel-2"},
                      root=melt.ROOT / "output" / "products")
    except Exception as e:   # products are additive; never fail a season over them
        print(f"  [products] WARNING: {type(e).__name__}: {str(e)[:60]}")
```

Add `output/products/` to `.gitignore` (regenerable, large).

- [ ] **Step 6: Full suite + commit**

Run: `python -m pytest tests -q` — all PASS.

```bash
git add src/core/products.py src/shelf.py tests/test_sources.py .gitignore
git commit -m "feat(m1): COG+STAC product writer wired into run_season (reserved M3/M4 fields)"
```

---

### Task 5: PC regression gate + decision log

**Files:**
- Create: `scripts/regression_pc.py`, `docs/decisions/0001-compute-backend.md`

**Interfaces:**
- Consumes: `PlanetaryComputerSource` (Task 2), `melt.set_source` (Task 3).

- [ ] **Step 1: Write the regression script**

```python
# scripts/regression_pc.py
"""M1 acceptance gate: George VI 2020-21 via Planetary Computer must match
the blind-validated Earth Search result within 3%.

Networked + slow (~10 min); run manually, not part of pytest:
    python scripts/regression_pc.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
import melt
import shelf
from core.pc import PlanetaryComputerSource

FROZEN = 142.48          # 2020-21, blind-validated (19DEA = 127.17)
TOL = 0.03

melt.set_source(PlanetaryComputerSource())
grid = shelf.build_fixed_grid()          # cached grid; footprints identical
r = shelf.run_season("2020-21", grid=grid)
got = r["shelf_km2"]
dev = abs(got - FROZEN) / FROZEN
print(f"\nPC 2020-21 = {got} km2 vs frozen {FROZEN} (dev {100*dev:.2f}%)")
print("REGRESSION:", "PASS" if dev <= TOL else "FAIL")
sys.exit(0 if dev <= TOL else 1)
```

- [ ] **Step 2: Run it**

Run: `python scripts/regression_pc.py`
Expected: `REGRESSION: PASS`. If FAIL: diff the per-tile lines against the Earth Search run (`output/shelf/history_run.log` 2020-21 block); the usual suspects are scene-inventory gaps in one catalog or an unsigned-href error (both fixable in the adapter, not the science).

- [ ] **Step 3: Write the decision log**

```markdown
# 0001 — First cloud backend: Microsoft Planetary Computer

**Date:** 2026-07-29 · **Status:** accepted

## Decision
Planetary Computer is the first cloud backend behind the SceneSource adapter.
Google Earth Engine remains available later through the same seam.

## Alternatives
- **Earth Engine:** free massive compute + easy hosted apps, but requires
  rewriting detection in EE's server-side idiom - abandons the blind-validated
  numpy code and its 51-test suite.
- **Stay on Earth Search only:** zero work, but no compute story and single
  provider risk.

## Reasoning
PC runs the EXISTING validated Python unchanged (same COGs, signed hrefs),
so the George VI record stays the regression gate. Open standards (STAC),
reproducible without an account-bound runtime. The adapter keeps this
reversible - the regret-proofing keystone from the design doc.

## Revisit when
M2 scale-out saturates laptop bandwidth (move detection loops into PC's Hub /
Dask), or M5 needs hosted interactive maps (EE apps become attractive).
```

- [ ] **Step 4: Commit**

```bash
git add scripts/regression_pc.py docs/decisions/0001-compute-backend.md
git commit -m "feat(m1): PC regression gate (PASS) + compute-backend decision log"
```

---

## Self-review notes

- Spec coverage: 6.1 adapter (Tasks 1–3), 6.2 data model (Task 4), 6.3 port + regression (Tasks 2, 5), 6.4 decision log (Task 5). The `detect/`-package split from 6.1 is deliberately deferred: `set_source` + logical band names already isolate the science from providers, and moving `melt.py` functions wholesale would churn 51 passing tests for zero behavior change — revisit at M2 when the PC Dask story needs importable pure modules (YAGNI).
- Type consistency: `SceneSource.search` returns a list everywhere (the `ids=` call site indexes `[0]`); `band_href(item, band)` signature identical in all three sources; `tile_query` consumed by shelf/shelves in Task 3.
- No placeholders: every step has runnable code or an exact command.

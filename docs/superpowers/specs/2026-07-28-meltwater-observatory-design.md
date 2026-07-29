# Antarctic Meltwater Observatory — Design & Roadmap

**Date:** 2026-07-28
**Status:** Approved (vision + roadmap). M1 ready to plan.
**Author:** Advait (with Claude)

---

## 1. Vision

Turn the validated single-shelf George VI meltwater detector into **the first open,
accessible, cloud-proof, volume-aware Antarctic surface-meltwater observatory** — a
living tool, not another static dataset.

One line: *pick any Antarctic ice shelf, any date, and see how much meltwater is on
it (area **and** volume), gap-free and year-round, from free satellites.*

### Why this, why now (the gap)
The field has the pieces but not the whole, and they live in separate papers as static
`.tif` dumps:

- Optical lake-**area** mapping is mature and now continent-wide — so an area dataset
  alone is no longer novel (Nature Climate Change 2025, s41558-025-02363-5).
- **SAR** melt detection (Sentinel-1; sees through cloud and polar night) is emerging
  but siloed (Frontiers 2025 George VI winter melt; ScienceDirect S1569843224002498).
- **Depth/volume** retrieval is being cracked but ICESat-2 at scale is "hundreds of TB
  of unstructured point cloud… challenging to use at scale" (TC 2024, 18/5173/2024).
- A **circum-Antarctic spatio-temporal record is "entirely lacking"** and Antarctic
  lakes "remain understudied" vs. crowded Greenland (Sci. Reports 2019, s41598-019-50343-5).
- NSIDC's near-real-time melt service was **defunded Oct 2025** — demand is real, and
  funded tools die.

**The unfilled niche = integration + accessibility + rigor:** fuse optical + SAR into a
gap-free year-round record, resolve **volume** not just area, and make it a *usable*
open tool. Nobody has built that. It is name-making because people cite tools they use.

## 2. Non-goals (YAGNI)
- Not another standalone lake-area dataset.
- Not a novel single-algorithm paper; the novelty is the integrated, accessible system
  (which can *host* a hero algorithm, e.g. the fusion or drainage detector).
- Not full ICESat-2 point-cloud processing at scale (use empirical optical depth
  calibrated on sparse ICESat-2 instead).
- Not funded infrastructure; must run at near-zero cost.

## 3. Target architecture (5 layers, bottom-up)

1. **Compute & data foundation** — validated detection on a free cloud-geospatial
   backend (Microsoft Planetary Computer or Earth Engine) hosting Sentinel-2/-1 and
   Landsat as analysis-ready data; incremental (process only new scenes).
2. **Detection & fusion core** — the cloud-robust optical detector **+** a Sentinel-1
   SAR detector, **fused** to a gap-free, year-round, cloud-independent per-pixel melt
   record. Fusion is designed as *"each sensor emits a melt observation with a
   confidence; combine observations"* — sensor-agnostic, not "optical product + SAR
   patch."
3. **Physical retrieval** — empirical pond **depth** (calibrated vs sparse ICESat-2) →
   **volume (km³)**; plus melt onset/duration and **drainage events** (a lake that
   vanishes between scenes → likely hydrofracture) → per-shelf vulnerability signal.
4. **Data products & open API** — standardized per-shelf, per-date area+volume on a
   fixed grid; STAC-cataloged, DOI'd open dataset; queryable API.
5. **Accessible frontend** — interactive Antarctica map: choose shelf + year, view
   area/volume trends, melt-year rankings, drainage alerts.

Cross-cutting: validation (blind labels, cross-sensor, uncertainty) and full
reproducibility (open source, containerized) at every layer.

## 4. Regret-proofing principles

Regret comes from a few irreversible foundational choices, not the reversible many.
Lock the few in correctly now; architect so the rest stays cheap to change.

**Irreversible — get right up front:**
1. **Backend-agnostic science core.** Detection stays portable Python (numpy/xarray)
   behind swappable compute + I/O adapters. Makes "Planetary Computer vs Earth Engine"
   reversible; never trapped on a platform. *This is the keystone.*
2. **Full-vision data model on day one** — multi-sensor, area + depth/volume + per-pixel
   uncertainty, continent-wide, open standards (STAC + COG/Zarr + xarray). Milestones
   fill it incrementally; no migration ever.
3. **Region-agnostic** — boundaries + tiles are parameters. Antarctica now, Greenland or
   any glacier later for free.
4. **Open-source + open-data + reproducible from commit #1.**
5. **Near-zero running cost** — serverless/static site, on-demand compute. Survives with
   no funding (the NSIDC lesson, turned into an advantage).

**Reversible — decide later, no agonizing:** exact backend; threshold-vs-ML detector
(pluggable); dashboard tech; depth algorithm; shelf order.

**Residual bets (can't eliminate) and hedges:**
- Solo-student sustainability → every milestone is independently shippable/citable; stop
  anywhere with a real result.
- Tool-vs-algorithm citability → the platform hosts a hero algorithm; get both.
- Speed vs. the field → ship the accessible tool early around the optical product;
  fusion/volume keep it ahead.

**Decision log:** record every foundational choice with its alternatives and reasoning
(`docs/decisions/NNNN-*.md`), so future revisiting is deliberate, not second-guessing.

## 5. Milestone roadmap (delivers everything, ordered by dependency)

- **M1 — Foundation.** Refactor detection into a backend-agnostic core + adapters; stand
  up the full-vision data model; port to the cloud backend; prove George VI reproduces
  the validated numbers in minutes. *Unblocks all scale.*
- **M2 — Optical scale-out.** Peninsula → all major shelves; publish open dataset v0.1
  (accessibility starts — build in public); per-region validation gate.
- **M3 — SAR fusion.** Add Sentinel-1 detection + fusion → gap-free/year-round; validate
  vs the 2025 George VI SAR paper. *The novel, cloud-proof core.*
- **M4 — Volume.** Empirical depth → volume; calibrate vs ICESat-2; drainage/vulnerability
  signals.
- **M5 — Platform.** Interactive dashboard + open API + DOI + short preprint. *The
  name-making deliverable.*

Each milestone is independently valuable; none is all-or-nothing.

## 6. M1 detailed design (first buildable sub-project)

**Goal:** same validated science, restructured so it is portable, testable, and
cloud-scalable — and pouring results into the day-one data model. No science change; the
George VI 9-season numbers must reproduce.

**6.1 Backend-agnostic core.** Split `src/melt.py`/`shelf.py` into:
- `detect/` — pure functions on arrays (NDWI, shadow test, hysteresis, cloud mask,
  clear-scene selection). No I/O, no globals. Already ~90% pure; formalize it.
- `io/` — an adapter interface (`SceneSource`) with methods to list scenes for a
  tile/date-range and read band windows. Implementations: `LocalSTAC`/`EarthSearch`
  (current), `PlanetaryComputer` (new). The core never imports an adapter.
- `grid/` — the fixed EPSG:3031 grid + shelf-boundary rasterization (exists).

**6.2 Data model (open standards).** Per (shelf, date, sensor): a record carrying
`area_km2`, `water_fraction` grid (COG/Zarr), and reserved fields for
`depth`/`volume`/`uncertainty`/`sensor_confidence` (null until M3/M4). Catalog as STAC
items; store arrays as cloud-optimized (COG or Zarr) so the frontend/API can read them
directly. Seasonal + per-date granularity.

**6.3 Backend port.** Add the `PlanetaryComputer` adapter; run George VI on it; assert
the 9-season history matches the local run within tolerance (regression test).

**6.4 Decision to make in M1 (logged):** Planetary Computer vs Earth Engine as the
first backend. Recommendation: **Planetary Computer** (Python/Dask ports existing code
directly; open/reproducible), with the adapter boundary keeping GEE available later,
especially for frontend hosting.

**M1 done =** backend-agnostic core with tests; PC adapter; George VI reproduces
validated numbers on the cloud; results written in the day-one data model; decision log
started.

## 7. Open decisions (deferred, not blocking M1)
- Depth-retrieval algorithm (M4).
- Dashboard/API stack (M5) — likely static site + serverless for zero cost.
- ML detector vs threshold (post-M2, pluggable behind the same core interface).

---

*This document is the strategic spec. M1 gets its own implementation plan
(writing-plans) next; M2–M5 each get their own spec → plan cycle when reached.*

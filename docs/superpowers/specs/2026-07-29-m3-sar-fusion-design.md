# M3 — SAR Fusion: the Melt-State Map

**Date:** 2026-07-29 · **Status:** approved (product choice: melt-state map)
**Depends on:** M1 (SceneSource adapter, products), the validated GVI optical record

## 1. Goal

Upgrade the observatory from a cloud-gated pond map to a **year-round, per-pixel
melt-state record**: every grid cell classified **DRY / WET / PONDED**, built by
fusing Sentinel-1 radar (wetness, sees through cloud and darkness) with the
blind-validated Sentinel-2 optical detector (ponds). George VI only in M3; other
shelves scale later via the same region-agnostic machinery.

**Why melt-state instead of "gap-filled pond area":** radar physics. Wet snow
collapses backscatter; radar therefore detects *melting*, not pond boundaries.
Pretending radar pixels are pond area would bake hidden assumptions into a
headline number. The three-state product is honest, richer (melt onset/duration
— quantities the field wants), and absorbs the Larsen C finding: slush is WET,
not PONDED, by construction.

## 2. Sentinel-1 detection (the new science)

- **Data:** `sentinel-1-rtc` on Planetary Computer (terrain-corrected gamma0,
  HH, IW mode; verified 34 scenes/month over GVI). Reached via a
  `Sentinel1Source(PlanetaryComputerSource)` with `COLLECTION="sentinel-1-rtc"`,
  `BAND_MAP={"hh": "hh"}` — the M1 seam, one small subclass.
- **Method (published standard, Nagler & Rott lineage):** per pixel, melt when
  backscatter drops **≥ 3 dB below that pixel's frozen-winter reference**.
- **Winter reference:** per-pixel mean gamma0 over austral winter (Jun–Aug),
  built **per relative orbit** (incidence-angle geometry differs by track;
  mixing tracks smears the baseline). Two winter baselines bracket each melt
  season; use the preceding winter.
- **Compositing:** per-pixel wet/dry time series at S1's natural ~6–12-day
  revisit, on the existing fixed 30 m EPSG:3031 GVI grid (RTC native ~10 m,
  averaged down like the optical path).

## 3. Fusion rule (simple, auditable)

For each grid cell and time step:

    PONDED  if the nearest clear optical scene within ±6 days shows pond
            AND radar says wet
    WET     if radar says wet
    DRY     otherwise

- Optical never overrides radar-dry to wet (cloud FP protection); radar never
  invents ponds (physics honesty). Conflicts (optical pond, radar dry) are
  flagged, counted, and reported — they are the diagnostic of both systems.
- Per-cell **confidence** = f(radar margin below baseline, optical clarity,
  time gap) → fills the reserved `sensor_confidence` field from M1's data
  model. Nothing migrates; the schema was built for this day.

## 4. Scope guard (v1 = proof on the anchor shelf)

- **Shelf:** George VI only.
- **Seasons:** 2019-20 (record), 2020-21 (blind-validated anchor), 2021-22
  (quiet year) + their preceding winters. Full 9-season + multi-shelf SAR
  scale-out is a follow-on batch, not M3.
- **Products:** melt-onset / melt-end / duration maps + wet-extent and ponded-
  extent series per season, written through the existing product writer.

## 5. Validation gates (all must pass before "M3 done")

1. **Containment:** ≥90% of blind-validated optical pond pixels are radar-WET
   on the nearest S1 pass (ponds are a subset of melt; if radar misses
   validated ponds, the baseline or threshold is wrong).
2. **Season ranking:** radar wet-extent ranks the seasons like the validated
   optical record (2019-20 ≫ 2020-21 > 2021-22).
3. **Winter sanity:** near-zero wet in deep winter, EXCEPT the buried/stored
   meltwater anomalies documented for GVI by the Frontiers 2025 S1 paper
   (feart.2025.1545009) — reproduce their qualitative winter-storage signal.
4. **Timing plausibility:** melt onset in Nov–Dec, end in Feb–Mar, spatially
   coherent (no salt-and-pepper).

## 6. Non-goals (YAGNI)

- No deep-learning SAR classifier (threshold-vs-baseline is published,
  auditable, and testable; ML is a later upgrade behind the same interface).
- No pond-depth from radar. No VV/HH dual-pol work (HH only available here).
- No attempt to make radar delineate pond edges.

## 7. Risks & honest unknowns

- **Winter baseline contamination** (early/late melt in "winter" months):
  mitigate with median instead of mean if contaminated; check histograms.
- **Orbit-track coverage unevenness** at -71°: some cells see 2 tracks, some 1;
  report per-cell revisit; never compare wet-DAYS across cells with different
  revisit without normalizing.
- **3 dB threshold transferability:** validated in literature for dry-snow
  zones; GVI's percolation regime may need 2–4 dB sensitivity check (sweep and
  report, pick by validation gate 1).

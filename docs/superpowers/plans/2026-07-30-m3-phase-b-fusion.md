# M3 Phase B — Fusion: the Melt-State Product

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fuse the gate-validated Sentinel-1 wet/dry record with the blind-validated
optical pond record into per-pixel DRY/WET/PONDED melt-state composites for the
three pilot George VI seasons, with per-cell sensor confidence filling M1's
reserved field.

**Architecture:** `src/fusion.py` holds pure, offline-tested state logic
(classify + confidence + conflict accounting). Optical pond masks are
reprojected per clear scene onto the 120 m SAR grid (generalizing the Gate-1
helper). A season builder walks the S1 scene series, attaches the nearest
optical evidence within ±6 days, and emits a state stack + summary + products.

**Tech Stack:** numpy, rasterio, existing sar.py/shelf.py/melt.py infrastructure.

## Global Constraints

- Suite green after every task (65 now); fusion logic gets synthetic tests.
- PONDED ⊆ WET **by construction** — the classifier cannot emit PONDED where
  radar says dry; optical-pond-but-radar-dry is a COUNTED CONFLICT, never a state.
- Threshold: use `gate1_best_thresh_db` from gates.json, **but final choice
  waits for Gate 3**: if winter wet-fraction at 2 dB exceeds the 2% bar while
  3 dB passes, adopt 3 dB (containment 0.999 vs 1.000 is an acceptable trade
  for specificity). Task 5 encodes this decision rule; do not hand-pick.
- Never modify optical outputs; fusion is additive.
- State encoding fixed: 0=UNOBSERVED, 1=DRY, 2=WET, 3=PONDED (uint8).

## File Structure

- `src/fusion.py` — pure logic + season builder
- `tests/test_fusion.py` — synthetic tests
- `scripts/fusion_pilot.py` — networked pilot runner (3 seasons) + validation report

### Task 1: classify_state + confidence (pure)

**Interfaces produced:**
- `classify_state(wet, obs, pond_evidence) -> ndarray uint8` where inputs are
  bool arrays; pond_evidence = optical pond within ±6 days on this cell.
- `confidence(margin_db, gap_days) -> ndarray f4 in [0,1]`:
  `clip(margin_db / 6, 0, 1) * exp(-gap_days / 6)` for PONDED cells (radar
  margin × optical freshness); WET cells use the margin term alone; DRY the
  margin below-threshold distance analogue; UNOBSERVED = 0.
- `conflict_mask(wet, obs, pond_evidence) -> ndarray bool` (pond & obs & ~wet).

- [ ] Step 1: failing tests (state table: all 8 combinations of wet/obs/pond;
  PONDED requires wet AND pond; conflict counted not classified)

```python
def test_state_table():
    import numpy as np, fusion
    wet  = np.array([[1,1,0,0],[1,0,1,0]], bool)
    obs  = np.array([[1,1,1,1],[1,1,0,0]], bool)
    pond = np.array([[1,0,1,0],[0,1,1,1]], bool)
    st = fusion.classify_state(wet, obs, pond)
    assert st.tolist() == [[3,2,1,1],[2,1,0,0]]
    cf = fusion.conflict_mask(wet, obs, pond)
    assert cf.tolist() == [[False,False,True,False],[False,True,False,False]]
```

- [ ] Step 2: run → FAIL. Step 3: implement (vectorized, ~15 lines).
- [ ] Step 4: PASS. Step 5: commit `feat(m3b): melt-state classifier`.

### Task 2: optical pond series on the SAR grid

**Interfaces produced:** `pond_series(season, sgrid) -> list[(date, mask)]` —
for each scene in the season's clean pool (select_pool over `_season_scenes`
for every GVI tile), detect ponds full-res, reproject to the SAR grid
(fraction ≥ 0.25 like Gate 1), OR per date across tiles. Cache npz per season.

- [ ] Step 1: failing test for the date-attach helper
  `nearest_pond_evidence(dates, masks, scene_date, window_days=6)` (synthetic:
  picks nearest within window, all-False outside window).
- [ ] Steps 2-4: implement (generalize `_validated_ponds_on_sargrid`), PASS.
- [ ] Step 5: commit `feat(m3b): optical pond series on SAR grid`.

### Task 3: season melt-state builder

**Interfaces produced:** `melt_state_season(season, grid, source, bbox, thresh)`
— walks S1 scenes date-ordered (per orbit, cached baselines), computes wet/obs,
attaches pond evidence, classifies, accumulates: per-cell final-state days
(days_wet, days_ponded), conflict count, mean confidence; writes
`output/sar/state_<season>.npz` + JSON summary (wet/ponded extent maxima,
conflict rate, onset/duration from state series) and a product entry filling
`sensor_confidence`.

- [ ] Step 1: failing test on the accumulator (synthetic 2-scene walk).
- [ ] Steps 2-4: implement, PASS. Step 5: commit.

### Task 4: pilot runs + fusion validation report (networked)

- [ ] `scripts/fusion_pilot.py`: run 2020-21, 2019-20, 2021-22; report must show
  (a) conflict rate < 5% per season (else diagnose before proceeding),
  (b) PONDED extent ≤ optical pond extent per date (construction check),
  (c) ponded-extent season ranking matches the validated optical ranking,
  (d) wet/ponded ratio per season (the new science number: how much melt
  ponds vs merely wets).
- [ ] Run detached, log to output/sar/fusion_pilot.log; verdicts to
  output/sar/fusion_report.json. Commit results.

### Task 5: threshold finalization (consumes Gate 3)

- [ ] Read gates.json gate3_winter: if 2 dB passed (<2% winter wet), keep 2 dB;
  else if 3 dB data needed, rebuild affected composites at 3 dB and rerun the
  pilot; record the decision + evidence in docs/decisions/0002-sar-threshold.md.
- [ ] Final commit: `feat(m3b): M3 complete - fused melt-state record on GVI`.

## Self-review
- Spec §3 fusion rule → T1; §2 grids/orbits reuse → T2/T3; §5-style gates → T4;
  threshold governance → T5. PONDED⊆WET enforced in classifier, checked again in
  T4(b) — belt and braces. No placeholders; all offline steps carry code.

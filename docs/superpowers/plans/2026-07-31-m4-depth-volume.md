# M4 Depth & Volume Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Physically-based dual-band pond depths on George VI, calibrated against
ICESat-2 crossovers, aggregated to per-season volume with Monte Carlo
uncertainty, plus the drainage-vs-refreeze event catalogue — all gated per spec
v1.1 (docs/superpowers/specs/2026-07-30-m4-depth-volume-design.md).

**Architecture:** `src/depth.py` holds pure attenuation math + Monte Carlo
(offline-tested). `src/icesat2.py` harvests ATL03 crossovers via earthaccess
(NASA creds already on disk). `scripts/depth_pilot.py` runs the two pilot
seasons end-to-end and writes gate verdicts. Validity mask v1 = GVI optical
ponds (blind-validated ×1.0 regime); swaps to fused PONDED when M3-B pilots
land — mask is an input, so the swap is one argument.

**Tech Stack:** numpy, rasterio, earthaccess (NASA ATL03), scipy, existing
melt/shelf infra.

## Global Constraints

- Suite stays green (68 now); all depth math offline-tested synthetically.
- Depth formula: `z = (ln(A_d - R_inf) - ln(R_w - R_inf)) / g`, per band;
  reflectance = DN / 10000 (melt.DN_PER_REFLECTANCE); z clipped to [0, 10 m].
- Dual-band rule: keep mean(z_red, z_green) when |z_red - z_green| <= 0.3 *
  max(z_red, z_green); else flag `depth_uncertain` (excluded from volume,
  counted in uncertainty).
- ICESat-2: ATL03 strong beams only; refraction factor 0.752; crossover =
  track ∩ pond polygon within ±3 days of the clear scene.
- Never touch existing detection outputs; depth/volume are additive product
  fields (they were RESERVED for this).
- Every gate result goes to output/depth/gates_m4.json before "done" claims.

## File Structure

- `src/depth.py` — attenuation math, dual-band merge, rim-albedo estimation,
  Monte Carlo volume
- `src/icesat2.py` — ATL03 search/download (earthaccess), photon-histogram
  depth extraction, crossover pairing
- `scripts/depth_pilot.py` — end-to-end pilot: calibrate g -> depth maps ->
  volumes -> gates
- `tests/test_depth.py` — synthetic tests (no network)
- `output/depth/` — crossovers, depth maps, gate verdicts (gitignored)

### Task 1: attenuation math + dual-band merge (pure)

**Interfaces produced:**
- `depth_single(R_w, A_d, R_inf, g) -> z` (vectorized, NaN-safe: R_w <= R_inf
  or R_w >= A_d -> NaN; clip 0..10)
- `dual_band_merge(z_red, z_green, tol=0.3) -> (z, uncertain_mask)`
- `rim_albedo(band_img, pond_mask, rim_px=2) -> A_d per pond` (label ponds
  with scipy.ndimage, dilate rim, mean of rim pixels per label)

- [ ] Step 1: failing tests

```python
# tests/test_depth.py
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
import numpy as np
import depth


def test_depth_single_analytic():
    # z = (ln(Ad-Rinf) - ln(Rw-Rinf)) / g ; with Ad=0.6, Rinf=0.05, g=0.8:
    # Rw chosen so z = 1.0 exactly: Rw = Rinf + (Ad-Rinf)*exp(-g*z)
    Rw = 0.05 + 0.55 * np.exp(-0.8)
    z = depth.depth_single(np.array([Rw]), 0.6, 0.05, 0.8)
    assert abs(z[0] - 1.0) < 1e-6


def test_depth_single_invalid_inputs_nan():
    z = depth.depth_single(np.array([0.04, 0.7]), 0.6, 0.05, 0.8)
    assert np.isnan(z).all()          # below Rinf / above Ad -> NaN


def test_dual_band_merge():
    zr = np.array([1.0, 1.0, np.nan]); zg = np.array([1.1, 2.0, 1.0])
    z, unc = depth.dual_band_merge(zr, zg)
    assert abs(z[0] - 1.05) < 1e-9 and not unc[0]
    assert unc[1] and np.isnan(z[1])
    assert unc[2] and np.isnan(z[2])  # one-band-only counts uncertain in v1


def test_rim_albedo_two_ponds():
    img = np.full((20, 20), 0.5, "f4")
    mask = np.zeros((20, 20), bool); mask[5:8, 5:8] = True; mask[12:15, 12:15] = True
    img[4:9, 4:9] = 0.5; img[11:16, 11:16] = 0.7   # different rim brightness
    labels, ad = depth.rim_albedo(img, mask, rim_px=1)
    vals = sorted(set(np.round(ad[mask], 2)))
    assert vals == [0.5, 0.7]
```

- [ ] Step 2: run -> FAIL (no module). Step 3: implement `src/depth.py`.
- [ ] Step 4: PASS. Step 5: commit `feat(m4): attenuation depths + dual-band merge + rim albedo`.

### Task 2: Monte Carlo volume (pure)

**Interfaces produced:** `mc_volume(z_fn, n=1000, seed=...) -> dict` — caller
passes a closure that recomputes the depth map from perturbed (g, A_d_shift,
R_inf); returns {volume_km3_median, p16, p84}. Perturbations: g ~ N(g_cal,
sigma_cal); A_d shift ~ N(0, rim_sigma); R_inf uniform over literature range.

- [ ] Step 1: failing test (deterministic seed; analytic sanity: symmetric
  perturbation of g around truth -> median volume within 5% of unperturbed).
- [ ] Steps 2-5: implement, PASS, commit.

### Task 3: ICESat-2 crossover harvest (networked, cached)

**Interfaces produced:**
- `search_atl03(bbox, t0, t1) -> granules` via earthaccess (login strategy
  netrc — already configured).
- `photon_depth(granule, pond_polygon) -> list[(lat, lon, depth_m)]`:
  strong-beam photons clipped to polygon; histogram elevations (0.1 m bins);
  surface peak = strongest, bottom peak = next local max >= 0.3 m below;
  depth = separation * 0.752; require >= 50 bottom-return photons else None.
- `harvest(season) -> crossovers.json` for the big validated GVI ponds
  (polygons from the season's detection, area >= 0.5 km2), ±3 days pairing.

- [ ] Step 1: offline test for the histogram peak-splitter with synthetic
  photon clouds (surface at 20.0 m, bottom at 18.5 m -> depth 1.128 after
  x0.752).
- [ ] Step 2-4: implement; run harvest for 2019-20 + 2020-21 detached
  (output/depth/crossovers_<season>.json); N total >= 30 or gate 1 reports
  calibration-starved per spec risk clause.
- [ ] Step 5: commit.

### Task 4: calibration + pilot depth maps + gates (networked)

- [ ] `scripts/depth_pilot.py`:
  1. Fit g_red (and g_green) by least squares: satellite reflectance at
     crossover pixels vs ICESat-2 depths; report RMSE -> **Gate 1: RMSE < 0.5 m,
     N >= 30**.
  2. Depth maps for the two pilot seasons over the validity mask (GVI optical
     ponds v1); **Gate 2**: centre-vs-rim gradient check (mean depth of pond
     interior > mean depth of pond edge for >= 90% of ponds >= 0.5 km2).
  3. Volumes + Monte Carlo; **Gate 3**: 2019-20 GVI volume within x2 of Corr
     et al. 2022.
  4. Write output/depth/gates_m4.json; products gain depth_mean_m, volume_km3,
     uncertainty via a small updater over output/products/george_vi/.
- [ ] Run detached; evaluate; commit results + decision record if any
  parameter choices were made.

### Task 5: drainage-vs-refreeze catalogue

- [ ] `drainage_events(season)` in src/depth.py or own module: consecutive
  clear-scene pond masks -> per-pond-polygon area series; candidates = >80%
  loss in <=14 days; classify DRAINAGE vs REFREEZE per spec v1.1 (loss speed,
  basin darkness, optional SAR wet->dry where composites exist); save pre/post
  chips; **Gate 4**: finds documented 2019-20 GVI drainage events, zero flags
  on stable ponds (spot-check chips).
- [ ] Offline test: synthetic area series classifier table.
- [ ] Run for pilots, evaluate, commit `feat(m4): drainage-event catalogue + gate 4`.

## Self-review
- Spec v1.1 coverage: §3 dual-band+validity -> T1/T4; ICESat-2 recipe -> T3;
  §4 MC -> T2/T4; §5 discrimination -> T5; §6 gates 1-4 -> T4/T5. Validity-mask
  swap (optical->PONDED) is a single argument by construction. No placeholders.

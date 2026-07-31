# M4 — Depth & Volume: from km² to km³

**Date:** 2026-07-30 · **Status:** draft for review (parallel-tracked with M3 Phase B)
**Depends on:** M1 (products, adapters), M2 (multi-shelf records). **Not** on M3 —
depth rides on optical pixels; fusion and depth meet later in the data model.

## 1. Goal

Give every detected pond a **depth**, and every shelf-season a **meltwater
volume (km³) with uncertainty** — the quantity hydrofracture models actually
consume. Plus the first **drainage-event detector** (rapid pond disappearance
= the vulnerability precursor).

## 2. Why this is the right next quantity

Area says "how much surface is wet"; volume says **how much water is pressing
on the shelf** — the physical driver of hydrofracture. Published work shows
satellite depths are retrievable to <0.3 m error; no open tool ships them
per-shelf, per-season, with uncertainty. That's our slot.

## 3. Depth method (published, physically based)

**Light attenuation:** water absorbs red light exponentially with depth. For a
pond pixel, depth `z = [ln(A_d − R_inf) − ln(R_w − R_inf)] / g` (Pope et al.
2016 lineage; Moussavi variants):
- `R_w` = observed red reflectance of the pond pixel
- `A_d` = bottom albedo, estimated per pond from its rim pixels (the ice just
  outside the pond edge)
- `R_inf` = optically-deep-water reflectance (constant, literature)
- `g` = red attenuation coefficient (literature starting value, then
  **calibrated against ICESat-2**)

**Dual-band retrieval:** depth is computed independently from red AND green
(red = primary: robust for shallow <3 m ponds; green = deeper penetration but
scattering-sensitive). Agreement within 30% -> high-confidence depth (mean);
disagreement -> pixel flagged `depth_uncertain` (turbidity/shadow indicator).
Two estimates beat one, and their disagreement maps where the physics breaks.

**Validity mask — the M3 synergy:** depth is retrieved ONLY where the fused
melt-state record classifies PONDED (radar-wet + optical pond). This
automatically excludes frozen lake lids (Amery) and slush (Larsen regimes),
where attenuation-model depths are meaningless — the corrections table told us
this qualitatively; the melt-state map enforces it per pixel.

**Calibration/validation — ICESat-2, concrete recipe:** ATL03 geolocated
photons, STRONG beams only (gt1l/gt2l/gt3l or their r-counterparts by
orientation). Per crossover segment over a pond polygon: histogram photon
elevations -> surface peak (air-water) and bottom peak; raw depth = peak
separation; true depth = raw x 0.752 (refraction + speed-of-light correction,
Parrish et al. 2019). Crossover = ICESat-2 track intersecting a validated pond
polygon within +/-3 days of a clear Sentinel-2 scene. We already hold NASA
Earthdata credentials (from the MEaSUReS download). v1 uses **sparse
crossovers** on the big 2019-20/2020-21 ponds - not bulk photon processing,
the published at-scale pain we explicitly avoid.

## 4. Volume & uncertainty

- `V = Σ depth × (10 m)²` over pond pixels, per shelf-season.
- Uncertainty via **Monte Carlo** (1,000 draws, seconds): jointly perturb
  g (calibration posterior), A_d (rim-pixel spread), and R_inf (literature
  range) -> per-season volume distribution; report median and 16-84th
  percentile band. Combined with the per-region area corrections
  (`reference/regional_corrections.json`) already in the data model.
- Products: `depth_mean_m`, `volume_km3`, `uncertainty` fill the fields M1
  **reserved on day one** — zero schema migration, as designed.

## 5. Drainage events (the vulnerability signal)

Between consecutive genuinely-clear scenes (same clean-scene gates as always):
any pond polygon losing **>80% of its area in <=14 days** flags a CANDIDATE.
Candidates are then classified:
  - **DRAINAGE** (the hydrofracture precursor): fast loss AND the emptied basin
    stays dark/deep-floored AND - radar corroboration - the cell flips wet->dry
    in the S1 record within the same window.
  - **REFREEZE** (benign): gradual area fade across multiple scenes and/or the
    surface brightens toward the pale-cyan lid signature while radar stays wet
    (liquid under lid).
Pre/post chips saved for visual confirmation either way. Output: per-shelf
event catalogue (date, location, class, lost area, est. lost volume). The
discrimination kills the refreeze false-alarm class that would otherwise
dominate any Antarctic drainage catalogue.

## 6. Scope guard (v1)

- **Shelf:** George VI only (the validated anchor; big documented ponds).
- **Seasons:** 2019-20 + 2020-21 (record + blind-validated).
- **Gates before "M4 done":**
  1. Depth RMSE vs ICESat-2 crossovers **< 0.5 m** (published bar ~0.3 m;
     student-honest margin), N ≥ 30 crossover points.
  2. Depth map sanity: deepest at pond centres, ~0 at rims (spatial-gradient
     check, automated).
  3. GVI volume within **x2 of Corr et al. (2022)** for the overlapping
     season - a quantitative bar against the published GVI depth/volume work,
     not order-of-magnitude hand-waving.
  4. Drainage detector finds the documented 2019-20 GVI drainage events (the
     literature records several) without flagging stable ponds.

## 7. Non-goals (YAGNI)

- No bulk ICESat-2 photon-cloud processing (sparse crossovers only).
- No turbidity modelling (dual-band disagreement flags it instead).
- No SAR in the depth path (fusion meets depth in the data model, not here).
- No pan-Antarctic depth in v1 — method must survive GVI gates first.

## 8. Risks, honestly

- **ICESat-2 crossover scarcity:** tracks may miss big ponds on clear days →
  mitigation: widen to ±3 days scene-to-track pairing, accept N≥30 rather than
  hundreds; if still starved, report calibration-limited uncertainty and say so.
- **Bottom-albedo estimation** is the known weak link (dark/dirty pond floors
  read too deep) → rim-based per-pond A_d + sensitivity band, never a global
  constant.
- **Slush regimes:** depth model assumes ponded water; on ×0.35–0.6 shelves it
  applies only to the open-water fraction — corrections table already encodes
  exactly this.

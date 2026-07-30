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

**Why red:** green penetrates deeper but is more scattering-sensitive; red is
the published compromise for shallow (<3 m) Antarctic ponds. Both bands already
stream through our pipeline at 10 m.

**Calibration/validation — ICESat-2:** ATL03 photon returns over GVI ponds give
direct depths (surface + bottom returns, refraction-corrected ×0.75). We
already hold NASA Earthdata credentials (MEaSUReS download). v1 uses **sparse
crossovers** (tracks × validated pond polygons, big 2019-20/2020-21 ponds), not
bulk photon processing — the published at-scale pain we explicitly avoid.

## 4. Volume & uncertainty

- `V = Σ depth × (10 m)²` over pond pixels, per shelf-season.
- Uncertainty stack: depth RMSE from the ICESat-2 comparison ⊕ area band from
  the per-region blind corrections (already published in
  `reference/regional_corrections.json`) ⊕ `A_d` rim-estimate sensitivity
  (report the depth spread from rim ±1σ).
- Products: `depth_mean_m`, `volume_km3`, `uncertainty` fill the fields M1
  **reserved on day one** — zero schema migration, as designed.

## 5. Drainage events (the vulnerability signal)

Between consecutive genuinely-clear scenes (same clean-scene gates as always):
any pond polygon losing **>80% of its area in ≤14 days** flags a candidate
drainage event, with pre/post chips saved for visual confirmation. Rapid
drainage = water forced into the shelf = the hydrofracture precursor. Output:
per-shelf event catalogue (date, location, lost area, est. lost volume).

## 6. Scope guard (v1)

- **Shelf:** George VI only (the validated anchor; big documented ponds).
- **Seasons:** 2019-20 + 2020-21 (record + blind-validated).
- **Gates before "M4 done":**
  1. Depth RMSE vs ICESat-2 crossovers **< 0.5 m** (published bar ~0.3 m;
     student-honest margin), N ≥ 30 crossover points.
  2. Depth map sanity: deepest at pond centres, ~0 at rims (spatial-gradient
     check, automated).
  3. GVI 2019-20 volume lands within the literature's reported range for that
     season (order-of-magnitude sanity vs Banwell/Corr).
  4. Drainage detector finds the documented 2019-20 GVI drainage events (the
     literature records several) without flagging stable ponds.

## 7. Non-goals (YAGNI)

- No bulk ICESat-2 photon-cloud processing (sparse crossovers only).
- No green-band dual-model fusion in v1; no turbidity modelling.
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

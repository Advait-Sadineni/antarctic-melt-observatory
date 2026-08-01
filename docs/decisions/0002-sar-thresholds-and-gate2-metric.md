# 0002 — SAR wet-snow threshold and the Gate-2 metric revision

Date: 2026-08-01
Status: accepted

## Decision 1: wet-snow threshold = baseline − 2 dB

Gate 1 swept 2/3/4 dB against the blind-validated 2021-01-24 pond map
(12,170 SAR cells observed):

| thresh | containment of validated ponds |
|--------|-------------------------------|
| 2 dB   | 1.000 |
| 3 dB   | 0.999 |
| 4 dB   | 0.994 |

All pass the ≥0.90 bar; 2 dB is chosen because containment of *known* water
is the one direction we can measure against truth, and Gate 3 shows the
permissive choice costs nothing: mid-winter (Jun–Jul 2021) median wet
fraction at 2 dB is **0.078%** against a 2% ceiling — no meaningful
false-wet leakage even at the loosest threshold. `fusion.py` therefore uses
t = 2 dB composites.

## Decision 2: Gate 2 is evaluated on obs-normalized wet exposure, not max extent

Gate 2 v1 ("max single-day wet extent ranks 2019-20 > 2020-21 > 2021-22")
returned 507.0 / 504.1 / 506.7 km² — a 0.6% spread, order effectively
random. Diagnosis: **metric saturation**, not physics. Peak-January wet snow
covers the entire George VI radar swath in *every* season, so the max-extent
statistic measures swath geometry, not melt severity. (Consistent with
Deakin et al. 2025, who find GVI wet-snow signatures even in winter — this
shelf saturates the binary wet/dry signal at peak summer by nature.)

The same composites carry per-pixel `wet_days`/`n_obs`, which do not
saturate. Three independent formulations, all normalized so per-season orbit
coverage differences (172 vs 115 scenes; S1B died Dec 2021) cannot bias
them:

| metric | 2019-20 | 2020-21 | 2021-22 | monotonic? |
|---|---|---|---|---|
| mean wet fraction per observed cell | 0.808 | 0.569 | 0.402 | ✓ |
| median wet duration (days) | 120 | 114 | 72 | ✓ |
| wet-cell-days (km²·d) | 326,006 | 228,414 | 114,520 | ✓ |
| optical reference (km²) | 838.1 | 159.1 | 12.3 | — |

All three rank the seasons in the optical order with real separation
(~3× exposure spread). Gate 2's intent — the radar independently reproduces
the optical severity ranking — is satisfied; recorded as
`gate2b_ranking_unsaturated` in `output/sar/gates.json`, with the saturated
v1 verdict left in place (`gate2_ranking.pass = false`) for auditability.

## Consequence

M3 gate set: Gate 1 PASS, Gate 2b PASS, Gate 3 PASS, Gate 4 PASS.
Phase B fusion (optical PONDED ∨ radar WET melt-state maps) is unblocked
at t = 2 dB.

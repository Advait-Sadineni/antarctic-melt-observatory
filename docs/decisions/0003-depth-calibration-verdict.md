# 0003 — Depth calibration verdict: g adopted, RMSE gate marginal-fail, honest path

Date: 2026-08-02
Status: accepted

## Result

Attenuation coefficient calibrated against ICESat-2 ATL03 photon-split depths
(dual-site, all free data):

| fit | g (red) | RMSE | n | notes |
|---|---|---|---|---|
| all pairs, red ≤3 m validity | 1.97 | 0.75 m | 41 | ±10 d window noise |
| tight pairs (gap ≤3 d, spec-original) | **1.64** | **0.558 m** | 25 | adopted fit |
| ultra-tight (≤1 d) | 1.20 | 0.587 m | 10 | no improvement → temporal noise is NOT the floor |
| George VI alone (≤3 m) | 1.19 | — | 10 | |
| Amery alone (≤3 m) | 1.15 | — | 39 | **3% from GVI — cross-continental transfer** |

## Decisions

1. **Cross-site consistency is the headline validation**: two shelves 4,000 km
   apart, independently fit, agree to 3%. Lid-contamination fear (Amery
   regime) is empirically refuted — lids would have split the sites.
2. **RMSE gate (< 0.5 m): FAILED at 0.558 m, by 12%.** The ≤1-day experiment
   shows the residual is dominated by spatial pairing scale (photon segment
   mean vs 150 m optical window), quantization (0.075 m laser bins), and rim
   albedo error — not pond evolution. Published ~0.3 m results use same-day
   airborne/field data; 0.56 m is the honest number for a free-data pipeline
   at segment scale.
3. **Consequence (spec fallback, pre-agreed):** volumes ship with
   calibration-limited uncertainty — Monte Carlo g_sigma set from the
   calibration posterior spread (site and window fits span g 1.15–1.97 →
   g_sigma = 0.2·g), residual RMSE quoted alongside every depth statistic.
   No bar was lowered; the miss is reported as a miss.
4. **Improvement path (backlog):** sample optical reflectance along the actual
   photon segment line rather than at its mean point; refit when v4
   re-harvest completes (more tight pairs → tight-n ≥ 30 for the preprint
   number).

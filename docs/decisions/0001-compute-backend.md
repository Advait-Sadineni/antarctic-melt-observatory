# 0001 — First cloud backend: Microsoft Planetary Computer

**Date:** 2026-07-29 · **Status:** accepted

## Decision
Planetary Computer is the first cloud backend behind the SceneSource adapter.
Google Earth Engine remains available later through the same seam.

## Alternatives
- **Earth Engine:** free massive compute + easy hosted apps, but requires
  rewriting detection in EE's server-side idiom - abandons the blind-validated
  numpy code and its test suite.
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

## Evidence
`scripts/regression_pc.py` — George VI 2020-21 through the PC adapter vs the
frozen blind-validated 142.48 km²; result recorded below after first run.

**First run (2026-07-29): PASS.** PC 2020-21 = 139.8 km² vs frozen 142.48
(deviation 1.88%, within the 3% gate). 19DEA picked the SAME blind-validated
scene date (24 Jan 2021) through PC's catalog; the small delta comes from minor
scene-inventory differences between catalogs, not the science.

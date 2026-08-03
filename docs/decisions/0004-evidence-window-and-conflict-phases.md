# 0004 — Melt-core evidence window; conflicts split by radar phase

Date: 2026-08-02
Status: accepted

## Problem

The 2020-21 melt-state pilot returned conflict_rate = 24.7% against the <5%
gate. Attribution analysis (per-evidence-date overlap with high-conflict
cells) showed the 2021-03-25 optical scene alone covered 100% of top-quartile
conflict cells (84% of its own cells conflicted), with early-December dates
contributing the rest. Radar phenology says the shelf refroze around Feb 23
(median last-wet day 115): the "ponds" optics sees in late March are frozen
lake lids — stored meltwater under ice, real physics, but not PONDING
evidence. The blind-validated 2021-01-24 scene was clean (4% own-conflict).

## Decisions

1. **Optical pond evidence is restricted to the melt core (Dec 1 – Feb 28)**
   (`fusion.pond_series`, caches renamed ponds_<season>_core.npz). Outside
   that window, detections are lid-dominated in every validated regime.
2. **Conflicts are phase-split against each cell's own radar wet season**
   (`fusion.split_conflict`, ±6-day margin, composite phenology):
   - wet-season conflicts = true sensor disagreement → `conflict_rate`,
     gated < 5%;
   - off-season conflicts = lid / diurnal-refreeze signal → reported as
     `conflict_days_offseason`, a product feature (it maps stored-water
     lids), never an error metric.
3. All three pilot seasons re-walked from zero under these definitions
   (pre-0004 outputs archived in output/{sar,depth}/pre0004/ for audit).

## Why not just drop the March scene

Ad-hoc scene exclusion is unauditable. A window stated in advance, physically
motivated, applied identically to every season and every future shelf, is.

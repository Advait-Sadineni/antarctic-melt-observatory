# An open, continuously updated multi-sensor observatory of surface meltwater
# on Antarctic ice shelves

**Advait Sadineni** — Pennsylvania State University
*Draft for EarthArXiv — numbers marked [F] freeze from the final pipeline run.*

## Abstract

Surface meltwater ponding preconditions Antarctic ice shelves for
hydrofracture-driven collapse, yet published meltwater records are typically
single-region or single-season snapshots frozen at publication. We present an
open-source observatory that maintains a continuously updated meltwater
record for 12 Antarctic ice shelves across nine melt seasons (2017-18 through
2025-26), built entirely from free data and compute. Sentinel-2 detection
(Moussavi-lineage thresholds with hysteresis growth and cloud-robust scene
selection) is validated by design: stratified blind point samples with hidden
answer keys yield per-regime precision and recall and Horvitz-Thompson
area corrections for four melt regimes. Sentinel-1 wet-snow cross-validation
passes four pre-registered gates, including 100% containment of
blind-validated ponds and mid-winter false-wet below 0.1%. A fused per-pixel
melt-state product (dry / wet / ponded, with sensor conflicts counted rather
than classified) separates melt that merely wets firn from melt that ponds:
in the record 2019-20 season on George VI Ice Shelf, radar-wetted area
exceeded concentrated ponding tenfold (peak 8,111 vs 797 km2) - and in the
quiet 2021-22 season, 8,100 km2 wetted while essentially nothing ponded:
ponding is strongly nonlinear in melt. Pond depths from a
dual-band attenuation model are calibrated against 111 ICESat-2 photon-split
crossover depths spanning George VI and Amery ice shelves; independently
fitted attenuation coefficients at the two sites agree to 3% (g = 1.19 vs
1.15 m^-1 red-band validity-domain fits), with segment-scale RMSE of 0.56 m
against a 0.5 m target, reported as a near-miss. All code, products,
validation artifacts, and decision records are public, and the record
updates each austral summer, including the current season.

## 1. Introduction

Ice shelves buttress the grounded ice behind them; their collapse accelerates
sea-level contribution. The 2002 Larsen B collapse followed a season of
exceptional surface ponding, and meltwater-driven hydrofracture remains the
mechanism by which ice-sheet models project the largest and most uncertain
future sea-level contributions (DeConto & Pollard lineage). Observations of
*where and how much* meltwater ponds - as opposed to where surfaces are
merely wet - are the constraint those projections lack.

Existing Antarctic meltwater datasets are one-off inventories (e.g. Corr et
al., 2022, for January 2017; Stokes et al., 2019 for East Antarctica) or
regional studies. Continuous multi-year records with validation that reports
its own failure modes are, to our knowledge, absent. This paper describes a
system engineered for that gap: a living observatory rather than a snapshot,
maintainable by a single person on free compute, with every scientific
decision recorded in versioned decision documents.

Three principles distinguish the design:
1. **Validation by design, not comparison.** Blind stratified point samples
   with hidden answer keys, scored once; per-regime corrections with
   confidence intervals rather than a single accuracy number.
2. **Failure as signal.** Sensor disagreements are counted, never silently
   classified; failed gates are published alongside passed ones (four of the
   project's decision records exist because a gate or metric failed).
3. **Free and reproducible end-to-end.** Sentinel-1/-2 via public STAC
   catalogs, ICESat-2 via Earthdata, consumer-laptop compute.

## 2. Data

Sentinel-2 L2A surface reflectance (AWS Earth Search; Microsoft Planetary
Computer as a validated alternate backend, cross-backend regression < 2%).
Sentinel-1 RTC gamma0 (Planetary Computer). ICESat-2 ATL03 geolocated
photons (NASA Earthdata). MEaSUReS Antarctic Boundaries v2 shelf polygons.

## 3. Methods

### 3.1 Optical detection
NDWI > 0.19 with green-red > 0.09 shadow rejection (Moussavi et al., 2020),
DN/10000 reflectance, plus hysteresis region growth to 0.14 from seed pixels
- capturing pond margins that fixed thresholds truncate. Per-shelf fixed
30 m EPSG:3031 grids make season-over-season areas denominator-identical.

### 3.2 Clear-scene selection
Cloud false positives, not detection thresholds, dominate error in Antarctic
optical melt mapping. Scene selection is two-branch: among scenes genuinely
clear over the shelf (halo fraction < 0.08) take the peak-melt scene; if no
scene clears the bar because melt itself depresses NDSI (record seasons are
self-obscuring), fall back to scenes with metadata cloud < 15% and halo
< 0.30. Selection windows are evaluated per month for evidence series
(decision 0005): a clean pre-melt month must not veto a hazy record-melt
month's fallback. Rejected alternatives (persistence compositing, clean-scene
union, peak-water selection) and their failure modes are documented in the
repository.

### 3.3 Blind validation and regional corrections
Per validated tile: 80 stratified points (40 detected, 26 near-boundary, 14
far), labeled blind against hidden answer keys, scored once. Horvitz-
Thompson estimation corrects area bias; Wilson intervals on precision.
Four regimes validated: George VI channel-fed ponds (precision 0.625
[0.47-0.76], area-weighted recall 0.51, correction x1.0), Larsen C slush
(0.60, recall 1.0, x0.60), Amery frozen lake lids (0.35, recall 0.70,
x0.50), Larsen B slush-dominated embayment (0.35, recall 1.0, x0.35).
Raw detections are meltwater-affected area; corrections convert to strict
open water. Unvalidated shelves carry no correction and say so.

### 3.4 Radar wet-snow record
Per-pixel, per-relative-orbit frozen-winter baselines (Jun-Aug median,
<= 12 scenes evenly sampled); wet where gamma0 falls >= 2 dB below baseline
(threshold selected by gate sweep, decision 0002). 120 m grid derived from
the optical grid (same origin, x4). Four pre-registered gates: (1)
containment of blind-validated ponds in radar wet masks: 100.0% at 2 dB
(12,170 cells); (2) season severity ranking matches the optical record on
three obs-normalized metrics (the naive max-extent metric saturates - every
George VI January wets the full swath - and its failure is retained in the
gate record); (3) mid-winter false-wet 0.08% against a 2% ceiling; (4) melt
onset/end timing within climatological bounds (Nov 2 / Feb 24).

### 3.5 Fused melt states
Each radar acquisition classifies cells UNOBSERVED / DRY / WET / PONDED.
PONDED requires radar-wet AND optical pond evidence within +/-6 days
(evidence restricted to the Dec-Feb melt core; decision 0004): PONDED is a
subset of WET by construction. Optical-pond-but-radar-dry is a counted
conflict, split by the cell's own radar phenology into wet-season conflicts
(true sensor disagreement, gated < 5% where ponding instances >= 500) and
off-season conflicts (the frozen-lid / diurnal-refreeze signal - reported as
a product feature mapping stored water under ice). Per-cell confidence =
clip(margin/6 dB) x exp(-gap/6 d).

### 3.6 Depth and volume
Dual-band (red primary <= 3 m validity; green secondary) Beer-Lambert
inversion with per-pond rim albedo (windowed per pond; never a global
constant) and R_inf = 0.05 [0.03-0.07]. The attenuation coefficient g is
never hand-picked: 111 ICESat-2 ATL03 crossover depths (strong beams,
histogram peak splitting, x0.752 refraction, QC 0.25-8 m, photon-in-pond-
mask selection) calibrate it. Tight (+/-3 d) window: g = 1.62-1.64 m^-1,
RMSE 0.557 m (n = 29); +/-1 d does not improve it, locating the residual in
segment-vs-pixel pairing scale, not pond evolution. Independently fitted
sites agree to 3% (George VI 1.19, Amery 1.15, validity-domain fits). The
0.5 m RMSE target is missed by 12% and reported as missed (decision 0003);
Monte Carlo volume uncertainty carries the calibration posterior
(sigma_g = 0.2 g), rim-albedo spread, and the R_inf range. Depth retrieval
is restricted to fused-PONDED cells, excluding lids and slush where
attenuation physics does not apply.

### 3.7 Drainage catalogue
Between consecutive evidence dates <= 14 days apart, pond components losing
>= 80% of cells - required to be >= 70% observed on the after-date
(coverage-blind detection produced hundreds of phantom events; the
constraint is part of the method) - flag candidates, classified drainage
(radar stays wet >= 3 d beyond the loss) vs refreeze (radar froze in the
window).

## 4. Results
- 12-shelf x 9-season record: 108 shelf-seasons, 3 flagged poorly observed.
- George VI fused states (all pilot checks pass): peak ponded 796.6 /
  186.2 / 0.5 km2 for 2019-20 / 2020-21 / 2021-22, reproducing the optical
  ranking (838 / 159 / 12 km2); wet-season conflict 1.09% / 1.25% (quiet
  season below the 500-instance floor); PONDED-outside-optical cells: zero
  in all seasons. Wet-to-ponded peak ratios 10.2 / 43.3 / ~16,000: the
  ponded (hydrofracture-relevant) fraction collapses far faster than wetted
  area as seasons cool - vulnerability is nonlinear in melt.
- Volumes (g = 1.62, PONDED-restricted): 2019-20 0.129 km3 [0.108-0.166];
  2020-21 0.050 km3 [0.042-0.064]. Centre-deeper-than-rim: 97% and 99% of
  ponds (n = 1,200 and 733).
- Volume-per-area vs Corr et al. (2022) assumed 0.716 m scaling: 2020-21
  0.42 m (1.7x, within the x2 gate); record 2019-20 0.24 m (3.0x outside) -
  a laser-calibrated measurement of the record season's shallow marginal
  ponding vs an assumed constant; we report the disagreement rather than
  tune to it.
- Drainage catalogue: zero events survive coverage screening at v1's
  clear-scene cadence (a coverage-blind detector produced 784 phantoms; the
  documented 2019-20 events fall between evidence dates >14 d apart) -
  detection at this cadence is a stated limitation, not a claim.

## 5. Limitations
Single-labeler validation (a three-labeler weighted protocol per Leeson,
pers. comm., is planned); depth RMSE 0.56 m at segment scale (0.5 m target
missed); red-band depths saturate beyond ~3 m; PONDED requires 25% cell
fill, insensitive to diffuse water thinner than the 120 m grid; seven
shelves carry no per-regime validation yet; 2017-18 optical coverage is
sparse (pre-S2B constellation completion).

## 6. Data and code availability
https://github.com/Advait-Sadineni/antarctic-melt-observatory
DOI: 10.5281/zenodo.21711608 · Dashboard: [PAGES_URL]

## Acknowledgments
A. Leeson (Lancaster) suggested the multi-labeler uncertainty protocol
adopted in the roadmap. This project was built with AI-assisted engineering
(Claude); methodology, validation design, and all scientific decisions were
reviewed and are documented in the repository's decision records.

## Key references
Moussavi et al. 2020 (Remote Sens.); Corr et al. 2022 (ESSD 14, 209);
Banwell et al. 2021 (GVI record ponding); Parrish et al. 2019 (ICESat-2
bathymetry); Pope et al. 2016 (attenuation depths); Deakin et al. 2025
(GVI winter SAR); Nagler & Rott lineage (wet-snow SAR); DeConto & Pollard
2016; Scambos et al. 2000 (Larsen B).

# Antarctic Surface Meltwater Detection

Detects surface meltwater ponds on Antarctic ice shelves from free Sentinel-2
satellite imagery. Meltwater matters because ponds drive *hydrofracture*: water
pooling on an ice shelf wedges open crevasses under its own weight and can
disintegrate the shelf, as happened to Larsen B in 2002. Once a shelf is gone,
the glaciers it buttressed accelerate into the ocean.

Current scope: single-scene detection, per-season pond-area time series across
nine melt seasons (2017-18 to 2025-26), a multi-year peak comparison, a
threshold sensitivity analysis, and two independent validations (cross-sensor
against Landsat 8, and against blind visual interpretation) — all over one
61 km study area on **George VI Ice Shelf** (Antarctic Peninsula). Roadmap
below covers full-shelf coverage, ICESat-2 elevation, and ERA5-driven
prediction.

## Run it

```
pip install rasterio matplotlib pystac-client numpy
```

| Command | Output |
|---|---|
| `python src/demo.py` | `output/george_vi_meltwater.png` — true colour beside the detected mask |
| `python src/tune_threshold.py` | `output/threshold_sensitivity.png` — how much the answer depends on the threshold |
| `python src/season.py` | `output/season_2020_21.png` + `.csv` — pond area across one melt season |
| `python src/season.py 2019-20` | same, for any season |
| `python src/multiyear.py` | `output/multiyear_trend.png` — peak meltwater by season |
| `python src/validate_landsat.py` | cross-sensor check against Landsat 8 |
| `python src/validate_points.py make` / `score` | blind visual accuracy check |
| `python src/retune.py eval` | score candidate detection methods vs blind labels |
| `python src/uncertainty.py` | `output/uncertainty_peaks.png` — peaks as ranges, not points |
| `python src/export_gis.py <scene>` | `output/gis/` — GeoTIFF mask + GeoJSON ponds for QGIS |

Install: `pip install -r requirements.txt`. Tests: `python -m pytest tests -q`
(46 tests, offline, ~1 s).

No credentials required. Imagery comes from the **AWS Open Data** Sentinel-2 L2A
archive as Cloud-Optimized GeoTIFFs, found via the Earth Search STAC API. Each
script reads only a 2048×2048 pixel window directly out of the remote COGs
rather than downloading the full 110 km tile, and caches it to `data/`.
`season.py` is resumable — it appends to the CSV and skips scenes already done.

## How it works

**NDWI** (Normalized Difference Water Index) = `(Green − NIR) / (Green + NIR)`,
from Sentinel-2 bands B03 and B08. Liquid water absorbs near-infrared very
strongly but still reflects green, so NDWI goes positive over ponds. Snow and
ice reflect both bands strongly and sit near or below zero. On a white ice
shelf the two separate cleanly, which is why a plain threshold works here at all.

Detection is the published **Moussavi et al. 2020** method for Antarctic
supraglacial lakes: `NDWI > 0.19`, a **shadow test** (green − red > 0.09
reflectance) that rejects crevasse and topographic shadow, a brightness floor
(green > 3000 DN), and SCL-based cloud rejection. Pond cores are then grown by
**hysteresis** into connected margin pixels above `NDWI > 0.14`, which recovers
under-drawn pond edges without flooding disconnected crevasse fields. This
configuration was chosen by validation, not taste — see
[Retuning against blind labels](#retuning-against-blind-labels).

### Two findings worth knowing

**SCL cannot be used to find the ponds — only to reject cloud.** Every L2A
product ships a Scene Classification Layer. Measured on this study area, SCL
labels **80% of visually-confirmed melt ponds as "snow/ice" rather than
"water"**, because its water class is tuned for open ocean and deep lakes, not
shallow meltwater lying on bright ice. Using SCL class 6 as the pond detector
would discard four fifths of the signal. It is used strictly as a reject mask
(cloud, cirrus, shadow, saturated, nodata).

**Cloud shadow is the dangerous false positive, and SCL's *water* class is what
catches it.** The first season run produced a 46.9 km² spike on 2020-11-04 —
before melt onset, against ~0.02 km² for every neighbouring November scene.
The cause was cloud shadow lying across the shelf, which defeats both guards
at once: shadow on ice is lit by blue skylight, raising green relative to NIR
until it mimics water's NDWI signature, while haze lifts green back above the
brightness floor. SCL's *cloud* classes under-flag shadow over bright ice, so
the scene looked clean.

The tell is SCL's *water* class. On an ice-shelf AOI it should be near zero:
clear scenes here measure 0.00–0.42%, and that upper value is itself real
ponds (SCL calls ~20% of pond pixels water, and peak ponds are ~2.1% of the
AOI). The contaminated scene reported **7.7% water**, implying an impossible
~38% pond cover. `season.py` now rejects any scene above 2%.

**Scene-level cloud metadata is the wrong filter.** `eo:cloud_cover` describes
the whole 110 km tile. One scene here reports 37% tile cloud while the study
area is 97% clear. `season.py` therefore computes cloud fraction *inside the
AOI* from SCL and filters on that, which keeps scenes that tile-level filtering
would have thrown away.

**SCL badly under-reports cloud over ice, so NDSI does that job instead.**
White cloud on a white shelf is genuinely hard for sen2cor. Three February
2021 scenes that SCL called 2–19% cloudy turned out, on inspection, to be
blanketed edge to edge, and two produced multi-km² phantom ponds along the
cloud edges — including a 5.04 km² reading whose mask sat squarely on bright
cloud rather than on any pond.

The fix is **NDSI** = `(Green − SWIR) / (Green + SWIR)`, from B03 and B11
(1610 nm). Snow, ice and liquid water all absorb SWIR strongly and sit high
(0.85–0.95 here); water cloud *reflects* SWIR and drops. Measured over this
AOI, clean scenes put 0.0–0.1% of pixels below NDSI 0.6 while cloud-covered
ones put 3.7–35.9% there, so `melt.ice_check()` gates scenes at 1%.

Thresholding NDSI per-pixel does *not* remove the phantoms — they sit beside
the cloud rather than under it and score high NDSI themselves. Dilating the
detected cloud by a **1 km halo** does remove them. On a confirmed cloudy
scene that cuts the phantom from 8.16 to 2.32 km² while costing a confirmed
clean scene 0.5% of its real detections.

The halo replaced an earlier whole-scene NDSI threshold, which turned out not
to survive a change of study area: a cloud bank that filled the old 20 km
window covers only ~11% of the current 61 km one, so a gate calibrated at 1%
silently stopped firing when the AOI grew. A local mask has no such
dependence.

**Grazing sunlight fakes meltwater.** Scenes at 9.5° and 11.7° solar elevation
returned the largest "melt" of their season — in March, at the end of the melt
season — with detections tracking the crevassed eastern margin rather than any
pond. At low sun, shadow from crevasses and surface topography lengthens, and
shadowed snow lit by blue skylight reproduces water's NDWI signature. Scenes
below **15° solar elevation** are rejected.

**A plausibility tripwire catches what the physics-based tests miss.** The
largest meltwater density on a visually confirmed clean scene is 0.49% of the
study area, so anything above 2% is treated as a detector failure rather than
a record melt event, and logged loudly. It has earned its place: 2020-01-29
reported 2.78% with a halo of only 0.23%, meaning the NDSI test saw almost no
cloud. Inspection showed a large dark cloud-shadow band with the mask sitting
squarely on it — shadow, not bright cloud, so NDSI was blind to it. The
tripwire was the only thing standing between that scene and the record books,
in a season that genuinely was a record-warm Antarctic Peninsula summer, which
is exactly when a wrong answer would have been most believable.

**Partial cloud no longer discards the whole scene.** Demanding a clear 61 km
square is far stricter than demanding a clear 20 km one — the 2026-02-10
record scene is 23.7% cloudy over the large area but only 5.8% over the small
one, so a whole-scene cloud gate was throwing away 76% of usable ground. A
scene is now kept if at least half the study area survives screening, and the
result is reported as a **density over the surviving ground**, scaled back to
the full area (`pond_km2_equiv`). That keeps scenes with different usable
footprints comparable.

### Threshold sensitivity

`tune_threshold.py` sweeps NDWI 0.04–0.40. The curve has three regimes:

- **Below ~0.08** — saturation. Bare ice itself goes slightly positive and
  essentially the entire 3775 km² study area is flagged as water. Unusable.
- **0.08–0.14** — a cliff, area dropping 30–90% per 0.02 step.
- **Above ~0.16** — a smooth decay easing from ~15% to ~6% per 0.02.

At the default 0.16, local sensitivity is **11% of area per 0.01 of
threshold**. Across the published range 0.10–0.25 the answer moves from 41.9
to ~6.3 km² — a factor of about seven.

That is markedly worse than the factor of 2.5 measured on the old hand-picked
window, and the reason is instructive: that window was almost entirely pond or
bare ice, so it sat well clear of the cliff. A representative area contains a
great deal of marginal, near-threshold surface, so the answer depends more on
where the line is drawn. The honest reading is that the *shape* of a season —
onset, peak timing, relative magnitude — is far more robust than any absolute
area, and the threshold must always be quoted alongside the number.

The practical consequence: **absolute pond area is threshold-dependent and
should be quoted with its threshold.** Comparisons across dates at a *fixed*
threshold are far more trustworthy than any single absolute number, which is
why the time series is the more defensible product.

### Melt season 2020-21 result

Of 94 Sentinel-2 scenes on the tile, **18 survived screening** (the rest were
cloud, off-swath nodata, or contamination-gated). The resulting curve is
physically coherent: near-zero through November and December, a sharp rise
through January, a **peak of 8.80 km² on 24 January 2021** at peak austral
summer, then decline through February and March.

One feature is *not* yet explained: area drops from 6.61 km² (13 Feb) to
0.23 km² (16 Feb) and rebounds to 5.04 km² (23 Feb). Fresh snowfall burying
ponds, or a brief refreeze, would both produce this and are well documented
in the literature — but it is equally consistent with undetected thin cloud.
It is left in the data rather than smoothed away. Treat it as an open question.

### Multi-year trend, 2017-18 to 2025-26

`multiyear.py` runs all nine Sentinel-2 melt seasons. Peak observed meltwater
area over the study area:

Areas are `pond_km2_equiv`: meltwater density over the ground that survived
screening, scaled to the full 3775 km² study area.

Numbers below use the adopted method (Moussavi shadow test + hysteresis).

| Season | Peak km² | Peak date | Usable scenes | |
|---|---|---|---|---|
| 2017-18 | 6.94 | 2018-02-22 | 2 | likely undercounted |
| 2018-19 | 4.37 | 2019-01-25 | 7 | |
| 2019-20 | 8.20 | 2020-01-19 | 4 | literature-confirmed peak day |
| 2020-21 | 9.90 | 2021-01-24 | 8 | peak scene visually verified |
| 2021-22 | 3.54 | 2022-01-12 | 9 | low-melt season |
| 2022-23 | 8.70 | 2023-02-22 | 6 | |
| 2023-24 | — | — | 0 | no scene passed screening |
| 2024-25 | 7.35 | 2025-02-05 | 3 | likely undercounted |
| 2025-26 | **16.98** | 2026-02-10 | 7 | largest in record, verified |

2025-26 is the largest in the record and its peak scene was inspected
directly rather than trusted — clean, with unmistakable dark meltwater
channels and the mask sitting exactly on them.

**The 2019-20 peak falls on 19 January 2020 — the exact day Banwell et al.
(2021) independently report as the record season's peak ponding** (see
[Literature comparison](#literature-comparison)). The pipeline reproduced that
date without ever using the paper.

**Do not read a trend off this chart.** Eight plotted seasons, three of them
undersampled and one with no usable data at all, cannot separate a trend from
year-to-year weather over a single 61 km window.

Note how much these numbers moved as the method matured. An earlier version of
this table had 2024-25 peaking at 34.2 km² (a cumulus band), 2019-20 at 8.24
from a cloud-shadow scene, and 2018-19 at 11.56 from a crevasse-shadow scene
that later failed cross-sensor validation — the shadow test now cuts that
2018-19 peak to a defensible 4.37. Every one looked reasonable at the time,
which is why the pipeline is now anchored to blind labels and a second sensor.

Three caveats that matter for reading the table:

- **There is a noise floor at 1.3 km².** Clean pre-melt November scenes report
  0.7–1.2 km² of residual detections; the floor sits just above that. With the
  shadow test the background dropped sharply (it was mostly crevasse shadow),
  and every plotted season now peaks well above the floor — no season reads
  "no melt detected" any more.
- **Peak observed ≠ seasonal peak.** If the true maximum fell on a cloudy day
  it is simply missed. Seasons with fewer than five clear scenes are greyed;
  seasons with fewer than two are dropped entirely, which is what happened to
  2017-18 (one scene, in March).
- **Sparsity is the weather, not the filter.** For 2023-24 (2 usable of 95),
  the rejections break down as 53 cloud, 35 off-swath nodata, 4 water, and
  only **1** from the NDSI gate. Screening is not what is throwing the data
  away. The 35 nodata scenes are partial-swath and are **not** recoverable by
  mosaicking: checked all five multi-scene dates in 2023-24, and in every case
  the same-date companion has the same footprint rather than a complementary
  one. When the AOI sits at an orbit edge, nobody imaged the missing ground
  that day.
- **Same-date scenes are true duplicates, not complementary swaths.** Checked
  before collapsing them: in every observed case both cover 376–419 km² of the
  ~419 km² AOI, so keeping the better one discards no coverage.

## Configuration

Everything lives in `src/melt.py`:

| Setting | Value | Why |
|---|---|---|
| Study area | tile `19CEV`, 6144×6144 px @ 10 m | ~61 km square, 3775 km², ~72.1°S 69°W |
| Screening resolution | 60 m | fractions only; a fraction of the bytes |
| NDWI threshold (core) | `0.19` | Moussavi et al. 2020; best F1 on blind labels |
| Shadow test | green − red > 0.09 | rejects crevasse/topographic shadow |
| Hysteresis grow | `NDWI > 0.14`, connected | recovers pond margins |
| Brightness floor | green > 3000 DN | rejects dark ocean/shadow |
| NDSI ice minimum | `0.60` | below this a pixel is not snow/ice/water |
| Cloud halo | 1 km | cloud contaminates ground beside it |
| Min solar elevation | 15° | grazing light fakes meltwater |
| Min usable fraction | 50% | keeps partly-cloudy scenes, normalised |
| Max halo fraction | 5% | heavy detected cloud implies undetected cloud |
| Max plausible density | 2% of usable | tripwire for unmodelled failures |
| Noise floor | 1.3 km² | above the pre-melt baseline |

Every one of these was set from a measurement, and most exist because a
specific scene produced a specific wrong answer. None is a default.

The study area is pinned to one Sentinel-2 tile on purpose. The S2 tiling grid
is fixed, so the same pixel window on tile 19CEV is the same ground on every
date — that is what makes areas comparable across a time series.

### The study area is chosen for coverage, not for ponds

The first version of this project used a 20 km window found by searching the
tile for the *densest ponding*. That makes every absolute number a best case,
and it is the kind of choice that quietly invalidates a result.

The current 61 km area is instead the largest window that is fully imaged and
100% snow/ice across four clear scenes drawn from different seasons (2019,
2020, 2021, 2026). Pond density played no part in selecting it, so it
includes a great deal of shelf that never ponds.

The difference is large enough to be worth stating plainly. On the same scene
(2021-01-24), the hand-picked window gave **2.10%** meltwater cover while the
representative area gives **0.302%** — a factor of seven. Any pond-fraction
number from the old window should be read as roughly 7× optimistic.

Screening runs at 60 m and only scenes that pass pay for full-resolution
reads, which is what keeps a nine-season run to about half an hour despite
the area being nine times larger.

## Validation

```
python src/validate_landsat.py
```

The pipeline was internally consistent but externally unanchored — every
threshold in it was calibrated against its own output. This compares it
against **Landsat 8**, an independent satellite with independent optics and an
independent processing chain, on the same ground on the same day.

The comparison uses the **sensor-transferable core** of the detector — NDWI +
hysteresis — applied identically to both. The shadow test is *excluded* here:
its 0.09 green−red cutoff is calibrated on Sentinel-2 bottom-of-atmosphere
reflectance, and Landsat is only available as top-of-atmosphere over ice,
where atmospheric haze compresses green−red and the same cutoff rejects nearly
all real water. The shadow test is validated instead by the blind reference
points. So this check answers one clean question: *do two independent
satellites compute the same NDWI water on the same ground?*

11 usable pairs, 2018–2021 (Landsat Collection-1 ends in 2021).

### Result

| Measure | All 11 pairs | Clean melt pairs |
|---|---|---|
| Pearson *r* on areas | **0.92** | |
| Spearman rank correlation | **0.85** | |
| Area agreement (Landsat ÷ Sentinel-2, 30 m) | 1.09 | 0.94 |
| IoU at matched 30 m | 0.33 | up to 0.58–0.62 |
| Agreement within 1 Landsat pixel | 0.38 | |

**How much meltwater there is, the pipeline gets right — now near-perfectly.**
Across all pairs the area correlation is *r* = 0.92 and the Landsat/Sentinel-2
area ratio is 1.09 (0.94 on clean pairs) — the two satellites agree on area to
within a few percent. Both improved over the pre-retune method (*r* 0.83,
ratio 0.90), because hysteresis made the two sensors draw pond margins more
alike.

**Which pixels, it gets right moderately** — IoU matched 0.33 median, reaching
**0.58–0.62 on the strongest melt scenes** (2019-01-25, 2020-01-19). Most of
the residual is resolution: a 10 m channel is a sub-pixel mixture to Landsat,
so its NDWI is diluted no matter how good either instrument is. The rest is
inter-sensor co-registration (~1 Landsat pixel) and real change between
overpasses.

### What validation caught

**2019-02-17 stays anomalous even in the core comparison.** This is the scene
that, at the old 0.16 threshold, gave a phantom 11.6 km² and made 2018-19 the
second-highest season. The shadow test now cuts its pipeline area to a
defensible **2.9 km²** — no longer the season peak. Cross-sensor agreement on
it remains the worst of any clean same-day pair (IoU matched **0.22**, against
0.58–0.62 for genuine melt scenes), independently confirming it as diffuse
marginal surface — slush or wet snow — rather than standing water. Two
independent lines (the shadow test and the second sensor) now agree on what a
single threshold got wrong.

### Against visual interpretation

```
python src/validate_points.py make     # build blind chips
python src/validate_points.py score    # score labels against the hidden key
```

Cross-sensor agreement shows two instruments measure the same thing. It cannot
show either is *right* — they can agree on a wrong answer. This compares the
mask against visual interpretation of the true-colour imagery.

96 points on the 2021-01-24 scene, stratified random (detected / near-detected
/ far), labelled **blind**: the chips are shuffled and anonymously numbered,
nothing reveals the detector's opinion, and the answer key is not opened until
scoring. Estimates are area-weighted, so the rarity of meltwater is respected.

| Metric | Value |
|---|---|
| Precision | **0.63** (95% CI 0.47–0.77) |
| Recall, area-weighted | **0.36** |
| F1 | 0.46 |
| Area bias | detector reports 0.57× the sampled truth |

**This is the least flattering result in the project, and the most useful.**
About a third of detected pixels are not meltwater, and roughly two thirds of
meltwater is missed.

**The false positives have a single cause: crevasse shadow.** All 14 of them
sit in heavily crevassed terrain, with the mask tracing the dark stripes of
the crevasse pattern. And they are numerically distinct from real ponds:

| Visual label | n | NDWI p25 | median | p75 | max |
|---|---|---|---|---|---|
| water | 28 | 0.239 | **0.382** | 0.505 | 0.686 |
| not water | 61 | 0.068 | **0.074** | 0.118 | **0.245** |

Every false positive falls between NDWI 0.168 and 0.245 — just above the 0.16
threshold — while genuine ponds sit far higher. The two classes separate
cleanly; the threshold is simply drawn inside the overlap.

Sweeping it against the labelled sample showed precision climbing from 0.63 at
0.16 to 1.00 by 0.25 as the crevasse-shadow band is excluded. Rather than
retune on 96 points from a single scene — the overfitting this project has
tried to avoid — this was recorded as a diagnosis and settled properly by
labelling two more scenes and comparing published methods, below.

### Retuning against blind labels

```
python src/retune.py chips     # blind chips for the extra scenes
python src/retune.py eval      # score every candidate method vs the labels
```

The diagnosis above pointed at the threshold; the literature pointed at the
same thing and named the fix. Nine candidate detection rules were scored
against **220 blind-labelled points across three scenes** (2021-01-24,
2020-01-19, 2026-02-10), area-weighted:

| Method | Precision | Recall (area) | F1 |
|---|---|---|---|
| current — NDWI > 0.16 | 0.55 | 0.40 | 0.46 |
| raise to 0.19 (Moussavi) | 0.68 | 0.39 | 0.50 |
| dual-NDWI (Corr) | 0.79 | 0.30 | 0.43 |
| 0.19 + shadow test | 0.71 | 0.39 | 0.50 |
| **0.19 + shadow + hysteresis (adopted)** | **0.68** | **0.45** | **0.54** |

The adopted method is the published **Moussavi et al. 2020** rule —
`NDWI > 0.19` with a `green − red > 0.09` shadow test — extended with
hysteresis margin recovery. It wins on F1 and lifts precision from 0.55 to
0.68 by removing exactly the crevasse-shadow false positives the diagnosis
identified. Adopting 0.19 also aligns the pipeline with the same papers whose
19 Jan 2020 peak date it reproduced.

**On recall, two numbers tell different truths.** By *count*, the detector
finds **0.89** of labelled water points — on par with published methods
(0.77–0.85). By *area*, recall is **0.45**, because every missed point is a
**pond margin** (stratum B), never a missed lake (stratum C had zero water)
and never a detection core. The detector finds ponds and under-draws their
edges; hysteresis recovered part of that (area recall 0.39 → 0.45), and the
residual is the sub-pixel limit of 10 m optical. Net area bias is **0.67×** —
the reported areas are more likely under- than over-estimates.

### Literature comparison

The pipeline was checked against the published George VI literature. All
figures below were verified against the source abstracts/papers; the peak-date
and method comparisons are the strongest anchors, absolute areas the weakest
(different study areas).

**Peak timing — an exact, independent match.** Banwell et al. (2021, *The
Cryosphere* 15:909) documented the record 2019-20 melt season on northern
George VI and dated peak ponding to **19 January 2020**. This pipeline,
built without reference to that paper, places the 2019-20 seasonal maximum on
**2020-01-19** — the same day. Dirscherl et al. (2021) likewise report the
regional peak in "late January", declining through February, which matches the
pipeline's per-season curves.

**Method — we now use the published detector.** The retune independently
converged on **Moussavi et al. (2020)**: NDWI > 0.19 for Sentinel-2 with a
green−red shadow test — the most-cited Antarctic supraglacial-lake method.
Corr et al. (2022) use a dual index (NDWI > 0.16 *and* blue−red > 0.18); that
was among the nine candidates tested and scored highest precision (0.79) but
lowest recall (0.30), so it was not adopted.

**Accuracy — in the published range.** The pipeline's count-recall of 0.89
sits at/above Corr's reported Sentinel-2 sensitivity of 85.3% and above their
Landsat 77.6%. Dirscherl's deep-learning F1 of ~0.95 remains ahead of this
threshold method — the motivation for the ML roadmap item.

**Absolute area — not yet comparable, by construction.** Banwell reports a
peak of **~1,200 km²** over a 7,850 km² northern-shelf AOI (≈15–23% cover);
Dirscherl **~805 km²** whole-shelf; Corr **29.4 km²** for a quiet early-January
2017 scene. This pipeline's 61 km window (3,775 km²) is a *subset* of the
northern shelf, so its 8–17 km² peaks cannot be compared to whole-shelf totals
until the AOI is expanded (roadmap item 2). What *is* comparable — peak timing,
method, per-pixel accuracy, and the ~40× dynamic range between quiet and record
seasons — all agree.

*Sources: Banwell et al. 2021 (tc-15-909-2021); Dirscherl et al. 2021
(tc-15-5205-2021); Corr et al. 2022 (essd-14-209-2022); Moussavi et al. 2020
(rs12010134). Absolute-area figures are from different study areas and are
context, not a pass/fail target.*

### Limits of this validation

- **Processing levels differ, so the shadow test is not cross-validated.**
  Sentinel-2 L2A is bottom-of-atmosphere, Landsat Level-1 top-of-atmosphere;
  the shadow test's absolute cutoff does not transfer, so the cross-sensor
  check covers only the NDWI + hysteresis core. The shadow test's validation
  comes from the reference points instead.
- **Small sample.** 11 pairs, 4 of them clean same-day; the strong *r* = 0.92
  rests on a handful of well-covered melt scenes.
- **Landsat Collection-1 ends in 2021**, so the 2022-26 seasons — including
  the record 2025-26 — have no cross-sensor check. Collection-2 Level-1 was
  not reachable without an account; Collection-2 *Level-2* was reachable and
  turned out to be unusable (below).
- **The visual reference is one interpreter.** Across all three labelled
  scenes, hard cases (water-filled crevasse vs crevasse shadow, slush vs pond)
  were labelled "unsure" and excluded, but a single interpreter remains a real
  limitation. The physical argument favours the labels — shadow dims all bands
  together and lifts NDWI only slightly, whereas water absorbs NIR strongly and
  pushes NDWI past 0.4.
- **Landsat Level-2 surface reflectance is invalid over ice.** Measured here,
  96.9% of its pixels over the study area fall outside the product's own
  documented valid range, reading green reflectance of 1.22 where reflectance
  cannot exceed 1.0. The atmospheric correction is not designed for bright
  cryospheric surfaces. Using it would have produced an official-looking
  validation built on unphysical numbers, so this module uses Level-1 and
  computes top-of-atmosphere reflectance itself.

### Peaks as ranges, not points

```
python src/uncertainty.py
```

Every headline figure rests on choices that were made, not measured.
`uncertainty.py` re-thresholds each peak scene across the NDWI band 0.16–0.25
— the range spanned by published methods — with the full shadow test and
hysteresis applied, and reports the span.

Before the retune this range was **98% of the point estimate** — threshold
choice alone swung each peak by ±half. **The shadow test collapsed it to a
median of 24%**, because the crevasse-shadow pixels that used to move with the
threshold are now removed regardless of it. The numbers are far more stable.

| Season | Point (0.19) | Range (0.16–0.25) | Coverage |
|---|---|---|---|
| 2017-18 | 6.94 | 5.9 – 7.6 | 72% |
| 2018-19 | 4.37 | 2.7 – 5.2 | 100% |
| 2019-20 | 8.20 | 7.0 – 8.9 | 100% |
| 2020-21 | 9.90 | 8.8 – 10.6 | 100% |
| 2021-22 | 3.54 | 1.5 – 4.7 | 97% |
| 2022-23 | 8.70 | 7.6 – 9.5 | 91% |
| 2024-25 | 7.35 | 5.7 – 8.4 | 53% |
| 2025-26 | **16.98** | 15.6 – 17.9 | 76% |

- **2025-26 is robustly the largest** — its whole range clears every other
  season's point estimate.
- **The rankings are now stable under threshold choice**, where before the
  retune 2018-19's range alone spanned 0.9–22.3 km². The one wide band left is
  2021-22 (a genuine low-melt season near the detection limit).

A second, independent source of error runs the other way. Blind reference-point
sampling gives an area bias of **0.67×** — the detector under-reports, because
the pond-margin area it still misses outweighs its false detections. That is
not folded into the ranges above, but it means true areas likely sit toward or
above the upper end of each range.

The takeaway the project now states plainly: **trust the ranking and the
shape of a season, quote the threshold with any number, and treat a single
absolute km² as indicative only.**

## Tests

```
python -m pytest tests -q
```

43 tests, no network — synthetic arrays and fake STAC items, so they run
offline in under a second. They cover the parts where a silent wrong answer
is possible: NDWI sign conventions, the shadow test and hysteresis margin
growth, nodata never becoming water, threshold monotonicity, the
pixel-count-to-km² conversion, division-by-zero on a fully masked scene, all
four branches of the baseline 04.00 offset, SCL fraction accounting, the
same-date dedup rule, and the cross-sensor reflectance conversions.

Several are regression guards for bugs found the hard way: that SCL class 6 is
never added to the reject set (it would silently delete 80% of real ponds),
that the noise floor stays above the pre-melt baseline, and that the adopted
threshold matches the published Moussavi value.

## Known limitations

- **Study area is one 61 km window, not the whole shelf.** Absolute numbers
  are not shelf-wide melt statistics, and cannot yet be compared to published
  whole-shelf areas until the AOI is expanded (roadmap item 1).
- **Area recall is ~0.45.** The detector finds 0.89 of water *features* but
  under-draws pond *margins*, so absolute area is under-reported by ~0.67×.
  Recovering the sub-pixel margin fraction (spectral unmixing) is the main
  remaining accuracy item.
- **Cloud reduces usable area**, which biases absolute pond area low on hazier
  days. The CSV records `usable_frac` and `pond_pct_of_usable` so this is
  visible rather than hidden.
- **No slope or projection-distortion correction.** Area is
  `pixel count × 100 m²`; at 72°S in UTM 19S the error is small but nonzero.
- **NDWI cannot fully separate shallow ponds from wet/slushy snow.** The
  shadow test removes crevasse shadow but not slush. 2019-02-17 is the clearest
  case: the shadow test cut its area from a phantom 11.6 km² to a defensible
  2.9 km², but cross-sensor agreement on it stays poor (IoU 0.22 vs ~0.58 for
  clean scenes), consistent with diffuse marginal surface. Left in and flagged.
- **No comparison against published George VI numbers yet.** Cross-sensor
  validation anchors the pipeline against another satellite; the literature
  peak *date* matches (19 Jan 2020), but absolute-area comparison needs the
  full-shelf AOI. This is the largest outstanding gap.

Handled correctly: the **Sentinel-2 baseline 04.00 offset**. From 2022-01-25
onward, L2A reflectance is shifted by −1000 DN, and since NDWI is a ratio of
differences an additive offset genuinely changes the result.

In this archive the correction turns out to be a no-op: Earth Search sets
`earthsearch:boa_offset_applied: true` and has already removed it, which
`melt.boa_offset()` honours. Verified rather than assumed — median snow DN is
~10000 on both sides of the 2022 boundary, so the pixel values really are
harmonised, and applying −1000 again would have inflated NDWI on every
post-2022 scene. Note the STAC `raster:bands` metadata still advertises
`offset: -0.1`, which contradicts the flag; the pixels agree with the flag.
The code is kept defensive so a non-harmonised source would still be correct.

## Accounts (still not needed)

Everything above runs anonymously. Register these for the roadmap:

| Service | For | Credentials go |
|---|---|---|
| [Copernicus Data Space](https://dataspace.copernicus.eu/) | Full S2 catalogue, newer than the AWS mirror | `.env` → `CDSE_USER` / `CDSE_PASSWORD` |
| [NASA Earthdata](https://urs.earthdata.nasa.gov/) | ICESat-2 ATL06/ATL11 via `icepyx` | `~/.netrc` |
| [Copernicus CDS](https://cds.climate.copernicus.eu/) | ERA5 reanalysis | `~/.cdsapirc` |

Google Earth Engine is deliberately skipped — it would move computation off the
local machine and cost the pipeline its "plain local script" property.

## Layout

```
src/melt.py             core: scene access, cloud screening, NDWI detection
src/demo.py             one scene -> side-by-side PNG
src/season.py           one melt season -> pond-area time series
src/multiyear.py        all seasons -> peak meltwater trend
src/tune_threshold.py   threshold sensitivity sweep
data/                   cached GeoTIFF windows (gitignored)
output/                 rendered PNGs and CSVs
```

## Roadmap

- **[done] Retune against labelled data.** Adopted the Moussavi shadow test +
  hysteresis after scoring nine methods on 220 blind labels. Precision 0.55 →
  0.68, area recall 0.39 → 0.45, threshold sensitivity 98% → 24%.
1. **Spectral unmixing (fractional water).** The main remaining accuracy
   ceiling: recover the sub-pixel pond-margin area the threshold misses,
   correcting the 0.67× under-report. Highest accuracy impact.
2. **Full-shelf coverage (multi-tile).** Expand from one 61 km window to the
   northern shelf, enabling direct comparison to Banwell's 1,200 km² and
   Dirscherl's 805 km². The main barrier to literature-level comparison.
3. **ML classifier** (random forest → U-Net) trained on the labels — the path
   to Dirscherl-level F1 (~0.95) on both precision and recall.
4. **ICESat-2 elevation** (`icepyx`) — pond depth and hence volume, which is
   what actually drives hydrofracture, rather than extent.
5. **ERA5 climate drivers** + prediction layer: seasonal weather → expected
   melt extent. (2019-20 was a record-warm Peninsula summer — a natural test.)
6. **Per-pixel cloud masking** (`s2cloudless` on L1C) to recover
   partially-clear scenes rejected whole today.

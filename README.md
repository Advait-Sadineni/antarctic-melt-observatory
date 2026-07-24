# Antarctic Surface Meltwater Detection

Detects surface meltwater ponds on Antarctic ice shelves from free Sentinel-2
satellite imagery. Meltwater matters because ponds drive *hydrofracture*: water
pooling on an ice shelf wedges open crevasses under its own weight and can
disintegrate the shelf, as happened to Larsen B in 2002. Once a shelf is gone,
the glaciers it buttressed accelerate into the ocean.

Current scope: single-scene detection, per-season pond-area time series across
nine melt seasons (2017-18 to 2025-26), a multi-year peak comparison, and a
threshold sensitivity analysis — all over one 20 km study area on **George VI
Ice Shelf** (Antarctic Peninsula). Roadmap below covers full-shelf coverage,
ICESat-2 elevation, and ERA5-driven prediction.

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

Detection is `NDWI > 0.16`, plus a brightness floor (green > 3000 DN) to reject
dark open ocean and deep shadow, plus SCL-based cloud rejection.

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

| Season | Peak km² | Peak date | Usable scenes | |
|---|---|---|---|---|
| 2017-18 | 8.90 | 2018-02-22 | 2 | likely undercounted |
| 2018-19 | *11.56* | 2019-02-17 | 7 | **peak fails cross-sensor validation** |
| 2019-20 | 10.03 | 2020-01-19 | 3 | likely undercounted |
| 2020-21 | 11.34 | 2021-01-24 | 8 | peak scene visually verified |
| 2021-22 | *(1.47)* | — | 8 | **no melt detected** |
| 2022-23 | 10.49 | 2023-02-26 | 6 | |
| 2023-24 | — | — | 0 | no scene passed screening |
| 2024-25 | 9.60 | 2025-02-05 | 3 | likely undercounted |
| 2025-26 | **18.66** | 2026-02-10 | 7 | largest in record, verified |

2025-26 is the largest in the record and its peak scene was inspected
directly rather than trusted — clean, with unmistakable dark meltwater
channels and the mask sitting exactly on them.

**Do not read a trend off this chart.** Eight plotted seasons, three of them
undersampled and one with no usable data at all, cannot separate a trend from
year-to-year weather over a single 61 km window.

Note how much these numbers moved as screening improved. An earlier version of
this table had 2024-25 peaking at 34.2 km², which turned out to be a cumulus
band; and 2019-20 at 8.24 km² from a scene later shown to be cloud shadow.
Both looked entirely reasonable at the time.

Four caveats that matter for reading the table:

- **There is a noise floor at ~0.7 km².** Two scenes that pass every screen
  still report meltwater *before melt onset* — 0.39 km² on 2021-11-09 and
  0.68 km² on 2022-11-25. Inspection shows both hugging the edges of a dark
  patch (bare or blue ice, or shadowed terrain) near the AOI boundary, not
  ponds. 0.68 km² is the largest confirmed non-melt detection, so anything
  below 0.7 km² is reported as "no melt detected" and drawn hollow in the
  chart. This is why 2021-22 and 2024-25 carry no number.
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
| NDWI threshold | `0.16` | tunable; see sensitivity above |
| Brightness floor | green > 3000 DN | rejects dark ocean/shadow |
| NDSI ice minimum | `0.60` | below this a pixel is not snow/ice/water |
| Cloud halo | 1 km | cloud contaminates ground beside it |
| Min solar elevation | 15° | grazing light fakes meltwater |
| Min usable fraction | 50% | keeps partly-cloudy scenes, normalised |
| Max halo fraction | 15% | heavy detected cloud implies undetected cloud |
| Max plausible density | 2% of usable | tripwire for unmodelled failures |
| Noise floor | 2.2 km² | above the pre-melt baseline |

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
independent processing chain, on the same ground on the same day. The same
detector is applied to both sensors, in reflectance units, so any difference
is attributable to the sensor rather than the method.

11 usable pairs, 2018–2021 (Landsat Collection-1 ends in 2021).

### Result

| Measure | All 11 pairs | Same-day, ≥80% shared (4) |
|---|---|---|
| Area agreement (Landsat ÷ Sentinel-2, at 30 m) | 1.34 | **0.90** |
| Pearson *r* on areas | **0.83** | |
| Spearman rank correlation | **0.86** | |
| IoU at native resolution | 0.12 | |
| IoU at matched 30 m | 0.30 | 0.35 |
| Agreement within 1 Landsat pixel | 0.39 | 0.43 |

**How much meltwater there is, the pipeline gets right.** On the cleanest
pairs the two satellites agree on area to within 10%, and across all pairs
the correlation is strong in both magnitude (*r* = 0.83) and ordering
(Spearman = 0.86). Since the ranking of dates and seasons is what this project
actually claims, that is the number that matters most.

**Exactly which pixels, it gets right only moderately** — IoU 0.30–0.35.
Most of that gap is resolution rather than error: IoU roughly **doubles**
(0.12 → 0.30) when Sentinel-2 is degraded to Landsat's 30 m grid before
detection, because a 10 m meltwater channel is a sub-pixel mixture to Landsat
and its NDWI is diluted below threshold no matter how good either instrument
is. The remainder is inter-sensor co-registration, which is itself of order
one Landsat pixel, and genuine surface change between overpasses.

### What validation caught

**2019-02-17 fails cross-sensor validation.** This is the 11.6 km² scene that
makes 2018-19 the second-highest season and that no internal screen could
explain — not cloud, not low sun, not nodata. It is a same-day pair with 98%
shared clear ground, so it should be one of the *best* comparisons available.
Instead it is the worst on every metric:

| | 2019-02-17 | 2019-01-25 (same season, same conditions) |
|---|---|---|
| IoU matched | **0.14** | 0.41 |
| Within 1 pixel | **0.28** | 0.77 |
| Landsat ÷ Sentinel-2 | **0.48** | 1.15 |

Landsat sees roughly half the area, in materially different places. That is
the signature of a diffuse, marginal, near-threshold surface — most likely
slush or wet snow, which NDWI cannot separate from ponded water — rather than
of standing meltwater, which reproduces well across sensors on the other
same-day pairs.

**Consequence: the 2018-19 peak of 11.56 km² should be treated as
unconfirmed**, and with it that season's rank. This is exactly the failure
mode that motivated doing validation at all, and it was invisible to every
internal check.

### Limits of this validation

- **Processing levels differ.** Sentinel-2 L2A is bottom-of-atmosphere;
  Landsat Level-1 is top-of-atmosphere. NDWI is a normalised ratio and far
  more robust to this than raw reflectance, but it is not immune, and the
  size of that effect has not yet been measured here.
- **Four clean pairs is a small sample.** The headline 0.90 area agreement
  rests on them.
- **Landsat Collection-1 ends in 2021**, so the 2022-26 seasons — including
  the record 2025-26 — have no cross-sensor check. Collection-2 Level-1 was
  not reachable without an account; Collection-2 *Level-2* was reachable and
  turned out to be unusable (below).
- **Landsat Level-2 surface reflectance is invalid over ice.** Measured here,
  96.9% of its pixels over the study area fall outside the product's own
  documented valid range, reading green reflectance of 1.22 where reflectance
  cannot exceed 1.0. The atmospheric correction is not designed for bright
  cryospheric surfaces. Using it would have produced an official-looking
  validation built on unphysical numbers, so this module uses Level-1 and
  computes top-of-atmosphere reflectance itself.

## Tests

```
python -m pytest tests -q
```

18 tests, no network — synthetic arrays and fake STAC items, so they run
offline in under a second. They cover the parts where a silent wrong answer
is possible: NDWI sign conventions, nodata never becoming water, threshold
monotonicity, the pixel-count-to-km² conversion, division-by-zero on a fully
masked scene, all four branches of the baseline 04.00 offset (including
unparseable metadata), SCL fraction accounting, and the same-date dedup rule.

Two are regression guards for bugs found the hard way: that SCL class 6 is
never added to the reject set (it would silently delete 80% of real ponds),
and that the noise floor stays above the largest confirmed artifact.

## Known limitations

- **Study area is one 20 km window, not the whole shelf.** It was picked for
  dense ponding. Absolute numbers are not shelf-wide melt statistics.
- **Cloud reduces usable area**, which biases absolute pond area low on hazier
  days. The CSV records `usable_frac` and `pond_pct_of_usable` so this is
  visible rather than hidden.
- **No slope or projection-distortion correction.** Area is
  `pixel count × 100 m²`; at 72°S in UTM 19S the error is small but nonzero.
- **NDWI cannot distinguish shallow ponds from wet/slushy snow**, a known
  ambiguity in the literature and a motivation for the planned ML classifier.
  The clearest case is 2019-02-17, which reports 11.6 km² over a large dark
  patch late in the season and survives every internal screen. Cross-sensor
  validation now gives an independent verdict on it: Landsat sees half the
  area in materially different places (IoU 0.14 against 0.41 for a comparable
  scene in the same season), which is the signature of diffuse marginal
  surface rather than standing water. The reading is left in the data and
  flagged, not tuned away — but the 2018-19 seasonal peak that rests on it is
  unconfirmed.
- **Blue ice and dark terrain can false-positive.** This is what sets the
  0.7 km² noise floor above. The brightness floor catches the worst of it,
  but it is not a physical discriminator.
- **No comparison against published George VI numbers yet.** Cross-sensor
  validation anchors the pipeline against another satellite, but not against
  the literature; that remains the largest outstanding gap.
- **Every threshold is calibrated on a handful of scenes.** NDWI 0.16,
  brightness floor 3000, NDSI gate 1%, water gate 2%, noise floor 0.7 km² —
  all set from one or two seasons. Since the gates decide which data survives,
  a bad calibration would quietly reshape the whole trend.

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

1. Expand the study area beyond one hand-picked window. This is the main
   barrier to comparing against published studies.
2. Per-pixel cloud masking (`s2cloudless` on L1C) instead of whole-scene
   rejection, recovering partially-clear scenes. Note that NDSI cannot do this
   job — the phantom pond pixels score high NDSI, which is why it only works
   as a scene gate.
4. ICESat-2 ATL06/ATL11 elevation overlay (`icepyx`) — pond depth and hence
   volume, which is what actually drives hydrofracture, rather than extent.
5. ERA5 climate drivers + ML classifier to replace NDWI thresholding, which
   also dissolves the threshold-sensitivity problem.
6. Prediction layer: seasonal weather → expected melt extent.

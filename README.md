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

It has to be a *scene* gate, not a pixel mask: the phantom pond pixels
themselves still score high NDSI, so masking per-pixel would not remove them.

The cost is temporal coverage. The 2020-21 season drops from 18 usable scenes
to 9. That is the right trade — the 9 are defensible and the discarded ones
included the phantoms — but it means gaps are wide, and a season's peak can
be missed entirely if it fell on a cloudy day.

### Threshold sensitivity

`tune_threshold.py` sweeps NDWI 0.04–0.40. The curve has three regimes:

- **Below ~0.08** — saturation. Bare ice itself goes slightly positive and the
  entire 419 km² AOI is flagged as water. Unusable.
- **0.08–0.12** — a cliff, area dropping ~80% per 0.02 step.
- **Above ~0.12** — a stable, smooth decay of roughly 5% per 0.02.

At the default 0.16, local sensitivity is **4% of area per 0.01 of threshold**.
Across the full published range 0.10–0.25 the answer moves from 17.0 to 6.8 km²
— a factor of 2.5.

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

| Season | Peak km² | Peak date | Usable scenes | |
|---|---|---|---|---|
| 2017-18 | — | — | 1 | omitted, not a season |
| 2018-19 | 1.02 | 2019-01-25 | 9 | marginal, just above floor |
| 2019-20 | 8.24 | 2020-01-29 | 9 | |
| 2020-21 | 8.80 | 2021-01-24 | 9 | |
| 2021-22 | *(0.39)* | — | 9 | **no melt detected** |
| 2022-23 | 5.28 | 2023-02-22 | 5 | |
| 2023-24 | 2.02 | 2024-02-08 | 2 | likely undercounted |
| 2024-25 | *(0.29)* | — | 4 | **no melt detected** |
| 2025-26 | **11.52** | 2026-02-10 | 7 | largest in record |

Inter-annual variability is large — from no detectable melt to 11.52 km² —
and 2025-26 is the largest in the record. **Do not read a trend off this
chart.** Eight plotted seasons over one 20 km window, several undersampled,
cannot separate a trend from year-to-year weather.

The 2025-26 peak was inspected directly rather than trusted: the scene is
clean, with unmistakable dark meltwater channels and the mask sitting exactly
on them.

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
  away. The 35 nodata scenes are partial-swath and could largely be recovered
  by mosaicking same-date scenes.
- **Same-date scenes are true duplicates, not complementary swaths.** Checked
  before collapsing them: in every observed case both cover 376–419 km² of the
  ~419 km² AOI, so keeping the better one discards no coverage.

## Configuration

Everything lives in `src/melt.py`:

| Setting | Value | Note |
|---|---|---|
| Study area | tile `19CEV`, 2048×2048 px @ 10 m | ~20.5 km square, ~72.1°S 69°W |
| NDWI threshold | `0.16` | tunable; see above |
| Brightness floor | green > 3000 DN | rejects dark ocean/shadow |
| AOI cloud limit | 20% | `season.py` scene rejection |

The study area is pinned to one Sentinel-2 tile on purpose. The S2 tiling grid
is fixed, so the same pixel window on tile 19CEV is the same ground on every
date — that is what makes areas comparable across a time series.

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
- **Blue ice and dark terrain can false-positive.** This is what sets the
  0.7 km² noise floor above. The brightness floor catches the worst of it,
  but it is not a physical discriminator.
- **No validation against published George VI numbers.** The pipeline is
  internally consistent but externally unanchored, which is the largest
  outstanding gap in the science.
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

1. Mosaic same-date partial-swath scenes. 35 of 95 scenes in 2023-24 were
   discarded for off-swath nodata; merging them would roughly double usable
   coverage in the worst seasons, which is the single biggest weakness above.
2. Expand the AOI from one window to the full shelf via a geographic polygon.
   This is the main barrier to comparing against published studies.
3. Per-pixel cloud masking (`s2cloudless` on L1C) instead of whole-scene
   rejection, recovering partially-clear scenes.
4. ICESat-2 ATL06/ATL11 elevation overlay (`icepyx`) — pond depth and hence
   volume, which is what actually drives hydrofracture, rather than extent.
5. ERA5 climate drivers + ML classifier to replace NDWI thresholding, which
   also dissolves the threshold-sensitivity problem.
6. Prediction layer: seasonal weather → expected melt extent.

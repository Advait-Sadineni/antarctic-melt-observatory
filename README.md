# Antarctic Surface Meltwater Detection

Detects surface meltwater ponds on Antarctic ice shelves from free Sentinel-2
satellite imagery. Meltwater matters because ponds drive *hydrofracture*: water
pooling on an ice shelf wedges open crevasses under its own weight and can
disintegrate the shelf, as happened to Larsen B in 2002. Once a shelf is gone,
the glaciers it buttressed accelerate into the ocean.

Current scope: single-scene detection, a full melt-season pond-area time series,
and a threshold sensitivity analysis, all over one 20 km study area on **George VI
Ice Shelf** (Antarctic Peninsula). Roadmap below covers multi-year trends,
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
- **Blue ice can false-positive.** The brightness floor catches the worst of it,
  but it is not a physical discriminator.

Handled correctly: the **Sentinel-2 baseline 04.00 offset**. From 2022-01-25
onward, L2A reflectance is shifted by −1000 DN. NDWI is a ratio of differences,
so an additive offset genuinely changes the result. `melt.boa_offset()` applies
it per-scene from the processing baseline, on both the cached and fresh paths.

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
src/tune_threshold.py   threshold sensitivity sweep
data/                   cached GeoTIFF windows (gitignored)
output/                 rendered PNGs and CSVs
```

## Roadmap

1. Multi-year trend — repeat the season run across melt years for a trend chart
2. Expand the AOI from one window to the full shelf via a geographic polygon
3. ICESat-2 ATL06/ATL11 elevation overlay (`icepyx`) — pond depth, not just extent
4. ERA5 climate drivers + ML classifier to replace NDWI thresholding
5. Prediction layer: seasonal weather → expected melt extent

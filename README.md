# Antarctic Surface Meltwater Detection

Detects surface meltwater ponds on Antarctic ice shelves from free Sentinel-2
satellite imagery. Meltwater matters because ponds drive *hydrofracture*: water
pooling on an ice shelf wedges open crevasses under its own weight and can
disintegrate the shelf, as happened to Larsen B in 2002. Once a shelf is gone,
the glaciers it buttressed accelerate into the ocean.

This is **v1 — a single-scene demo**, not the full pipeline. It pulls one
cloud-free January scene over George VI Ice Shelf (Antarctic Peninsula), computes
NDWI, thresholds it into a meltwater mask, and renders a side-by-side PNG with
pond pixel count and area. The eventual goal is multi-season pond-area time
series, ICESat-2 elevation overlays, and an ERA5-driven prediction layer.

## Run it

```
pip install rasterio matplotlib pystac-client numpy
python src/demo.py
```

Output: `output/george_vi_meltwater.png`, plus pond count and km² on stdout.

First run streams ~16 MB from AWS (about a minute) and caches it to `data/`.
Later runs read the cache and finish in seconds. Delete `data/*.tif` to refetch.

## How it works

**NDWI** (Normalized Difference Water Index) = `(Green − NIR) / (Green + NIR)`,
using Sentinel-2 bands B03 and B08. Liquid water absorbs near-infrared very
strongly but still reflects green, so NDWI goes positive over ponds. Snow and
ice reflect both bands strongly and land near zero or negative. On a white ice
shelf the two separate cleanly, which is why a plain threshold works at all here.

Imagery comes from the **AWS Open Data** Sentinel-2 L2A archive as
Cloud-Optimized GeoTIFFs, found via the Earth Search STAC API. The script reads
only a 2048×2048 pixel window directly out of the remote COGs rather than
downloading the full 110 km tile.

## What's hardcoded / stubbed

Everything below is a deliberate v1 shortcut, not a design decision.

| Thing | Current value | Note |
|---|---|---|
| Scene | `S2B_19CEV_20210124_0_L2A` | One fixed scene, 2021-01-24, 0% cloud |
| Region | George VI Ice Shelf, ~72.1°S 69°W | Window hand-picked for dense ponding |
| Window | 2048×2048 px @ 10 m (~20.5 km) | Fixed offset into the tile |
| NDWI threshold | `0.16` | **Tunable.** Literature spans 0.10–0.25 |
| Brightness floor | Green > 3000 DN | Crude reject of dark ocean/shadow |
| Cloud masking | **None** | Relies on picking a 0%-cloud scene |
| Reflectance offset | 0 | Correct pre-2022 only; see below |

Two known correctness caveats:

- **No cloud mask.** Cloud shadow on ice can mimic a pond's NDWI signature. The
  SCL band ships with every L2A product and isn't used yet.
- **Reflectance offset.** Sentinel-2 processing baseline 04.00 (2022-01-25
  onward) subtracts 1000 DN from L2A reflectance. Our scene predates it so
  `BOA_OFFSET = 0` is correct, but any post-2022 scene needs −1000 applied or
  NDWI will be biased.

Pond area is `pixel count × 100 m²`. No slope or projection-distortion
correction; at 72°S in UTM 19S that error is small but nonzero.

## Accounts (not needed yet)

The demo runs anonymously — AWS Open Data requires no credentials. Register these
for the roadmap items:

| Service | For | Where credentials go |
|---|---|---|
| [Copernicus Data Space](https://dataspace.copernicus.eu/) | Full S2 catalogue, newer scenes than AWS mirrors | `.env` → `CDSE_USER` / `CDSE_PASSWORD` |
| [NASA Earthdata](https://urs.earthdata.nasa.gov/) | ICESat-2 ATL06/ATL11 via `icepyx` | `~/.netrc` (icepyx reads this) |
| [Copernicus CDS](https://cds.climate.copernicus.eu/) | ERA5 climate reanalysis | `~/.cdsapirc` |

Google Earth Engine is deliberately skipped — it would move computation off the
local machine and make the pipeline harder to show as a plain script.

## Layout

```
src/demo.py   the entire pipeline
data/         cached GeoTIFF windows (gitignored)
output/       rendered PNGs
```

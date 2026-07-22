"""
Core meltwater detection: scene access, cloud screening, NDWI thresholding.

Shared by demo.py (one scene -> PNG), season.py (one melt season -> time
series) and tune_threshold.py (threshold sensitivity sweep).
"""

import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.windows import Window

os.environ.setdefault("GDAL_DISABLE_READDIR_ON_OPEN", "EMPTY_DIR")
os.environ.setdefault("AWS_NO_SIGN_REQUEST", "YES")
os.environ.setdefault("GDAL_HTTP_MAX_RETRY", "3")
os.environ.setdefault("GDAL_HTTP_RETRY_DELAY", "2")

STAC_API = "https://earth-search.aws.element84.com/v1"
COLLECTION = "sentinel-2-l2a"

# --- Area of interest --------------------------------------------------------
# George VI Ice Shelf, Antarctic Peninsula. Pinned to a single Sentinel-2 tile
# so that every scene shares one pixel grid: the S2 tiling grid is fixed, so
# the same pixel window on tile 19CEV is the same patch of ground on every
# date. That is what makes areas comparable across a time series.
TILE = "19CEV"
WIN_ROW, WIN_COL, WIN_SIZE = 1280, 800, 2048  # 2048 px * 10 m = ~20.5 km square
PIXEL_M = 10.0
PIXEL_KM2 = (PIXEL_M**2) / 1e6

# --- Detection parameters ----------------------------------------------------
NDWI_THRESHOLD = 0.16  # TUNABLE - see tune_threshold.py
BRIGHTNESS_FLOOR = 3000  # rejects dark open ocean and deep shadow
NDSI_ICE_MIN = 0.60  # below this a pixel is not snow/ice/water - see ice_check()

# Detections at or below this are not distinguishable from artifacts.
#
# Set from evidence, not taste. Two scenes that pass every screen still report
# meltwater before melt onset: 2021-11-09 gives 0.39 km2 and 2022-11-25 gives
# 0.68 km2, both in early November. Inspecting them shows the detections
# hugging the edges of a dark patch - bare or blue ice, or shadowed terrain -
# near the AOI boundary, not ponds. 0.68 km2 is the largest confirmed
# non-melt detection, so anything under ~0.7 km2 is reported as "no melt
# detected" rather than as a measured area.
NOISE_FLOOR_KM2 = 0.7

# --- SCL (Scene Classification Layer) ----------------------------------------
# sen2cor ships a per-pixel class map with every L2A product.
#
# We use SCL ONLY to reject pixels, never to confirm water. Measured on this
# AOI: SCL labels 80% of visually-confirmed melt ponds as "snow/ice" rather
# than "water", because its water class is tuned for open ocean and deep
# lakes, not shallow meltwater lying on bright ice. Trusting SCL class 6 as
# the pond detector would throw away four fifths of the signal.
#
# The same class map is useful in the other direction as a scene-quality
# gate. Our AOI sits on an ice shelf, so SCL "open water" (class 6) should be
# near zero. Measured: clear scenes run 0.00-0.42% water, and that 0.42%
# tracks real ponds (SCL calls ~20% of pond pixels water, and peak ponds are
# ~2.1% of the AOI). A November scene reporting 7.7% water - which would
# imply ~38% pond cover - turned out to be cloud shadow across the shelf
# front, lit by blue skylight, which mimics water's NDWI signature almost
# exactly. See MAX_WATER_IN_AOI in season.py.
SCL_REJECT = {
    0,  # no data
    1,  # saturated / defective
    3,  # cloud shadow
    8,  # cloud, medium probability
    9,  # cloud, high probability
    10,  # thin cirrus
}

ROOT = Path(__file__).resolve().parent.parent
CACHE_DIR = ROOT / "data"


def _window(scale=1):
    """Read window at 10 m (scale=1) or 20 m (scale=2, for SCL)."""
    return Window(WIN_COL // scale, WIN_ROW // scale, WIN_SIZE // scale, WIN_SIZE // scale)


def search_scenes(start, end, max_cloud_tile=100):
    """All scenes on our tile between two dates, newest processing first."""
    from pystac_client import Client

    items = list(
        Client.open(STAC_API)
        .search(
            collections=[COLLECTION],
            bbox=[-70.0, -72.5, -66.0, -70.5],
            datetime=f"{start}/{end}",
            limit=100,
        )
        .items()
    )
    items = [i for i in items if TILE in i.id]
    if max_cloud_tile < 100:
        items = [i for i in items if (i.properties.get("eo:cloud_cover") or 0) <= max_cloud_tile]
    return sorted(items, key=lambda x: x.datetime)


def get_item(item_id):
    from pystac_client import Client

    return next(Client.open(STAC_API).search(collections=[COLLECTION], ids=[item_id]).items())


def boa_offset(item):
    """DN offset to ADD before computing reflectance ratios.

    Sentinel-2 processing baseline 04.00 (2022-01-25 onward) shifted L2A
    reflectance by -1000 DN. NDWI is a ratio of differences, so an additive
    offset does change the result - this is not cosmetic.
    """
    if item.properties.get("earthsearch:boa_offset_applied"):
        return 0  # the archive already corrected it
    baseline = item.properties.get("s2:processing_baseline")
    if baseline is not None:
        try:
            return -1000 if float(baseline) >= 4.0 else 0
        except (TypeError, ValueError):
            pass
    return -1000 if item.datetime >= datetime(2022, 1, 25, tzinfo=timezone.utc) else 0


def read_coarse(item, asset):
    """Read an asset over the AOI at 20 m.

    For a 20 m asset that is a native read. For a 10 m asset, asking for a
    half-size output makes GDAL pull the COG's first overview instead of the
    full-resolution pixels, so screening costs a quarter of the bytes.
    """
    with rasterio.open(item.assets[asset].href) as src:
        native10 = abs(src.res[0] - 10.0) < 1e-6
        return src.read(
            1,
            window=_window(scale=1 if native10 else 2),
            out_shape=(WIN_SIZE // 2, WIN_SIZE // 2),
            resampling=Resampling.average if native10 else Resampling.nearest,
        ).astype("f4")


def ice_check(item):
    """Fraction of the AOI that does not look like snow, ice, or water.

    NDSI = (green - SWIR) / (green + SWIR), using B03 and B11 (1610 nm).
    Snow, ice and liquid water all absorb SWIR strongly, so they sit high
    (0.85-0.95 here). Water cloud *reflects* SWIR, so it drops.

    This exists because SCL badly under-reports cloud over ice - white cloud
    on a white shelf is genuinely hard for it. Three February scenes that SCL
    called 2-19% cloudy were, on inspection, blanketed edge to edge, and two
    of them produced multi-km2 phantom ponds along the cloud edges.

    Measured over this AOI: clean scenes put 0.0-0.1% of pixels below NDSI
    0.6, while cloud-covered ones put 3.7-35.9% there. Used as a scene-level
    gate, not a pixel mask - the phantom pond pixels themselves still score
    high NDSI, so masking per-pixel would not remove them.
    """
    offset = boa_offset(item)
    green = read_coarse(item, "green")
    swir = read_coarse(item, "swir16")
    valid = (green > 0) & (swir > 0)
    green = np.where(valid, green + offset, 0)
    swir = np.where(valid, swir + offset, 0)

    ndsi = np.where(valid, (green - swir) / np.maximum(green + swir, 1e-6), np.nan)
    with np.errstate(invalid="ignore"):
        return float(np.nanmean(ndsi < NDSI_ICE_MIN))


def read_scl(item):
    """SCL is 20 m; nearest-neighbour it up to the 10 m grid."""
    with rasterio.open(item.assets["scl"].href) as src:
        return src.read(
            1,
            window=_window(scale=2),
            out_shape=(WIN_SIZE, WIN_SIZE),
            resampling=Resampling.nearest,
        )


def read_band(item, name):
    with rasterio.open(item.assets[name].href) as src:
        return src.read(1, window=_window()).astype("f4")


def screen(scl):
    """Cloud/shadow statistics inside the AOI.

    Deliberately computed from SCL over our window rather than from the
    scene's `eo:cloud_cover` property: that property describes the whole
    110 km tile. Measured here, a scene reporting 37% tile cloud was 97%
    clear over this AOI. Filtering on tile metadata discards good data.
    """
    rejected = np.isin(scl, list(SCL_REJECT))
    nodata = scl == 0
    return {
        "cloud_frac": float((rejected & ~nodata).mean()),
        "nodata_frac": float(nodata.mean()),
        "usable_frac": float((~rejected).mean()),
        "water_frac": float((scl == 6).mean()),
        "snow_frac": float((scl == 11).mean()),
        "reject_mask": rejected,
    }


def detect(green, nir, reject_mask=None, threshold=NDWI_THRESHOLD):
    """NDWI threshold -> meltwater mask.

    NDWI = (green - NIR) / (green + NIR). Liquid water absorbs near-infrared
    strongly but still reflects green, so ponds go positive; snow and ice
    reflect both bands and sit near or below zero.
    """
    valid = (green > 0) & (nir > 0)
    ndwi = np.where(valid, (green - nir) / np.maximum(green + nir, 1e-6), np.nan)

    ponds = (ndwi > threshold) & valid & (green > BRIGHTNESS_FLOOR)
    if reject_mask is not None:
        ponds &= ~reject_mask

    return ndwi, ponds, valid


def pond_stats(ponds, valid):
    n = int(ponds.sum())
    return {
        "pond_px": n,
        "pond_km2": n * PIXEL_KM2,
        "valid_km2": int(valid.sum()) * PIXEL_KM2,
        "pond_pct": 100.0 * n / max(int(valid.sum()), 1),
    }


def load_scene(item, bands=("red", "green", "blue", "nir"), use_cache=True):
    """Windowed multi-band read, cached to data/ as a GeoTIFF.

    Reads only our 2048x2048 window out of the remote Cloud-Optimized
    GeoTIFFs rather than pulling the whole 110 km tile.
    """
    tag = "-".join(bands)
    cache = CACHE_DIR / f"{item.id}_{WIN_ROW}_{WIN_COL}_{WIN_SIZE}_{tag}.tif"

    if use_cache and cache.exists():
        with rasterio.open(cache) as src:
            arrays = [src.read(i + 1) for i in range(len(bands))]
    else:
        arrays, profile = [], None
        for name in bands:
            with rasterio.open(item.assets[name].href) as src:
                arrays.append(src.read(1, window=_window()))
                if profile is None:
                    profile = src.profile.copy()
                    profile.update(
                        count=len(bands),
                        height=WIN_SIZE,
                        width=WIN_SIZE,
                        transform=src.window_transform(_window()),
                        compress="deflate",
                    )
        CACHE_DIR.mkdir(parents=True, exist_ok=True)
        with rasterio.open(cache, "w", **profile) as dst:
            for i, arr in enumerate(arrays, start=1):
                dst.write(arr, i)
            dst.descriptions = tuple(bands)

    # The cache holds raw archive DN. The baseline offset is applied here, on
    # both paths, so a cached rerun cannot disagree with a fresh fetch.
    offset = boa_offset(item)
    return {
        b: np.where(a > 0, a.astype("f4") + offset, 0.0) for b, a in zip(bands, arrays)
    }

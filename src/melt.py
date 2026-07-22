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
from scipy import ndimage

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

# The study area is chosen for coverage and representativeness, NOT for pond
# density. An earlier 20 km window was picked by searching for the densest
# ponding, which made every absolute number a best case. This one is the
# largest window that is fully imaged and 100% snow/ice across four clear
# scenes drawn from different seasons (2019, 2020, 2021, 2026), so it
# includes plenty of shelf that never ponds.
WIN_ROW, WIN_COL, WIN_SIZE = 1180, 920, 6144  # 6144 px * 10 m = ~61 km square
PIXEL_M = 10.0
PIXEL_KM2 = (PIXEL_M**2) / 1e6

# Screening only needs fractions, not detail, so it runs on heavily decimated
# reads. At 60 m a 61 km AOI is ~1024 px square, which costs a fraction of a
# percent of the full-resolution bytes.
SCREEN_PIXEL_M = 60.0

# --- Detection parameters ----------------------------------------------------
NDWI_THRESHOLD = 0.16  # TUNABLE - see tune_threshold.py
BRIGHTNESS_FLOOR = 3000  # rejects dark open ocean and deep shadow
NDSI_ICE_MIN = 0.60  # below this a pixel is not snow/ice/water - see cloud_mask()

# Cloud contaminates the ground around it, not just under it. Detected cloud
# is dilated by this distance before being excluded. Measured on a confirmed
# cloudy scene, a 1 km halo removes 72% of the phantom ponds while costing a
# confirmed clean scene 0.5% of its real detections.
CLOUD_HALO_KM = 1.0

# Below this solar elevation the surface is lit at a grazing angle, shadows
# from crevasses and surface topography lengthen, and shadowed snow lit by
# blue skylight mimics water's NDWI signature. Measured: scenes at 9.5 deg
# and 11.7 deg return the largest "melt" of their season, in March, with the
# detections tracking the crevassed eastern margin rather than any pond.
MIN_SUN_ELEVATION_DEG = 15.0

# Detections at or below this are not distinguishable from a persistent
# background and are reported as "no melt detected" rather than as an area.
#
# Set from evidence, and re-derived whenever the study area changes, because
# it does not scale with area in any predictable way. The calibration is the
# pre-melt baseline: early-November scenes, before any melt onset, that pass
# every screen and still report meltwater.
#
# On the 20 km window the baseline was 0.39-0.68 km2, and inspection showed
# the detections hugging the edges of a dark patch - bare or blue ice, or
# shadowed terrain - rather than ponds.
#
# On the current 61 km study area the November baseline is 1.21-1.75 km2 and
# strikingly stable, which says it is a persistent surface feature rather
# than random noise. The floor sits just above the observed maximum.
NOV_BASELINE_MAX_KM2 = 1.95  # largest pre-melt reading observed on this AOI
NOISE_FLOOR_KM2 = 2.2

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


def _window(scale=1, row=None, col=None, size=None):
    """Read window in the source's own pixel grid.

    scale is the source pixel size divided by 10 m: 1 for a 10 m band, 2 for
    a 20 m band such as SCL or SWIR.
    """
    row = WIN_ROW if row is None else row
    col = WIN_COL if col is None else col
    size = WIN_SIZE if size is None else size
    return Window(col // scale, row // scale, size // scale, size // scale)


def _screen_shape(size=None):
    size = WIN_SIZE if size is None else size
    return int(size * PIXEL_M / SCREEN_PIXEL_M)


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
    """Read an asset over the AOI at SCREEN_PIXEL_M.

    Asking for a much smaller output makes GDAL pull a COG overview rather
    than full-resolution pixels, so screening reads a tiny fraction of the
    bytes a detection read would.
    """
    with rasterio.open(item.assets[asset].href) as src:
        native10 = abs(src.res[0] - 10.0) < 1e-6
        n = _screen_shape()
        return src.read(
            1,
            window=_window(scale=1 if native10 else 2),
            out_shape=(n, n),
            resampling=Resampling.average,
        ).astype("f4")


def cloud_mask(item, halo_km=CLOUD_HALO_KM):
    """Locate cloud over ice, and mask the ground it contaminates.

    NDSI = (green - SWIR) / (green + SWIR), using B03 and B11 (1610 nm).
    Snow, ice and liquid water all absorb SWIR strongly, so they sit high
    (0.85-0.95 here). Water cloud *reflects* SWIR, so it drops.

    This exists because SCL badly under-reports cloud over ice - white cloud
    on a white shelf is genuinely hard for it. Three February scenes that SCL
    called 2-19% cloudy were, on inspection, blanketed edge to edge, and two
    of them produced multi-km2 phantom ponds along the cloud edges.

    The phantom pixels themselves score *high* NDSI, so thresholding NDSI
    per-pixel does not remove them - they sit beside the cloud rather than
    under it. Dilating the low-NDSI cells by CLOUD_HALO_KM does remove them,
    and unlike a whole-scene NDSI gate it does not depend on the size of the
    study area. That mattered: a cloud bank that filled the old 20 km window
    covers only ~11% of the current 61 km one, so the old scene-level
    threshold silently stopped firing when the AOI grew.

    Returns the full-resolution contaminated mask plus coarse diagnostics.
    """
    offset = boa_offset(item)
    green = read_coarse(item, "green")
    swir = read_coarse(item, "swir16")
    valid = (green > 0) & (swir > 0)
    g = np.where(valid, green + offset, 0)
    s = np.where(valid, swir + offset, 0)

    ndsi = np.where(valid, (g - s) / np.maximum(g + s, 1e-6), np.nan)
    low = np.nan_to_num(ndsi, nan=1.0) < NDSI_ICE_MIN

    cells = int(round(halo_km * 1000 / SCREEN_PIXEL_M))
    halo = ndimage.binary_dilation(low, iterations=cells) if cells else low

    scale = WIN_SIZE // low.shape[0]
    return {
        "nonice_frac": float(low.mean()),
        "halo_frac": float(halo.mean()),
        "mask": np.kron(halo, np.ones((scale, scale), bool)),
    }


def sun_elevation(item):
    return float(item.properties.get("view:sun_elevation") or 90.0)


def read_scl(item, full_res=False):
    """SCL is 20 m. Screening wants it decimated; masking wants the 10 m grid."""
    with rasterio.open(item.assets["scl"].href) as src:
        n = WIN_SIZE if full_res else _screen_shape()
        return src.read(
            1,
            window=_window(scale=2),
            out_shape=(n, n),
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


def reject_mask(item):
    """Full-resolution cloud/shadow mask, for scenes that passed screening."""
    return np.isin(read_scl(item, full_res=True), list(SCL_REJECT))


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


def load_scene(item, bands=("red", "green", "blue", "nir"), use_cache=False,
               row=None, col=None, size=None):
    """Windowed multi-band read, optionally cached to data/ as a GeoTIFF.

    Reads only the requested window out of the remote Cloud-Optimized GeoTIFFs
    rather than pulling the whole 110 km tile.

    Caching defaults off. At the 61 km study area one scene is ~150 MB, so
    caching a multi-year run would run to tens of gigabytes for no benefit -
    season results are already cached as CSV, which makes reruns free. It is
    worth enabling for the small detail window that demo.py re-reads.
    """
    row = WIN_ROW if row is None else row
    col = WIN_COL if col is None else col
    size = WIN_SIZE if size is None else size
    win = _window(row=row, col=col, size=size)

    tag = "-".join(bands)
    cache = CACHE_DIR / f"{item.id}_{row}_{col}_{size}_{tag}.tif"

    if use_cache and cache.exists():
        with rasterio.open(cache) as src:
            arrays = [src.read(i + 1) for i in range(len(bands))]
    else:
        arrays, profile = [], None
        for name in bands:
            with rasterio.open(item.assets[name].href) as src:
                arrays.append(src.read(1, window=win))
                if profile is None:
                    profile = src.profile.copy()
                    profile.update(
                        count=len(bands),
                        height=size,
                        width=size,
                        transform=src.window_transform(win),
                        compress="deflate",
                    )
        if not use_cache:
            offset = boa_offset(item)
            return {b: np.where(a > 0, a.astype("f4") + offset, 0.0)
                    for b, a in zip(bands, arrays)}
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

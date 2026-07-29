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
#
# This is the published Moussavi et al. 2020 method for Antarctic supraglacial
# lakes (Remote Sensing 12:134): NDWI (green-NIR) > 0.19 for Sentinel-2, plus a
# green-minus-red shadow test. It replaced a hand-tuned 0.16 threshold with no
# shadow test, which validation showed was too permissive - every false
# positive was crevasse shadow sitting at NDWI 0.17-0.25, just above 0.16.
#
# The choice is not taste. On 220 blind-labelled points across three scenes,
# this configuration scored the best F1 (0.50 vs 0.46) of nine candidates and
# lifted precision from 0.55 to 0.71, while the literature independently uses
# exactly this threshold - the same papers whose 19 Jan 2020 peak date this
# pipeline reproduced. See retune.py for the full comparison.
NDWI_THRESHOLD = 0.19
BRIGHTNESS_FLOOR = 3000  # rejects dark open ocean and deep shadow

# Moussavi shadow / crevasse test: real meltwater keeps green well above red,
# whereas shaded snow and crevasse shadow - which dim all visible bands roughly
# together - do not. Expressed in reflectance; the imagery is DN, so it is
# scaled at the point of use. (green - red) is invariant to the baseline offset
# because both bands carry it, so no offset handling is needed here.
SHADOW_GREEN_MINUS_RED = 0.09
DN_PER_REFLECTANCE = 10000.0

# Hysteresis. A single threshold detects pond cores well but under-draws their
# margins: the fringe pixels are water/ice mixtures whose NDWI is diluted below
# 0.19. Blind reference points confirmed every missed water point was a margin
# adjacent to a detection - never a missed lake and never far-field. So cores
# are found at 0.19 and then grown into connected pixels above this lower
# threshold, which recovers margins without flooding crevasse fields, because
# those are not connected to a genuine core. The shadow test is kept on the
# grown pixels too - dropping it there re-admits crevasse shadow and collapses
# precision. Measured over 220 labelled points, this lifts area recall from
# 0.39 to 0.45 and F1 from 0.50 to 0.54 at a 0.03 cost in precision.
GROW_THRESHOLD = 0.14

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
# With the shadow test and hysteresis, the pre-melt background collapsed: the
# clean mid-November readings now cluster at 0.66-1.23 km2, down from 1.2-1.95
# before, because the old baseline was mostly crevasse-shadow false positives
# that the shadow test now removes. One early-November scene (2022-11-01) reads
# 3.2 km2, but inspection shows genuine dark margin melt, not an artifact, so
# it is treated as real early melt rather than used to set the floor. The floor
# sits just above the clean background.
NOV_BASELINE_MAX_KM2 = 1.23  # largest clean pre-melt reading (excl. real early melt)
NOISE_FLOOR_KM2 = 1.3

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


# Sentinel-2 tile 19CEV sits in UTM zone 19S. The S2 tiling grid is fixed, so
# the tile's georeferencing is a constant of the tile rather than of any
# particular scene - which is what lets a pixel window mean the same ground on
# every date, and lets another sensor be warped onto exactly this grid.
#
# Read from a real asset rather than hardcoded: an origin wrong by a couple of
# pixels would misalign every cross-sensor comparison silently, and the error
# would look like genuine disagreement between satellites.
TILE_CRS = "EPSG:32719"  # default tile 19CEV; other tiles carry their own
_GEOREF_CACHE = {}


_SOURCE = None


def get_source():
    """Active SceneSource. Defaults to Earth Search, constructed lazily and
    offline (tests never touch the network just by importing melt)."""
    global _SOURCE
    if _SOURCE is None:
        from core.sources import EarthSearchSource
        _SOURCE = EarthSearchSource()
    return _SOURCE


def set_source(source):
    """Swap the scene provider. Clears per-tile georefs: they were read from
    the old provider's files and must not leak across sources."""
    global _SOURCE
    _SOURCE = source
    _GEOREF_CACHE.clear()


def _tile_of(item):
    return get_source().tile_of(item)


def tile_georeference(item):
    """(crs, transform) of a tile's 10 m grid, cached per tile.

    Keyed by MGRS tile, not a single slot: the shelf spans several tiles across
    three UTM zones, and handing them all the first tile's transform would
    misplace every one silently.
    """
    key = _tile_of(item)
    if key not in _GEOREF_CACHE:
        with rasterio.open(get_source().band_href(item, "green")) as src:
            if abs(src.res[0] - PIXEL_M) > 1e-6:
                raise ValueError(f"expected a {PIXEL_M} m asset, got {src.res[0]} m")
            _GEOREF_CACHE[key] = (src.crs, src.transform)
    return _GEOREF_CACHE[key]


def aoi_transform(item):
    """Affine transform of the study area, in the tile's own CRS."""
    _, tile = tile_georeference(item)
    return tile * tile.identity().translation(WIN_COL, WIN_ROW)


# --- multi-tile support ------------------------------------------------------
# The single-tile pipeline is driven by the module globals above. For full-shelf
# processing, set_aoi() repoints them at another tile/window plus an optional
# shelf mask (which pixels are ice shelf vs ocean/mountain). It is used one tile
# at a time, so the existing per-scene threading inside a tile stays safe; every
# function keeps its signature, so single-tile callers and the tests are
# unaffected (the defaults are exactly the original 19CEV window).
AOI_SHELF_MASK = None  # full-res bool, True = ice shelf; None = whole window


def set_aoi(tile, win_row, win_col, win_size, shelf_mask=None):
    """Point the pipeline at a tile/window (and optional shelf mask)."""
    global TILE, WIN_ROW, WIN_COL, WIN_SIZE, AOI_SHELF_MASK
    TILE, WIN_ROW, WIN_COL, WIN_SIZE = tile, win_row, win_col, win_size
    AOI_SHELF_MASK = shelf_mask


def search_scenes(start, end, max_cloud_tile=100):
    """All scenes on our tile between two dates, newest processing first."""
    items = get_source().search(
        bbox=[-70.0, -72.5, -66.0, -70.5],
        datetime=f"{start}/{end}",
        limit=100,
    )
    items = [i for i in items if get_source().tile_of(i) == TILE]
    if max_cloud_tile < 100:
        items = [i for i in items if (i.properties.get("eo:cloud_cover") or 0) <= max_cloud_tile]
    return sorted(items, key=lambda x: x.datetime)


def get_item(item_id):
    return get_source().search(ids=[item_id])[0]


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
    with rasterio.open(get_source().band_href(item, asset)) as src:
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

    return {
        "nonice_frac": float(low.mean()),
        "halo_frac": float(halo.mean()),
        "mask": _upsample_to_window(halo),
    }


def _upsample_to_window(coarse):
    """Nearest-neighbour upsample a coarse mask to exactly (WIN_SIZE, WIN_SIZE).

    ceil-then-crop, so it is correct for any window size, not only those that
    are a clean multiple of the 60 m screen grid (they are not, once the AOI
    is a shelf tile of arbitrary extent).
    """
    scale = int(np.ceil(WIN_SIZE / coarse.shape[0]))
    return np.kron(coarse, np.ones((scale, scale), bool))[:WIN_SIZE, :WIN_SIZE]


def sun_elevation(item):
    return float(item.properties.get("view:sun_elevation") or 90.0)


def read_scl(item, full_res=False):
    """SCL is 20 m. Screening wants it decimated; masking wants the 10 m grid."""
    with rasterio.open(get_source().band_href(item, "scl")) as src:
        n = WIN_SIZE if full_res else _screen_shape()
        return src.read(
            1,
            window=_window(scale=2),
            out_shape=(n, n),
            resampling=Resampling.nearest,
        )


def read_band(item, name):
    with rasterio.open(get_source().band_href(item, name)) as src:
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
    """Full-resolution reject mask: cloud/shadow, plus non-shelf if an AOI shelf
    mask is set. Everything downstream (detection, usable fraction) then counts
    only ice-shelf ground."""
    rej = np.isin(read_scl(item, full_res=True), list(SCL_REJECT))
    if AOI_SHELF_MASK is not None:
        rej = rej | ~AOI_SHELF_MASK
    return rej


def detect(green, nir, reject_mask=None, threshold=NDWI_THRESHOLD, red=None,
           grow_threshold=GROW_THRESHOLD):
    """NDWI threshold, a shadow test, and hysteresis -> meltwater mask.

    NDWI = (green - NIR) / (green + NIR). Liquid water absorbs near-infrared
    strongly but still reflects green, so ponds go positive; snow and ice
    reflect both bands and sit near or below zero.

    When the red band is supplied the full published method runs: the Moussavi
    shadow test (green - red > 0.09 reflectance) rejects crevasse and
    topographic shadow, and hysteresis grows pond cores into connected margin
    pixels above ``grow_threshold``. Red is optional only so a pure-NDWI
    sensitivity sweep can isolate the threshold; real detection always passes
    it. Pass grow_threshold=None to get cores without margin growth.
    """
    valid = (green > 0) & (nir > 0)
    ndwi = np.where(valid, (green - nir) / np.maximum(green + nir, 1e-6), np.nan)
    ndwi0 = np.nan_to_num(ndwi, nan=-9)

    ok = valid & (green > BRIGHTNESS_FLOOR)
    if red is not None:
        ok = ok & ((green - red) > SHADOW_GREEN_MINUS_RED * DN_PER_REFLECTANCE)
    if reject_mask is not None:
        ok = ok & ~reject_mask

    core = ok & (ndwi0 > threshold)
    if red is not None and grow_threshold is not None:
        grow = ok & (ndwi0 > grow_threshold)
        labels, n = ndimage.label(grow)
        keep = np.zeros(n + 1, bool)
        keep[np.unique(labels[core])] = True
        keep[0] = False  # background label never counts
        ponds = keep[labels]
    else:
        ponds = core

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
            with rasterio.open(get_source().band_href(item, name)) as src:
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

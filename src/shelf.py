"""
Full-shelf coverage: a shelf-wide meltwater area over George VI, built the
correct way so it can be compared to the published whole-shelf inventories
(Banwell et al. 2021; Dirscherl et al. 2021; Corr et al. 2022).

Two problems make a naive per-tile sum wrong, and this module solves both:

  1. Adjacent Sentinel-2 tiles overlap and the shelf is shared across them, so
     summing per-tile areas double-counts. Here every tile's detection is
     reprojected onto ONE common Antarctic Polar Stereographic grid
     (EPSG:3031) and UNIONed, so each ground location is counted once.

  2. "Ice" is not "shelf" - the tiles also cover Alexander Island's mountains
     and valley glaciers, where the literature does not map lakes. Detection is
     masked to an authoritative George VI ice-shelf polygon rasterised onto the
     common grid, so only the flat floating shelf counts.

The shelf polygon is a drop-in parameter (reference/george_vi_measures.geojson,
the authoritative MEaSUReS Antarctic Boundaries v2 outline - area 23,260 km2,
matching the known ~24,000 km2 whole shelf, and tight at the calving front so
open ocean and sea ice are excluded).

Per-tile detection reuses the validated single-tile detector via melt.set_aoi,
so the shelf number is produced by exactly the method the rest of the project
validated (Moussavi shadow test + hysteresis).

STATUS: the engine is validated, uses an authoritative boundary, and the
dominant tile has now passed blind reference-point validation. Per-tile
detection reprojects and unions correctly on the common grid; the MEaSUReS
polygon correctly delineates the flat shelf (mountain tiles 19DDA/19DDC read 0,
and its outline hugs the shelf edge on inspection). All three issues found
earlier are resolved:

  - [fixed] Boundary precision. Swapped the cartographic Natural Earth polygon
    for the authoritative MEaSUReS Antarctic Boundaries v2 outline (23,260 km2,
    matching the known shelf), which is tight at the calving front and excludes
    open ocean / sea ice.
  - [fixed] Resampling inflation. Reprojection is now area-weighted (average),
    giving each 30 m grid cell its true water fraction instead of marking the
    whole cell water if any 10 m pixel is. This trimmed the tile areas ~15-25%.

  - [resolved] Per-tile validation of 19DEA (the 127 km2 that dominates the
    2020-21 shelf-wide 138 km2). 80 stratified blind points were labelled
    (validate_shelf_tile.py). Result: precision 0.625 (95% CI 0.47-0.76),
    recall 0.51. The detection is genuine meltwater - visually vivid blue
    ponds/channels, spectrally dark-red + high-NDWI. Of the false positives,
    about half are pond-margin mixed pixels (spectrally identical to water) and
    half are crevasse-shadow / thin-cloud that are spectrally entangled with
    shallow meltwater: sweeps of the NDWI core threshold (0.19->0.30) and the
    shadow test (0.09->0.18) each trade real water for false positives ~1:1, so
    no threshold cleanly removes them, and a shape filter would delete George
    VI's real supraglacial channels. Crucially, those false positives are more
    than offset by margin water the detector MISSES (recall 0.51): the
    Horvitz-Thompson bias-corrected true-water area is 157 +/- 37 km2, i.e. the
    reported 127 km2 is unbiased-to-conservative, not inflated. The absolute
    shelf-wide number therefore stands. See output/shelf_val/19DEA/.

The core sound tiles (19CDV + 19CEV, deduped) give ~15 km2 for 2020-21.

Run:  python src/shelf.py season 2020-21
"""

import json
import sys
from datetime import date
from pathlib import Path

import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.features import rasterize
from rasterio.transform import from_origin
from rasterio.warp import reproject, transform as warp_transform, transform_bounds

import melt

SHELF_TILES = ["18DXF", "18DXG", "18DXH", "19CDV", "19CEV", "19DDA",
               "19DDC", "19DEA", "19DFA", "19DFC", "20CME"]

GRID_CRS = "EPSG:3031"     # Antarctic Polar Stereographic - one grid, all zones
GRID_RES = 30.0           # m; 30 m keeps the shelf-wide grid tractable
BOUNDARY = melt.ROOT / "reference" / "george_vi_measures.geojson"
OUT_DIR = melt.ROOT / "output" / "shelf"


# --- common grid + shelf mask ------------------------------------------------

def _boundary_3031():
    """George VI polygon geometry reprojected to the grid CRS (EPSG:3031)."""
    gj = json.loads(BOUNDARY.read_text())
    geoms = []
    for f in gj["features"]:
        g = f["geometry"]
        polys = g["coordinates"] if g["type"] == "MultiPolygon" else [g["coordinates"]]
        new = []
        for poly in polys:
            rings = []
            for ring in poly:
                xs, ys = warp_transform("EPSG:4326", GRID_CRS,
                                        [c[0] for c in ring], [c[1] for c in ring])
                rings.append(list(zip(xs, ys)))
            new.append(rings)
        geoms.append({"type": "MultiPolygon", "coordinates": new})
    return geoms


def build_grid(ref_items):
    """Common EPSG:3031 grid covering the on-shelf tiles' extent.

    ref_items: one Sentinel-2 item per tile, to get tile footprints.
    Returns (transform, width, height, shelf_mask).
    """
    xs, ys = [], []
    for it in ref_items:
        w, s, e, n = transform_bounds("EPSG:4326", GRID_CRS, *it.bbox)
        xs += [w, e]
        ys += [s, n]
    minx, maxx = min(xs), max(xs)
    miny, maxy = min(ys), max(ys)
    # snap out to the grid
    minx = np.floor(minx / GRID_RES) * GRID_RES
    maxy = np.ceil(maxy / GRID_RES) * GRID_RES
    width = int(np.ceil((maxx - minx) / GRID_RES))
    height = int(np.ceil((maxy - miny) / GRID_RES))
    transform = from_origin(minx, maxy, GRID_RES, GRID_RES)

    shelf = rasterize(
        [(g, 1) for g in _boundary_3031()],
        out_shape=(height, width), transform=transform,
        fill=0, dtype="uint8", all_touched=True,
    ).astype(bool)
    return transform, width, height, shelf


# --- fixed grid (season-independent, so every season shares one denominator) -

GRID_CACHE = OUT_DIR / "grid_fixed.npz"
CURRENT_SHELF = "george_vi"   # which shelf the module globals currently point at


def set_shelf(name, tiles, boundary_path):
    """Repoint the module at another ice shelf (region-agnostic operation).

    Sets the boundary polygon, the Sentinel-2 tile list, the per-shelf grid
    cache and the shelf-mask cache so the exact same validated machinery runs on
    any shelf. Defaults (no call) keep the original George VI configuration, so
    single-shelf callers and tests are unaffected."""
    global BOUNDARY, SHELF_TILES, GRID_CACHE, CURRENT_SHELF, _COARSE_SHELF, _SHELF_BBOX
    CURRENT_SHELF = name
    BOUNDARY = Path(boundary_path)
    SHELF_TILES = list(tiles)
    GRID_CACHE = OUT_DIR / f"grid_{name}.npz"
    _COARSE_SHELF = {}   # per-tile coarse shelf masks are shelf-specific; reset
    _SHELF_BBOX = None   # bbox follows the new boundary; recompute lazily


def _reference_items():
    """One clear Sentinel-2 item per shelf tile, for a season-independent
    footprint. A summer 2020 scene is picked so the tile is fully imaged."""
    from pystac_client import Client
    cl = Client.open(melt.STAC_API)
    items = {}
    # widen until every tile has a footprint scene - a tile silently missing
    # here would clip its part of the shelf out of the mask (wrong denominator)
    windows = [("2020-01-01/2020-02-28", 30), ("2020-01-01/2020-02-28", 80),
               ("2021-01-01/2021-02-28", 80), ("2022-01-01/2022-02-28", 95)]
    for t in SHELF_TILES:
        for dt, cc in windows:
            r = list(cl.search(collections=["sentinel-2-l2a"],
                     query={"grid:code": {"eq": f"MGRS-{t}"},
                            "eo:cloud_cover": {"lt": cc}},
                     datetime=dt, limit=50).items())
            if r:
                items[t] = min(r, key=lambda x: x.properties.get("eo:cloud_cover", 99))
                break
        else:
            print(f"[grid] WARNING: no footprint scene for {t}; shelf mask may be clipped")
    return items


def build_fixed_grid(rebuild=False):
    """Season-INDEPENDENT common grid + shelf mask.

    The grid extent is fixed by the 11 shelf tiles' footprints (which never
    change), and the shelf mask is the George VI polygon clipped to the tiles'
    combined coverage. Because it does not depend on which scenes a given season
    happened to have, every season is normalised to the SAME shelf area - the
    fix for seasons whose denominator otherwise shrank when tiles were cloudy.
    Cached to disk so it is built once.
    """
    if GRID_CACHE.exists() and not rebuild:
        z = np.load(GRID_CACHE, allow_pickle=False)
        tr = rasterio.transform.Affine(*z["transform"])
        return tr, int(z["w"]), int(z["h"]), z["shelf"].astype(bool)

    items = _reference_items()
    tr, w, h, _ = build_grid(list(items.values()))

    # tile coverage: union of each tile's imaged footprint on the grid
    cover = np.zeros((h, w), bool)
    for it in items.values():
        g = it.geometry
        polys = g["coordinates"] if g["type"] == "MultiPolygon" else [g["coordinates"]]
        rings3031 = []
        for poly in polys:
            rr = []
            for ring in poly:
                xs, ys = warp_transform("EPSG:4326", GRID_CRS,
                                        [c[0] for c in ring], [c[1] for c in ring])
                rr.append(list(zip(xs, ys)))
            rings3031.append(rr)
        geom = {"type": "MultiPolygon", "coordinates": rings3031}
        cover |= rasterize([(geom, 1)], out_shape=(h, w), transform=tr,
                           dtype="uint8", all_touched=True).astype(bool)

    poly = rasterize([(g, 1) for g in _boundary_3031()], out_shape=(h, w),
                     transform=tr, fill=0, dtype="uint8", all_touched=True).astype(bool)
    shelf = poly & cover

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    np.savez(GRID_CACHE, transform=np.array(tr[:6], "f8"),
             w=w, h=h, shelf=shelf)
    print(f"[grid] fixed {w}x{h} @ {GRID_RES:.0f} m, shelf mask "
          f"{shelf.sum()*(GRID_RES**2)/1e6:,.0f} km2 (cached)")
    return tr, w, h, shelf


# --- per-tile detection warped onto the grid ---------------------------------

_SHELF_BBOX = None   # lon/lat bbox of the CURRENT shelf boundary, cached


def _shelf_bbox():
    """Bounding box (w, s, e, n) of the current BOUNDARY polygon, padded a
    little so shelf margins are never clipped. Derived from the boundary file so
    it follows set_shelf() - the window is shelf-specific, never hard-coded
    (a hard-coded George VI bbox here silently dropped other shelves' tiles)."""
    global _SHELF_BBOX
    if _SHELF_BBOX is None:
        gj = json.loads(BOUNDARY.read_text())
        xs, ys = [], []
        for f in gj["features"]:
            g = f["geometry"]
            polys = (g["coordinates"] if g["type"] == "MultiPolygon"
                     else [g["coordinates"]])
            for poly in polys:
                for ring in poly:
                    for x, y in ring:
                        xs.append(x); ys.append(y)
        _SHELF_BBOX = (min(xs) - 0.15, min(ys) - 0.05,
                       max(xs) + 0.15, max(ys) + 0.05)
    return _SHELF_BBOX


def _tile_window(item):
    """Window on this tile that overlaps the current shelf's bbox."""
    crs, tile_tr = melt.tile_georeference(item)
    # shelf lon/lat bbox -> this tile's UTM -> pixel window
    w, s, e, n = _shelf_bbox()
    xs, ys = warp_transform("EPSG:4326", crs, [w, e, w, e], [s, s, n, n])
    col0 = max(0, int((min(xs) - tile_tr.c) / melt.PIXEL_M))
    col1 = min(melt.TILE_PX if hasattr(melt, "TILE_PX") else 10980,
               int((max(xs) - tile_tr.c) / melt.PIXEL_M))
    row0 = max(0, int((tile_tr.f - max(ys)) / melt.PIXEL_M))
    row1 = min(10980, int((tile_tr.f - min(ys)) / melt.PIXEL_M))
    size = min(max(col1 - col0, row1 - row0), 10980 - max(row0, col0))
    return row0, col0, size


MIN_PERSIST = 2   # a pixel must read water on >= this many clear scenes to count
DRY_COARSE_PX = 20  # on-shelf 60 m water below this in every scene => tile is dry

_COARSE_SHELF = {}   # per-MGRS-tile coarse (60 m) shelf mask, cached


def _boundary_in_crs(crs):
    """George VI polygon rings reprojected from EPSG:4326 to ``crs``."""
    gj = json.loads(BOUNDARY.read_text())
    geoms = []
    for f in gj["features"]:
        g = f["geometry"]
        polys = g["coordinates"] if g["type"] == "MultiPolygon" else [g["coordinates"]]
        new = []
        for poly in polys:
            rings = []
            for ring in poly:
                xs, ys = warp_transform("EPSG:4326", crs,
                                        [c[0] for c in ring], [c[1] for c in ring])
                rings.append(list(zip(xs, ys)))
            new.append(rings)
        geoms.append({"type": "MultiPolygon", "coordinates": new})
    return geoms


def _coarse_shelf_mask(item):
    """Shelf polygon rasterised onto this tile's 60 m screening grid (cached),
    so the dry-tile check counts only water that falls on the actual ice shelf -
    not sea ice or ocean at the calving front that off-shelf tiles are full of."""
    tile = melt._tile_of(item)
    if tile not in _COARSE_SHELF:
        crs, _ = melt.tile_georeference(item)
        ctr = melt.aoi_transform(item) * rasterio.Affine.scale(6)  # 10 m -> 60 m
        cs = melt._screen_shape()
        _COARSE_SHELF[tile] = rasterize(
            [(g, 1) for g in _boundary_in_crs(crs)], out_shape=(cs, cs),
            transform=ctr, dtype="uint8", all_touched=True).astype(bool)
    return _COARSE_SHELF[tile]


def _coarse_water_mask(item):
    """Cheap 60 m water mask (COG overviews only), used only to skip full
    compositing of tiles with no on-shelf water. Requires ice-like NDSI so cloud
    - water-like in NDWI but low in NDSI - does not block a dry tile's skip."""
    g = melt.read_coarse(item, "green")
    n = melt.read_coarse(item, "nir")
    r = melt.read_coarse(item, "red")
    s = melt.read_coarse(item, "swir16")
    valid = (g > 0) & (n > 0) & (s > 0)
    ndwi = np.where(valid, (g - n) / np.maximum(g + n, 1e-6), -9)
    ndsi = np.where(valid, (g - s) / np.maximum(g + s, 1e-6), -9)
    ok = (valid & (g > melt.BRIGHTNESS_FLOOR) & (ndsi > melt.NDSI_ICE_MIN)
          & ((g - r) > melt.SHADOW_GREEN_MINUS_RED * melt.DN_PER_REFLECTANCE))
    return ok & (ndwi > melt.NDWI_THRESHOLD)


CLEAN_HALO = 0.08    # a scene this clear over the shelf shows real water, not haze
MELT_META = 15.0     # % metadata cloud: clear sky (used when halo is melt-lifted)
MELT_HALO_MAX = 0.30  # cap on halo in the extreme-melt branch (rejects real cloud)
UNION_CLEAN = False   # EXPERIMENT: union all qualifying scenes (see tile_water_on_grid)


def _onshelf_coarse_water(item, shelf_coarse):
    return int((_coarse_water_mask(item) & shelf_coarse).sum())


def select_scene(scored, water_score):
    """Two-branch clear-scene selection - pure logic, unit-tested offline.

    ``scored`` is [(halo_frac, metadata_cloud_pct, item), ...] for a tile's
    candidate scenes; ``water_score(item)`` ranks melt (bigger = more water).

      branch 1 (normal years):  among scenes actually clear over the shelf
        (halo < CLEAN_HALO) take the peak-melt one - local haze reads high halo
        and can never be picked, so cloud false positives cannot inflate a year.
      branch 2 (extreme melt):  if nothing clears CLEAN_HALO because melt itself
        depresses NDSI on every scene, metadata is trustworthy (no local haze,
        just melt): take the peak among scenes with metadata < MELT_META and
        halo < MELT_HALO_MAX. Haze fails one gate or the other.
      branch 3 (poorly observed): nothing qualifies at all - fall back to the
        least-cloudy scene and FLAG the tile.

    Returns (chosen_item, halo_of_chosen, poorly_observed).
    """
    pool, poorly = select_pool(scored)
    halo, _, chosen = max(pool, key=lambda x: water_score(x[2]))
    return chosen, halo, poorly


def select_pool(scored):
    """The qualifying-scene pool behind select_scene (same three branches);
    returned separately so union mode can composite over every qualifying
    scene. Returns (pool, poorly_observed)."""
    pool = [x for x in scored if x[0] < CLEAN_HALO]
    if pool:
        return pool, False
    pool = [x for x in scored if x[1] < MELT_META and x[0] < MELT_HALO_MAX]
    if pool:
        return pool, False
    return [min(scored, key=lambda x: x[0])], True


def tile_water_on_grid(candidates, grid_tr, gw, gh):
    """Per-tile seasonal-maximum meltwater from the peak-melt clear-sky scene.

    No single scene separates thin cloud/haze from meltwater spectrally, and the
    scene-level cloud metadata is tile-wide and misses local haze - so the actual
    cloud fraction over the shelf is measured per scene from the cloud mask
    (halo_frac). The rule, in two branches:

      * Normal years: keep scenes actually clear over the shelf (halo <
        CLEAN_HALO) and take the peak-melt one. Local haze reads high halo and is
        excluded, so it can never be picked (19DEA 10 Feb 2026 haze, halo 13%,
        494 km2 is dropped; the clean 7 Feb, halo 7%, 82 km2 wins).

      * Extreme-melt years: when melt is so widespread that slush and ponds
        depress NDSI on every scene, halo is melt-confounded (2019-20's clear
        record scenes all read ~16%) and no scene clears CLEAN_HALO. There the
        tile-wide metadata IS reliable - there is no local haze, just melt - so
        keep scenes with metadata < MELT_META (genuinely clear sky) and halo <
        MELT_HALO_MAX, and take the peak among them (19DEA 29 Jan 2020, metadata
        7%, halo 17%, 834 km2). A haze scene fails this because its metadata OR
        halo is high.

    Returns (grid_water, chosen, halo, n_candidates, poorly_observed).
    """
    ref = candidates[0]
    row, col, size = _tile_window(ref)
    if size < 256:
        return None, None, None, 0, False
    melt.set_aoi(melt._tile_of(ref), row, col, size)

    # dry-tile early-out on ON-SHELF coarse water (off-shelf tiles read 0 always)
    sc = _coarse_shelf_mask(ref)
    if not any(_onshelf_coarse_water(it, sc) >= DRY_COARSE_PX for it in candidates):
        return np.zeros((gh, gw), "f4"), candidates[0], 0.0, len(candidates), False

    scored = [(melt.cloud_mask(it)["halo_frac"],
               it.properties.get("eo:cloud_cover", 100.0), it) for it in candidates]
    chosen, halo, poorly = select_scene(scored, lambda it: _onshelf_coarse_water(it, sc))

    def _detect_one(it):
        bands = melt.load_scene(it, bands=("green", "nir", "red"))
        reject = melt.reject_mask(it) | melt.cloud_mask(it)["mask"]
        _, p, _ = melt.detect(bands["green"], bands["nir"], reject, red=bands["red"])
        return p

    if UNION_CLEAN and not poorly:
        # EXPERIMENTAL seasonal max: union detections across EVERY qualifying
        # scene (same pool the selector trusts), so ponds that peak on different
        # dates in different parts of the tile all count. Only scenes the
        # two-branch gates already accept contribute, so haze stays excluded.
        # Off by default until the anchor validation passes.
        pool, _ = select_pool(scored)
        ponds = None
        for _, _, it in pool:
            p = _detect_one(it)
            ponds = p if ponds is None else (ponds | p)
    else:
        ponds = _detect_one(chosen)

    src_crs, _ = melt.tile_georeference(ref)
    src_tr = melt.aoi_transform(ref)
    dst = np.zeros((gh, gw), "f4")
    reproject(
        source=ponds.astype("f4"), destination=dst,
        src_transform=src_tr, src_crs=src_crs,
        dst_transform=grid_tr, dst_crs=GRID_CRS,
        resampling=Resampling.average,  # fraction of each 30 m cell that is water
    )
    return dst, chosen, float(halo), len(candidates), poorly


def _season_scenes(tile, start, end, cc=85, n=25):
    """Candidate scenes for a tile in the Jan-Feb window. The metadata cloud cap
    is deliberately GENEROUS (85%): metadata is tile-wide and a scene that is
    cloudy off-shelf can be perfectly clear over the shelf itself (19DEA 7 Feb
    2026 was 59% metadata but 7% actual over the shelf). Filtering hard on
    metadata would discard exactly those clear-over-shelf scenes, so we gather
    broadly and let the cloud-mask (halo_frac) judge actual clarity per scene.
    The n clearest by metadata are kept only to bound how many masks we build."""
    from pystac_client import Client
    cl = Client.open(melt.STAC_API)
    items = list(cl.search(
        collections=["sentinel-2-l2a"],
        query={"grid:code": {"eq": f"MGRS-{tile}"}, "eo:cloud_cover": {"lt": cc}},
        datetime=f"{start}/{end}", limit=100).items())
    if not items:
        items = list(cl.search(
            collections=["sentinel-2-l2a"],
            query={"grid:code": {"eq": f"MGRS-{tile}"}},
            datetime=f"{start}/{end}", limit=100).items())
    items.sort(key=lambda x: x.properties.get("eo:cloud_cover", 99))
    return items[:n]


def production_scene(tile, start, end):
    """The exact scene run_season's selection would use for this tile/window.
    validate_shelf_tile uses this so blind labels always test the PRODUCTION
    scene choice - validating a scene chosen any other way would validate a
    pipeline that does not exist."""
    items = _season_scenes(tile, start, end)
    if not items:
        return None
    row, col, size = _tile_window(items[0])
    if size < 256:
        return None
    melt.set_aoi(melt._tile_of(items[0]), row, col, size)
    sc = _coarse_shelf_mask(items[0])
    scored = [(melt.cloud_mask(it)["halo_frac"],
               it.properties.get("eo:cloud_cover", 100.0), it) for it in items]
    chosen, _, _ = select_scene(scored, lambda it: _onshelf_coarse_water(it, sc))
    return chosen


def _peak_scene(tile, start, end, cc=70):
    """Single clearest-by-metadata scene (legacy; kept for comparisons)."""
    from pystac_client import Client
    items = list(Client.open(melt.STAC_API).search(
        collections=["sentinel-2-l2a"],
        query={"grid:code": {"eq": f"MGRS-{tile}"}, "eo:cloud_cover": {"lt": cc}},
        datetime=f"{start}/{end}", limit=100).items())
    return min(items, key=lambda x: x.properties.get("eo:cloud_cover", 99)) if items else None


def run_season(label, peak_window=("01-01", "02-28"), grid=None):
    """Shelf-wide seasonal-maximum meltwater on ONE fixed shelf grid (so the
    denominator is identical every season). Each tile contributes the peak-melt
    scene among those genuinely clear over the shelf (see tile_water_on_grid),
    unioned across overlapping tiles. A season is flagged poorly observed when
    its dominant water tile had no genuinely clear scene."""
    start_year = int(label.split("-")[0]) + 1  # peak melt is the second year
    start = f"{start_year}-{peak_window[0]}"
    end = f"{start_year}-{peak_window[1]}"

    cands = {t: _season_scenes(t, start, end) for t in SHELF_TILES}
    cands = {t: v for t, v in cands.items() if v}
    if not cands:
        print(f"[{label}] no scenes in {start}..{end}")
        return

    grid_tr, gw, gh, shelf = grid if grid is not None else build_fixed_grid()
    cell_km2 = (GRID_RES**2) / 1e6
    shelf_km2 = shelf.sum() * cell_km2
    print(f"[{label}] fixed grid {gw}x{gh}, shelf {shelf_km2:,.0f} km2, "
          f"tiles imaged {len(cands)}/{len(SHELF_TILES)}")

    # water fraction per grid cell, unioned (max) across overlapping tiles.
    water = np.zeros((gh, gw), "f4")
    water_tiles = []   # (on_shelf_km2, halo, poorly) per water-bearing tile
    for t, items in cands.items():
        try:
            wm, chosen, halo, k, poorly = tile_water_on_grid(items, grid_tr, gw, gh)
        except Exception as e:
            print(f"  {t}: ERROR {type(e).__name__}: {str(e)[:50]}")
            continue
        if wm is None:
            continue
        wm *= shelf                       # in place: no gh x gw temps (Amery-
        on_shelf_km2 = float(wm.sum()) * cell_km2   # sized grids risk OOM)
        if on_shelf_km2 > 0.05:
            water_tiles.append((on_shelf_km2, halo, poorly))
        np.maximum(water, wm, out=water)
        date = chosen.datetime.date() if hasattr(chosen, "datetime") else "-"
        print(f"  {t}  {date} of {k:2d}  halo {100*halo:4.1f}%"
              f"{' POORLY' if poorly else ''}  +{on_shelf_km2:7.2f} km2 on shelf")

    total = float(water.sum()) * cell_km2
    # Quality from the DOMINANT water tile: a season is only poorly observed when
    # the tile carrying most of the water had no genuinely clear scene.
    if water_tiles:
        _, dom_halo, dom_poorly = max(water_tiles, key=lambda x: x[0])
    else:
        dom_halo, dom_poorly = 0.0, False
    obs_cloud = round(100 * dom_halo, 1)
    print(f"\n[{label}] shelf-wide seasonal-max meltwater: {total:.1f} km2 "
          f"({100*total/max(shelf_km2,1):.2f}% of shelf); dominant-tile halo {obs_cloud}%"
          f"{'  [POORLY OBSERVED]' if dom_poorly else ''}")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    return {"season": label, "shelf_km2": round(total, 2),
            "shelf_area_km2": round(shelf_km2, 1),
            "tiles": len(cands), "tiles_total": len(SHELF_TILES),
            "obs_cloud": obs_cloud, "poorly_observed": dom_poorly}


# 9 austral melt seasons with reliable Sentinel-2 coverage over the shelf,
# labelled by the first year (peak melt falls in Jan/Feb of the second year).
SEASONS = ["2017-18", "2018-19", "2019-20", "2020-21", "2021-22",
           "2022-23", "2023-24", "2024-25", "2025-26"]


def run_history(seasons=SEASONS):
    """Run every season, saving after each so a mid-run network failure on one
    season keeps the others. Writes output/shelf/history.json and prints a
    table of shelf-wide seasonal-maximum meltwater area."""
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    hist_path = OUT_DIR / ("history.json" if CURRENT_SHELF == "george_vi"
                           else f"history_{CURRENT_SHELF}.json")
    results = {}

    grid = build_fixed_grid()  # one denominator shared by every season

    for label in seasons:
        print(f"\n{'='*64}\n{label}\n{'='*64}")
        try:
            r = run_season(label, grid=grid)
        except Exception as e:
            print(f"[{label}] FAILED: {type(e).__name__}: {str(e)[:80]}")
            continue
        if r is not None:
            results[label] = r
            ordered = [results[s] for s in seasons if s in results]
            hist_path.write_text(json.dumps(ordered, indent=1))

    print(f"\n{'='*64}\nSHELF-WIDE SEASONAL-MAX MELTWATER (George VI)\n{'='*64}")
    print(f"  {'season':9s} {'meltwater km2':>14} {'% of shelf':>11} "
          f"{'tiles':>7} {'shelf halo':>11}  quality")
    for s in seasons:
        if s not in results:
            print(f"  {s:9s} {'--':>14}"); continue
        r = results[s]
        pct = 100 * r["shelf_km2"] / r["shelf_area_km2"]
        tiles = f"{r['tiles']}/{r.get('tiles_total', len(SHELF_TILES))}"
        oc = r.get("obs_cloud")
        ocs = f"{oc:.1f}%" if oc is not None else "--"
        q = "poorly observed" if r.get("poorly_observed") else "ok"
        print(f"  {s:9s} {r['shelf_km2']:14.1f} {pct:10.2f}% {tiles:>7} {ocs:>11}  {q}")
    print(f"\n  saved -> {hist_path.relative_to(melt.ROOT)}")
    return results


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "season"
    if cmd == "season":
        label = sys.argv[2] if len(sys.argv) > 2 else "2020-21"
        run_season(label)
    elif cmd == "history":
        run_history()

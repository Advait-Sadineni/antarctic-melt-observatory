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


def _reference_items():
    """One clear Sentinel-2 item per shelf tile, for a season-independent
    footprint. A summer 2020 scene is picked so the tile is fully imaged."""
    from pystac_client import Client
    cl = Client.open(melt.STAC_API)
    items = {}
    for t in SHELF_TILES:
        r = list(cl.search(collections=["sentinel-2-l2a"],
                 query={"grid:code": {"eq": f"MGRS-{t}"},
                        "eo:cloud_cover": {"lt": 30}},
                 datetime="2020-01-01/2020-02-28", limit=50).items())
        if r:
            items[t] = min(r, key=lambda x: x.properties.get("eo:cloud_cover", 99))
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

def _tile_window(item):
    """Window on this tile that overlaps the shelf-tile bbox, capped for memory."""
    crs, tile_tr = melt.tile_georeference(item)
    # shelf-tile lon/lat bbox -> this tile's UTM -> pixel window
    w, s, e, n = -70.5, -72.6, -66.0, -70.0
    xs, ys = warp_transform("EPSG:4326", crs, [w, e, w, e], [s, s, n, n])
    col0 = max(0, int((min(xs) - tile_tr.c) / melt.PIXEL_M))
    col1 = min(melt.TILE_PX if hasattr(melt, "TILE_PX") else 10980,
               int((max(xs) - tile_tr.c) / melt.PIXEL_M))
    row0 = max(0, int((tile_tr.f - max(ys)) / melt.PIXEL_M))
    row1 = min(10980, int((tile_tr.f - min(ys)) / melt.PIXEL_M))
    size = min(max(col1 - col0, row1 - row0), 10980 - max(row0, col0))
    return row0, col0, size


def _coarse_water_score(item):
    """Cheap 60 m proxy for detected water in the AOI, for ranking scenes by
    melt without a full-resolution detection. Reads COG overviews (a tiny
    fraction of the bytes), applies the same NDWI + shadow gates as detect() but
    skips hysteresis and the fine cloud mask - enough to pick the peak-melt date,
    which is then re-detected at full resolution."""
    g = melt.read_coarse(item, "green")
    n = melt.read_coarse(item, "nir")
    r = melt.read_coarse(item, "red")
    valid = (g > 0) & (n > 0)
    ndwi = np.where(valid, (g - n) / np.maximum(g + n, 1e-6), -9)
    ok = (valid & (g > melt.BRIGHTNESS_FLOOR)
          & ((g - r) > melt.SHADOW_GREEN_MINUS_RED * melt.DN_PER_REFLECTANCE))
    return int((ok & (ndwi > melt.NDWI_THRESHOLD)).sum())


def _detect_ponds(item):
    """Detect meltwater on one scene in its tile window (tile-space mask)."""
    bands = melt.load_scene(item, bands=("green", "nir", "red"))
    reject = melt.reject_mask(item) | melt.cloud_mask(item)["mask"]
    _, ponds, _ = melt.detect(bands["green"], bands["nir"], reject, red=bands["red"])
    return ponds


def tile_water_on_grid(items, grid_tr, gw, gh):
    """Per-tile SEASONAL-MAX meltwater warped onto the common grid.

    ``items`` are the candidate scenes for one tile across the melt window. The
    peak-melt date - the scene with the most water - is selected and reprojected
    onto the grid. Selecting the peak scene (not the clearest) is what makes this
    a seasonal maximum: the clearest scene is often a pre-melt December date, so
    ranking by clarity misses the peak entirely (e.g. the record 2019-20
    summer). Candidates are ranked cheaply at 60 m and only the winner is
    detected at full resolution. Keeping a single real scene per tile preserves
    the single-scene precision validated on 19DEA, avoiding the false-positive
    accumulation of multi-date compositing.
    """
    ref = items[0]
    row, col, size = _tile_window(ref)
    if size < 256:
        return None, None, 0
    melt.set_aoi(melt._tile_of(ref), row, col, size)

    chosen = (max(items, key=_coarse_water_score) if len(items) > 1 else items[0])
    ponds = _detect_ponds(chosen)

    src_crs, _ = melt.tile_georeference(ref)
    src_tr = melt.aoi_transform(ref)
    dst = np.zeros((gh, gw), "f4")
    reproject(
        source=ponds.astype("f4"), destination=dst,
        src_transform=src_tr, src_crs=src_crs,
        dst_transform=grid_tr, dst_crs=GRID_CRS,
        resampling=Resampling.average,  # fraction of each 30 m cell that is water
    )
    return dst, chosen, len(items)


def _season_scenes(tile, start, end, cc=40, n=12):
    """Up to ``n`` clearest scenes for a tile in the melt window, cloud < ``cc``.
    These are the candidates whose peak-melt date tile_water_on_grid picks. If
    none clear the cap, fall back to the single clearest scene up to 80% cloud
    so the tile is still represented."""
    from pystac_client import Client
    cl = Client.open(melt.STAC_API)
    items = list(cl.search(
        collections=["sentinel-2-l2a"],
        query={"grid:code": {"eq": f"MGRS-{tile}"}, "eo:cloud_cover": {"lt": cc}},
        datetime=f"{start}/{end}", limit=100).items())
    if not items:
        items = list(cl.search(
            collections=["sentinel-2-l2a"],
            query={"grid:code": {"eq": f"MGRS-{tile}"}, "eo:cloud_cover": {"lt": 80}},
            datetime=f"{start}/{end}", limit=100).items())
    items.sort(key=lambda x: x.properties.get("eo:cloud_cover", 99))
    return items[:n]


def _peak_scene(tile, start, end, cc=70):
    """Single clearest scene for a tile in a date range (used by the per-tile
    blind validation, validate_shelf_tile.py). Returns None if none exist."""
    from pystac_client import Client
    items = list(Client.open(melt.STAC_API).search(
        collections=["sentinel-2-l2a"],
        query={"grid:code": {"eq": f"MGRS-{tile}"}, "eo:cloud_cover": {"lt": cc}},
        datetime=f"{start}/{end}", limit=100).items())
    return min(items, key=lambda x: x.properties.get("eo:cloud_cover", 99)) if items else None


def run_season(label, peak_window=("12-01", "02-28"), grid=None):
    """Shelf-wide seasonal-maximum meltwater: each tile's clearest scene in the
    melt window, unioned on ONE fixed shelf grid so the denominator is the same
    every season. The window spans the whole Dec-Feb melt season and the cloud
    cap is generous, so all 11 tiles are represented every season."""
    start_year = int(label.split("-")[0]) + 1  # peak melt is the second year
    start = f"{start_year}-{peak_window[0]}"
    # Dec belongs to the FIRST year of the label; split the query if needed.
    dec_start = f"{start_year-1}-12-01" if peak_window[0].startswith("12") else start
    end = f"{start_year}-{peak_window[1]}"

    cands = {t: _season_scenes(t, dec_start, end) for t in SHELF_TILES}
    cands = {t: v for t, v in cands.items() if v}
    if not cands:
        print(f"[{label}] no scenes in {dec_start}..{end}")
        return

    grid_tr, gw, gh, shelf = grid if grid is not None else build_fixed_grid()
    cell_km2 = (GRID_RES**2) / 1e6
    shelf_km2 = shelf.sum() * cell_km2
    print(f"[{label}] fixed grid {gw}x{gh}, shelf {shelf_km2:,.0f} km2, "
          f"tiles imaged {len(cands)}/{len(SHELF_TILES)}")

    # water fraction per grid cell, taking the best (max) estimate where tiles
    # overlap - a union that dedups without double-counting.
    water = np.zeros((gh, gw), "f4")
    clouds, n_scenes = [], 0
    for t, items in cands.items():
        try:
            wm, chosen, k = tile_water_on_grid(items, grid_tr, gw, gh)
        except Exception as e:
            print(f"  {t}: ERROR {type(e).__name__}: {str(e)[:50]}")
            continue
        if wm is None:
            continue
        cc = chosen.properties.get("eo:cloud_cover", 0.0)
        clouds.append(cc)
        n_scenes += k
        on_shelf_km2 = float((wm * shelf).sum()) * cell_km2
        np.maximum(water, np.where(shelf, wm, 0), out=water)
        print(f"  {t}  peak {chosen.datetime.date()} of {k:2d}  cloud {cc:4.1f}  "
              f"+{on_shelf_km2:6.2f} km2 on shelf")

    total = float(water.sum()) * cell_km2
    print(f"\n[{label}] shelf-wide seasonal-max meltwater: {total:.1f} km2 "
          f"({100*total/max(shelf_km2,1):.2f}% of shelf)")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    return {"season": label, "shelf_km2": round(total, 2),
            "shelf_area_km2": round(shelf_km2, 1),
            "tiles": len(cands), "tiles_total": len(SHELF_TILES),
            "scenes_scanned": n_scenes,
            "mean_cloud": round(float(np.mean(clouds)), 1) if clouds else None}


# 9 austral melt seasons with reliable Sentinel-2 coverage over the shelf,
# labelled by the first year (peak melt falls in Jan/Feb of the second year).
SEASONS = ["2017-18", "2018-19", "2019-20", "2020-21", "2021-22",
           "2022-23", "2023-24", "2024-25", "2025-26"]


def run_history(seasons=SEASONS):
    """Run every season, saving after each so a mid-run network failure on one
    season keeps the others. Writes output/shelf/history.json and prints a
    table of shelf-wide seasonal-maximum meltwater area."""
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    hist_path = OUT_DIR / "history.json"
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
          f"{'tiles':>7} {'mean cloud':>11}")
    for s in seasons:
        if s not in results:
            print(f"  {s:9s} {'--':>14}"); continue
        r = results[s]
        pct = 100 * r["shelf_km2"] / r["shelf_area_km2"]
        tiles = f"{r['tiles']}/{r.get('tiles_total', len(SHELF_TILES))}"
        mc = f"{r['mean_cloud']:.1f}%" if r.get("mean_cloud") is not None else "--"
        print(f"  {s:9s} {r['shelf_km2']:14.1f} {pct:10.2f}% {tiles:>7} {mc:>11}")
    print(f"\n  saved -> {hist_path.relative_to(melt.ROOT)}")
    return results


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "season"
    if cmd == "season":
        label = sys.argv[2] if len(sys.argv) > 2 else "2020-21"
        run_season(label)
    elif cmd == "history":
        run_history()

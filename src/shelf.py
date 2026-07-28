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

STATUS: the engine is validated and now uses an authoritative boundary. Per-
tile detection reprojects and unions correctly on the common grid; the MEaSUReS
polygon correctly delineates the flat shelf (mountain tiles 19DDA/19DDC read 0,
and its outline hugs the shelf edge on inspection). Two of the three issues
found earlier are resolved:

  - [fixed] Boundary precision. Swapped the cartographic Natural Earth polygon
    for the authoritative MEaSUReS Antarctic Boundaries v2 outline (23,260 km2,
    matching the known shelf), which is tight at the calving front and excludes
    open ocean / sea ice.
  - [fixed] Resampling inflation. Reprojection is now area-weighted (average),
    giving each 30 m grid cell its true water fraction instead of marking the
    whole cell water if any 10 m pixel is. This trimmed the tile areas ~15-25%.

  - [OPEN] Per-tile validation. For 2020-21 the shelf-wide seasonal-max is
    138 km2 (0.77% of an 17,878 km2 shelf), but it is dominated by tile 19DEA
    (127 km2). Inspection shows that detection is spectrally water-like (green
    ~7100, the confirmed-pond range) and spatially an organised field of
    channel-like streaks on the shelf - consistent with the extensive
    supraglacial meltwater channels George VI is documented to have, but the
    regularity could also be flow-stripe features. It is real water, not slush
    or cloud, but whether all of it is meltwater needs the same blind-label
    validation used elsewhere in this project, applied per shelf tile. Until
    that is done the absolute shelf-wide number is provisional.

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


def tile_water_on_grid(item, grid_tr, gw, gh):
    """Detect meltwater on one tile and reproject the mask onto the common grid."""
    row, col, size = _tile_window(item)
    if size < 256:
        return None
    melt.set_aoi(melt._tile_of(item), row, col, size)
    bands = melt.load_scene(item, bands=("green", "nir", "red"))
    reject = melt.reject_mask(item) | melt.cloud_mask(item)["mask"]
    _, ponds, _ = melt.detect(bands["green"], bands["nir"], reject, red=bands["red"])

    src_crs, _ = melt.tile_georeference(item)
    src_tr = melt.aoi_transform(item)
    dst = np.zeros((gh, gw), "f4")
    reproject(
        source=ponds.astype("f4"), destination=dst,
        src_transform=src_tr, src_crs=src_crs,
        dst_transform=grid_tr, dst_crs=GRID_CRS,
        resampling=Resampling.average,  # fraction of each 30 m cell that is water
    )
    return dst  # 0..1 water fraction per grid cell (no coarsening inflation)


def _peak_scene(tile, start, end):
    """Clearest scene for a tile in a date range (crude peak-melt pick)."""
    from pystac_client import Client
    items = list(Client.open(melt.STAC_API).search(
        collections=["sentinel-2-l2a"],
        query={"grid:code": {"eq": f"MGRS-{tile}"}, "eo:cloud_cover": {"lt": 15}},
        datetime=f"{start}/{end}", limit=100).items())
    return min(items, key=lambda x: x.properties.get("eo:cloud_cover", 99)) if items else None


def run_season(label, peak_window=("01-05", "02-20")):
    """Shelf-wide seasonal-maximum meltwater: each tile's clearest peak-melt
    scene, unioned on the common grid, masked to the shelf polygon."""
    start_year = int(label.split("-")[0]) + 1  # Jan/Feb of the second year
    start = f"{start_year}-{peak_window[0]}"
    end = f"{start_year}-{peak_window[1]}"

    scenes = {t: _peak_scene(t, start, end) for t in SHELF_TILES}
    scenes = {t: it for t, it in scenes.items() if it is not None}
    if not scenes:
        print(f"[{label}] no scenes in {start}..{end}")
        return

    grid_tr, gw, gh, shelf = build_grid(list(scenes.values()))
    print(f"[{label}] grid {gw}x{gh} @ {GRID_RES:.0f} m, "
          f"shelf mask {shelf.sum() * (GRID_RES**2) / 1e6:,.0f} km2")

    cell_km2 = (GRID_RES**2) / 1e6
    # water fraction per grid cell, taking the best (max) estimate where tiles
    # overlap - a union that dedups without double-counting.
    water = np.zeros((gh, gw), "f4")
    for t, it in scenes.items():
        try:
            wm = tile_water_on_grid(it, grid_tr, gw, gh)
        except Exception as e:
            print(f"  {t}: ERROR {type(e).__name__}: {str(e)[:50]}")
            continue
        if wm is None:
            continue
        on_shelf_km2 = float((wm * shelf).sum()) * cell_km2
        np.maximum(water, np.where(shelf, wm, 0), out=water)
        print(f"  {t}  {it.datetime.date()}  cloud {it.properties.get('eo:cloud_cover'):4.1f}  "
              f"+{on_shelf_km2:6.2f} km2 on shelf")

    total = float(water.sum()) * cell_km2
    print(f"\n[{label}] shelf-wide seasonal-max meltwater: {total:.1f} km2 "
          f"({100*total/(shelf.sum()*cell_km2):.2f}% of shelf)")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    return {"season": label, "shelf_km2": round(total, 2),
            "shelf_area_km2": round(shelf.sum() * (GRID_RES**2) / 1e6, 1),
            "tiles": len(scenes)}


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "season"
    if cmd == "season":
        label = sys.argv[2] if len(sys.argv) > 2 else "2020-21"
        run_season(label)

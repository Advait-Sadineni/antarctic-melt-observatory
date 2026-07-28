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

The shelf polygon is a drop-in parameter (reference/george_vi_boundary.geojson,
currently Natural Earth 10 m - area 24,163 km2, matching the known ~24,000 km2
whole shelf). Swap in MEaSUREs/SCAR for final precision without touching code.

Per-tile detection reuses the validated single-tile detector via melt.set_aoi,
so the shelf number is produced by exactly the method the rest of the project
validated (Moussavi shadow test + hysteresis).

STATUS: the engine is validated (per-tile detection reprojects and unions
correctly on the common grid; the shelf polygon correctly excludes the
Alexander Island mountain tiles - 19DDA/19DDC read 0). The absolute shelf-wide
number is NOT yet trustworthy, for three reasons found by inspection and
recorded here rather than hidden:

  1. Boundary precision at the shelf front. The Natural Earth polygon is
     cartographic-grade; near the calving front it can include a strip of open
     ocean / sea ice that reads as high-NDWI water. A precise MEaSUREs/SCAR
     boundary (a polygon swap) is the fix.
  2. Slush / thin cloud at shelf scale. Tile 19DEA on 2021-01-24 contributes a
     large bright detection (~140 km2) that inspection shows is mostly extensive
     slush or thin cloud, not deep ponds - the known NDWI/slush ambiguity, now
     visible at shelf scale where it dominates. Needs slush handling or a
     per-tile plausibility gate.
  3. Max-resampling inflation. Reprojecting the 10 m mask onto the 30 m grid
     with Resampling.max marks a whole 30 m cell water if any 10 m pixel is,
     inflating area ~20%. A fractional (area-weighted) reprojection removes it.

The core sound tiles (19CDV + 19CEV, deduped) give ~15 km2 for 2020-21, which
is the sensible part of the number; 19DEA's 176 km2 is the questionable part.

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
BOUNDARY = melt.ROOT / "reference" / "george_vi_boundary.geojson"
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
    dst = np.zeros((gh, gw), "uint8")
    reproject(
        source=ponds.astype("uint8"), destination=dst,
        src_transform=src_tr, src_crs=src_crs,
        dst_transform=grid_tr, dst_crs=GRID_CRS,
        resampling=Resampling.max,  # a pond present in any source pixel -> present
    )
    return dst.astype(bool)


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

    water = np.zeros((gh, gw), bool)
    for t, it in scenes.items():
        try:
            wm = tile_water_on_grid(it, grid_tr, gw, gh)
        except Exception as e:
            print(f"  {t}: ERROR {type(e).__name__}: {str(e)[:50]}")
            continue
        if wm is None:
            continue
        on_shelf = wm & shelf
        water |= on_shelf
        print(f"  {t}  {it.datetime.date()}  cloud {it.properties.get('eo:cloud_cover'):4.1f}  "
              f"+{on_shelf.sum() * (GRID_RES**2)/1e6:6.2f} km2 on shelf")

    total = water.sum() * (GRID_RES**2) / 1e6
    print(f"\n[{label}] shelf-wide seasonal-max meltwater: {total:.1f} km2 "
          f"({100*total/(shelf.sum()*(GRID_RES**2)/1e6):.2f}% of shelf)")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    return {"season": label, "shelf_km2": round(total, 2),
            "shelf_area_km2": round(shelf.sum() * (GRID_RES**2) / 1e6, 1),
            "tiles": len(scenes)}


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "season"
    if cmd == "season":
        label = sys.argv[2] if len(sys.argv) > 2 else "2020-21"
        run_season(label)

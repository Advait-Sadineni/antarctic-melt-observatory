"""Region-agnostic multi-shelf driver.

Runs the exact validated George VI machinery (shelf.py) on ANY Antarctic ice
shelf, straight from its MEaSUReS boundary:

  1. extract_boundary(name) - pull the shelf polygon from the MEaSUReS Antarctic
     Boundaries v2 shapefile and write it as GeoJSON (EPSG:4326).
  2. discover_tiles(boundary) - find the Sentinel-2 MGRS tiles that cover it.
  3. run_shelf(name) - set_shelf() + build the fixed grid + run the 9-season
     history, saving to output/shelf/history_<name>.json.

This is the region-agnostic keystone from the observatory design: a shelf is a
parameter, not hard-coded. New shelves are first-pass (not yet blind-validated
per region); George VI remains the validated anchor.

Run:  python src/shelves.py Bach Stange Wilkins LarsenC
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import melt
import shelf

MEASURES = melt.ROOT / "reference" / "measures" / "IceShelf_Antarctica_v02.shp"
BND_DIR = melt.ROOT / "reference" / "boundaries"


def extract_boundary(name):
    """Write reference/boundaries/<name>.geojson (EPSG:4326) from MEaSUReS."""
    BND_DIR.mkdir(parents=True, exist_ok=True)
    out = BND_DIR / f"{name}.geojson"
    if out.exists():
        return out
    import geopandas as gpd
    gdf = gpd.read_file(MEASURES)
    sub = gdf[gdf["NAME"] == name]
    if sub.empty:
        raise ValueError(f"shelf {name!r} not found in MEaSUReS")
    sub = sub.to_crs(4326).dissolve()          # one (multi)polygon in lon/lat
    out.write_text(sub.to_json())
    return out


def _bbox(boundary_path):
    gj = json.loads(Path(boundary_path).read_text())
    xs, ys = [], []
    for f in gj["features"]:
        g = f["geometry"]
        polys = g["coordinates"] if g["type"] == "MultiPolygon" else [g["coordinates"]]
        for poly in polys:
            for ring in poly:
                for x, y in ring:
                    xs.append(x); ys.append(y)
    return [min(xs), min(ys), max(xs), max(ys)]


def discover_tiles(boundary_path):
    """Sentinel-2 MGRS tiles intersecting the shelf's bounding box. Tiles that
    do not actually overlap the polygon read zero (on-shelf mask + dry-skip), so
    a bbox-level list is safe; it just may include a few empty tiles."""
    src = melt.get_source()
    bbox = _bbox(boundary_path)
    tiles = set()
    for yr in ("2020", "2021", "2022"):        # a few summers to enumerate tiles
        items = src.search(bbox=bbox,
                           datetime=f"{yr}-12-15/{int(yr)+1}-02-15",
                           query={"eo:cloud_cover": {"lt": 50}}, limit=100)
        for it in items:
            try:
                tiles.add(src.tile_of(it))
            except (KeyError, IndexError):
                pass
    return sorted(tiles), bbox


def prepare(name):
    bnd = extract_boundary(name)
    tiles, bbox = discover_tiles(bnd)
    print(f"[{name}] boundary {bnd.name}, bbox {[round(v,1) for v in bbox]}, "
          f"{len(tiles)} tiles: {tiles}")
    return {"name": name, "boundary": str(bnd), "tiles": tiles}


def run_shelf(name, seasons=None):
    cfg = prepare(name)
    if not cfg["tiles"]:
        print(f"[{name}] no tiles found - skipping")
        return
    shelf.set_shelf(cfg["name"], cfg["tiles"], cfg["boundary"])
    shelf.build_fixed_grid(rebuild=True)
    shelf.run_history(seasons=seasons or shelf.SEASONS)


if __name__ == "__main__":
    targets = sys.argv[1:] or ["Bach", "Stange", "Wilkins", "LarsenC"]
    for t in targets:
        print(f"\n{'#'*68}\n#### SHELF: {t}\n{'#'*68}")
        try:
            run_shelf(t)
        except Exception as e:
            print(f"#### SHELF {t} FAILED: {type(e).__name__}: {str(e)[:120]}")

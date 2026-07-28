"""
Export detected meltwater as GIS-ready layers: a georeferenced GeoTIFF mask
and a GeoJSON of pond polygons with per-pond areas.

PNGs are for looking at; glaciologists work in QGIS/ArcGIS. This writes the
detection in the study area's real projected coordinates (UTM 19S), so it drops
straight onto a basemap and overlays with ICESat-2, velocity, or DEM layers.

Run:  python src/export_gis.py S2B_19CEV_20210124_0_L2A
      python src/export_gis.py                 # default peak scene
"""

import json
import sys

import numpy as np
import rasterio
from rasterio.features import shapes

import melt

DEFAULT_ITEM = "S2B_19CEV_20210124_0_L2A"
OUT_DIR = melt.ROOT / "output" / "gis"

MASK_NODATA = 255  # cloud/shadow/off-swath, distinct from 0 (dry) and 1 (water)


def export(item_id, min_pond_px=4):
    item = melt.get_item(item_id)
    bands = melt.load_scene(item, bands=("green", "nir", "red"))
    reject = melt.reject_mask(item) | melt.cloud_mask(item)["mask"]
    _, ponds, valid = melt.detect(bands["green"], bands["nir"], reject, red=bands["red"])

    usable = valid & ~reject
    mask = np.full(ponds.shape, MASK_NODATA, dtype="uint8")
    mask[usable] = 0
    mask[ponds & usable] = 1

    crs, _ = melt.tile_georeference(item)
    transform = melt.aoi_transform(item)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    tif = OUT_DIR / f"{item_id}_meltwater.tif"
    profile = {
        "driver": "GTiff", "height": mask.shape[0], "width": mask.shape[1],
        "count": 1, "dtype": "uint8", "crs": crs, "transform": transform,
        "nodata": MASK_NODATA, "compress": "deflate", "tiled": True,
    }
    with rasterio.open(tif, "w", **profile) as dst:
        dst.write(mask, 1)
        dst.update_tags(1, meltwater="1", dry="0", nodata=str(MASK_NODATA))
        dst.update_tags(
            scene=item_id, date=item.datetime.date().isoformat(),
            method="NDWI>0.19 + green-red>0.09 shadow test + hysteresis(0.14)",
            ndwi_threshold=str(melt.NDWI_THRESHOLD),
        )

    # Vectorise the water pixels into pond polygons with areas.
    features = []
    for geom, val in shapes(mask, mask=(mask == 1), transform=transform):
        area_m2 = _ring_area(geom)
        if area_m2 < min_pond_px * melt.PIXEL_M**2:
            continue
        features.append({
            "type": "Feature",
            "geometry": geom,
            "properties": {"area_km2": round(area_m2 / 1e6, 5),
                           "area_m2": round(area_m2, 1)},
        })
    features.sort(key=lambda f: -f["properties"]["area_m2"])

    gj = OUT_DIR / f"{item_id}_ponds.geojson"
    fc = {
        "type": "FeatureCollection",
        "name": f"{item_id}_meltwater_ponds",
        "crs": {"type": "name", "properties": {"name": str(crs)}},
        "metadata": {
            "scene": item_id, "date": item.datetime.date().isoformat(),
            "n_ponds": len(features),
            "total_km2": round(sum(f["properties"]["area_km2"] for f in features), 4),
        },
        "features": features,
    }
    gj.write_text(json.dumps(fc))

    total = fc["metadata"]["total_km2"]
    print(f"  scene     {item_id}")
    print(f"  water     {int((mask == 1).sum()):,} px = {total:.3f} km2")
    print(f"  ponds     {len(features)} (>= {min_pond_px} px), largest "
          f"{features[0]['properties']['area_km2']:.3f} km2" if features else "  ponds     none")
    print(f"[write] {tif.relative_to(melt.ROOT)}")
    print(f"[write] {gj.relative_to(melt.ROOT)}")
    return tif, gj


def _ring_area(geom):
    """Planar polygon area (m²) via the shoelace formula, minus holes.

    Coordinates are already projected metres (UTM), so shoelace gives true
    area without any geographic correction.
    """
    def ring(coords):
        x = np.array([c[0] for c in coords])
        y = np.array([c[1] for c in coords])
        return 0.5 * abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))

    rings = geom["coordinates"]
    if not rings:
        return 0.0
    area = ring(rings[0])
    for hole in rings[1:]:
        area -= ring(hole)
    return area


def main():
    item_id = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_ITEM
    export(item_id)


if __name__ == "__main__":
    main()

"""The observatory's data model, v0: one product per (shelf, season).

A cloud-optimized GeoTIFF of per-cell water fraction on the fixed EPSG:3031
grid, plus a STAC-style item.json. Reserved fields (null for now) are the
full-vision schema from the design doc - M3 (sensor_confidence) and M4
(depth/volume/uncertainty) fill them in place; nothing ever migrates.
"""
import json
from pathlib import Path

import rasterio

RESERVED = ("volume_km3", "depth_mean_m", "uncertainty_km2", "sensor_confidence")


def write_product(shelf_name, season, water, transform, crs, meta, root):
    """Write water_fraction.tif (tiled + overviews) and item.json; returns dir."""
    d = Path(root) / shelf_name / season
    d.mkdir(parents=True, exist_ok=True)

    profile = dict(driver="GTiff", height=water.shape[0], width=water.shape[1],
                   count=1, dtype="float32", crs=crs, transform=transform,
                   tiled=True, blockxsize=512, blockysize=512,
                   compress="deflate", predictor=2)
    tif = d / "water_fraction.tif"
    with rasterio.open(tif, "w", **profile) as dst:
        dst.write(water.astype("float32"), 1)
        dst.build_overviews([2, 4, 8, 16], rasterio.enums.Resampling.average)

    props = {"shelf": shelf_name, "season": season, **meta}
    for k in RESERVED:
        props.setdefault(k, None)
    item = {"type": "Feature", "stac_version": "1.0.0",
            "id": f"{shelf_name}-{season}",
            "properties": props,
            "assets": {"water_fraction": {
                "href": "./water_fraction.tif",
                "type": "image/tiff; application=geotiff",
                "roles": ["data"]}}}
    (d / "item.json").write_text(json.dumps(item, indent=1))
    return d

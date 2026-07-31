"""M4 Task 3 runner: pond polygons -> ICESat-2 crossover harvest (one season).

ATL03 granules are multi-GB, so: only the biggest ponds (they're where the
laser statistics live anyway), max 2 granules per pond-date, and each granule
is deleted right after its photons are extracted.

Run:  python scripts/depth_harvest.py 2020-21
"""
import json
import sys
from pathlib import Path

import numpy as np
from rasterio.warp import transform as warp_transform
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
import icesat2
import melt
import shelf

MIN_POND_KM2 = 0.5
MAX_PONDS = 8
MAX_GRANULES_PER_POND = 2
WATER_TILES = ("19DEA", "19CEV")


def pond_polygons(season):
    """Biggest ponds of the season's chosen scenes as lon/lat bboxes."""
    year = int(season.split("-")[0]) + 1
    feats = []
    for tile in WATER_TILES:
        it = shelf.production_scene(tile, f"{year}-01-01", f"{year}-02-28")
        if it is None:
            continue
        bands = melt.load_scene(it, bands=("green", "nir", "red"))
        reject = melt.reject_mask(it) | melt.cloud_mask(it)["mask"]
        _, ponds, _ = melt.detect(bands["green"], bands["nir"], reject,
                                  red=bands["red"])
        labels, n = ndimage.label(ponds)
        sizes = ndimage.sum(ponds, labels, range(1, n + 1))
        crs, _ = melt.tile_georeference(it)
        tr = melt.aoi_transform(it)
        big = np.argsort(sizes)[::-1][:MAX_PONDS]
        for k in big:
            if sizes[k] * melt.PIXEL_KM2 < MIN_POND_KM2:
                break
            rows, cols = np.where(labels == (k + 1))
            x0, y0 = tr * (cols.min() - 2, rows.max() + 2)
            x1, y1 = tr * (cols.max() + 2, rows.min() - 2)
            lons, lats = warp_transform(crs, "EPSG:4326",
                                        [x0, x1], [y0, y1])
            feats.append({
                "type": "Feature",
                "bbox": [min(lons), min(lats), max(lons), max(lats)],
                "properties": {
                    "pond_id": f"{tile}-{k + 1}",
                    "area_km2": round(float(sizes[k] * melt.PIXEL_KM2), 2),
                    "scene_date": it.datetime.date().isoformat()},
                "geometry": None})
    return {"type": "FeatureCollection", "features": feats}


def main(season):
    icesat2.OUT.mkdir(parents=True, exist_ok=True)
    gj = pond_polygons(season)
    p = icesat2.OUT / f"ponds_{season}.geojson"
    p.write_text(json.dumps(gj, indent=1))
    print(f"[{season}] {len(gj['features'])} big ponds", flush=True)

    import earthaccess
    earthaccess.login(strategy="netrc")
    rows = []
    from datetime import date
    for feat in gj["features"]:
        w, s, e, n = feat["bbox"]
        day = date.fromisoformat(feat["properties"]["scene_date"])
        try:
            grans = icesat2.search_atl03([w, s, e, n], day)[:MAX_GRANULES_PER_POND]
            if not grans:
                continue
            files = earthaccess.download(grans, str(icesat2.OUT / "atl03"))
        except Exception as ex:
            print(f"  [skip {feat['properties']['pond_id']}] {type(ex).__name__}", flush=True)
            continue
        for fp in files:
            fp = Path(fp)
            try:
                for _, la, lo, d in icesat2.granule_pond_depths(fp, [(w, s, e, n)]):
                    rows.append({"pond": feat["properties"]["pond_id"],
                                 "scene_date": feat["properties"]["scene_date"],
                                 "lat": la, "lon": lo, "depth_m": round(d, 3)})
            except Exception as ex:
                print(f"  [granule error] {type(ex).__name__}", flush=True)
            finally:
                fp.unlink(missing_ok=True)          # multi-GB files: never keep
        print(f"  {feat['properties']['pond_id']}: running total {len(rows)} depths", flush=True)

    out = icesat2.OUT / f"crossovers_{season}.json"
    out.write_text(json.dumps(rows, indent=1))
    print(f"[{season}] DONE: {len(rows)} crossover depths -> {out.name}", flush=True)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "2020-21")

"""M4 Task 3 runner: pond masks -> ICESat-2 crossover harvest (one season).

v2: harvests over the full clean-scene POOL (up to 3 scenes per tile — each
clean scene date is an independent ICESat-2 crossover window), and photons are
tested against the actual pond MASK, not just its bbox (a sprawling drainage
network's bbox is mostly crevassed ice; v1's 667 m "depth" came from exactly
that). ATL03 granules are multi-GB, so each is deleted right after extraction.

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
MAX_PONDS = 20
MAX_GRANULES_PER_POND = 3
MAX_POOL_SCENES = 3
WATER_TILES = ("19DEA", "19CEV")


def pond_targets(season):
    """Biggest ponds across every pool scene of the season: bbox features for
    the granule search plus (mask, transform, crs) for photon-level selection.
    The same pond on two dates is two legitimate crossover windows."""
    year = int(season.split("-")[0]) + 1
    feats, masks = [], []
    for tile in WATER_TILES:
        pool = shelf.production_pool(tile, f"{year}-01-01", f"{year}-02-28",
                                     max_scenes=MAX_POOL_SCENES)
        for it in pool:
            bands = melt.load_scene(it, bands=("green", "nir", "red"))
            reject = melt.reject_mask(it) | melt.cloud_mask(it)["mask"]
            _, ponds, _ = melt.detect(bands["green"], bands["nir"], reject,
                                      red=bands["red"])
            labels, n = ndimage.label(ponds)
            if not n:
                continue
            sizes = ndimage.sum(ponds, labels, range(1, n + 1))
            crs, _ = melt.tile_georeference(it)
            tr = melt.aoi_transform(it)
            scene_date = it.datetime.date().isoformat()
            big = np.argsort(sizes)[::-1][:MAX_PONDS]
            for k in big:
                if sizes[k] * melt.PIXEL_KM2 < MIN_POND_KM2:
                    break
                comp = labels == (k + 1)
                rows, cols = np.where(comp)
                r0, r1 = rows.min(), rows.max()
                c0, c1 = cols.min(), cols.max()
                x0, y0 = tr * (c0 - 2, r1 + 2)
                x1, y1 = tr * (c1 + 2, r0 - 2)
                lons, lats = warp_transform(crs, "EPSG:4326",
                                            [x0, x1], [y0, y1])
                # mask crop (+1 px dilation absorbs georef jitter) with its
                # own window transform, so photon tests stay cheap
                crop = ndimage.binary_dilation(comp[r0:r1 + 1, c0:c1 + 1])
                crop_tr = tr * tr.translation(c0, r0)
                feats.append({
                    "type": "Feature",
                    "bbox": [min(lons), min(lats), max(lons), max(lats)],
                    "properties": {
                        "pond_id": f"{tile}-{scene_date}-{k + 1}",
                        "area_km2": round(float(sizes[k] * melt.PIXEL_KM2), 2),
                        "scene_date": scene_date},
                    "geometry": None})
                masks.append((crop, crop_tr, crs))
    return {"type": "FeatureCollection", "features": feats}, masks


def main(season):
    icesat2.OUT.mkdir(parents=True, exist_ok=True)
    gj, masks = pond_targets(season)
    p = icesat2.OUT / f"ponds_{season}.geojson"
    p.write_text(json.dumps(gj, indent=1))
    print(f"[{season}] {len(gj['features'])} pond-date targets", flush=True)

    import earthaccess
    earthaccess.login(strategy="netrc")
    rows = []
    from datetime import date
    for feat, maskinfo in zip(gj["features"], masks):
        w, s, e, n = feat["bbox"]
        day = date.fromisoformat(feat["properties"]["scene_date"])
        try:
            grans = icesat2.search_atl03([w, s, e, n], day, pad_days=7)[:MAX_GRANULES_PER_POND]
            if not grans:
                continue
            files = earthaccess.download(grans, str(icesat2.OUT / "atl03"))
        except Exception as ex:
            print(f"  [skip {feat['properties']['pond_id']}] {type(ex).__name__}", flush=True)
            continue
        for fp in files:
            fp = Path(fp)
            try:
                for _, la, lo, d in icesat2.granule_pond_depths(
                        fp, [(w, s, e, n)], masks=[maskinfo]):
                    rows.append({"pond": feat["properties"]["pond_id"],
                                 "scene_date": feat["properties"]["scene_date"],
                                 "lat": la, "lon": lo, "depth_m": round(d, 3),
                                 "qc_pass": icesat2.qc_pass(d)})
            except Exception as ex:
                print(f"  [granule error] {type(ex).__name__}", flush=True)
            finally:
                fp.unlink(missing_ok=True)          # multi-GB files: never keep
        print(f"  {feat['properties']['pond_id']}: running total {len(rows)} depths", flush=True)

    out = icesat2.OUT / f"crossovers_{season}.json"
    out.write_text(json.dumps(rows, indent=1))
    nqc = sum(r["qc_pass"] for r in rows)
    print(f"[{season}] DONE: {len(rows)} depths ({nqc} QC-pass) -> {out.name}", flush=True)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "2020-21")

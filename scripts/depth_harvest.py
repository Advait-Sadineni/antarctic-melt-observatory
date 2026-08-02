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
MAX_GRANULES_PER_POND = 4
MAX_POOL_SCENES = 3
MAX_MASK_DIM = 2000          # decimate huge (Amery-scale) pond masks

# Harvest sites: GVI is the validated anchor; Amery's far larger lakes are
# added for crossover statistics (peak melt there is Dec-Jan, not Jan-Feb).
# Window months-days carry a year offset relative to the season's first year.
SITES = (
    {"shelf": None,                                   # default George VI state
     "tiles": ("19DEA", "19CEV"),
     "window": ((1, "01-01"), (1, "02-28"))},
    {"shelf": ("Amery", ["41CPV"], "reference/boundaries/Amery.geojson"),
     "tiles": ("41CPV",),
     "window": ((0, "12-01"), (1, "01-31"))},
)


def _site_pool_scenes(site, season):
    """(tile, item) pairs for a site's clean pool in its own melt window."""
    y0 = int(season.split("-")[0])
    (soff, smd), (eoff, emd) = site["window"]
    start, end = f"{y0 + soff}-{smd}", f"{y0 + eoff}-{emd}"
    if site["shelf"] is not None:
        name, tiles, bnd = site["shelf"]
        shelf.set_shelf(name, tiles, str(Path(bnd).resolve()))
    out = []
    for tile in site["tiles"]:
        for it in shelf.production_pool(tile, start, end,
                                        max_scenes=MAX_POOL_SCENES):
            out.append((tile, it))
    return out


def pond_targets(season):
    """Biggest ponds across every pool scene of the season at every site:
    bbox features for the granule search plus (mask, transform, crs) for
    photon-level selection. The same pond on two dates is two legitimate
    crossover windows."""
    feats, masks = [], []
    for site in SITES:
        pool = _site_pool_scenes(site, season)
        for tile, it in pool:
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
                step = max(1, int(np.ceil(max(crop.shape) / MAX_MASK_DIM)))
                if step > 1:      # Amery-scale lakes: bound mask memory
                    crop = crop[::step, ::step]
                    crop_tr = crop_tr * crop_tr.scale(step)
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
            grans = icesat2.search_atl03([w, s, e, n], day, pad_days=10)[:MAX_GRANULES_PER_POND]
            if not grans:
                continue
            files = earthaccess.download(grans, str(icesat2.OUT / "atl03"))
        except Exception as ex:
            print(f"  [skip {feat['properties']['pond_id']}] {type(ex).__name__}", flush=True)
            continue
        for fp in files:
            fp = Path(fp)
            try:
                # ATL03_YYYYMMDDHHMMSS_... -> laser acquisition date, so the
                # calibration can restrict to tight temporal crossovers
                gdate = date(int(fp.name[6:10]), int(fp.name[10:12]),
                             int(fp.name[12:14]))
                gap = abs((gdate - day).days)
                for _, la, lo, d in icesat2.granule_pond_depths(
                        fp, [(w, s, e, n)], masks=[maskinfo]):
                    rows.append({"pond": feat["properties"]["pond_id"],
                                 "scene_date": feat["properties"]["scene_date"],
                                 "laser_date": gdate.isoformat(),
                                 "gap_days": gap,
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

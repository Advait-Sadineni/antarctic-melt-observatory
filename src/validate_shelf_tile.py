"""
Blind reference-point validation for a single shelf tile, so a tile that
dominates the shelf-wide total (e.g. 19DEA) is not trusted on spectral
plausibility alone. Same method as validate_points/retune: stratified random
points, blind chips, a hidden key, area-weighted precision/recall.

Run:  python src/validate_shelf_tile.py 19DEA make    # build blind chips
      python src/validate_shelf_tile.py 19DEA score   # score labels vs key
"""

import csv
import json
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from rasterio.features import rasterize
from rasterio.warp import transform as warp_transform
from scipy import ndimage

import melt
import shelf

N_PER_STRATUM = {"A": 40, "B": 26, "C": 14}   # A=detected, B=near, C=far shelf
NEAR_M = 300
CHIP_PX = 61
GRID = (3, 4)
GAMMA = 0.55
SEED = 20260728
PEAK = ("01-05", "02-20")

OUT = melt.ROOT / "output" / "shelf_val"


def _setup(tile):
    it = shelf._peak_scene(tile, f"2021-{PEAK[0]}", f"2021-{PEAK[1]}")
    row, col, size = shelf._tile_window(it)
    melt.set_aoi(tile, row, col, size)
    bands = melt.load_scene(it, bands=("red", "green", "blue", "nir"))
    reject = melt.reject_mask(it) | melt.cloud_mask(it)["mask"]
    _, ponds, valid = melt.detect(bands["green"], bands["nir"], reject, red=bands["red"])

    # shelf polygon rasterised onto this tile's grid
    crs, _ = melt.tile_georeference(it)
    tr = melt.aoi_transform(it)
    geoms = []
    for g in shelf._boundary_3031():
        parts = []
        for poly in g["coordinates"]:
            rings = []
            for ring in poly:
                xs, ys = warp_transform(shelf.GRID_CRS, crs,
                                        [c[0] for c in ring], [c[1] for c in ring])
                rings.append(list(zip(xs, ys)))
            parts.append(rings)
        geoms.append({"type": "MultiPolygon", "coordinates": parts})
    shelfmask = rasterize([(g, 1) for g in geoms], out_shape=ponds.shape,
                          transform=tr, dtype="uint8").astype(bool)
    return it, bands, ponds, valid, shelfmask


def make(tile):
    it, bands, ponds, valid, shelfmask = _setup(tile)
    on = valid & shelfmask
    det = ponds & on
    near = ndimage.binary_dilation(det, iterations=int(NEAR_M / melt.PIXEL_M))
    strata = {"A": det, "B": near & ~det & on, "C": ~near & on}
    print(f"[{tile}] shelf {on.sum()*melt.PIXEL_KM2:.0f} km2, "
          f"detected {det.sum()*melt.PIXEL_KM2:.1f} km2 "
          f"({100*det.sum()/max(on.sum(),1):.2f}% of shelf-in-tile)")

    rng = np.random.default_rng(SEED)
    half = CHIP_PX // 2
    picks = []
    for k, n in N_PER_STRATUM.items():
        idx = np.argwhere(strata[k])
        inside = ((idx[:, 0] >= half) & (idx[:, 0] < ponds.shape[0] - half)
                  & (idx[:, 1] >= half) & (idx[:, 1] < ponds.shape[1] - half))
        idx = idx[inside]
        if len(idx) == 0:
            continue
        take = rng.choice(len(idx), size=min(n, len(idx)), replace=False)
        for r, c in idx[take]:
            picks.append({"stratum": k, "row": int(r), "col": int(c)})
    rng.shuffle(picks)
    for i, p in enumerate(picks, 1):
        p["id"] = i

    rgb = np.dstack([bands["red"], bands["green"], bands["blue"]])
    s = rgb[on]
    lo, hi = np.percentile(s[s > 0], [1, 99])
    disp = np.clip((rgb - lo) / max(hi - lo, 1e-6), 0, 1) ** GAMMA

    d = OUT / tile
    d.mkdir(parents=True, exist_ok=True)
    per = GRID[0] * GRID[1]
    for f0 in range(0, len(picks), per):
        chunk = picks[f0:f0 + per]
        fig, axes = plt.subplots(*GRID, figsize=(GRID[1] * 4.0, GRID[0] * 4.3))
        for ax, p in zip(axes.ravel(), chunk):
            r, c = p["row"], p["col"]
            ax.imshow(disp[r-half:r+half+1, c-half:c+half+1], interpolation="nearest")
            for a, b in ((-12, -4), (4, 12)):
                ax.plot([half, half], [half+a, half+b], color="#ff2d55", lw=1.3)
                ax.plot([half+a, half+b], [half, half], color="#ff2d55", lw=1.3)
            ax.set_title(f"#{p['id']}", fontsize=11)
            ax.set_xticks([]); ax.set_yticks([])
        for ax in axes.ravel()[len(chunk):]:
            ax.axis("off")
        fig.suptitle(f"{tile} — label the centre pixel: water / not / unsure", fontsize=11)
        fig.tight_layout()
        fig.savefig(d / f"chips_{f0//per+1:02d}.png", dpi=110, bbox_inches="tight")
        plt.close(fig)

    (d / "answer_key.json").write_text(json.dumps(
        {"scene": it.id, "tile": tile,
         "stratum_px": {k: int(m.sum()) for k, m in strata.items()},
         "shelf_px": int(on.sum()), "picks": picks}, indent=1))
    lp = d / "labels.csv"
    if not lp.exists():
        with open(lp, "w", newline="") as f:
            w = csv.writer(f); w.writerow(["id", "label"])
            for p in picks:
                w.writerow([p["id"], ""])
    print(f"[{tile}] wrote {len(picks)} chips to {d.relative_to(melt.ROOT)}")


def score(tile):
    d = OUT / tile
    key = json.loads((d / "answer_key.json").read_text())
    truth = {p["id"]: p for p in key["picks"]}
    labels = {int(r["id"]): r["label"].strip().lower()
              for r in csv.DictReader(open(d / "labels.csv")) if r["label"].strip()}
    if not labels:
        print("no labels yet"); return

    weights = {k: v / key["shelf_px"] for k, v in key["stratum_px"].items()}
    cnt = {k: {"n": 0, "w": 0} for k in key["stratum_px"]}
    for i, lab in labels.items():
        if lab not in ("water", "not"):
            continue
        s = truth[i]["stratum"]
        cnt[s]["n"] += 1
        cnt[s]["w"] += (lab == "water")

    print(f"[{tile}] {len(labels)} labelled, scene {key['scene']}\n")
    print(f"  {'stratum':8s} {'shelf share':>12} {'labelled':>9} {'water':>6} {'rate':>7}")
    for k in ("A", "B", "C"):
        c = cnt[k]
        rate = c["w"] / c["n"] if c["n"] else float("nan")
        print(f"  {k:8s} {weights[k]*100:11.2f}% {c['n']:9d} {c['w']:6d} {rate:7.2f}")

    a = cnt["A"]
    prec = a["w"] / a["n"] if a["n"] else float("nan")
    from validate_points import wilson
    lo, hi = wilson(a["w"], a["n"])
    dens = {k: (cnt[k]["w"]/cnt[k]["n"]) if cnt[k]["n"] else 0.0 for k in weights}
    true_frac = sum(weights[k]*dens[k] for k in weights)
    rec = (weights["A"]*dens["A"]) / true_frac if true_frac else float("nan")
    print(f"\n  precision (detected streaks that are water): {prec:.3f}  95% CI [{lo:.3f},{hi:.3f}]")
    print(f"  recall (area-weighted):                      {rec:.3f}")
    print(f"\n  VERDICT: {'detections are real meltwater' if prec>=0.7 else 'many detections are NOT water - artifact'}")


if __name__ == "__main__":
    tile = sys.argv[1] if len(sys.argv) > 1 else "19DEA"
    cmd = sys.argv[2] if len(sys.argv) > 2 else "make"
    {"make": make, "score": score}[cmd](tile)

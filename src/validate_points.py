"""
Reference-point validation: does the mask agree with what a human sees?

Cross-sensor validation shows Sentinel-2 and Landsat measure the same thing.
It cannot show that either one is *right* - two instruments can agree on a
wrong answer. This compares the detector against visual interpretation of the
true-colour imagery, which is the closest thing to ground truth available for
a place nobody is standing in.

Method: stratified random point sampling, following standard remote-sensing
accuracy practice rather than hand-drawn polygons. Delineating hundreds of
small ponds by hand is slow and its errors are hard to quantify; a random
sample gives unbiased estimates with real confidence intervals.

Meltwater covers ~0.3% of the study area, so a simple random sample would be
almost entirely blank snow and would pin down precision while saying nothing
useful about recall. Three strata fix that:

  A  detected as meltwater                     - estimates precision
  B  not detected, but within 200 m of a       - where misses concentrate
     detection
  C  everything else                           - confirms misses are rare
                                                 in the vast empty majority

Estimates are area-weighted back to the whole study area, so the rarity of
strata A and B is accounted for rather than ignored.

The labelling is blind. `make` writes shuffled, anonymously numbered chips
plus a separate answer key that is not consulted until `score` runs. Nothing
in a chip reveals its stratum or the detector's opinion of it.

Run:  python src/validate_points.py make     # build chips + hidden key
      python src/validate_points.py score    # score labels against the key
"""

import csv
import json
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage

import melt

SCENE = "S2B_19CEV_20210124_0_L2A"  # visually confirmed clean, peak melt

N_PER_STRATUM = {"A": 40, "B": 32, "C": 24}
NEAR_DISTANCE_M = 200
CHIP_PX = 61  # ~610 m across at 10 m; odd so there is a true centre pixel
GRID = (3, 4)  # chips per figure

# Snow is so bright that a plain percentile stretch crushes ponds and crevasse
# shadows into the same black. Gamma lifts the dark end so the two can be told
# apart by shape and tone, which is the whole point of looking.
DISPLAY_GAMMA = 0.55
SEED = 20260722

OUT_DIR = melt.ROOT / "output" / "pointval"
KEY_PATH = OUT_DIR / "answer_key.json"  # not read until scoring
LABELS_PATH = OUT_DIR / "labels.csv"  # filled in by the interpreter
RESULT_PATH = melt.ROOT / "output" / "validation_points.csv"


def build_strata(ponds, valid):
    near = ndimage.binary_dilation(
        ponds, iterations=int(round(NEAR_DISTANCE_M / melt.PIXEL_M))
    )
    return {
        "A": ponds & valid,
        "B": near & ~ponds & valid,
        "C": ~near & valid,
    }


def make():
    item = melt.get_item(SCENE)
    bands = melt.load_scene(item, use_cache=False)
    reject = melt.reject_mask(item) | melt.cloud_mask(item)["mask"]
    _, ponds, valid = melt.detect(bands["green"], bands["nir"], reject)
    valid = valid & ~reject

    strata = build_strata(ponds, valid)
    total_valid = int(valid.sum())
    print(f"[strata] valid {total_valid:,} px "
          f"({total_valid * melt.PIXEL_KM2:.0f} km2)")
    for k, m in strata.items():
        print(f"   {k}: {int(m.sum()):>10,} px  {m.sum()/total_valid*100:6.3f}% of valid")

    rng = np.random.default_rng(SEED)
    picks = []
    half = CHIP_PX // 2
    for k, n in N_PER_STRATUM.items():
        idx = np.argwhere(strata[k])
        # keep chips fully inside the array
        inside = ((idx[:, 0] >= half) & (idx[:, 0] < melt.WIN_SIZE - half)
                  & (idx[:, 1] >= half) & (idx[:, 1] < melt.WIN_SIZE - half))
        idx = idx[inside]
        take = rng.choice(len(idx), size=min(n, len(idx)), replace=False)
        for r, c in idx[take]:
            picks.append({"stratum": k, "row": int(r), "col": int(c)})

    rng.shuffle(picks)
    for i, p in enumerate(picks, 1):
        p["id"] = i

    # A single stretch for every chip, so brightness is comparable and I am
    # not misled into calling a dark-looking chip water because it was
    # stretched differently from its neighbours.
    rgb = np.dstack([bands["red"], bands["green"], bands["blue"]])
    sample = rgb[valid]
    lo, hi = np.percentile(sample[sample > 0], [1, 99])
    rgb_disp = np.clip((rgb - lo) / max(hi - lo, 1e-6), 0, 1) ** DISPLAY_GAMMA

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    per_fig = GRID[0] * GRID[1]
    for f0 in range(0, len(picks), per_fig):
        chunk = picks[f0:f0 + per_fig]
        fig, axes = plt.subplots(*GRID, figsize=(GRID[1] * 4.0, GRID[0] * 4.3))
        for ax, p in zip(axes.ravel(), chunk):
            r, c = p["row"], p["col"]
            ax.imshow(rgb_disp[r - half:r + half + 1, c - half:c + half + 1],
                      interpolation="nearest")
            # crosshair marks the pixel being judged, without covering it
            ax.plot([half, half], [half - 12, half - 4], color="#ff2d55", lw=1.3)
            ax.plot([half, half], [half + 4, half + 12], color="#ff2d55", lw=1.3)
            ax.plot([half - 12, half - 4], [half, half], color="#ff2d55", lw=1.3)
            ax.plot([half + 4, half + 12], [half, half], color="#ff2d55", lw=1.3)
            ax.set_title(f"#{p['id']}", fontsize=11)
            ax.set_xticks([])
            ax.set_yticks([])
        for ax in axes.ravel()[len(chunk):]:
            ax.axis("off")
        fig.suptitle("Label the centre pixel only:  water / not water",
                     fontsize=11)
        fig.tight_layout()
        page = f0 // per_fig + 1
        fig.savefig(OUT_DIR / f"chips_{page:02d}.png", dpi=110,
                    bbox_inches="tight")
        plt.close(fig)
        print(f"[write] {(OUT_DIR / f'chips_{page:02d}.png').relative_to(melt.ROOT)}"
              f"  ({len(chunk)} chips)")

    key = {
        "scene": SCENE,
        "seed": SEED,
        "stratum_px": {k: int(m.sum()) for k, m in strata.items()},
        "valid_px": total_valid,
        "picks": picks,
    }
    KEY_PATH.write_text(json.dumps(key, indent=1))
    print(f"[write] {KEY_PATH.relative_to(melt.ROOT)}  (do not read before labelling)")

    if not LABELS_PATH.exists():
        with open(LABELS_PATH, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["id", "label"])
            for p in picks:
                w.writerow([p["id"], ""])
        print(f"[write] {LABELS_PATH.relative_to(melt.ROOT)}  (fill in water/not)")


def wilson(k, n, z=1.96):
    """Wilson score interval - behaves sensibly at small n and near 0 or 1."""
    if n == 0:
        return float("nan"), float("nan")
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return max(c - h, 0.0), min(c + h, 1.0)


def score():
    key = json.loads(KEY_PATH.read_text())
    truth = {p["id"]: p for p in key["picks"]}

    with open(LABELS_PATH, newline="") as f:
        labels = {int(r["id"]): r["label"].strip().lower()
                  for r in csv.DictReader(f) if r["label"].strip()}
    if not labels:
        print("No labels yet. Look at output/pointval/chips_*.png and fill in "
              f"{LABELS_PATH.relative_to(melt.ROOT)}.")
        return

    valid_px = key["valid_px"]
    weights = {k: v / valid_px for k, v in key["stratum_px"].items()}

    counts = {k: {"n": 0, "water": 0} for k in key["stratum_px"]}
    for i, lab in labels.items():
        if lab not in ("water", "not"):
            continue
        s = truth[i]["stratum"]
        counts[s]["n"] += 1
        counts[s]["water"] += (lab == "water")

    print(f"scene {key['scene']}   {len(labels)} labelled points\n")
    print(f"  {'stratum':8s} {'area share':>11} {'labelled':>9} {'water':>7} {'rate':>8}")
    for k in ("A", "B", "C"):
        c = counts[k]
        rate = c["water"] / c["n"] if c["n"] else float("nan")
        print(f"  {k:8s} {weights[k]*100:10.3f}% {c['n']:9d} {c['water']:7d} {rate:8.3f}")

    # Precision: of pixels the detector called water, how many are water.
    a = counts["A"]
    prec = a["water"] / a["n"] if a["n"] else float("nan")
    plo, phi = wilson(a["water"], a["n"])

    # Recall: detected water as a share of all water, area-weighted so the
    # rarity of each stratum is respected.
    water_density = {k: (counts[k]["water"] / counts[k]["n"]) if counts[k]["n"] else 0.0
                     for k in weights}
    total_water = sum(weights[k] * water_density[k] for k in weights)
    detected_water = weights["A"] * water_density["A"]
    recall = detected_water / total_water if total_water else float("nan")

    f1 = 2 * prec * recall / (prec + recall) if (prec + recall) else float("nan")

    det_frac = weights["A"]
    true_frac = total_water
    aoi_km2 = valid_px * melt.PIXEL_KM2

    print(f"\n  precision            {prec:.3f}   95% CI [{plo:.3f}, {phi:.3f}]")
    print(f"  recall (area-wtd)    {recall:.3f}")
    print(f"  F1                   {f1:.3f}")
    print(f"\n  detector says water  {det_frac*100:.3f}% of AOI "
          f"= {det_frac*aoi_km2:.2f} km2")
    print(f"  sample says water    {true_frac*100:.3f}% of AOI "
          f"= {true_frac*aoi_km2:.2f} km2")
    if true_frac:
        print(f"  bias                 {det_frac/true_frac:.2f}x")

    # What the sample says about the threshold itself. This is the first
    # evidence about NDWI 0.16 that does not come from the detector's own
    # output, so it is worth printing even though 96 points on one scene is
    # too thin a basis for actually retuning.
    item = melt.get_item(key["scene"])
    bands = melt.load_scene(item, bands=("green", "nir"))
    reject = melt.reject_mask(item) | melt.cloud_mask(item)["mask"]
    ndwi, _, _ = melt.detect(bands["green"], bands["nir"], reject)

    pts = []
    for i, lab in labels.items():
        if lab in ("water", "not"):
            p = truth[i]
            pts.append((p["stratum"], lab, float(ndwi[p["row"], p["col"]])))

    print("\n  NDWI by visual label")
    for lab in ("water", "not"):
        v = np.array([p[2] for p in pts if p[1] == lab])
        if len(v):
            print(f"    {lab:6s} n={len(v):3d}  p25 {np.percentile(v, 25):.3f}  "
                  f"median {np.median(v):.3f}  p75 {np.percentile(v, 75):.3f}  "
                  f"max {v.max():.3f}")

    print(f"\n  {'thresh':>7} {'precision':>10} {'recall':>8} {'F1':>7}")
    sweep = []
    for t in (0.10, 0.14, 0.16, 0.18, 0.20, 0.22, 0.25, 0.30, 0.40):
        tp = sum(1 for s, l, v in pts if v > t and l == "water")
        fp = sum(1 for s, l, v in pts if v > t and l == "not")
        p_ = tp / (tp + fp) if (tp + fp) else float("nan")
        true_f = sum(weights[s] * np.mean([1.0 if l == "water" else 0.0
                                           for ss, l, _ in pts if ss == s])
                     for s in weights if any(ss == s for ss, _, _ in pts))
        det_f = sum(weights[s] * np.mean([1.0 if (l == "water" and v > t) else 0.0
                                          for ss, l, v in pts if ss == s])
                    for s in weights if any(ss == s for ss, _, _ in pts))
        r_ = det_f / true_f if true_f else float("nan")
        f_ = 2 * p_ * r_ / (p_ + r_) if (p_ and r_ and p_ + r_) else float("nan")
        mark = "  <-- current" if abs(t - melt.NDWI_THRESHOLD) < 1e-9 else ""
        print(f"  {t:7.2f} {p_:10.3f} {r_:8.3f} {f_:7.3f}{mark}")
        sweep.append((t, p_, r_, f_))

    with open(melt.ROOT / "output" / "validation_points_threshold.csv", "w",
              newline="") as f:
        w = csv.writer(f)
        w.writerow(["ndwi_threshold", "precision", "recall", "f1"])
        for t, p_, r_, f_ in sweep:
            w.writerow([t, round(p_, 4), round(r_, 4), round(f_, 4)])

    RESULT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(RESULT_PATH, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["metric", "value"])
        for name, v in (("scene", key["scene"]), ("n_labelled", len(labels)),
                        ("precision", round(prec, 4)),
                        ("precision_ci_lo", round(plo, 4)),
                        ("precision_ci_hi", round(phi, 4)),
                        ("recall_area_weighted", round(recall, 4)),
                        ("f1", round(f1, 4)),
                        ("detected_frac_of_aoi", round(det_frac, 6)),
                        ("sampled_true_frac_of_aoi", round(true_frac, 6)),
                        ("area_bias", round(det_frac / true_frac, 4) if true_frac else "")):
            w.writerow([name, v])
    print(f"\n[write] {RESULT_PATH.relative_to(melt.ROOT)}")


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "make"
    {"make": make, "score": score}[cmd]()

"""
Retune the detector against blind labels from several scenes, using the
methods the published literature actually uses for Antarctic supraglacial
lakes rather than a single hand-tuned threshold.

Reference-point validation on one scene (validate_points.py) showed precision
0.63 / recall 0.36, with every false positive being crevasse shadow at NDWI
0.17-0.25 while real ponds sit at 0.38+. The literature agrees the threshold
is too low and names the fixes:

  - Moussavi et al. 2020 (Remote Sensing 12:134): NDWI (green-NIR) > ~0.18-0.19
    for Sentinel-2, plus a shadow test (green - red) > 0.09 in reflectance to
    reject shaded snow and crevasse shadow, which is spectrally shadow.
  - Corr et al. 2022 (ESSD 14:209): a dual index - NDWI_GNIR > 0.16 AND
    NDWI_BR = (blue-red)/(blue+red) > 0.18 - reaching 85% recall on Sentinel-2.

This module labels points across several scenes, then scores each candidate
configuration against the pooled labels with Horvitz-Thompson weights, so the
comparison is unbiased despite the stratified sampling. The labels are pixel
truth and are reused across every configuration - only the decision rule
changes.

Run:  python src/retune.py chips     # build blind chips for the new scenes
      python src/retune.py eval      # score every candidate config vs labels
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

# Scenes spanning seasons and conditions. 2021-01-24 reuses the labels already
# made in output/pointval; the other two are labelled fresh here. 2020-01-19
# is the literature-verified peak day of the record 2019-20 season; 2026-02-10
# is the current record season.
SCENES = [
    ("2021-01-24", "S2B_19CEV_20210124_0_L2A", "reuse"),
    ("2020-01-19", "S2B_19CEV_20200119_0_L2A", "new"),
    ("2026-02-10", "S2B_19CEV_20260210_0_L2A", "new"),
]

N_PER_STRATUM = {"A": 30, "B": 22, "C": 14}  # per new scene
NEAR_DISTANCE_M = 200
CHIP_PX = 61
GRID = (3, 4)
DISPLAY_GAMMA = 0.55
SEED = 20260725

S2_SCALE = 1.0 / 10000.0
FLOOR_REFL = melt.BRIGHTNESS_FLOOR * S2_SCALE

RETUNE_DIR = melt.ROOT / "output" / "retune"
POINTVAL_DIR = melt.ROOT / "output" / "pointval"  # existing 2021-01-24 labels
OUT_CSV = melt.ROOT / "output" / "retune_configs.csv"
OUT_PNG = melt.ROOT / "output" / "retune_configs.png"


# --- candidate detection rules, all in reflectance ---------------------------

def _indices(g, n, r, b):
    with np.errstate(invalid="ignore", divide="ignore"):
        ndwi = (g - n) / np.maximum(g + n, 1e-6)
        ndwi_br = (b - r) / np.maximum(b + r, 1e-6)
    return ndwi, ndwi_br


def config_masks(g, n, r, b):
    """Every candidate decision rule, as boolean arrays over the same pixels.

    Point-sampled configs are pixel-wise. The final adopted method also grows
    cores into connected margins (hysteresis), which is spatial - it cannot be
    evaluated from isolated points here, so it is measured on full masks in
    the hysteresis probe and reported separately. This function still lets the
    pixel-wise candidates be compared apples-to-apples.
    """
    from scipy import ndimage

    ndwi, ndwi_br = _indices(g, n, r, b)
    bright = g > FLOOR_REFL
    shadow_ok = (g - r) > 0.09  # Moussavi shadow / crevasse test
    br_ok = ndwi_br > 0.18      # Corr blue-red index

    moussavi = (ndwi > 0.19) & bright & shadow_ok
    ok = bright & shadow_ok
    grow = ok & (ndwi > 0.14)
    labels, ncc = ndimage.label(grow)
    keep = np.zeros(ncc + 1, bool)
    keep[np.unique(labels[moussavi])] = True
    keep[0] = False
    hysteresis = keep[labels]

    return {
        "current (NDWI>0.16)": (ndwi > 0.16) & bright,
        "raise 0.19 (Moussavi)": (ndwi > 0.19) & bright,
        "0.16 + shadow test": (ndwi > 0.16) & bright & shadow_ok,
        "dual-NDWI (Corr)": (ndwi > 0.16) & bright & br_ok,
        "dual + shadow": (ndwi > 0.16) & bright & br_ok & shadow_ok,
        "0.19 + shadow (Moussavi full)": moussavi,
        "0.19 + shadow + hysteresis (adopted)": hysteresis,
    }


# --- labelling ---------------------------------------------------------------

def build_strata(ponds, valid):
    near = ndimage.binary_dilation(ponds, iterations=int(round(NEAR_DISTANCE_M / melt.PIXEL_M)))
    return {"A": ponds & valid, "B": near & ~ponds & valid, "C": ~near & valid}


def make_chips():
    for i, (tag, item_id, mode) in enumerate(SCENES):
        if mode != "new":
            print(f"[{tag}] reusing existing labels in output/pointval")
            continue
        _make_one(tag, item_id, seed=SEED + i)


def _make_one(tag, item_id, seed):
    item = melt.get_item(item_id)
    bands = melt.load_scene(item, use_cache=False)
    reject = melt.reject_mask(item) | melt.cloud_mask(item)["mask"]
    _, ponds, valid = melt.detect(bands["green"], bands["nir"], reject)
    valid = valid & ~reject
    strata = build_strata(ponds, valid)

    print(f"[{tag}] valid {int(valid.sum()):,} px; "
          f"A {int(strata['A'].sum()):,}  B {int(strata['B'].sum()):,}  C {int(strata['C'].sum()):,}")

    rng = np.random.default_rng(seed)
    half = CHIP_PX // 2
    picks = []
    for k, n in N_PER_STRATUM.items():
        idx = np.argwhere(strata[k])
        inside = ((idx[:, 0] >= half) & (idx[:, 0] < melt.WIN_SIZE - half)
                  & (idx[:, 1] >= half) & (idx[:, 1] < melt.WIN_SIZE - half))
        idx = idx[inside]
        take = rng.choice(len(idx), size=min(n, len(idx)), replace=False)
        for r, c in idx[take]:
            picks.append({"stratum": k, "row": int(r), "col": int(c)})
    rng.shuffle(picks)
    for i, p in enumerate(picks, 1):
        p["id"] = i

    rgb = np.dstack([bands["red"], bands["green"], bands["blue"]])
    s = rgb[valid]
    lo, hi = np.percentile(s[s > 0], [1, 99])
    disp = np.clip((rgb - lo) / max(hi - lo, 1e-6), 0, 1) ** DISPLAY_GAMMA

    scene_dir = RETUNE_DIR / tag
    scene_dir.mkdir(parents=True, exist_ok=True)
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
        fig.suptitle(f"{tag}  —  label the centre pixel:  water / not / unsure", fontsize=11)
        fig.tight_layout()
        fig.savefig(scene_dir / f"chips_{f0//per+1:02d}.png", dpi=110, bbox_inches="tight")
        plt.close(fig)

    (scene_dir / "answer_key.json").write_text(json.dumps(
        {"scene": item_id, "tag": tag, "seed": seed,
         "stratum_px": {k: int(m.sum()) for k, m in strata.items()},
         "valid_px": int(valid.sum()), "picks": picks}, indent=1))
    lp = scene_dir / "labels.csv"
    if not lp.exists():
        with open(lp, "w", newline="") as f:
            w = csv.writer(f); w.writerow(["id", "label"])
            for p in picks:
                w.writerow([p["id"], ""])
    print(f"[{tag}] wrote {len(picks)} chips to {scene_dir.relative_to(melt.ROOT)}")


# --- evaluation --------------------------------------------------------------

def _load_scene_points(tag, item_id, key_path, labels_path):
    key = json.loads(key_path.read_text())
    truth = {p["id"]: p for p in key["picks"]}
    with open(labels_path, newline="") as f:
        labels = {int(r["id"]): r["label"].strip().lower()
                  for r in csv.DictReader(f) if r["label"].strip()}

    item = melt.get_item(item_id)
    bands = melt.load_scene(item, bands=("green", "nir", "red", "blue"))
    g = bands["green"] * S2_SCALE
    n = bands["nir"] * S2_SCALE
    r = bands["red"] * S2_SCALE
    b = bands["blue"] * S2_SCALE
    masks = config_masks(g, n, r, b)

    weights = {k: v / key["valid_px"] for k, v in key["stratum_px"].items()}
    n_sampled = {k: sum(1 for i, l in labels.items()
                        if truth[i]["stratum"] == k and l in ("water", "not"))
                 for k in key["stratum_px"]}

    pts = []
    for i, lab in labels.items():
        if lab not in ("water", "not"):
            continue
        p = truth[i]
        s = p["stratum"]
        # Horvitz-Thompson weight: stratum area share / points drawn there.
        w = weights[s] / max(n_sampled[s], 1)
        decisions = {name: bool(m[p["row"], p["col"]]) for name, m in masks.items()}
        pts.append({"w": w, "water": lab == "water", "cfg": decisions})
    return pts, sum(1 for l in labels.values() if l == "unsure")


def evaluate():
    all_pts, unsure = [], 0
    for tag, item_id, mode in SCENES:
        d = (POINTVAL_DIR if mode == "reuse" else RETUNE_DIR / tag)
        kp, lp = d / "answer_key.json", d / "labels.csv"
        if not lp.exists() or not any(csv.DictReader(open(lp)).__iter__()):
            print(f"[{tag}] no labels yet at {lp.relative_to(melt.ROOT)}")
            continue
        pts, u = _load_scene_points(tag, item_id, kp, lp)
        all_pts += pts
        unsure += u
        print(f"[{tag}] {len(pts)} labelled points ({u} unsure excluded)")

    if not all_pts:
        print("nothing labelled yet")
        return

    names = list(all_pts[0]["cfg"].keys())
    rows = []
    for name in names:
        tp = sum(p["w"] for p in all_pts if p["cfg"][name] and p["water"])
        fp = sum(p["w"] for p in all_pts if p["cfg"][name] and not p["water"])
        fn = sum(p["w"] for p in all_pts if not p["cfg"][name] and p["water"])
        prec = tp / (tp + fp) if (tp + fp) else float("nan")
        rec = tp / (tp + fn) if (tp + fn) else float("nan")
        f1 = 2 * prec * rec / (prec + rec) if prec and rec and (prec + rec) else float("nan")
        rows.append({"config": name, "precision": prec, "recall": rec, "f1": f1})

    rows.sort(key=lambda r: (-r["f1"] if r["f1"] == r["f1"] else 1))
    print(f"\n  pooled points: {len(all_pts)} across {sum(1 for _,_,m in SCENES)} scenes\n")
    print(f"  {'config':32s} {'prec':>6} {'recall':>7} {'F1':>6}")
    for r in rows:
        star = "  *" if r is rows[0] else ""
        print(f"  {r['config']:32s} {r['precision']:6.3f} {r['recall']:7.3f} {r['f1']:6.3f}{star}")

    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["config", "precision", "recall", "f1"])
        w.writeheader()
        for r in rows:
            w.writerow({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()})

    _plot(rows)
    print(f"\n  best by F1: {rows[0]['config']}")
    print(f"[write] {OUT_CSV.relative_to(melt.ROOT)}")


def _plot(rows):
    rows = [r for r in rows if r["f1"] == r["f1"]]
    fig, ax = plt.subplots(figsize=(9, 5.5))
    for r in rows:
        cur = r["config"].startswith("current")
        ax.scatter(r["recall"], r["precision"], s=90,
                   color="#c1121f" if cur else "#0077b6",
                   zorder=3, edgecolor="white")
        ax.annotate(r["config"], (r["recall"], r["precision"]),
                    textcoords="offset points", xytext=(7, 3), fontsize=7.5,
                    color="#c1121f" if cur else "#333")
    ax.set_xlabel("recall (area-weighted)")
    ax.set_ylabel("precision")
    ax.set_title("Detection configurations vs blind labels\n"
                 "literature methods (shadow test, dual-NDWI) against the current single threshold",
                 fontsize=11)
    ax.grid(alpha=0.3)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.02)
    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=130, bbox_inches="tight")
    print(f"[write] {OUT_PNG.relative_to(melt.ROOT)}")


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "eval"
    {"chips": make_chips, "eval": evaluate}[cmd]()

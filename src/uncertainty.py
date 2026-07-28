"""
Turn each season's peak into a defensible range instead of a single number.

Every headline figure this project reports is a point estimate resting on
choices that were made, not measured. Validation has now quantified the two
that dominate, so the honest form of the result is an interval.

Sources of uncertainty, in order of size:

1. THRESHOLD CHOICE. Local sensitivity at NDWI 0.16 is ~11% of area per 0.01,
   and the reference-point sample shows the classes actually separate around
   0.20-0.25 rather than 0.16. The band here spans 0.14 to 0.25 - from just
   below the current value to the point where sampled precision reaches 1.00.
   This is computed directly by re-thresholding each peak scene, not modelled.

2. PARTIAL COVERAGE. Scenes are kept when at least half the study area
   survives screening, and the density is scaled back to the full area. That
   extrapolation assumes the obscured part behaves like the visible part. The
   penalty grows as coverage falls and is reported per scene.

3. DETECTOR BIAS. Blind reference-point sampling on one scene gives precision
   0.63 and area bias 0.57x - the detector under-reports overall, because
   missed pond margins outweigh false crevasse-shadow detections. One scene is
   too thin to correct with, so this is reported alongside rather than applied.

The output is deliberately conservative: a range, its width, and the reasons.

Run:  python src/uncertainty.py
"""

import csv
import glob

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import melt
import season as S

# Just below the current value, to where sampled precision reaches 1.00.
THRESHOLD_BAND = (0.14, 0.25)
THRESHOLDS = (0.14, 0.16, 0.20, 0.25)

# From validate_points.py on the 2021-01-24 scene.
MEASURED_AREA_BIAS = 0.57
MEASURED_PRECISION = 0.63

OUT_CSV = melt.ROOT / "output" / "uncertainty_peaks.csv"
OUT_PNG = melt.ROOT / "output" / "uncertainty_peaks.png"


def season_peaks():
    """The peak scene of each season, from the committed season CSVs."""
    peaks = []
    for path in sorted(glob.glob(str(melt.ROOT / "output" / "season_*.csv"))):
        with open(path, newline="") as f:
            rows = [r for r in csv.DictReader(f)]
        if not rows:
            continue
        label = path.rsplit("season_", 1)[1][:7].replace("_", "-")
        best = max(rows, key=lambda r: float(r["pond_km2_equiv"]))
        peaks.append((label, best))
    return peaks


def scene_areas(item_id, thresholds=THRESHOLDS):
    """Meltwater area at several thresholds, from one pass over the imagery."""
    item = melt.get_item(item_id)
    bands = melt.load_scene(item, bands=("green", "nir"))
    reject = melt.reject_mask(item) | melt.cloud_mask(item)["mask"]

    green, nir = bands["green"], bands["nir"]
    valid = (green > 0) & (nir > 0) & ~reject
    with np.errstate(invalid="ignore", divide="ignore"):
        ndwi = np.where(valid, (green - nir) / np.maximum(green + nir, 1e-6), np.nan)
    bright = green > melt.BRIGHTNESS_FLOOR

    usable = int(valid.sum())
    aoi_px = melt.WIN_SIZE ** 2
    out = {}
    for t in thresholds:
        n = int(((np.nan_to_num(ndwi, nan=-9) > t) & valid & bright).sum())
        # density over usable ground, scaled to the whole study area
        out[t] = n / max(usable, 1) * aoi_px * melt.PIXEL_KM2
    return out, usable / aoi_px


def main():
    peaks = season_peaks()
    print(f"[peaks] {len(peaks)} seasons\n")

    rows = []
    for label, r in peaks:
        try:
            areas, cover = scene_areas(r["item_id"])
        except Exception as e:
            print(f"  {label}  ERROR {type(e).__name__}: {str(e)[:60]}")
            continue

        lo = min(areas[t] for t in THRESHOLDS)
        hi = max(areas[t] for t in THRESHOLDS)
        mid = areas[melt.NDWI_THRESHOLD]
        below_floor = hi <= melt.NOISE_FLOOR_KM2

        rows.append({
            "season": label,
            "date": r["date"],
            "item_id": r["item_id"],
            "coverage": round(cover, 4),
            "km2_at_0.14": round(areas[0.14], 3),
            "km2_at_0.16": round(areas[0.16], 3),
            "km2_at_0.20": round(areas[0.20], 3),
            "km2_at_0.25": round(areas[0.25], 3),
            "range_lo": round(lo, 3),
            "range_hi": round(hi, 3),
            "range_width_pct": round((hi - lo) / max(mid, 1e-9) * 100, 1),
            "below_noise_floor": below_floor,
        })
        flag = "  (entirely below noise floor)" if below_floor else ""
        print(f"  {label}  {r['date']}  coverage {cover*100:5.1f}%   "
              f"{lo:6.2f} - {hi:6.2f} km2   "
              f"(point estimate {mid:6.2f}, range is {(hi-lo)/max(mid,1e-9)*100:5.0f}% of it){flag}")

    if not rows:
        print("nothing to report")
        return

    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    widths = [r["range_width_pct"] for r in rows]
    print(f"\n  median range width: {np.median(widths):.0f}% of the point estimate")
    print(f"  threshold band:     {THRESHOLD_BAND[0]} to {THRESHOLD_BAND[1]}")
    print(f"  measured area bias: {MEASURED_AREA_BIAS:.2f}x "
          f"(detector under-reports; not applied)")

    plot(rows)


def plot(rows):
    labels = [r["season"] for r in rows]
    lo = np.array([r["range_lo"] for r in rows])
    hi = np.array([r["range_hi"] for r in rows])
    mid = np.array([r["km2_at_0.16"] for r in rows])
    cover = np.array([r["coverage"] for r in rows])
    x = np.arange(len(rows))

    fig, ax = plt.subplots(figsize=(11.5, 6.2))
    ax.vlines(x, lo, hi, color="#0077b6", lw=9, alpha=0.32,
              label=f"NDWI {THRESHOLD_BAND[0]}-{THRESHOLD_BAND[1]}")
    ax.plot(x, mid, "o", color="#0077b6", ms=8, zorder=3,
            label=f"point estimate (NDWI {melt.NDWI_THRESHOLD})")

    for i, r in enumerate(rows):
        if r["coverage"] < 0.85:
            ax.annotate(f"{r['coverage']*100:.0f}% cover",
                        (x[i], hi[i]), textcoords="offset points",
                        xytext=(0, 8), ha="center", fontsize=7.5, color="#a15c00")

    ax.axhline(melt.NOISE_FLOOR_KM2, color="#c1121f", ls="--", lw=1.2,
               label=f"noise floor {melt.NOISE_FLOOR_KM2} km$^2$")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_ylabel("Peak meltwater area (km$^2$)")
    ax.set_title("George VI peak meltwater, with threshold uncertainty\n"
                 "range spans the plausible NDWI threshold; "
                 "detector separately under-reports by ~0.57x",
                 fontsize=12)
    ax.grid(alpha=0.3, axis="y")
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=130, bbox_inches="tight")
    print(f"\n[write] {OUT_CSV.relative_to(melt.ROOT)}")
    print(f"[write] {OUT_PNG.relative_to(melt.ROOT)}")


if __name__ == "__main__":
    main()

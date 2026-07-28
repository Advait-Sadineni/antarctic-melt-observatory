"""
How much does the headline pond area depend on the NDWI threshold?

0.16 is an arbitrary constant inherited from the literature. This sweeps it
across the published range and plots detected area against it. A plateau
means the number is robust; a steep slope means it is mostly a choice.

Run:  python src/tune_threshold.py
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import melt

ITEM_ID = "S2B_19CEV_20210124_0_L2A"
THRESHOLDS = np.round(np.arange(0.04, 0.42, 0.02), 2)
OUT_PNG = melt.ROOT / "output" / "threshold_sensitivity.png"


def main():
    item = melt.get_item(ITEM_ID)
    bands = melt.load_scene(item, bands=("green", "nir", "red"))
    reject = melt.reject_mask(item)  # full-res; screen() returns a decimated one

    areas = []
    print(f"[sweep] {ITEM_ID}\n")
    print("  thresh   pond km2   change vs prev")
    for t in THRESHOLDS:
        # Core sweep: shadow test on, hysteresis off, so the curve isolates
        # sensitivity to the NDWI cutoff itself. The full method's sensitivity
        # (with hysteresis) is measured in uncertainty.py.
        _, ponds, valid = melt.detect(bands["green"], bands["nir"], reject,
                                      threshold=float(t), red=bands["red"],
                                      grow_threshold=None)
        km2 = melt.pond_stats(ponds, valid)["pond_km2"]
        delta = "" if not areas else f"{(km2 - areas[-1]) / max(areas[-1], 1e-9) * 100:+7.1f}%"
        marker = "  <-- current default" if abs(t - melt.NDWI_THRESHOLD) < 1e-9 else ""
        print(f"  {t:5.2f}   {km2:8.3f}   {delta:>8}{marker}")
        areas.append(km2)

    areas = np.array(areas)

    # Local sensitivity around the default: how much area moves per 0.01 of
    # threshold, expressed as a percentage of the value at the default.
    i = int(np.argmin(np.abs(THRESHOLDS - melt.NDWI_THRESHOLD)))
    lo, hi = areas[max(i - 1, 0)], areas[min(i + 1, len(areas) - 1)]
    span = (THRESHOLDS[min(i + 1, len(areas) - 1)] - THRESHOLDS[max(i - 1, 0)]) * 100
    sens = (lo - hi) / max(areas[i], 1e-9) * 100 / max(span, 1e-9)
    print(f"\n  at threshold {melt.NDWI_THRESHOLD}: {areas[i]:.3f} km2")
    print(f"  local sensitivity: {sens:.1f}% of area per 0.01 of threshold")

    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.plot(THRESHOLDS, areas, "o-", color="#0077b6", lw=2, ms=5)
    ax.axvline(melt.NDWI_THRESHOLD, color="#d00000", ls="--", lw=1.4,
               label=f"default = {melt.NDWI_THRESHOLD}")
    ax.axvspan(0.10, 0.25, color="#0077b6", alpha=0.08, label="published range 0.10-0.25")
    ax.axhline(melt.WIN_SIZE**2 * melt.PIXEL_KM2, color="#888", ls=":", lw=1.2,
               label="whole AOI (total saturation)")
    # Log scale: the saturation branch below 0.08 is ~50x the plateau, and on a
    # linear axis it flattens the plateau - the part we actually care about -
    # into an unreadable line along the bottom.
    ax.set_yscale("log")
    ax.set_xlabel("NDWI threshold")
    ax.set_ylabel("Detected meltwater area (km$^2$, log scale)")
    ax.set_title(f"Threshold sensitivity - George VI, {ITEM_ID.split('_')[2]}")
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()
    OUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, dpi=130)
    print(f"\n[write] {OUT_PNG.relative_to(melt.ROOT)}")


if __name__ == "__main__":
    main()

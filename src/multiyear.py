"""
Multi-year melt trend for George VI Ice Shelf.

Runs the season pipeline across every Sentinel-2 melt season and charts peak
observed meltwater area per season.

Read the output carefully. This is *peak observed* area, not true seasonal
peak: if the real maximum fell on a cloudy day, that season is undercounted.
The number of usable scenes is plotted alongside for exactly that reason - a
season with three clear scenes is not comparable to one with twelve.

Run:  python src/multiyear.py
"""

import csv

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import melt
import season as S

# Sentinel-2A launched 2015, 2B in 2017. Antarctic coverage is thin early on.
SEASONS = ["2017-18", "2018-19", "2019-20", "2020-21", "2021-22",
           "2022-23", "2023-24", "2024-25", "2025-26"]

# Below this many clear scenes, the seasonal peak is likely simply missed.
MIN_SCENES_TRUSTED = 5

# Below this many, the season is not characterised at all and is dropped from
# the chart rather than drawn as a bar. 2017-18 had exactly one usable scene,
# in March, which is not a melt season by any reading.
MIN_SCENES_PLOTTED = 2

OUT_CSV = melt.ROOT / "output" / "multiyear_summary.csv"
OUT_PNG = melt.ROOT / "output" / "multiyear_trend.png"


def main():
    summary = []
    for label in SEASONS:
        print(f"\n{'=' * 70}")
        rows = S.run_season(label, quiet=True)
        S.plot(label)

        if not rows:
            print(f"[{label}] no usable scenes")
            summary.append({"season": label, "n_usable": 0, "peak_km2": float("nan"),
                            "peak_date": "", "mean_top3_km2": float("nan")})
            continue

        areas = sorted((float(r["pond_km2"]) for r in rows), reverse=True)
        peak_row = max(rows, key=lambda r: float(r["pond_km2"]))
        summary.append({
            "season": label,
            "n_usable": len(rows),
            "peak_km2": round(float(peak_row["pond_km2"]), 4),
            "peak_date": peak_row["date"],
            "mean_top3_km2": round(float(np.mean(areas[:3])), 4),
        })
        print(f"[{label}] {len(rows)} usable, peak {summary[-1]['peak_km2']:.2f} km2"
              f" on {summary[-1]['peak_date']}")

    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(summary[0].keys()))
        w.writeheader()
        w.writerows(summary)

    plot(summary)


def plot(summary):
    dropped = [s for s in summary if s["n_usable"] < MIN_SCENES_PLOTTED]
    for s in dropped:
        print(f"[plot] {s['season']} omitted - only {s['n_usable']} usable scene(s)")
    summary = [s for s in summary if s["n_usable"] >= MIN_SCENES_PLOTTED]

    labels = [s["season"] for s in summary]
    peaks = [s["peak_km2"] for s in summary]
    counts = [s["n_usable"] for s in summary]
    x = np.arange(len(labels))

    fig, (ax, ax2) = plt.subplots(
        2, 1, figsize=(11, 7.6), sharex=True,
        gridspec_kw={"height_ratios": [3, 1], "hspace": 0.12})

    floor = melt.NOISE_FLOOR_KM2
    bars = ax.bar(x, peaks, color="#0077b6", width=0.62)

    for b, xi, p, c in zip(bars, x, peaks, counts):
        below = (not np.isnan(p)) and p <= floor
        if below:
            # Not a measurement. Draw it hollow so it cannot be read as one.
            b.set_facecolor("none")
            b.set_edgecolor("#aab4bc")
            b.set_hatch("///")
            ax.text(xi, floor * 1.06, "no melt\ndetected", ha="center", va="bottom",
                    fontsize=7.5, color="#77828a", linespacing=1.15)
        elif c < MIN_SCENES_TRUSTED:
            b.set_color("#b8c4cc")
            ax.text(xi, p, f"{p:.1f}", ha="center", va="bottom", fontsize=9)
        else:
            ax.text(xi, p, f"{p:.1f}", ha="center", va="bottom", fontsize=9)

    ax.axhline(floor, color="#c1121f", ls="--", lw=1.2, zorder=3)

    ax.set_ylabel("Peak observed meltwater area (km$^2$)")
    ax.set_title("George VI Ice Shelf - peak meltwater by melt season\n"
                 f"{melt.WIN_SIZE * int(melt.PIXEL_M) / 1000:.0f} km AOI, "
                 f"NDWI > {melt.NDWI_THRESHOLD}, Sentinel-2 L2A",
                 fontsize=12)
    ax.grid(alpha=0.3, axis="y")
    ax.set_ylim(0, max(1.0, np.nanmax(peaks)) * 1.18)  # headroom for labels
    ax.legend(handles=[
        plt.Rectangle((0, 0), 1, 1, color="#0077b6", label="measured peak"),
        plt.Rectangle((0, 0), 1, 1, color="#b8c4cc",
                      label=f"<{MIN_SCENES_TRUSTED} clear scenes - likely undercounted"),
        plt.Rectangle((0, 0), 1, 1, facecolor="none", edgecolor="#aab4bc", hatch="///",
                      label="at or below noise floor - not a measurement"),
        plt.Line2D([0], [0], color="#c1121f", ls="--", lw=1.2,
                   label=f"noise floor, {floor:.1f} km$^2$"),
    ], loc="upper left", fontsize=8.5, framealpha=0.9)

    ax2.bar(x, counts, color="#8d99ae", width=0.62)
    ax2.set_ylabel("usable\nscenes", fontsize=9)
    ax2.set_xticks(x)
    ax2.set_xticklabels(labels, rotation=45, ha="right")
    ax2.grid(alpha=0.3, axis="y")

    fig.savefig(OUT_PNG, dpi=130, bbox_inches="tight")
    print(f"\n[write] {OUT_CSV.relative_to(melt.ROOT)}")
    print(f"[write] {OUT_PNG.relative_to(melt.ROOT)}")


if __name__ == "__main__":
    main()

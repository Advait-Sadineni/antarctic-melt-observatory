"""Multi-shelf comparison figure: the observatory's first cross-shelf view.

Six shelves x nine seasons, small multiples (scales differ 40x, so no shared or
dual axis). Single-hue magnitude bars; each shelf's record season highlighted;
poorly-observed seasons hatched; first-pass shelves labelled as such.

Run:  python src/multishelf_figure.py  ->  output/multishelf_comparison.png
"""

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

import melt

OUT = melt.ROOT / "output" / "multishelf_comparison.png"
SEASONS = ["2017-18", "2018-19", "2019-20", "2020-21", "2021-22",
           "2022-23", "2023-24", "2024-25", "2025-26"]

# (title, history file, status note)
SHELVES = [
    ("George VI", "history.json", "validated"),
    ("Wilkins", "history_Wilkins.json", "first pass"),
    ("Larsen C", "history_LarsenC.json", "first pass"),
    ("Bach", "history_Bach.json", "first pass"),
    ("Stange", "history_Stange.json", "first pass"),
    ("Amery", "history_Amery.json", "first pass - in progress"),
]

INK, DIM, FAINT = "#122031", "#4c5f74", "#8fa3b8"
BAR, HOT, WARN = "#1f93cf", "#0c5f8f", "#b4801f"
GROUND, PANEL = "#f4f8fb", "#ffffff"


def main():
    fig, axes = plt.subplots(2, 3, figsize=(15.5, 8.2), facecolor=GROUND)
    fig.suptitle("Antarctic ice-shelf meltwater, season by season",
                 fontsize=17, fontweight="bold", color=INK, y=0.985)
    fig.text(0.5, 0.938, "seasonal-maximum meltwater area (km²) · same validated "
             "detector on every shelf · note the different y-scales",
             ha="center", fontsize=10.5, color=DIM)

    for ax, (name, fname, status) in zip(axes.ravel(), SHELVES):
        ax.set_facecolor(PANEL)
        p = melt.ROOT / "output" / "shelf" / fname
        rows = {r["season"]: r for r in json.loads(p.read_text())} if p.exists() else {}
        vals = [rows.get(s, {}).get("shelf_km2") for s in SEASONS]
        present = [v for v in vals if v is not None]
        peak = max(present) if present else 0.0
        area = next((r["shelf_area_km2"] for r in rows.values()), None)

        for i, s in enumerate(SEASONS):
            r = rows.get(s)
            if r is None:
                continue
            v = r["shelf_km2"]
            is_peak = v == peak and peak > 0
            ax.bar(i, v, width=0.66, color=HOT if is_peak else BAR,
                   hatch="///" if r.get("poorly_observed") else None,
                   edgecolor=WARN if r.get("poorly_observed") else "none",
                   linewidth=0.8, zorder=3)
            if is_peak:
                ax.annotate(f"{v:,.0f}", (i, v), textcoords="offset points",
                            xytext=(0, 3), ha="center", fontsize=9.5,
                            fontweight="bold", color=INK, zorder=4)

        pct = f"  ·  peak = {100*peak/area:.1f}% of shelf" if area and peak else ""
        ax.set_title(f"{name}", fontsize=12.5, fontweight="bold",
                     color=INK, loc="left", pad=10)
        ax.text(0.0, 1.012, f"\n{status}{pct}", transform=ax.transAxes,
                fontsize=8.5, color=WARN if "first" in status else DIM, va="bottom")
        ax.set_xticks(range(len(SEASONS)))
        ax.set_xticklabels([s.replace("20", "'", 1)[:6] for s in SEASONS],
                           fontsize=8, color=DIM, rotation=45, ha="right")
        ax.tick_params(axis="y", labelsize=8.5, colors=DIM, length=0)
        ax.grid(axis="y", color="#dbe4ec", linewidth=0.7, zorder=0)
        for sp in ax.spines.values():
            sp.set_visible(False)
        ax.margins(y=0.16)

    fig.legend(handles=[
        Patch(color=BAR, label="seasonal max (km²)"),
        Patch(color=HOT, label="shelf's record season"),
        Patch(facecolor=BAR, hatch="///", edgecolor=WARN, label="poorly observed (flagged)"),
    ], loc="lower center", ncol=3, frameon=False, fontsize=9.5,
        bbox_to_anchor=(0.5, -0.005))
    fig.tight_layout(rect=(0, 0.03, 1, 0.925))
    fig.savefig(OUT, dpi=125, facecolor=GROUND, bbox_inches="tight")
    print("wrote", OUT)


if __name__ == "__main__":
    main()

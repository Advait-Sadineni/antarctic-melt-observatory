"""
Seasonal evolution of meltwater across all melt seasons.

The multi-year chart shows each season's single peak. This shows the whole
within-season curve — the rise from November, the January/February maximum,
and the autumn decline — as small multiples aligned on a common calendar, so
seasons can be compared directly. This is the standard way the supraglacial-
lake literature presents seasonal evolution (e.g. Dirscherl et al. 2021).

Uses the committed season CSVs; no new detection, no network.

Run:  python src/evolution.py
"""

import csv
import glob
from datetime import date

import matplotlib

matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np

import melt
import season as S

OUT_PNG = melt.ROOT / "output" / "seasonal_evolution.png"


def _season_series(label):
    """One point per date (deduped), sorted, as (date, km2_equiv, cloud)."""
    rows = S.series(label)
    out = []
    for r in rows:
        out.append((date.fromisoformat(r["date"]),
                    float(r["pond_km2_equiv"]),
                    100 - float(r["usable_frac"]) * 100))
    return sorted(out)


def _day_of_season(d):
    """Days since 1 November of that melt season (which spans the new year)."""
    start_year = d.year if d.month >= 7 else d.year - 1
    return (d - date(start_year, 11, 1)).days


def main():
    labels = sorted(p.rsplit("season_", 1)[1][:7].replace("_", "-")
                    for p in glob.glob(str(melt.ROOT / "output" / "season_*.csv")))
    series = {lab: _season_series(lab) for lab in labels}
    series = {k: v for k, v in series.items() if v}

    ncol = 3
    nrow = int(np.ceil(len(series) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(ncol * 4.2, nrow * 3.0),
                             sharex=True, sharey=True)
    axes = np.atleast_1d(axes).ravel()

    ymax = max(v for s in series.values() for _, v, _ in s) * 1.1
    floor = melt.NOISE_FLOOR_KM2

    # Month ticks in day-of-season units (1 Nov = 0).
    month_starts = [(date(2019, m, 1) if m >= 11 else date(2020, m, 1))
                    for m in (11, 12, 1, 2, 3)]
    ticks = [_day_of_season(d) for d in month_starts]
    ticklabels = ["Nov", "Dec", "Jan", "Feb", "Mar"]

    for ax, (lab, s) in zip(axes, series.items()):
        xs = [_day_of_season(d) for d, _, _ in s]
        ys = [v for _, v, _ in s]
        cl = [c for _, _, c in s]
        ax.axhline(floor, color="#c1121f", ls="--", lw=0.9, alpha=0.6)
        ax.plot(xs, ys, "-", color="#0077b6", lw=1.2, alpha=0.5, zorder=1)
        ax.scatter(xs, ys, c=cl, cmap="viridis_r", vmin=0, vmax=50, s=34,
                   zorder=2, edgecolor="white", linewidth=0.5)
        peak = max(range(len(ys)), key=lambda i: ys[i])
        ax.annotate(f"{ys[peak]:.1f}", (xs[peak], ys[peak]),
                    textcoords="offset points", xytext=(0, 5), ha="center",
                    fontsize=8, color="#0a3d5c")
        ax.set_title(lab, fontsize=11, loc="left")
        ax.set_xticks(ticks)
        ax.set_xticklabels(ticklabels, fontsize=8)
        ax.grid(alpha=0.25)

    for ax in axes[len(series):]:
        ax.axis("off")
    axes[0].set_ylim(0, ymax)

    fig.suptitle("George VI Ice Shelf — seasonal meltwater evolution\n"
                 "area over the 61 km study area; point colour = % screened out; "
                 "dashed line = noise floor",
                 fontsize=12)
    fig.supylabel("Meltwater area (km$^2$)", fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=130, bbox_inches="tight")

    # A one-line summary of peak timing, to compare with the literature's
    # documented mid-to-late-January maximum.
    print("season   peak date    day-of-season   peak km2")
    for lab, s in series.items():
        d, v, _ = max(s, key=lambda t: t[1])
        print(f"{lab}   {d.isoformat()}   {_day_of_season(d):>3d}            {v:.2f}")
    print(f"\n[write] {OUT_PNG.relative_to(melt.ROOT)}")


if __name__ == "__main__":
    main()

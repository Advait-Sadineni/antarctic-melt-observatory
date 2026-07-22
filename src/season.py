"""
One melt season -> pond-area time series.

Screens every Sentinel-2 scene on the tile by cloud fraction *inside the AOI*,
runs NDWI detection on the survivors, and writes a CSV plus a chart of
meltwater area against date.

Resumable: results append to the CSV, and a rerun skips scenes already done.

Run:  python src/season.py
"""

import csv
from datetime import date

import matplotlib

matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt

import melt

# Austral melt season: spans the new year, so November -> March.
SEASON_START, SEASON_END = "2020-11-01", "2021-03-31"
SEASON_LABEL = "2020-21"

MAX_CLOUD_IN_AOI = 0.20  # reject scene if >20% of the AOI is cloud/shadow
MAX_NODATA = 0.05

# Contamination gate. Cloud shadow on an ice shelf is lit by blue skylight,
# which raises green relative to NIR and mimics water's NDWI signature - it
# passes both the cloud test (SCL under-flags shadow over ice) and the
# brightness floor (haze lifts green above it). SCL's water class does catch
# it, so an implausible amount of "open water" over an ice shelf is the
# tell. Clear scenes here run <=0.42%; the one contaminated scene hit 7.7%.
MAX_WATER_IN_AOI = 0.02

OUT_CSV = melt.ROOT / "output" / f"season_{SEASON_LABEL.replace('-', '_')}.csv"
OUT_PNG = melt.ROOT / "output" / f"season_{SEASON_LABEL.replace('-', '_')}.png"

FIELDS = ["item_id", "date", "cloud_frac_aoi", "water_frac_aoi", "usable_frac",
          "tile_cloud", "pond_px", "pond_km2", "usable_km2", "pond_pct_of_usable"]


def load_done():
    if not OUT_CSV.exists():
        return {}
    with open(OUT_CSV, newline="") as f:
        return {r["item_id"]: r for r in csv.DictReader(f)}


def append_row(row):
    exists = OUT_CSV.exists()
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_CSV, "a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        if not exists:
            w.writeheader()
        w.writerow(row)


def process(item):
    """Screen on SCL first; only fetch the 10 m bands if the AOI is clear."""
    scr = melt.screen(melt.read_scl(item))

    if (scr["nodata_frac"] > MAX_NODATA
            or scr["cloud_frac"] > MAX_CLOUD_IN_AOI
            or scr["water_frac"] > MAX_WATER_IN_AOI):
        return None, scr

    bands = melt.load_scene(item, bands=("green", "nir"))
    _, ponds, valid = melt.detect(bands["green"], bands["nir"], scr["reject_mask"])

    # Cloud removes pixels from consideration, which would drag absolute area
    # down on hazier days. Track usable area so the two can be compared.
    usable = valid & ~scr["reject_mask"]
    st = melt.pond_stats(ponds, usable)
    return st, scr


def main():
    print(f"[season] {SEASON_LABEL}  tile {melt.TILE}  {SEASON_START} -> {SEASON_END}")
    items = melt.search_scenes(SEASON_START, SEASON_END)
    print(f"[season] {len(items)} scenes on tile\n")

    done = load_done()
    for n, item in enumerate(items, 1):
        if item.id in done:
            continue
        tag = f"[{n:3d}/{len(items)}] {item.datetime.date()} {item.id[-14:]}"
        try:
            st, scr = process(item)
        except Exception as e:
            print(f"{tag}  ERROR {type(e).__name__}: {str(e)[:60]}")
            continue

        if st is None:
            why = ("water" if scr["water_frac"] > MAX_WATER_IN_AOI
                   else "nodata" if scr["nodata_frac"] > MAX_NODATA else "cloud")
            print(f"{tag}  skip [{why:6s}] (cloud {scr['cloud_frac']*100:5.1f}%, "
                  f"nodata {scr['nodata_frac']*100:5.1f}%, "
                  f"water {scr['water_frac']*100:5.2f}%)")
            continue

        print(f"{tag}  cloud {scr['cloud_frac']*100:5.1f}%  "
              f"pond {st['pond_km2']:7.3f} km2  ({st['pond_pct']:5.2f}% of usable)")
        append_row({
            "item_id": item.id,
            "date": item.datetime.date().isoformat(),
            "cloud_frac_aoi": round(scr["cloud_frac"], 5),
            "water_frac_aoi": round(scr["water_frac"], 5),
            "usable_frac": round(scr["usable_frac"], 5),
            "tile_cloud": round(item.properties.get("eo:cloud_cover") or 0, 3),
            "pond_px": st["pond_px"],
            "pond_km2": round(st["pond_km2"], 4),
            "usable_km2": round(st["valid_km2"], 2),
            "pond_pct_of_usable": round(st["pond_pct"], 4),
        })

    plot()


def plot():
    rows = list(load_done().values())
    if not rows:
        print("[plot] nothing to plot")
        return

    # A date can yield more than one usable scene (overlapping orbits, or the
    # archive holding several processings). Keep the clearest one per date so
    # the series is one measurement per day. All of them stay in the CSV.
    best = {}
    for r in rows:
        cur = best.get(r["date"])
        if cur is None or float(r["cloud_frac_aoi"]) < float(cur["cloud_frac_aoi"]):
            best[r["date"]] = r
    dropped = len(rows) - len(best)
    rows = sorted(best.values(), key=lambda r: r["date"])
    if dropped:
        print(f"[plot] {dropped} same-date duplicate(s) collapsed to clearest scene")

    dates = [date.fromisoformat(r["date"]) for r in rows]
    areas = [float(r["pond_km2"]) for r in rows]
    clouds = [float(r["cloud_frac_aoi"]) * 100 for r in rows]

    fig, ax = plt.subplots(figsize=(11, 5.8))
    ax.plot(dates, areas, "-", color="#0077b6", lw=1.4, alpha=0.7, zorder=1)
    sc = ax.scatter(dates, areas, c=clouds, cmap="viridis_r", vmin=0,
                    vmax=MAX_CLOUD_IN_AOI * 100, s=55, zorder=2,
                    edgecolor="white", linewidth=0.6)
    fig.colorbar(sc, ax=ax, label="cloud in AOI (%)")

    peak = max(range(len(areas)), key=lambda i: areas[i])
    ax.annotate(f"peak {areas[peak]:.1f} km$^2$\n{dates[peak]:%d %b %Y}",
                xy=(dates[peak], areas[peak]),
                xytext=(12, -28), textcoords="offset points", fontsize=9,
                arrowprops=dict(arrowstyle="->", color="#444", lw=1))

    ax.set_ylabel("Meltwater area (km$^2$)")
    ax.set_title(f"George VI Ice Shelf meltwater, {SEASON_LABEL} melt season\n"
                 f"{melt.WIN_SIZE * int(melt.PIXEL_M) / 1000:.0f} km AOI, "
                 f"NDWI > {melt.NDWI_THRESHOLD}, {len(rows)} usable scenes",
                 fontsize=12)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%d %b"))
    ax.grid(alpha=0.3)
    fig.autofmt_xdate()
    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=130)
    print(f"\n[write] {OUT_CSV.relative_to(melt.ROOT)}")
    print(f"[write] {OUT_PNG.relative_to(melt.ROOT)}")


if __name__ == "__main__":
    main()

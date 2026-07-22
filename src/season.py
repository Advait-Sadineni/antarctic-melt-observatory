"""
One melt season -> pond-area time series.

Screens every Sentinel-2 scene on the tile, runs NDWI detection on the
survivors, and writes a CSV plus a chart of meltwater area against date.

Screening is layered cheapest-first, and every layer earned its place by
catching a specific false positive - see melt.py for the measurements:

  1. nodata      off-swath scenes
  2. SCL cloud   the easy cases
  3. SCL water   cloud shadow mimicking ponds (a 46.9 km2 phantom in November)
  4. NDSI        cloud over ice, which SCL badly under-reports

Resumable: results append to the CSV, and a rerun skips scenes already done.

Run:  python src/season.py            (default season)
      python src/season.py 2019-20    (any season)
"""

import csv
import sys
from concurrent.futures import ThreadPoolExecutor
from datetime import date

import matplotlib

matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt

import melt

DEFAULT_SEASON = "2020-21"

MAX_CLOUD_IN_AOI = 0.20  # reject scene if >20% of the AOI is cloud/shadow
MAX_NODATA = 0.05

# Contamination gate. Cloud shadow on an ice shelf is lit by blue skylight,
# which raises green relative to NIR and mimics water's NDWI signature - it
# passes both the cloud test (SCL under-flags shadow over ice) and the
# brightness floor (haze lifts green above it). SCL's water class does catch
# it, so an implausible amount of "open water" over an ice shelf is the
# tell. Clear scenes here run <=0.42%; the one contaminated scene hit 7.7%.
MAX_WATER_IN_AOI = 0.02

# Cloud gate that works over ice, where SCL does not. See melt.ice_check().
# Clean scenes measure 0.0-0.1% non-ice; cloud-covered ones 3.7-35.9%.
MAX_NONICE_IN_AOI = 0.01

FIELDS = ["item_id", "date", "cloud_frac_aoi", "water_frac_aoi", "nonice_frac_aoi",
          "usable_frac", "tile_cloud", "pond_px", "pond_km2", "usable_km2",
          "pond_pct_of_usable"]


def season_dates(label):
    """'2020-21' -> ('2020-11-01', '2021-03-31'). Austral melt spans the year."""
    start_year = int(label.split("-")[0])
    return f"{start_year}-11-01", f"{start_year + 1}-03-31"


def csv_path(label):
    return melt.ROOT / "output" / f"season_{label.replace('-', '_')}.csv"


def png_path(label):
    return melt.ROOT / "output" / f"season_{label.replace('-', '_')}.png"


def load_done(label):
    path = csv_path(label)
    if not path.exists():
        return {}
    with open(path, newline="") as f:
        return {r["item_id"]: r for r in csv.DictReader(f)}


def append_row(label, row):
    path = csv_path(label)
    exists = path.exists()
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        if not exists:
            w.writeheader()
        w.writerow(row)


def process(item):
    """Screen cheapest-first; only fetch 10 m bands once a scene has passed."""
    scr = melt.screen(melt.read_scl(item))
    scr["nonice_frac"] = float("nan")

    if (scr["nodata_frac"] > MAX_NODATA
            or scr["cloud_frac"] > MAX_CLOUD_IN_AOI
            or scr["water_frac"] > MAX_WATER_IN_AOI):
        return None, scr

    # NDSI gate last: it costs two extra reads, so only scenes that already
    # passed the cheap SCL checks pay for it.
    scr["nonice_frac"] = melt.ice_check(item)
    if scr["nonice_frac"] > MAX_NONICE_IN_AOI:
        return None, scr

    bands = melt.load_scene(item, bands=("green", "nir"))
    _, ponds, valid = melt.detect(bands["green"], bands["nir"], scr["reject_mask"])

    # Cloud removes pixels from consideration, which would drag absolute area
    # down on hazier days. Track usable area so the two can be compared.
    usable = valid & ~scr["reject_mask"]
    return melt.pond_stats(ponds, usable), scr


def _safe_process(item):
    try:
        st, scr = process(item)
        return st, scr, None
    except Exception as e:
        return None, None, e


def run_season(label, quiet=False, workers=6):
    start, end = season_dates(label)
    print(f"[season] {label}  tile {melt.TILE}  {start} -> {end}")
    items = melt.search_scenes(start, end)
    print(f"[season] {len(items)} scenes on tile")

    done = load_done(label)
    todo = [i for i in items if i.id not in done]
    if not todo:
        print("[season] all scenes already done")
        return series(label)

    # Screening is almost entirely waiting on HTTP range requests, so threads
    # help a lot here even under the GIL. ex.map yields in submission order,
    # which keeps the log readable and the CSV chronological.
    with ThreadPoolExecutor(max_workers=workers) as ex:
        results = ex.map(_safe_process, todo)

        for n, (item, (st, scr, err)) in enumerate(zip(todo, results), 1):
            tag = f"[{n:3d}/{len(todo)}] {item.datetime.date()} {item.id[-14:]}"
            if err is not None:
                print(f"{tag}  ERROR {type(err).__name__}: {str(err)[:60]}")
                continue
            _record(label, item, st, scr, tag, quiet)

    return series(label)


def _record(label, item, st, scr, tag, quiet):
    if st is None:
        why = ("water" if scr["water_frac"] > MAX_WATER_IN_AOI
               else "nodata" if scr["nodata_frac"] > MAX_NODATA
               else "cloud" if scr["cloud_frac"] > MAX_CLOUD_IN_AOI else "nonice")
        if not quiet:
            print(f"{tag}  skip [{why:6s}] (cloud {scr['cloud_frac']*100:5.1f}%, "
                  f"nodata {scr['nodata_frac']*100:5.1f}%, "
                  f"water {scr['water_frac']*100:5.2f}%, "
                  f"nonice {scr['nonice_frac']*100:6.2f}%)")
        return

    print(f"{tag}  cloud {scr['cloud_frac']*100:5.1f}%  nonice {scr['nonice_frac']*100:4.2f}%  "
          f"pond {st['pond_km2']:7.3f} km2  ({st['pond_pct']:5.2f}% of usable)")
    append_row(label, {
        "item_id": item.id,
        "date": item.datetime.date().isoformat(),
        "cloud_frac_aoi": round(scr["cloud_frac"], 5),
        "water_frac_aoi": round(scr["water_frac"], 5),
        "nonice_frac_aoi": round(scr["nonice_frac"], 5),
        "usable_frac": round(scr["usable_frac"], 5),
        "tile_cloud": round(item.properties.get("eo:cloud_cover") or 0, 3),
        "pond_px": st["pond_px"],
        "pond_km2": round(st["pond_km2"], 4),
        "usable_km2": round(st["valid_km2"], 2),
        "pond_pct_of_usable": round(st["pond_pct"], 4),
    })


def series(label):
    """One measurement per date: the clearest scene. All rows stay in the CSV."""
    rows = list(load_done(label).values())
    if not rows:
        return []
    best = {}
    for r in rows:
        cur = best.get(r["date"])
        if cur is None or float(r["cloud_frac_aoi"]) < float(cur["cloud_frac_aoi"]):
            best[r["date"]] = r
    return sorted(best.values(), key=lambda r: r["date"])


def plot(label):
    rows = series(label)
    if not rows:
        print("[plot] nothing to plot")
        return

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
    ax.set_title(f"George VI Ice Shelf meltwater, {label} melt season\n"
                 f"{melt.WIN_SIZE * int(melt.PIXEL_M) / 1000:.0f} km AOI, "
                 f"NDWI > {melt.NDWI_THRESHOLD}, {len(rows)} usable scenes",
                 fontsize=12)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%d %b"))
    ax.grid(alpha=0.3)
    fig.autofmt_xdate()
    fig.tight_layout()
    fig.savefig(png_path(label), dpi=130)
    print(f"\n[write] {csv_path(label).relative_to(melt.ROOT)}")
    print(f"[write] {png_path(label).relative_to(melt.ROOT)}")


def main():
    label = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_SEASON
    run_season(label)
    plot(label)


if __name__ == "__main__":
    main()

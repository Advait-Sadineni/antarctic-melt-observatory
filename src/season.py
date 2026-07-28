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

# A scene is kept if enough of the study area survives screening, and the
# result is reported as a *density* over that surviving ground rather than as
# a raw area. Demanding a clear 61 km square is far stricter than demanding a
# clear 20 km one: the 2026-02-10 record scene has 23.7% cloud over the large
# AOI and only 5.8% over the small one, so a whole-scene cloud gate threw away
# 76% of perfectly good ground. Requiring half the AOI and normalising keeps
# it, and keeps scenes comparable when their usable areas differ.
MIN_USABLE_FRAC = 0.50

# Usable fraction alone is not enough. A scene with a third of the AOI under
# detected cloud still had two thirds "usable", passed, and returned 365 km2
# of meltwater - thirty times anything plausible. Heavy detected cloud implies
# heavy *undetected* cloud in what is left, so the halo fraction gates the
# scene as well.
#
# First set at 15%, which was still too loose: 2024-12-20 passed at 10% halo
# and returned 34 km2, and inspection showed a cumulus band across the scene
# with the mask lighting up all over it. Scenes at 11% and 13% were the same
# story. Observed halo: clean scenes 0.03-0.06%, a salvageable partly-cloudy
# one 1.99%, everything at or above 10% contaminated. 5% cuts between them.
MAX_HALO_FRAC = 0.05

# Last-resort tripwire for failure modes not yet modelled. The largest
# meltwater density measured on a visually confirmed clean scene is 0.49% of
# the study area, so anything above 2% is treated as a detector failure
# rather than a record melt event. Rejections are logged, never silent - if
# this ever fires on a real scene, that is a finding, not a nuisance.
MAX_POND_FRAC_PLAUSIBLE = 0.02

# Contamination gate. Cloud shadow on an ice shelf is lit by blue skylight,
# which raises green relative to NIR and mimics water's NDWI signature - it
# passes both the cloud test (SCL under-flags shadow over ice) and the
# brightness floor (haze lifts green above it). SCL's water class does catch
# it, so an implausible amount of "open water" over an ice shelf is the
# tell. Clear scenes here run <=0.42%; the one contaminated scene hit 7.7%.
MAX_WATER_IN_AOI = 0.02

FIELDS = ["item_id", "date", "sun_elev", "cloud_frac_aoi", "water_frac_aoi",
          "nonice_frac_aoi", "halo_frac_aoi", "usable_frac", "tile_cloud",
          "pond_px", "pond_km2", "usable_km2", "pond_pct_of_usable",
          "pond_km2_equiv"]


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


AOI_KM2 = melt.WIN_SIZE**2 * melt.PIXEL_KM2


def process(item):
    """Screen cheapest-first; only fetch 10 m bands once a scene has passed."""
    scr = melt.screen(melt.read_scl(item))
    scr["sun_elev"] = melt.sun_elevation(item)
    scr["nonice_frac"] = float("nan")
    scr["halo_frac"] = float("nan")

    # Grazing illumination makes crevasse and topographic shadow indistinguishable
    # from meltwater, so those scenes are unusable rather than merely noisy.
    if scr["sun_elev"] < melt.MIN_SUN_ELEVATION_DEG:
        return None, scr
    if scr["water_frac"] > MAX_WATER_IN_AOI:
        return None, scr

    # NDSI costs two extra reads, so only scenes that passed the cheap SCL
    # checks pay for it.
    cm = melt.cloud_mask(item)
    scr["nonice_frac"] = cm["nonice_frac"]
    scr["halo_frac"] = cm["halo_frac"]

    if cm["halo_frac"] > MAX_HALO_FRAC:
        scr["usable_frac"] = float("nan")
        return None, scr

    # Cheap upper bound on usable ground: SCL-clear minus the cloud halo.
    # Bail before paying for full-resolution reads if it cannot pass.
    approx_usable = float((~scr["reject_mask"]).mean()) - cm["halo_frac"]
    if approx_usable < MIN_USABLE_FRAC:
        scr["usable_frac"] = max(approx_usable, 0.0)
        return None, scr

    # Only now pay for full-resolution reads. Red is needed for the shadow test.
    bands = melt.load_scene(item, bands=("green", "nir", "red"))
    reject = melt.reject_mask(item) | cm["mask"]
    _, ponds, valid = melt.detect(bands["green"], bands["nir"], reject, red=bands["red"])

    usable = valid & ~reject
    scr["usable_frac"] = float(usable.mean())
    if scr["usable_frac"] < MIN_USABLE_FRAC:
        return None, scr

    st = melt.pond_stats(ponds, usable)
    if st["pond_pct"] / 100.0 > MAX_POND_FRAC_PLAUSIBLE:
        scr["implausible_pct"] = st["pond_pct"]
        return None, scr

    # Density scaled back to the whole study area, so scenes with different
    # usable footprints stay comparable.
    st["pond_km2_equiv"] = st["pond_pct"] / 100.0 * AOI_KM2
    return st, scr


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
        if "implausible_pct" in scr:
            # Never silent: this means the detector failed in a new way.
            print(f"{tag}  REJECT [implausible] {scr['implausible_pct']:.2f}% of usable "
                  f"is meltwater (cap {MAX_POND_FRAC_PLAUSIBLE*100:.0f}%), "
                  f"halo {scr['halo_frac']*100:.1f}%")
            return
        why = ("sun" if scr["sun_elev"] < melt.MIN_SUN_ELEVATION_DEG
               else "water" if scr["water_frac"] > MAX_WATER_IN_AOI
               else "halo" if scr["halo_frac"] > MAX_HALO_FRAC else "usable")
        if not quiet:
            print(f"{tag}  skip [{why:6s}] (sun {scr['sun_elev']:4.1f}, "
                  f"cloud {scr['cloud_frac']*100:5.1f}%, "
                  f"nodata {scr['nodata_frac']*100:5.1f}%, "
                  f"water {scr['water_frac']*100:5.2f}%, "
                  f"halo {scr['halo_frac']*100:6.2f}%, "
                  f"usable {scr['usable_frac']*100:5.1f}%)")
        return

    print(f"{tag}  sun {scr['sun_elev']:4.1f}  usable {scr['usable_frac']*100:5.1f}%  "
          f"halo {scr['halo_frac']*100:5.2f}%  "
          f"pond {st['pond_km2_equiv']:7.3f} km2-eq  ({st['pond_pct']:5.3f}% of usable)")
    append_row(label, {
        "item_id": item.id,
        "date": item.datetime.date().isoformat(),
        "sun_elev": round(scr["sun_elev"], 2),
        "cloud_frac_aoi": round(scr["cloud_frac"], 5),
        "water_frac_aoi": round(scr["water_frac"], 5),
        "nonice_frac_aoi": round(scr["nonice_frac"], 5),
        "halo_frac_aoi": round(scr["halo_frac"], 5),
        "usable_frac": round(scr["usable_frac"], 5),
        "tile_cloud": round(item.properties.get("eo:cloud_cover") or 0, 3),
        "pond_px": st["pond_px"],
        "pond_km2": round(st["pond_km2"], 4),
        "usable_km2": round(st["valid_km2"], 2),
        "pond_pct_of_usable": round(st["pond_pct"], 4),
        "pond_km2_equiv": round(st["pond_km2_equiv"], 4),
    })


def _quality(row):
    """Rank same-date scenes: most usable ground first, then least cloud."""
    return (float(row["usable_km2"]), -float(row["cloud_frac_aoi"]))


def series(label):
    """One measurement per date: the best scene. All rows stay in the CSV.

    A date can yield more than one scene from overlapping orbits. Checked
    that these really are the same ground rather than complementary swaths -
    in every observed case both cover 376-419 km2 of the ~419 km2 AOI - so
    picking one is correct and does not discard coverage.
    """
    rows = list(load_done(label).values())
    if not rows:
        return []
    best = {}
    for r in rows:
        cur = best.get(r["date"])
        if cur is None or _quality(r) > _quality(cur):
            best[r["date"]] = r
    return sorted(best.values(), key=lambda r: r["date"])


def plot(label):
    rows = series(label)
    if not rows:
        print("[plot] nothing to plot")
        return

    dates = [date.fromisoformat(r["date"]) for r in rows]
    areas = [float(r["pond_km2_equiv"]) for r in rows]
    clouds = [100 - float(r["usable_frac"]) * 100 for r in rows]

    fig, ax = plt.subplots(figsize=(11, 5.8))
    ax.plot(dates, areas, "-", color="#0077b6", lw=1.4, alpha=0.7, zorder=1)
    sc = ax.scatter(dates, areas, c=clouds, cmap="viridis_r", vmin=0,
                    vmax=(1 - MIN_USABLE_FRAC) * 100, s=55, zorder=2,
                    edgecolor="white", linewidth=0.6)
    fig.colorbar(sc, ax=ax, label="screened out of AOI (%)")

    peak = max(range(len(areas)), key=lambda i: areas[i])
    ax.annotate(f"peak {areas[peak]:.1f} km$^2$\n{dates[peak]:%d %b %Y}",
                xy=(dates[peak], areas[peak]),
                xytext=(12, -28), textcoords="offset points", fontsize=9,
                arrowprops=dict(arrowstyle="->", color="#444", lw=1))

    ax.set_ylabel("Meltwater area (km$^2$, scaled to full AOI)")
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

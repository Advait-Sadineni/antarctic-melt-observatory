"""M4 Task 5: drainage-vs-refreeze event catalogue (spec section 5), offline.

    python scripts/drainage_catalogue.py 2019-20

Candidates: pond components losing >=80% of cells between consecutive optical
evidence dates <=14 days apart (fusion.drainage_candidates, on the cached
melt-core evidence series). Classification against the radar record:

  DRAINAGE  - the radar kept the site wet after the loss (median last_wet >=
              event end + 3 d): the water left while melt continued - the
              hydrofracture-relevant signature.
  REFREEZE  - the radar froze in/around the loss window: the pond iced over
              (benign; the pale-lid class decision 0004 made a signal).

Writes output/sar/drainage_<season>.json.
"""
import json
import sys
from datetime import date
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
import fusion
import sar

DRAIN_MARGIN_D = 3


def run(season):
    z = np.load(sar.OUT / f"ponds_{season}_corem.npz", allow_pickle=False)
    dates = [date.fromordinal(int(o)) for o in z["ordinals"]]
    masks = list(z["masks"])
    comp = np.load(sar.OUT / f"season_{season}_t2.npz")
    lw = comp["last_wet"]
    y0 = int(season.split("-")[0])
    season_zero = date(y0, 11, 1).toordinal()

    events = fusion.drainage_candidates(dates, masks)
    for e in events:
        # the component's radar phenology, sampled at the event centroid's
        # neighbourhood (5x5) to ride out speckle
        r, c = e["row"], e["col"]
        win = lw[max(0, r - 2):r + 3, max(0, c - 2):c + 3]
        win = win[win >= 0]
        end_day = date.fromisoformat(e["date_after"]).toordinal() - season_zero
        med_lw = float(np.median(win)) if win.size else None
        e["median_last_wet_day"] = med_lw
        e["event_end_day"] = end_day
        e["class"] = ("drainage" if med_lw is not None
                      and med_lw >= end_day + DRAIN_MARGIN_D else "refreeze")
    out = sar.OUT / f"drainage_{season}.json"
    out.write_text(json.dumps({"season": season, "events": events,
                               "n_drainage": sum(e["class"] == "drainage" for e in events),
                               "n_refreeze": sum(e["class"] == "refreeze" for e in events)},
                              indent=1))
    print(f"[{season}] {len(events)} candidates -> "
          f"{sum(e['class']=='drainage' for e in events)} drainage, "
          f"{sum(e['class']=='refreeze' for e in events)} refreeze")


if __name__ == "__main__":
    run(sys.argv[1] if len(sys.argv) > 1 else "2019-20")

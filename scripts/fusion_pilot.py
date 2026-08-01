"""M3-B Task 4: melt-state pilot for one GVI season + report assembly.

    python scripts/fusion_pilot.py 2020-21     # build one season's state product
    python scripts/fusion_pilot.py report      # assemble fusion_report.json

Pilot checks (plan 2026-07-30, Task 4):
  (a) conflict rate < 5% per season
  (b) PONDED extent <= optical pond evidence extent (construction check)
  (c) ponded-extent season ranking matches the validated optical ranking
  (d) wet/ponded ratio per season - the new science number
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
import fusion
import sar
import shelf

GVI_BBOX = [-70.5, -72.6, -66.0, -70.0]
SEASONS = ("2019-20", "2020-21", "2021-22")
OPTICAL_KM2 = {"2019-20": 838.1, "2020-21": 159.1, "2021-22": 12.3}


def thresh():
    g = sar.OUT / "gates.json"
    return (json.loads(g.read_text()).get("gate1_best_thresh_db", 3.0)
            if g.exists() else 3.0)


def run_season(season):
    grid = shelf.build_fixed_grid()
    src = sar.s1_source()
    fusion.melt_state_season(season, grid, src, GVI_BBOX, thresh=thresh())


def report():
    t = thresh()
    out = {"thresh_db": t, "seasons": {}, "checks": {}}
    ponded_max = {}
    for season in SEASONS:
        s = json.loads((sar.OUT / f"state_{season}_t{t:g}.json").read_text())
        out["seasons"][season] = s
        ponded_max[season] = s["ponded_extent_max_km2"]

        # (b) PONDED <= optical evidence, cell-level, from the stored stacks
        z = np.load(sar.OUT / f"state_{season}_t{t:g}.npz", allow_pickle=False)
        pz = np.load(sar.OUT / f"ponds_{season}.npz", allow_pickle=False)
        ever_pond_optical = pz["masks"].any(axis=0)
        ponded_radar = z["days_ponded"] > 0
        outside = int((ponded_radar & ~ever_pond_optical).sum())
        out["seasons"][season]["ponded_cells_outside_optical"] = outside

    conflicts = {s: out["seasons"][s]["conflict_rate"] for s in SEASONS}
    out["checks"]["a_conflict_lt_5pct"] = {
        "rates": conflicts,
        "pass": all(c is not None and c < 0.05 for c in conflicts.values())}
    out["checks"]["b_ponded_subset_optical"] = {
        "cells_outside": {s: out["seasons"][s]["ponded_cells_outside_optical"]
                          for s in SEASONS},
        "pass": all(out["seasons"][s]["ponded_cells_outside_optical"] == 0
                    for s in SEASONS)}
    optical_order = sorted(SEASONS, key=lambda s: -OPTICAL_KM2[s])
    fused_order = sorted(SEASONS, key=lambda s: -ponded_max[s])
    out["checks"]["c_ranking"] = {
        "optical_km2": OPTICAL_KM2, "ponded_max_km2": ponded_max,
        "optical_order": optical_order, "fused_order": fused_order,
        "pass": optical_order == fused_order}
    out["checks"]["d_wet_to_ponded"] = {
        s: round(out["seasons"][s]["wet_extent_max_km2"]
                 / max(ponded_max[s], 0.1), 1) for s in SEASONS}
    (sar.OUT / "fusion_report.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out["checks"], indent=1))


if __name__ == "__main__":
    arg = sys.argv[1] if len(sys.argv) > 1 else "report"
    report() if arg == "report" else run_season(arg)

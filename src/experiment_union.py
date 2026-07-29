"""EXPERIMENT: is clean-scene UNION a better seasonal maximum?

Hypothesis: unioning detections across every qualifying (genuinely clear) scene
captures ponds that peak on different dates in different parts of a tile, which
a single peak scene misses - a recall gain with no haze risk, because only
scenes the validated two-branch selector already trusts contribute.

Acceptance gates (checked against the anchors before this ever becomes default):
  2020-21  single=142.5, blind-validated 19DEA=127, HT true-water band 120-194
           -> union must stay inside ~145-195 (recall gain, not inflation)
  2024-25  single=98 (cloud-robust)   -> union should stay ~98-140
  2019-20  single=~840 (record)       -> union should stay ~840-1000; a huge
           jump would mean the metadata-gated pool is admitting haze after all

Writes a verdict per season; adopting the mode stays a HUMAN decision.

Run:  python src/experiment_union.py
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import melt
import shelf

ANCHORS = {
    "2020-21": (145.0, 195.0),
    "2024-25": (95.0, 140.0),
    "2019-20": (830.0, 1000.0),
}

OUT = melt.ROOT / "output" / "shelf" / "union_experiment.json"


def main():
    shelf.UNION_CLEAN = True
    grid = shelf.build_fixed_grid()
    results = {}
    for season, (lo, hi) in ANCHORS.items():
        print(f"\n=== union experiment: {season} (accept {lo}-{hi} km2) ===")
        r = shelf.run_season(season, grid=grid)
        km2 = r["shelf_km2"]
        verdict = "PASS" if lo <= km2 <= hi else "FAIL"
        results[season] = {"union_km2": km2, "accept": [lo, hi], "verdict": verdict}
        print(f">>> {season} union={km2} km2  [{verdict}]")
        OUT.write_text(json.dumps(results, indent=1))

    n_pass = sum(1 for v in results.values() if v["verdict"] == "PASS")
    print(f"\nUNION EXPERIMENT: {n_pass}/{len(results)} anchors pass -> "
          f"{'candidate for adoption (human review)' if n_pass == len(results) else 'DO NOT ADOPT'}")


if __name__ == "__main__":
    main()

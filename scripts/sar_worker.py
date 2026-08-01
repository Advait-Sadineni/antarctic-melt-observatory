"""One-composite SAR worker so independent seasons can run in parallel.

    python scripts/sar_worker.py 2021-22

Writes the same npz/json cache season_composite always writes (checkpointed,
resumable); sar_gates.py gate23 then cache-hits it. Never run two workers on
the SAME season - they would race on one checkpoint file.
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
import sar
import shelf

GVI_BBOX = [-70.5, -72.6, -66.0, -70.0]

season = sys.argv[1]
grid = shelf.build_fixed_grid()
sgrid = sar.sar_grid(grid)
src = sar.s1_source()
g = sar.OUT / "gates.json"
thresh = (json.loads(g.read_text()).get("gate1_best_thresh_db", 3.0)
          if g.exists() else 3.0)
print(f"[worker] {season} at t={thresh}", flush=True)
sar.season_composite(season, sgrid, src, GVI_BBOX, thresh=thresh)
print(f"[worker done] {season}", flush=True)

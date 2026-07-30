"""Daytime autonomous batch (2026-07-30, user back ~6pm).

Sequential, checkpointed, detached-safe:
  1. Larsen B: find dominant 2024-25 tile, generate blind-validation chips
     (fast first, so labelling can happen mid-day).
  2. George VI: rerun 9-season history -> data products now exist for it.
  3. Bach Stange Wilkins LarsenC Amery: rerun -> products for every shelf.
     (Numbers should reproduce the stored records - doubles as a
     reproducibility check; per-season saves mean partial progress survives.)

Run:  python scripts/day_batch.py   (detached, logs to output/day_batch.log)
"""
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))


def step1_larsenb_chips():
    import melt
    import shelf
    import shelves
    print("=== STEP 1: LarsenB dominant tile + chips ===", flush=True)
    cfg = shelves.prepare("LarsenB")
    shelf.set_shelf(cfg["name"], cfg["tiles"], cfg["boundary"])
    best_tile, best_water = None, -1
    for t in cfg["tiles"]:
        try:
            it = shelf.production_scene(t, "2025-01-01", "2025-02-28")
            if it is None:
                continue
            sc = shelf._coarse_shelf_mask(it)
            w = shelf._onshelf_coarse_water(it, sc)
            print(f"  {t}: coarse water {w}", flush=True)
            if w > best_water:
                best_tile, best_water = t, w
        except Exception as e:
            print(f"  {t}: skip {type(e).__name__}", flush=True)
    print(f"  dominant: {best_tile}", flush=True)
    subprocess.run([sys.executable, "-u", "src/validate_shelf_tile.py",
                    best_tile, "make", "LarsenB", "2024-25"],
                   cwd=ROOT, check=False)
    (ROOT / "output" / "day_batch_state.json").write_text(
        json.dumps({"larsenb_tile": best_tile}))


def step2_gvi_products():
    print("=== STEP 2: George VI product rerun ===", flush=True)
    subprocess.run([sys.executable, "-u", "src/shelf.py", "history"],
                   cwd=ROOT, check=False)


def step3_remaining_products():
    print("=== STEP 3: remaining shelves product rerun ===", flush=True)
    subprocess.run([sys.executable, "-u", "src/shelves.py",
                    "Bach", "Stange", "Wilkins", "LarsenC", "Amery"],
                   cwd=ROOT, check=False)


if __name__ == "__main__":
    step1_larsenb_chips()
    step2_gvi_products()
    step3_remaining_products()
    print("=== DAY BATCH COMPLETE ===", flush=True)

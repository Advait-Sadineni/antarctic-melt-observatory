"""M1 acceptance gate: George VI 2020-21 via Planetary Computer must match
the blind-validated Earth Search result within 3%.

Networked + slow (~10 min); run manually, not part of pytest:
    python scripts/regression_pc.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
import melt
import shelf
from core.pc import PlanetaryComputerSource

FROZEN = 142.48          # 2020-21, blind-validated (19DEA = 127.17)
TOL = 0.03

melt.set_source(PlanetaryComputerSource())
grid = shelf.build_fixed_grid()          # cached grid; footprints identical
r = shelf.run_season("2020-21", grid=grid)
got = r["shelf_km2"]
dev = abs(got - FROZEN) / FROZEN
print(f"\nPC 2020-21 = {got} km2 vs frozen {FROZEN} (dev {100*dev:.2f}%)")
print("REGRESSION:", "PASS" if dev <= TOL else "FAIL")
sys.exit(0 if dev <= TOL else 1)

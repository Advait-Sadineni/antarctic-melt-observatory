"""Compile observatory outputs into dashboard/data.js (static, no server).

Run after any history/gates/validation update:
    python scripts/build_dashboard.py
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "dashboard"

SHELVES = {  # data file stem -> display name (history.json is George VI)
    "history": "George VI",
    "history_Amery": "Amery",
    "history_Bach": "Bach",
    "history_Baudouin": "Baudouin",
    "history_Fimbul": "Fimbul",
    "history_LarsenB": "Larsen B",
    "history_LarsenC": "Larsen C",
    "history_Nivl": "Nivl",
    "history_Riiser-Larsen": "Riiser-Larsen",
    "history_Shackleton": "Shackleton",
    "history_Stange": "Stange",
    "history_Wilkins": "Wilkins",
}


def main():
    shelves = {}
    for stem, name in SHELVES.items():
        p = ROOT / "output" / "shelf" / f"{stem}.json"
        if not p.exists():
            print(f"  [skip] {p.name} missing")
            continue
        rows = json.loads(p.read_text())
        shelves[name] = {
            "area_km2": rows[0]["shelf_area_km2"] if rows else None,
            "seasons": [{
                "season": r["season"],
                "km2": r["shelf_km2"],
                "pct": round(100 * r["shelf_km2"] / r["shelf_area_km2"], 3),
                "obs_cloud": r.get("obs_cloud"),
                "poorly": bool(r.get("poorly_observed")),
            } for r in rows],
        }

    corr = json.loads((ROOT / "reference" / "regional_corrections.json").read_text())
    gates_p = ROOT / "output" / "sar" / "gates.json"
    gates = json.loads(gates_p.read_text()) if gates_p.exists() else {}
    fusion_p = ROOT / "output" / "sar" / "fusion_report.json"
    fusion = json.loads(fusion_p.read_text()) if fusion_p.exists() else None

    # per-season fused state summaries (wet vs ponded - the M3 science number)
    states = {}
    for p in sorted((ROOT / "output" / "sar").glob("state_*_t2.json")):
        s = json.loads(p.read_text())
        states[s["season"]] = s

    calib_p = ROOT / "output" / "depth" / "calibration.json"
    calib = json.loads(calib_p.read_text()) if calib_p.exists() else None
    depth_products = {}
    for p in sorted((ROOT / "output" / "depth").glob("product_*.json")):
        d = json.loads(p.read_text())
        depth_products[d["season"]] = d

    data = {"shelves": shelves, "corrections": corr, "gates": gates,
            "fusion": fusion, "states": states,
            "calibration": calib, "depth_products": depth_products,
            "doi": "10.5281/zenodo.21711608",
            "repo": "https://github.com/Advait-Sadineni/antarctic-melt-observatory"}
    OUT.mkdir(exist_ok=True)
    js = "window.OBS = " + json.dumps(data) + ";\n"
    (OUT / "data.js").write_text(js)
    n = sum(len(s["seasons"]) for s in shelves.values())
    print(f"dashboard/data.js: {len(shelves)} shelves, {n} shelf-seasons")


if __name__ == "__main__":
    sys.exit(main())

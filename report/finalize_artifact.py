"""Inject the finished 9-season history into the artifact's DATA array.

Reads output/shelf/history.json and rewrites the `const DATA = [...]` block in
report/george_vi_meltwater.html so the published chart always matches the run.
Confidence flag: a season is 'prov' (cloud-limited, lower confidence) when the
worst water-bearing tile's actual cloud fraction (obs_cloud) exceeds 15%.
"""
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
HIST = ROOT / "output" / "shelf" / "history.json"
HTML = ROOT / "report" / "george_vi_meltwater.html"

SHELF = 14390.6
LABELS = {  # nice en-dash season labels
    s: s.replace("-", "–") for s in
    ["2017-18", "2018-19", "2019-20", "2020-21", "2021-22",
     "2022-23", "2023-24", "2024-25", "2025-26"]
}


def main():
    hist = {r["season"]: r for r in json.loads(HIST.read_text())}
    rows = [hist[s] for s in LABELS if s in hist]
    record_km2 = max(r["shelf_km2"] for r in rows)

    lines = []
    for r in rows:
        s = r["season"]
        cloud = r.get("obs_meta", r.get("obs_cloud", 0.0)) or 0.0
        prov = bool(r.get("poorly_observed", False))
        flags = f'prov: {"true " if prov else "false"}'
        if r["shelf_km2"] == record_km2:
            flags += ", record: true"
        if s == "2020-21":
            flags += ", validated: true"
        lines.append(
            f'    {{ season: "{LABELS[s]}", km2: {r["shelf_km2"]:.1f}, '
            f'cloud: {cloud:.1f}, tiles: {r["tiles"]}, {flags} }},'
        )
    block = "const DATA = [\n" + "\n".join(lines) + "\n  ];"

    html = HTML.read_text(encoding="utf-8")
    html = re.sub(r"const DATA = \[.*?\];", block, html, count=1, flags=re.S)
    html = re.sub(r"const SHELF = [\d.]+;", f"const SHELF = {SHELF};", html)
    HTML.write_text(html, encoding="utf-8")
    print("updated DATA with", len(rows), "seasons; record =", record_km2, "km2")
    for l in lines:
        print(l)


if __name__ == "__main__":
    main()

"""
Antarctic surface meltwater detection - v1 demo.

Downloads (or loads from cache) one Sentinel-2 L2A scene over George VI Ice
Shelf, computes NDWI, thresholds it into a meltwater mask, and writes a
side-by-side PNG of true colour vs. detected ponds.

Everything here is hardcoded on purpose. See README for what's stubbed.

Run:  python src/demo.py
"""

import os
from pathlib import Path

import numpy as np
import rasterio
from rasterio.windows import Window

os.environ.setdefault("GDAL_DISABLE_READDIR_ON_OPEN", "EMPTY_DIR")
os.environ.setdefault("AWS_NO_SIGN_REQUEST", "YES")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

# --- HARDCODED SETTINGS (v1) -------------------------------------------------
# Scene: Sentinel-2B, tile 19CEV, 2021-01-24, 0% reported cloud.
# Peak austral summer melt over George VI Ice Shelf, Antarctic Peninsula.
# Hosted as public Cloud-Optimized GeoTIFFs on AWS Open Data - no account needed.
ITEM_ID = "S2B_19CEV_20210124_0_L2A"
STAC_API = "https://earth-search.aws.element84.com/v1"

# Sub-window of the 110 km tile, chosen because it has the densest ponding.
# Origin is in full-resolution (10 m) pixels from the tile's top-left corner.
WIN_ROW, WIN_COL, WIN_SIZE = 1280, 800, 2048  # 2048 px * 10 m = ~20.5 km square

NDWI_THRESHOLD = 0.16  # TUNABLE. Published Antarctic pond studies use 0.10-0.25.
BRIGHTNESS_FLOOR = 3000  # TUNABLE. Rejects dark open ocean / deep shadow.
PIXEL_M = 10.0  # Sentinel-2 B03/B08 native resolution

ROOT = Path(__file__).resolve().parent.parent
CACHE = ROOT / "data" / f"{ITEM_ID}_{WIN_ROW}_{WIN_COL}_{WIN_SIZE}.tif"
OUT_PNG = ROOT / "output" / "george_vi_meltwater.png"

# Sentinel-2 processing baseline 04.00 (2022-01-25 onward) added a -1000 DN
# offset to L2A reflectance. Our scene predates it, so no correction needed.
BOA_OFFSET = 0

BANDS = ["red", "green", "blue", "nir"]  # B04, B03, B02, B08 - all 10 m


def fetch_window():
    """Read a 4-band window straight out of the remote COGs and cache it."""
    from pystac_client import Client

    print(f"[fetch] querying {STAC_API} for {ITEM_ID}")
    item = next(Client.open(STAC_API).search(collections=["sentinel-2-l2a"], ids=[ITEM_ID]).items())

    win = Window(WIN_COL, WIN_ROW, WIN_SIZE, WIN_SIZE)
    stack, profile = [], None
    for name in BANDS:
        href = item.assets[name].href
        print(f"[fetch] {name:5s} <- {href.rsplit('/', 1)[-1]}")
        with rasterio.open(href) as src:
            stack.append(src.read(1, window=win))
            if profile is None:
                profile = src.profile.copy()
                profile.update(
                    count=len(BANDS),
                    height=WIN_SIZE,
                    width=WIN_SIZE,
                    transform=src.window_transform(win),
                    compress="deflate",
                )

    CACHE.parent.mkdir(parents=True, exist_ok=True)
    with rasterio.open(CACHE, "w", **profile) as dst:
        for i, arr in enumerate(stack, start=1):
            dst.write(arr, i)
        dst.descriptions = tuple(BANDS)
    print(f"[fetch] cached -> {CACHE.relative_to(ROOT)}")
    return np.stack(stack)


def load():
    if CACHE.exists():
        print(f"[cache] loading {CACHE.relative_to(ROOT)}")
        with rasterio.open(CACHE) as src:
            return src.read()
    return fetch_window()


def stretch_rgb(red, green, blue, lo=2, hi=98):
    """Percentile stretch to 0-1 for display.

    The three bands share one set of limits on purpose. Stretching each band
    independently would rebalance them against each other and turn the snow
    a false peach colour; a common scale keeps the band ratios intact, so
    ice stays white and ponds stay blue.
    """
    rgb = np.dstack([red, green, blue])
    sample = rgb[rgb > 0]
    a, b = np.percentile(sample, [lo, hi]) if sample.size else (0, 1)
    return np.clip((rgb - a) / max(b - a, 1e-6), 0, 1)


def main():
    red, green, blue, nir = (load().astype("f4") - BOA_OFFSET)

    valid = (green > 0) & (nir > 0)

    # NDWI: water absorbs near-infrared strongly but reflects green, so
    # (green - NIR) / (green + NIR) goes positive over water and negative
    # over snow and ice. Ponds on a white ice shelf separate cleanly.
    ndwi = np.where(valid, (green - nir) / np.maximum(green + nir, 1e-6), np.nan)

    ponds = (ndwi > NDWI_THRESHOLD) & valid & (green > BRIGHTNESS_FLOOR)

    n_px = int(ponds.sum())
    area_km2 = n_px * (PIXEL_M**2) / 1e6
    scene_km2 = int(valid.sum()) * (PIXEL_M**2) / 1e6

    print()
    print(f"  scene        {ITEM_ID}")
    print(f"  window       {WIN_SIZE}x{WIN_SIZE} px @ {PIXEL_M:.0f} m  ({scene_km2:.1f} km2 valid)")
    print(f"  NDWI thresh  {NDWI_THRESHOLD}")
    print(f"  pond pixels  {n_px:,}")
    print(f"  pond area    {area_km2:.3f} km2  ({100 * n_px / max(int(valid.sum()), 1):.2f}% of scene)")
    print()

    rgb = stretch_rgb(red, green, blue)

    fig, axes = plt.subplots(1, 2, figsize=(16, 8.6))
    axes[0].imshow(rgb)
    axes[0].set_title("Sentinel-2 true colour (B4/B3/B2)")

    axes[1].imshow(rgb * 0.35 + 0.15)  # dimmed backdrop for context
    axes[1].imshow(
        np.ma.masked_where(~ponds, ponds),
        cmap=ListedColormap(["#00b3ff"]),
        interpolation="nearest",
    )
    axes[1].set_title(f"Meltwater mask (NDWI > {NDWI_THRESHOLD})")

    for ax in axes:
        ax.set_xticks([])
        ax.set_yticks([])

    fig.suptitle(
        f"George VI Ice Shelf  -  {ITEM_ID.split('_')[2]}  -  "
        f"{n_px:,} pond pixels  =  {area_km2:.2f} km$^2$",
        fontsize=14,
    )
    fig.tight_layout()
    OUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, dpi=130, bbox_inches="tight")
    print(f"[write] {OUT_PNG.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

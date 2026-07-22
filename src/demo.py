"""
Single-scene meltwater demo: one Sentinel-2 scene -> side-by-side PNG.

Run:  python src/demo.py
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap

import melt

ITEM_ID = "S2B_19CEV_20210124_0_L2A"  # 2021-01-24, cloud-free, peak melt
OUT_PNG = melt.ROOT / "output" / "george_vi_meltwater.png"


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
    print(f"[scene] {ITEM_ID}")
    item = melt.get_item(ITEM_ID)
    bands = melt.load_scene(item)
    scl = melt.read_scl(item)
    scr = melt.screen(scl)

    ndwi, ponds, valid = melt.detect(bands["green"], bands["nir"], scr["reject_mask"])
    st = melt.pond_stats(ponds, valid)

    print()
    print(f"  window       {melt.WIN_SIZE}x{melt.WIN_SIZE} px @ {melt.PIXEL_M:.0f} m"
          f"  ({st['valid_km2']:.1f} km2 valid)")
    print(f"  BOA offset   {melt.boa_offset(item)} DN")
    print(f"  cloud in AOI {scr['cloud_frac'] * 100:.2f}%  (SCL-derived)")
    print(f"  NDWI thresh  {melt.NDWI_THRESHOLD}")
    print(f"  pond pixels  {st['pond_px']:,}")
    print(f"  pond area    {st['pond_km2']:.3f} km2  ({st['pond_pct']:.2f}% of scene)")
    print()

    rgb = stretch_rgb(bands["red"], bands["green"], bands["blue"])

    fig, axes = plt.subplots(1, 2, figsize=(16, 8.6))
    axes[0].imshow(rgb)
    axes[0].set_title("Sentinel-2 true colour (B4/B3/B2)")

    axes[1].imshow(rgb * 0.35 + 0.15)  # dimmed backdrop for context
    axes[1].imshow(
        np.ma.masked_where(~ponds, ponds),
        cmap=ListedColormap(["#00b3ff"]),
        interpolation="nearest",
    )
    axes[1].set_title(f"Meltwater mask (NDWI > {melt.NDWI_THRESHOLD})")

    for ax in axes:
        ax.set_xticks([])
        ax.set_yticks([])

    fig.suptitle(
        f"George VI Ice Shelf  -  {ITEM_ID.split('_')[2]}  -  "
        f"{st['pond_px']:,} pond pixels  =  {st['pond_km2']:.2f} km$^2$",
        fontsize=14,
    )
    fig.tight_layout()
    OUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, dpi=130, bbox_inches="tight")
    print(f"[write] {OUT_PNG.relative_to(melt.ROOT)}")


if __name__ == "__main__":
    main()

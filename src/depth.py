"""M4: physically-based pond depths and meltwater volume (spec v1.1).

Light attenuates exponentially with water depth; inverting Beer-Lambert over a
pond gives per-pixel depth:  z = [ln(A_d - R_inf) - ln(R_w - R_inf)] / g
(Pope et al. 2016 lineage). Retrieved independently for red and green;
agreement is the confidence signal. g is CALIBRATED against ICESat-2 photon
depths (src/icesat2.py) - never hand-picked.
"""
import numpy as np
from scipy import ndimage

Z_MAX = 10.0          # m; physical clip for shallow supraglacial ponds
DUAL_TOL = 0.30       # bands must agree within 30% of the larger estimate


def depth_single(R_w, A_d, R_inf, g):
    """Beer-Lambert inversion, NaN-safe. R_w outside (R_inf, A_d) -> NaN
    (physically impossible: brighter than dry bottom or darker than deep
    water means the model does not apply at that pixel)."""
    R_w = np.asarray(R_w, "f4")
    with np.errstate(invalid="ignore", divide="ignore"):
        z = (np.log(A_d - R_inf) - np.log(R_w - R_inf)) / g
    bad = ~np.isfinite(z) | (R_w <= R_inf) | (R_w >= A_d) | (z < 0)
    z = np.where(bad, np.nan, np.clip(z, 0.0, Z_MAX))
    return z.astype("f4")


def dual_band_merge(z_red, z_green, tol=DUAL_TOL):
    """Mean of the two estimates where they agree within tol of the larger;
    disagreement or a missing band -> NaN + uncertain flag. Two estimates beat
    one, and their disagreement maps where the physics breaks."""
    z_red = np.asarray(z_red, "f4")
    z_green = np.asarray(z_green, "f4")
    both = np.isfinite(z_red) & np.isfinite(z_green)
    larger = np.maximum(z_red, z_green)
    with np.errstate(invalid="ignore", divide="ignore"):
        agree = both & (np.abs(z_red - z_green) <= tol * larger)
    z = np.where(agree, (z_red + z_green) / 2.0, np.nan).astype("f4")
    uncertain = ~agree
    return z, uncertain


def rim_albedo(band_img, pond_mask, rim_px=2):
    """Per-pond bottom-albedo estimate from the ice ring just outside each
    pond (never a global constant - dirty/dark floors are the known weak
    link). Returns (labels, A_d image aligned to ponds)."""
    labels, n = ndimage.label(pond_mask)
    ad = np.full(band_img.shape, np.nan, "f4")
    struct = ndimage.generate_binary_structure(2, 2)
    for i in range(1, n + 1):
        pond = labels == i
        rim = ndimage.binary_dilation(pond, struct, iterations=rim_px) & ~pond_mask
        if rim.any():
            ad[pond] = float(np.nanmean(band_img[rim]))
    return labels, ad


def mc_volume(z_fn, pixel_area_m2, n=1000, seed=20260731,
              g_sigma=0.1, ad_sigma=0.02, rinf_range=(0.03, 0.07)):
    """Monte Carlo volume: jointly perturb calibration (g), rim albedo shift,
    and R_inf; z_fn(g_shift, ad_shift, rinf) recomputes the depth map. Returns
    median and 16-84th percentile volume in km3."""
    rng = np.random.default_rng(seed)
    vols = np.empty(n, "f8")
    for k in range(n):
        z = z_fn(rng.normal(0.0, g_sigma), rng.normal(0.0, ad_sigma),
                 rng.uniform(*rinf_range))
        vols[k] = np.nansum(z) * pixel_area_m2 / 1e9
    return {"volume_km3_median": float(np.median(vols)),
            "p16": float(np.percentile(vols, 16)),
            "p84": float(np.percentile(vols, 84)),
            "n_draws": n}

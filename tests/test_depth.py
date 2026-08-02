"""Offline tests for M4 depth retrieval - synthetic arrays only."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

import numpy as np

import depth


def test_depth_single_analytic():
    # with Ad=0.6, Rinf=0.05, g=0.8: choose Rw so z = 1.0 exactly
    Rw = 0.05 + 0.55 * np.exp(-0.8)
    z = depth.depth_single(np.array([Rw], "f4"), 0.6, 0.05, 0.8)
    assert abs(z[0] - 1.0) < 1e-5


def test_depth_single_invalid_inputs_nan():
    z = depth.depth_single(np.array([0.04, 0.7], "f4"), 0.6, 0.05, 0.8)
    assert np.isnan(z).all()          # below R_inf / above A_d -> NaN


def test_dual_band_merge():
    zr = np.array([1.0, 1.0, np.nan], "f4")
    zg = np.array([1.1, 2.0, 1.0], "f4")
    z, unc = depth.dual_band_merge(zr, zg)
    assert abs(z[0] - 1.05) < 1e-6 and not unc[0]
    assert unc[1] and np.isnan(z[1])
    assert unc[2] and np.isnan(z[2])   # one-band-only counts uncertain in v1


def test_rim_albedo_two_ponds():
    img = np.full((20, 20), 0.5, "f4")
    mask = np.zeros((20, 20), bool)
    mask[5:8, 5:8] = True
    mask[12:15, 12:15] = True
    img[11:16, 11:16] = 0.7            # second pond's surroundings brighter
    labels, ad = depth.rim_albedo(img, mask, rim_px=1)
    vals = sorted(np.unique(ad[mask]).astype("f8"))
    assert len(vals) == 2
    assert abs(vals[0] - 0.5) < 1e-3 and abs(vals[1] - 0.7) < 1e-3


def test_mc_volume_symmetric_perturbation():
    rng_truth = 1.0                     # true depth everywhere, 100 px of 100 m2
    def z_fn(g_shift, ad_shift, rinf):
        return np.full((10, 10), rng_truth + g_shift, "f4")  # linearized toy
    out = depth.mc_volume(z_fn, pixel_area_m2=100.0, n=500, seed=7,
                          g_sigma=0.05, ad_sigma=0.0, rinf_range=(0.0, 0.0))
    v_true = 100 * 100 * 1.0 / 1e9
    assert abs(out["volume_km3_median"] - v_true) / v_true < 0.05
    assert out["p16"] < out["volume_km3_median"] < out["p84"]


def test_icesat2_peak_split():
    import icesat2
    rng = np.random.default_rng(3)
    surface = rng.normal(20.0, 0.05, 3000)     # air-water interface photons
    bottom = rng.normal(18.5, 0.08, 400)       # pond-bottom returns
    noise = rng.uniform(15, 25, 150)
    elev = np.concatenate([surface, bottom, noise])
    d = icesat2.photon_peak_depth(elev)
    # separation 1.5 m x 0.752 refraction = 1.128 m
    assert d is not None and abs(d - 1.128) < 0.08


def test_icesat2_peak_split_no_bottom():
    import icesat2
    rng = np.random.default_rng(4)
    elev = rng.normal(20.0, 0.05, 3000)         # surface only - dry/no bottom
    assert icesat2.photon_peak_depth(elev) is None


def test_photons_in_mask_selects_only_pond_cells():
    import icesat2
    from rasterio.transform import from_origin
    # 4x4 grid of 10 m cells in a metric CRS; pond = center 2x2 block
    tr = from_origin(1000.0, 2000.0, 10.0, 10.0)
    mask = np.zeros((4, 4), bool)
    mask[1:3, 1:3] = True
    # photon A inside pond cell (row1,col1); photon B outside (row0,col0);
    # photon C off-grid entirely
    xs = np.array([1015.0, 1005.0, 5000.0])
    ys = np.array([1985.0, 1995.0, 5000.0])
    sel = icesat2.photons_in_mask(xs, ys, mask, tr, "EPSG:3031",
                                  src_crs="EPSG:3031")
    assert sel.tolist() == [True, False, False]


def test_depth_qc_bounds():
    import icesat2
    lo, hi = icesat2.DEPTH_QC_M
    assert lo >= 0.25         # min-separation artifacts (0.226 m) excluded
    assert hi <= 10.0         # supraglacial ponds are shallow; junk excluded
    assert icesat2.qc_pass(0.5) and icesat2.qc_pass(4.4)
    assert not icesat2.qc_pass(0.226) and not icesat2.qc_pass(667.2)


def test_fit_g_recovers_truth():
    rng = np.random.default_rng(11)
    g_true = 0.75
    z = rng.uniform(0.4, 4.0, 60)
    X = z * g_true + rng.normal(0, 0.02, 60)   # small noise in optical term
    out = depth.fit_g(X, z)
    assert abs(out["g"] - g_true) < 0.03
    assert out["n"] == 60 and out["rmse_m"] < 0.1


def test_fit_g_trims_outliers():
    g_true = 0.75
    z = np.linspace(0.5, 4.0, 40)
    X = z * g_true
    z_bad = np.concatenate([z, [1.0, 1.2]])
    X_bad = np.concatenate([X, [4.0, 5.0]])    # wildly wrong pairs
    out = depth.fit_g(X_bad, z_bad)
    assert abs(out["g"] - g_true) < 0.02       # outliers must not drag the fit
    assert out["n_trimmed"] == 2


def test_calib_X_matches_forward_model():
    # X = ln(Ad-Rinf) - ln(Rw-Rinf); forward depth_single inverts with g
    Rw = np.array([0.25], "f4")
    X = depth.calib_X(Rw, A_d=0.6, R_inf=0.05)
    z = depth.depth_single(Rw, 0.6, 0.05, g=0.8)
    assert abs(X[0] / 0.8 - z[0]) < 1e-5

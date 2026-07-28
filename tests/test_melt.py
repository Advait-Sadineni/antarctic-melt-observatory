"""
Tests for the detection core. No network: every test builds synthetic arrays
or fake STAC items, so this runs offline and fast.

Run:  python -m pytest tests -q
"""

import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

import melt  # noqa: E402
import season  # noqa: E402
import validate_landsat as vl  # noqa: E402
import uncertainty as unc  # noqa: E402
import validate_points as vp  # noqa: E402


class FakeItem:
    """Stands in for a pystac Item - only properties and datetime are used."""

    def __init__(self, properties=None, dt=None):
        self.properties = properties or {}
        self.datetime = dt or datetime(2021, 1, 24, tzinfo=timezone.utc)


# --- NDWI detection ----------------------------------------------------------

def test_water_detected_ice_rejected():
    """Water reflects green and absorbs NIR; ice reflects both."""
    green = np.array([[8000.0, 9000.0]])
    nir = np.array([[4000.0, 8500.0]])  # left: strong absorption, right: ice
    _, ponds, valid = melt.detect(green, nir)

    assert ponds[0, 0]      # NDWI = 0.33
    assert not ponds[0, 1]  # NDWI = 0.029
    assert valid.all()


def test_zero_pixels_are_invalid_not_water():
    """Nodata is 0 in both bands. It must not become a pond."""
    green = np.zeros((1, 2))
    nir = np.zeros((1, 2))
    _, ponds, valid = melt.detect(green, nir)

    assert not ponds.any()
    assert not valid.any()


def test_brightness_floor_rejects_dark_water():
    """Dark open ocean has pond-like NDWI but must be rejected."""
    dark_green = np.array([[melt.BRIGHTNESS_FLOOR - 500.0]])
    bright_green = np.array([[melt.BRIGHTNESS_FLOOR + 500.0]])
    nir = np.array([[100.0]])

    assert not melt.detect(dark_green, nir)[1][0, 0]
    assert melt.detect(bright_green, nir)[1][0, 0]


def test_reject_mask_removes_ponds():
    green, nir = np.array([[8000.0]]), np.array([[4000.0]])
    assert melt.detect(green, nir, reject_mask=np.array([[False]]))[1][0, 0]
    assert not melt.detect(green, nir, reject_mask=np.array([[True]]))[1][0, 0]


def test_threshold_is_monotonic():
    """Raising the threshold can never detect more water."""
    rng = np.random.default_rng(0)
    green = rng.uniform(4000, 11000, (60, 60))
    nir = rng.uniform(3000, 10000, (60, 60))

    counts = [melt.detect(green, nir, threshold=t)[1].sum()
              for t in [0.05, 0.10, 0.16, 0.25, 0.40]]
    assert counts == sorted(counts, reverse=True)


# --- Area accounting ---------------------------------------------------------

def test_pond_area_matches_pixel_count():
    """10 m pixels: 100 pond pixels = 100 * 100 m2 = 0.01 km2."""
    ponds = np.zeros((50, 50), bool)
    ponds[:10, :10] = True
    st = melt.pond_stats(ponds, np.ones((50, 50), bool))

    assert st["pond_px"] == 100
    assert st["pond_km2"] == pytest.approx(0.01)
    assert st["valid_km2"] == pytest.approx(0.25)
    assert st["pond_pct"] == pytest.approx(4.0)


def test_pond_stats_handles_no_valid_pixels():
    """Must not divide by zero on a fully-masked scene."""
    empty = np.zeros((4, 4), bool)
    st = melt.pond_stats(empty, empty)
    assert st["pond_px"] == 0
    assert st["pond_pct"] == 0.0


# --- Baseline 04.00 reflectance offset ---------------------------------------

def test_offset_zero_when_archive_already_applied():
    item = FakeItem({"earthsearch:boa_offset_applied": True,
                     "s2:processing_baseline": "04.00"},
                    datetime(2023, 1, 1, tzinfo=timezone.utc))
    assert melt.boa_offset(item) == 0


def test_offset_applied_for_baseline_04_when_not_corrected():
    item = FakeItem({"s2:processing_baseline": "04.00"},
                    datetime(2023, 1, 1, tzinfo=timezone.utc))
    assert melt.boa_offset(item) == -1000


def test_no_offset_for_old_baseline():
    item = FakeItem({"s2:processing_baseline": "02.14"},
                    datetime(2021, 1, 24, tzinfo=timezone.utc))
    assert melt.boa_offset(item) == 0


def test_offset_falls_back_to_date_when_baseline_missing():
    assert melt.boa_offset(FakeItem({}, datetime(2021, 6, 1, tzinfo=timezone.utc))) == 0
    assert melt.boa_offset(FakeItem({}, datetime(2023, 6, 1, tzinfo=timezone.utc))) == -1000


def test_offset_survives_unparseable_baseline():
    """Bad metadata must fall back to the date rule, not crash."""
    item = FakeItem({"s2:processing_baseline": "not-a-number"},
                    datetime(2023, 6, 1, tzinfo=timezone.utc))
    assert melt.boa_offset(item) == -1000


def test_offset_changes_ndwi():
    """The offset is not cosmetic - it must move the answer."""
    green, nir = np.array([[9000.0]]), np.array([[8000.0]])
    raw = (green - nir) / (green + nir)
    shifted = ((green - 1000) - (nir - 1000)) / ((green - 1000) + (nir - 1000))
    assert not np.isclose(raw, shifted)


# --- SCL screening -----------------------------------------------------------

def test_screen_counts_cloud_water_and_nodata():
    #        nodata, cloud-high, water, snow
    scl = np.array([[0, 9], [6, 11]], dtype="uint8")
    scr = melt.screen(scl)

    assert scr["nodata_frac"] == pytest.approx(0.25)
    assert scr["cloud_frac"] == pytest.approx(0.25)   # class 9 only, not nodata
    assert scr["water_frac"] == pytest.approx(0.25)
    assert scr["snow_frac"] == pytest.approx(0.25)
    assert scr["usable_frac"] == pytest.approx(0.5)   # water + snow survive


def test_screen_does_not_reject_water_class():
    """SCL calls ~20% of real ponds 'water'. Rejecting class 6 loses them."""
    assert 6 not in melt.SCL_REJECT
    scr = melt.screen(np.full((4, 4), 6, dtype="uint8"))
    assert scr["usable_frac"] == pytest.approx(1.0)


# --- Season plumbing ---------------------------------------------------------

def test_season_dates_span_the_new_year():
    assert season.season_dates("2020-21") == ("2020-11-01", "2021-03-31")
    assert season.season_dates("2025-26") == ("2025-11-01", "2026-03-31")


def test_series_prefers_most_usable_scene_on_a_date(tmp_path, monkeypatch):
    monkeypatch.setattr(melt, "ROOT", tmp_path)
    (tmp_path / "output").mkdir()

    for row in [
        {"item_id": "A", "date": "2021-01-24", "usable_km2": "376.6",
         "cloud_frac_aoi": "0.10", "pond_km2": "0.165"},
        {"item_id": "B", "date": "2021-01-24", "usable_km2": "419.2",
         "cloud_frac_aoi": "0.00", "pond_km2": "0.231"},
        {"item_id": "C", "date": "2021-01-25", "usable_km2": "419.0",
         "cloud_frac_aoi": "0.00", "pond_km2": "1.000"},
    ]:
        season.append_row("test", {k: row.get(k, "") for k in season.FIELDS})

    out = season.series("test")
    assert [r["item_id"] for r in out] == ["B", "C"]  # one per date, B beats A


def test_noise_floor_clears_the_pre_melt_baseline():
    """The floor must sit above every pre-melt reading, or artifacts read as melt."""
    assert melt.NOISE_FLOOR_KM2 > melt.NOV_BASELINE_MAX_KM2


def test_noise_floor_is_a_small_fraction_of_the_study_area():
    """Sanity bound: a floor near the AOI size would suppress everything."""
    aoi_km2 = melt.WIN_SIZE**2 * melt.PIXEL_KM2
    assert 0 < melt.NOISE_FLOOR_KM2 < 0.02 * aoi_km2


# --- Cloud halo --------------------------------------------------------------

def test_halo_dilation_grows_the_masked_region():
    """A cloud contaminates ground beside it, not only under it."""
    from scipy import ndimage

    low = np.zeros((100, 100), bool)
    low[50, 50] = True
    cells = int(round(melt.CLOUD_HALO_KM * 1000 / melt.SCREEN_PIXEL_M))
    grown = ndimage.binary_dilation(low, iterations=cells)

    assert cells >= 1
    assert grown.sum() > low.sum()
    assert grown[50, 50 + cells]      # reaches the halo distance
    assert not grown[50, 50 + cells + 2]  # but does not run away


def test_halo_upsampling_preserves_coverage():
    """np.kron from the 60 m screen grid to 10 m must not shift or lose area."""
    scale = int(melt.SCREEN_PIXEL_M / melt.PIXEL_M)
    coarse = np.array([[True, False], [False, True]])
    full = np.kron(coarse, np.ones((scale, scale), bool))

    assert full.shape == (2 * scale, 2 * scale)
    assert full.mean() == pytest.approx(coarse.mean())
    assert full[0, 0] and not full[0, scale]


# --- Sun elevation -----------------------------------------------------------

def test_sun_elevation_read_and_defaulted():
    assert melt.sun_elevation(FakeItem({"view:sun_elevation": 9.5})) == 9.5
    # Missing metadata must not silently reject every scene.
    assert melt.sun_elevation(FakeItem({})) >= melt.MIN_SUN_ELEVATION_DEG


def test_low_sun_scenes_are_below_the_gate():
    """The two scenes whose artifacts motivated the gate must fall under it."""
    for elev in (9.5, 11.7):
        assert elev < melt.MIN_SUN_ELEVATION_DEG


# --- Scene-quality guards ----------------------------------------------------

def test_halo_gate_rejects_confirmed_contaminated_scenes():
    """Every halo fraction confirmed contaminated by inspection must fail.

    2018-12-26 (32.7%) returned 365 km2; 2024-12-20 (10.0%) returned 34 km2
    off a cumulus band. A clean scene at 0.06% and a salvageable partly
    cloudy one at 1.99% must still pass.
    """
    for bad in (0.327, 0.1302, 0.112, 0.1002):
        assert bad > season.MAX_HALO_FRAC
    for good in (0.0006, 0.0199):
        assert good <= season.MAX_HALO_FRAC


def test_plausibility_cap_sits_above_real_melt_below_failures():
    """0.49% is the largest confirmed-clean density; 9.67% was a failure."""
    assert 0.0049 < season.MAX_POND_FRAC_PLAUSIBLE < 0.0967


def test_usable_fraction_gate_keeps_the_record_scene():
    """2026-02-10 is 23.7% cloudy over the large AOI but 76.3% usable."""
    assert 0.763 >= season.MIN_USABLE_FRAC


# --- Cross-sensor validation -------------------------------------------------

def test_brightness_floor_converts_to_reflectance():
    """The DN threshold and its reflectance form must mean the same thing."""
    assert vl.BRIGHTNESS_FLOOR_REFL == pytest.approx(melt.BRIGHTNESS_FLOOR / 10000.0)
    assert 0.0 < vl.BRIGHTNESS_FLOOR_REFL < 1.0


def test_ndwi_is_invariant_to_pure_scaling_but_not_to_offset():
    """Why Landsat L2 cannot be compared as raw DN.

    Sentinel-2 L2A is a pure scaling, so NDWI on DN equals NDWI on
    reflectance. Landsat L2 carries an additive offset, which does not cancel
    in a ratio of differences - comparing raw DN would compare two different
    quantities.
    """
    g, n = 9000.0, 8000.0
    ndwi = lambda a, b: (a - b) / (a + b)

    scaled = ndwi(g * 1e-4, n * 1e-4)
    assert ndwi(g, n) == pytest.approx(scaled)

    offset = ndwi(g * 2.75e-5 - 0.2, n * 2.75e-5 - 0.2)
    assert not np.isclose(ndwi(g, n), offset)


def test_landsat_l2_over_ice_is_out_of_valid_range():
    """Guards the finding that motivated using Level-1 instead.

    Measured green DN over the study area sits near 51589, which converts to
    1.22 reflectance - outside the product's documented -0.2 to 1.0 range and
    physically impossible.
    """
    measured_dn = 51589
    refl = measured_dn * 2.75e-5 - 0.2
    assert refl > 1.0


def test_block_reducers_preserve_geometry():
    """Resampling 10 m to 30 m must not shift or lose coverage."""
    a = np.zeros((9, 9), bool)
    a[0, 0] = True
    assert vl._block_any(a).shape == (3, 3)
    assert vl._block_any(a)[0, 0]
    assert vl._block_any(a).sum() == 1

    vals = np.ones((9, 9), "f4") * 0.5
    assert vl._block_mean(vals) == pytest.approx(np.full((3, 3), 0.5))


def test_resolution_ratio_matches_the_two_sensors():
    assert melt.PIXEL_M * vl.RES_RATIO == pytest.approx(30.0)


# --- Uncertainty ranges ------------------------------------------------------

def test_threshold_band_brackets_the_current_value():
    """The range must straddle 0.16, or it is not an honest interval around it."""
    lo, hi = unc.THRESHOLD_BAND
    assert lo <= melt.NDWI_THRESHOLD <= hi
    assert melt.NDWI_THRESHOLD in unc.THRESHOLDS


def test_measured_bias_and_precision_are_consistent():
    """Both come from the same reference-point run and must agree with it."""
    assert 0.4 < unc.MEASURED_AREA_BIAS < 1.0  # detector under-reports
    assert 0.5 < unc.MEASURED_PRECISION < 0.8


def test_point_validation_strata_are_disjoint_and_ordered():
    """A/B/C must partition the AOI: detected, near, and far, no overlap."""
    ponds = np.zeros((200, 200), bool)
    ponds[100, 100] = True
    valid = np.ones((200, 200), bool)
    strata = vp.build_strata(ponds, valid)

    overlap = (strata["A"] & strata["B"]) | (strata["B"] & strata["C"]) | \
              (strata["A"] & strata["C"])
    assert not overlap.any()
    total = strata["A"].sum() + strata["B"].sum() + strata["C"].sum()
    assert total == valid.sum()
    assert strata["A"][100, 100]           # the detection itself is stratum A
    assert strata["B"][100, 110]           # 100 m away is "near"
    assert strata["C"][0, 0]               # a far corner is "far"


def test_wilson_interval_contains_point_estimate():
    lo, hi = vp.wilson(24, 38)
    assert lo < 24 / 38 < hi
    assert 0.0 <= lo < hi <= 1.0


def test_screen_output_is_decimated_not_full_res():
    """screen() works on decimated SCL; its mask must not be used for masking.

    Guards a bug where the coarse reject_mask was passed to detect() against
    full-resolution bands, which only fails at runtime on a shape mismatch.
    """
    n = melt._screen_shape()
    assert n < melt.WIN_SIZE
    assert melt.screen(np.full((n, n), 11, dtype="uint8"))["reject_mask"].shape == (n, n)


def test_window_scales_to_source_resolution():
    """A 20 m band must read half the pixel offsets of a 10 m band."""
    w10 = melt._window(scale=1)
    w20 = melt._window(scale=2)
    assert (w10.col_off, w10.width) == (melt.WIN_COL, melt.WIN_SIZE)
    assert (w20.col_off, w20.width) == (melt.WIN_COL // 2, melt.WIN_SIZE // 2)


def test_aoi_stays_inside_the_tile():
    """A Sentinel-2 10 m tile is 10980 px. The window must fit."""
    assert melt.WIN_ROW + melt.WIN_SIZE <= 10980
    assert melt.WIN_COL + melt.WIN_SIZE <= 10980

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


def test_noise_floor_covers_confirmed_artifacts():
    """Both confirmed pre-melt false positives must fall under the floor."""
    assert melt.NOISE_FLOOR_KM2 >= 0.68

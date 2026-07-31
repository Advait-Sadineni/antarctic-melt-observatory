"""Offline tests for the M3 SAR melt-state core - synthetic arrays only."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

import numpy as np

from core.pc import Sentinel1Source


class FakeS1Item:
    def __init__(self, orbit=65):
        self.id = "S1B_IW_GRDH_1SSH_20210131T073912_rtc"
        self.properties = {"sat:relative_orbit": orbit}

        class A:
            href = "https://pc/hh.tif"
        self.assets = {"hh": A()}


def test_s1_source_band_and_orbit():
    src = Sentinel1Source.__new__(Sentinel1Source)
    it = FakeS1Item(orbit=65)
    assert src.band_href(it, "hh") == "https://pc/hh.tif"
    assert src.relative_orbit(it) == 65
    assert Sentinel1Source.COLLECTION == "sentinel-1-rtc"


def test_to_db_handles_zeros_and_scales():
    from sar import to_db
    g = np.array([[1.0, 0.1], [0.0, -1.0]], "f4")
    db = to_db(g)
    assert db[0, 0] == 0.0 and abs(db[0, 1] + 10.0) < 1e-4
    assert np.isnan(db[1, 0]) and np.isnan(db[1, 1])


def test_baseline_median_ignores_nan():
    import sar
    stack = [np.array([[1.0, np.nan]], "f4"),
             np.array([[3.0, 5.0]], "f4"),
             np.array([[2.0, np.nan]], "f4")]
    out = sar._median_stack(stack)
    assert out[0, 0] == 2.0 and out[0, 1] == 5.0


def test_wet_mask_nan_safe():
    import sar
    db = np.array([[-20.0, -10.0], [np.nan, -14.0]], "f4")
    base = np.array([[-15.0, -12.0], [-13.0, np.nan]], "f4")
    wet, obs = sar.wet_mask(db, base)
    assert wet[0, 0] and not wet[0, 1]
    assert not wet[1, 0] and not wet[1, 1]
    assert obs[0, 0] and obs[0, 1] and not obs[1, 0] and not obs[1, 1]


def test_composite_accumulator_first_last_wet():
    import sar
    acc = sar.CompositeAccumulator(shape=(1, 2))
    w1 = np.array([[True, False]]); o1 = np.array([[True, True]])
    w2 = np.array([[True, True]]); o2 = np.array([[True, True]])
    w3 = np.array([[False, True]]); o3 = np.array([[True, True]])
    acc.add(w1, o1, day=10); acc.add(w2, o2, day=22); acc.add(w3, o3, day=40)
    assert acc.wet_days[0, 0] == 2 and acc.wet_days[0, 1] == 2
    assert acc.n_obs[0, 0] == 3
    assert acc.first_wet[0, 0] == 10 and acc.first_wet[0, 1] == 22
    assert acc.last_wet[0, 0] == 22 and acc.last_wet[0, 1] == 40


def test_checkpoint_roundtrip(tmp_path):
    import sar
    acc = sar.CompositeAccumulator((4, 5))
    wet = np.zeros((4, 5), bool); wet[1, 2] = True
    obs = np.ones((4, 5), bool)
    acc.add(wet, obs, day=10)
    p = tmp_path / "ckpt.npz"
    sar.save_checkpoint(p, acc, {"a", "b"}, used=2, skipped=1, wet_max=3.5)
    acc2, done, used, skipped, wet_max = sar.load_checkpoint(p, (4, 5))
    assert done == {"a", "b"} and used == 2 and skipped == 1
    assert abs(wet_max - 3.5) < 1e-9
    assert np.array_equal(acc2.wet_days, acc.wet_days)
    assert np.array_equal(acc2.first_wet, acc.first_wet)
    assert np.array_equal(acc2.n_obs, acc.n_obs)


def test_checkpoint_missing_returns_fresh(tmp_path):
    import sar
    acc, done, used, skipped, wet_max = sar.load_checkpoint(
        tmp_path / "nope.npz", (3, 3))
    assert done == set() and used == 0 and skipped == 0 and wet_max == 0.0
    assert acc.wet_days.shape == (3, 3) and not acc.wet_days.any()

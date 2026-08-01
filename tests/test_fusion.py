"""Offline tests for the M3 Phase B melt-state fusion - synthetic arrays only."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

import numpy as np

import fusion


def test_state_table():
    wet = np.array([[1, 1, 0, 0], [1, 0, 1, 0]], bool)
    obs = np.array([[1, 1, 1, 1], [1, 1, 0, 0]], bool)
    pond = np.array([[1, 0, 1, 0], [0, 1, 1, 1]], bool)
    st = fusion.classify_state(wet, obs, pond)
    # PONDED(3) needs wet+obs+pond; WET(2) wet+obs; DRY(1) obs only;
    # UNOBSERVED(0) no obs. Radar-dry + pond evidence is a CONFLICT, not a state.
    assert st.tolist() == [[3, 2, 1, 1], [2, 1, 0, 0]]
    cf = fusion.conflict_mask(wet, obs, pond)
    assert cf.tolist() == [[False, False, True, False], [False, True, False, False]]


def test_confidence_bounds_and_decay():
    margin = np.array([[0.0, 3.0, 12.0]], "f4")
    gap = np.array([[0.0, 6.0, 0.0]], "f4")
    c = fusion.confidence(margin, gap)
    assert c[0, 0] == 0.0
    assert 0.17 < c[0, 1] < 0.19          # (3/6) * exp(-1) = 0.5*0.3679
    assert c[0, 2] == 1.0                  # capped at 1, gap 0
    assert c.min() >= 0.0 and c.max() <= 1.0


def test_nearest_pond_evidence_window():
    from datetime import date
    dates = [date(2021, 1, 10), date(2021, 1, 24)]
    m1 = np.array([[True, False]])
    m2 = np.array([[False, True]])
    ev = fusion.nearest_pond_evidence(dates, [m1, m2], date(2021, 1, 22), window_days=6)
    assert ev.tolist() == [[False, True]]          # nearest = Jan 24 (2 days)
    ev2 = fusion.nearest_pond_evidence(dates, [m1, m2], date(2021, 2, 15), window_days=6)
    assert not ev2.any()                            # nothing within window


def test_nearest_pond_with_gap():
    from datetime import date
    import fusion
    m1 = np.zeros((2, 2), bool); m1[0, 0] = True
    m2 = np.zeros((2, 2), bool); m2[1, 1] = True
    dates = [date(2021, 1, 10), date(2021, 1, 20)]
    m, gap = fusion.nearest_pond_with_gap(dates, [m1, m2], date(2021, 1, 12))
    assert m[0, 0] and not m[1, 1] and gap == 2
    m, gap = fusion.nearest_pond_with_gap(dates, [m1, m2], date(2021, 3, 1))
    assert not m.any() and gap is None


def test_state_accumulator_two_scene_walk():
    import fusion
    acc = fusion.StateAccumulator((2, 2))
    wet = np.array([[1, 0], [1, 0]], bool)
    obs = np.array([[1, 1], [1, 0]], bool)
    pond = np.array([[1, 0], [0, 0]], bool)
    st = fusion.classify_state(wet, obs, pond)
    conf = np.full((2, 2), 0.5, "f4")
    acc.add(st, fusion.conflict_mask(wet, obs, pond), conf, day=10)
    # scene 2: pond cell goes radar-dry -> conflict, not PONDED
    wet2 = np.array([[0, 0], [1, 0]], bool)
    st2 = fusion.classify_state(wet2, obs, pond)
    acc.add(st2, fusion.conflict_mask(wet2, obs, pond), conf, day=20)
    assert acc.days_ponded[0, 0] == 1 and acc.days_wet[0, 0] == 1
    assert acc.days_wet[1, 0] == 2 and acc.days_ponded[1, 0] == 0
    assert acc.conflict_days[0, 0] == 1          # scene-2 pond-but-dry
    assert acc.n_obs[0, 1] == 2 and acc.n_obs[1, 1] == 0
    assert acc.first_ponded[0, 0] == 10 and acc.last_ponded[0, 0] == 10


def test_state_checkpoint_roundtrip(tmp_path):
    import fusion
    acc = fusion.StateAccumulator((2, 3))
    st = np.array([[3, 2, 1], [0, 1, 2]], "u1")
    acc.add(st, np.zeros((2, 3), bool), np.full((2, 3), 0.25, "f4"), day=5)
    p = tmp_path / "st.npz"
    fusion.save_state_checkpoint(p, acc, {"x"}, 1, 0, wet_max=9.5, pond_max=2.5)
    acc2, done, used, skipped, wm, pm = fusion.load_state_checkpoint(p, (2, 3))
    assert done == {"x"} and used == 1
    assert wm == 9.5 and pm == 2.5
    assert np.array_equal(acc2.days_ponded, acc.days_ponded)
    assert np.array_equal(acc2.conf_sum, acc.conf_sum)

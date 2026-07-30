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

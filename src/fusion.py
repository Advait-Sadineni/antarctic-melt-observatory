"""M3 Phase B: fuse radar wetness with optical ponds into melt states.

State encoding (uint8): 0=UNOBSERVED, 1=DRY, 2=WET, 3=PONDED.
PONDED requires radar-wet AND optical pond evidence - by construction the
classifier cannot emit PONDED where radar says dry; optical-pond-but-radar-dry
is a counted CONFLICT (the diagnostic of both systems), never a state.
"""
import numpy as np

UNOBSERVED, DRY, WET, PONDED = 0, 1, 2, 3


def classify_state(wet, obs, pond_evidence):
    """Vectorized state table; see module docstring."""
    st = np.zeros(wet.shape, "u1")
    st[obs] = DRY
    st[obs & wet] = WET
    st[obs & wet & pond_evidence] = PONDED
    return st


def conflict_mask(wet, obs, pond_evidence):
    """Optical says pond, radar (observed) says dry - count, don't classify."""
    return pond_evidence & obs & ~wet


def confidence(margin_db, gap_days):
    """Per-cell confidence in [0,1]: radar margin below baseline (saturating at
    6 dB) times optical-evidence freshness (e-folding 6 days)."""
    m = np.clip(np.asarray(margin_db, "f4") / 6.0, 0.0, 1.0)
    g = np.exp(-np.asarray(gap_days, "f4") / 6.0)
    return (m * g).astype("f4")


def nearest_pond_evidence(dates, masks, scene_date, window_days=6):
    """The pond mask nearest in time to scene_date, all-False if none within
    the window. dates: list[date]; masks: list[bool ndarray] (same shapes)."""
    best, best_gap = None, window_days + 1
    for d, m in zip(dates, masks):
        gap = abs((d - scene_date).days)
        if gap <= window_days and gap < best_gap:
            best, best_gap = m, gap
    if best is None:
        return np.zeros(masks[0].shape, bool)
    return best

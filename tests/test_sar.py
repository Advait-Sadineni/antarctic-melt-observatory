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

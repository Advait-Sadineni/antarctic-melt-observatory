"""Adapter tests - offline, FakeItem pattern from test_melt.py.

The SceneSource adapters are the ONLY code that knows a provider's asset keys
or item-id conventions; these tests pin that mapping per provider.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from core.sources import EarthSearchSource  # noqa: E402


class FakeAsset:
    def __init__(self, href):
        self.href = href


class FakeESItem:
    """Earth Search style: id S2B_19DEA_20210124_0_L2A, friendly asset keys."""
    def __init__(self):
        self.id = "S2B_19DEA_20210124_0_L2A"
        self.properties = {}
        self.assets = {k: FakeAsset(f"https://es/{k}.tif")
                       for k in ("red", "green", "blue", "nir", "swir16", "scl")}


def test_earthsearch_band_href_uses_friendly_names():
    src = EarthSearchSource.__new__(EarthSearchSource)  # no network in __init__
    it = FakeESItem()
    assert src.band_href(it, "green") == "https://es/green.tif"
    assert src.band_href(it, "swir16") == "https://es/swir16.tif"


def test_earthsearch_tile_of_parses_id():
    src = EarthSearchSource.__new__(EarthSearchSource)
    assert src.tile_of(FakeESItem()) == "19DEA"


def test_earthsearch_tile_query_shape():
    src = EarthSearchSource.__new__(EarthSearchSource)
    assert src.tile_query("19DEA") == {"grid:code": {"eq": "MGRS-19DEA"}}

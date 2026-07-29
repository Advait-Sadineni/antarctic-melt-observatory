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


from core.pc import PlanetaryComputerSource  # noqa: E402


class FakePCItem:
    """PC style: id S2B_MSIL2A_..._T19DEA_..., B-number asset keys."""
    def __init__(self):
        self.id = "S2B_MSIL2A_20210124T120000_R000_T19DEA_20210124T150000"
        self.properties = {"s2:mgrs_tile": "19DEA"}
        self.assets = {k: FakeAsset(f"https://pc/{k}.tif")
                       for k in ("B02", "B03", "B04", "B08", "B11", "SCL")}


def test_pc_band_href_maps_logical_to_b_numbers():
    src = PlanetaryComputerSource.__new__(PlanetaryComputerSource)
    it = FakePCItem()
    assert src.band_href(it, "green") == "https://pc/B03.tif"
    assert src.band_href(it, "scl") == "https://pc/SCL.tif"


def test_pc_tile_of_uses_property_not_id():
    src = PlanetaryComputerSource.__new__(PlanetaryComputerSource)
    assert src.tile_of(FakePCItem()) == "19DEA"


def test_pc_tile_query_shape():
    src = PlanetaryComputerSource.__new__(PlanetaryComputerSource)
    assert src.tile_query("19DEA") == {"s2:mgrs_tile": {"eq": "19DEA"}}


import melt  # noqa: E402


class RecordingSource:
    BAND_MAP = {b: b for b in ("red", "green", "blue", "nir", "swir16", "scl")}
    def __init__(self):
        self.calls = []
    def band_href(self, item, band):
        self.calls.append(("band_href", band))
        return item.assets[band].href
    def tile_of(self, item):
        return "19DEA"
    def search(self, **kw):
        self.calls.append(("search", kw)); return []
    def tile_query(self, tile):
        return {"grid:code": {"eq": f"MGRS-{tile}"}}


def test_set_source_swaps_and_tile_of_delegates():
    rec = RecordingSource()
    old = melt.get_source()
    try:
        melt.set_source(rec)
        assert melt.get_source() is rec
        assert melt._tile_of(FakeESItem()) == "19DEA"   # delegated, not id-parsed
    finally:
        melt.set_source(old)


import json as _json
import numpy as np
import rasterio
from rasterio.transform import from_origin
from core.products import write_product  # noqa: E402


def test_write_product_cog_and_item(tmp_path):
    water = np.zeros((40, 50), "f4"); water[10:14, 20:26] = 0.5
    tr = from_origin(-2_000_000, 1_000_000, 30, 30)
    d = write_product("george_vi", "2020-21", water, tr, "EPSG:3031",
                      {"area_km2": 142.48, "shelf_area_km2": 14390.6,
                       "obs_cloud": 1.4, "poorly_observed": False,
                       "sensor": "sentinel-2"}, root=tmp_path)
    with rasterio.open(d / "water_fraction.tif") as src:
        assert src.crs.to_string() == "EPSG:3031"
        assert src.transform == tr
        assert float(src.read(1)[12, 22]) == 0.5
    item = _json.loads((d / "item.json").read_text())
    p = item["properties"]
    assert p["area_km2"] == 142.48 and p["sensor"] == "sentinel-2"
    # reserved fields exist and are null until M3/M4 fill them
    for k in ("volume_km3", "depth_mean_m", "uncertainty_km2", "sensor_confidence"):
        assert k in p and p[k] is None

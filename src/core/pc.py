"""Microsoft Planetary Computer adapter. Same Sentinel-2 archive, different
asset keys and item ids, and hrefs must be SAS-signed before rasterio reads
them - sign_inplace on the client handles that transparently."""
from .sources import SceneSource


class PlanetaryComputerSource(SceneSource):
    STAC_API = "https://planetarycomputer.microsoft.com/api/stac/v1"
    BAND_MAP = {"red": "B04", "green": "B03", "blue": "B02",
                "nir": "B08", "swir16": "B11", "scl": "SCL"}

    def band_href(self, item, band):
        """Re-sign at ACCESS time, not search time. PC SAS tokens live ~1 h;
        long runs (winter baselines, season composites) outlive them and every
        later read 403s. planetary_computer.sign() keeps a token cache and
        refreshes expired ones, so signing the bare URL here makes reads
        expiry-proof for any run length."""
        import planetary_computer
        href = item.assets[self.BAND_MAP[band]].href.split("?")[0]
        return planetary_computer.sign(href)

    def _make_client(self):
        import planetary_computer
        from pystac_client import Client
        return Client.open(self.STAC_API,
                           modifier=planetary_computer.sign_inplace)

    def tile_of(self, item):
        return item.properties["s2:mgrs_tile"]

    def tile_query(self, tile):
        return {"s2:mgrs_tile": {"eq": tile}}


class Sentinel1Source(PlanetaryComputerSource):
    """Sentinel-1 RTC (terrain-corrected gamma0). Not MGRS-tiled: search by
    bbox and group by relative orbit - radar geometry differs per track, so
    baselines and detections must never mix orbits."""
    COLLECTION = "sentinel-1-rtc"
    BAND_MAP = {"hh": "hh"}

    def relative_orbit(self, item):
        return int(item.properties["sat:relative_orbit"])

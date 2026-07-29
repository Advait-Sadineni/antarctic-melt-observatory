"""Microsoft Planetary Computer adapter. Same Sentinel-2 archive, different
asset keys and item ids, and hrefs must be SAS-signed before rasterio reads
them - sign_inplace on the client handles that transparently."""
from .sources import SceneSource


class PlanetaryComputerSource(SceneSource):
    STAC_API = "https://planetarycomputer.microsoft.com/api/stac/v1"
    BAND_MAP = {"red": "B04", "green": "B03", "blue": "B02",
                "nir": "B08", "swir16": "B11", "scl": "SCL"}

    def _make_client(self):
        import planetary_computer
        from pystac_client import Client
        return Client.open(self.STAC_API,
                           modifier=planetary_computer.sign_inplace)

    def tile_of(self, item):
        return item.properties["s2:mgrs_tile"]

    def tile_query(self, tile):
        return {"s2:mgrs_tile": {"eq": tile}}

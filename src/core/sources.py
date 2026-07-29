"""Scene-source adapters: the ONLY place that knows a provider's STAC API,
asset keys, or item-id conventions. Core code speaks logical band names
(red/green/blue/nir/swir16/scl) and calls search/band_href/tile_of/tile_query.

This seam is the regret-proofing keystone from the observatory design: the
science never imports a provider, so the compute backend stays swappable.
"""


class SceneSource:
    STAC_API = None          # subclass sets
    COLLECTION = "sentinel-2-l2a"
    BAND_MAP = {}            # logical band -> provider asset key

    def __init__(self):
        from pystac_client import Client
        self.client = Client.open(self.STAC_API)

    def search(self, collections=None, query=None, datetime=None,
               bbox=None, ids=None, limit=100):
        """Thin pystac-client wrapper; returns a LIST of items (limit is the
        page size - pystac paginates through everything)."""
        kw = {"collections": collections or [self.COLLECTION], "limit": limit}
        if query is not None:
            kw["query"] = query
        if datetime is not None:
            kw["datetime"] = datetime
        if bbox is not None:
            kw["bbox"] = bbox
        if ids is not None:
            kw["ids"] = ids
        return list(self.client.search(**kw).items())

    def band_href(self, item, band):
        return item.assets[self.BAND_MAP[band]].href

    def tile_of(self, item):
        """MGRS tile of an item, e.g. '19DEA'."""
        raise NotImplementedError

    def tile_query(self, tile):
        """Provider-specific STAC query fragment selecting one MGRS tile."""
        raise NotImplementedError


class EarthSearchSource(SceneSource):
    STAC_API = "https://earth-search.aws.element84.com/v1"
    BAND_MAP = {b: b for b in ("red", "green", "blue", "nir", "swir16", "scl")}

    def tile_of(self, item):
        # id like S2B_19DEA_20210124_0_L2A
        return item.id.split("_")[1]

    def tile_query(self, tile):
        return {"grid:code": {"eq": f"MGRS-{tile}"}}

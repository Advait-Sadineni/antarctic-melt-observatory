# Validation plan

Written 2026-07-21. **Done so far: Tier 1 (both items), Tier 4.** Outstanding:
Tier 2 (literature + published methods), Tier 3 (ERA5). See the README's
Validation section for full results.

Status by item:

- **[DONE] Landsat cross-sensor** — area agreement 0.90 on clean same-day
  pairs, *r* = 0.83, Spearman = 0.86; per-pixel IoU 0.30–0.35, dominated by
  the 30 m vs 10 m resolution gap. 2019-02-17 failed outright.
- **[DONE] Reference-point accuracy** — precision 0.63, recall 0.36, F1 0.46
  on 96 blind-labelled points. Every false positive is crevasse shadow at
  NDWI 0.17–0.25; classes separate cleanly and the threshold sits in the
  overlap. Precision reaches 1.00 by NDWI 0.25.
- **[DONE] Uncertainty as ranges** (`src/uncertainty.py`) — median range is
  98% of the point estimate. 2018-19's apparent peak collapses; 2025-26
  survives as the robust maximum.
- **[TODO] Literature comparison** — the one external anchor still missing.
- **[TODO] Published detection methods** (Moussavi et al. thresholds).
- **[TODO] ERA5 physical consistency** — needs a Copernicus CDS account.

The single most actionable finding: the NDWI threshold is set too low.
Reference points put every crevasse-shadow false positive at 0.17–0.25 and
real ponds at 0.38+. Retuning is roadmap item 0, to be done on labels from
*several* scenes, not the one done so far.

## Why this comes before new features

The pipeline is internally consistent and externally unanchored. Its headline
number moved four times during development — 8.80 → 11.52 → 34.23 → 9.60 km²
for various seasons — every time because a scene that looked fine turned out
to be cloud, shadow, or a cherry-picked window. Each of those was caught by
inspecting imagery, which does not scale and is not evidence anyone else can
check.

Until the detector is measured against something outside itself, no number it
produces should be quoted. That makes validation a prerequisite, not a
polish step.

Success is **not** matching published absolute areas — different study areas
and thresholds guarantee divergence. Success is agreement on **timing,
ranking, and order of magnitude**, plus a defensible IoU against hand labels.

---

## Tier 1 — Does the detector find the right pixels?

### 1. Hand-labelled reference scenes → precision / recall / IoU

Every number is downstream of "does the mask match reality", and that has
only ever been eyeballed.

- Digitize ponds manually at full zoom on 4 scenes:
  - peak melt (`S2B_19CEV_20210124_0_L2A`, visually confirmed clean)
  - shoulder season
  - pre-melt November (expect near-zero — tests the false-positive rate)
  - a known-contaminated scene (`S2B_19CEV_20210123_0_L2A`) as a negative control
- Score precision, recall, F1, IoU per scene
- Report the false-positive rate on the pre-melt scene separately; that is
  the number the 2.2 km² noise floor is standing in for

No accounts, no new data. Effort: about half a day, mostly careful clicking.

### 2. Cross-sensor agreement with Landsat 8 — **DONE**

Implemented in `src/validate_landsat.py`. Results and caveats in the README.

Notes for anyone repeating this, since the data routing was the hard part:

- `landsat-c2-l2` in the Earth Search STAC is **unusable over ice** — 96.9% of
  its pixels here fall outside the product's own valid range. Level-2 surface
  reflectance is not designed for bright cryospheric surfaces. Use Level-1 and
  compute TOA reflectance.
- Level-1 pixel access took four attempts: AWS `usgs-landsat` is Requester
  Pays, USGS LandsatLook returns an auth page, Planetary Computer's
  `landsat-c2-l1` holds only Landsat 1–5 MSS. Google's public
  `gcp-public-data-landsat` bucket works anonymously and covers 2013–2021.
- Compare at **matched resolution**. Sentinel-2 must be averaged to 30 m
  *before* detection, or the comparison measures the resolution gap rather
  than the detector.
- Do not hardcode the tile origin. Read georeferencing from the imagery; a
  two-pixel error looks exactly like sensor disagreement.

Still open from this item: Landsat Collection-1 ends in 2021, so seasons
2022-26 — including the record 2025-26 — remain unvalidated. Collection-2
Level-1 needs a USGS ERS account, which would close that gap.

---

## Tier 2 — Do the numbers match published science?

### 3. Literature comparison, George VI specifically

Candidate sources (**verify these citations before relying on them** — they
are recalled, not checked):

- Arthur et al. — supraglacial lakes on George VI across multiple seasons;
  the most directly comparable work
- Stokes et al. 2019 — continent-wide Antarctic supraglacial lake inventory
- Banwell et al. — George VI surface meltwater and ice-shelf flexure
- Dell et al. — slush vs ponded water on ice shelves, relevant to the
  unresolved 2019-02-17 scene

Compare pond-area *fractions* and season *rankings*, not raw km². This is
what settles whether 2025-26 being the record season is real.

### 4. Compare against a published detection method

- Moussavi et al. 2020 published NDWI thresholds plus supporting band tests
  for Antarctic supraglacial lakes
- Run their method beside ours on the same scenes and compare masks
- Matters because local threshold sensitivity is **11% of area per 0.01**, so
  0.16 resting only on our own sweep is a real weakness

---

## Tier 3 — Is it physically consistent?

### 5. ERA5 climate cross-check

- Needs a [Copernicus CDS](https://cds.climate.copernicus.eu/) account →
  credentials in `~/.cdsapirc`
- Pull summer air temperature and positive-degree-days over George VI per season
- Correlate against our peak areas
- 2019-20 was a record-warm Antarctic Peninsula summer (18.3 °C at Esperanza,
  Feb 2020). A correct detector should reflect that
- Doubles as groundwork for the ERA5 + ML roadmap item

### 6. Resolve 2019-02-17

The 11.6 km² scene that makes 2018-19 the second-highest season and that no
screen explains — not cloud (halo catches 1%), not low sun (21.4°), not
nodata. Probably slush, which NDWI cannot separate from ponded water.

Tier 1 labelling plus a Landsat pair on that date should settle it.

---

## Tier 4 — Report uncertainty honestly

### 7. Replace point estimates with ranges

- Quote every area as a band across the plausible threshold range rather than
  a single figure
- Propagate partial-coverage sampling bias into the interval — a 56%-coverage
  scene on 2020-01-19 extrapolates to 15.8 km² where full coverage gives 10.0
- Turns "18.66 km²" into something defensible

---

## Order of work

1. **Landsat cross-check** — cheapest, no accounts, validates or breaks the
   pipeline in an afternoon
2. **Hand labels** — hard metrics
3. **Literature** — external anchor
4. Then Tier 3 and 4

Register the Copernicus CDS account whenever convenient; Tiers 1–2 need nothing.

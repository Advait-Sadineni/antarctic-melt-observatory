# Validation plan

Written 2026-07-21. Nothing here is done yet.

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

### 2. Cross-sensor agreement with Landsat 8/9

Independent satellite, independent band passes, same ground. Agreement is
strong evidence; disagreement localizes the fault.

- `landsat-c2-l2` is already in the Earth Search STAC this project uses, so
  **no new credentials**
- Find Landsat/Sentinel-2 pairs ≤1 day apart over the AOI
- Run NDWI on both, compare total area and per-pixel overlap
- Landsat is 30 m vs our 10 m, so resample and expect Landsat to under-detect
  narrow channels — quantify that bias rather than treating it as error

Highest evidence per unit effort. Do this first.

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

# UNDERSTANDING.md — the whole observatory, at three depths

*This document exists so that Advait can defend every layer of this system
without notes — to a professor, a poster visitor, or a hostile reviewer.
Read it top to bottom once, then use it as the syllabus for the teaching
sessions. Every claim here traces to a file in this repo.*

---

# DEPTH 1 — The elevator (30 seconds)

> "Meltwater ponding is what kills Antarctic ice shelves — it wedges cracks
> open until the shelf shatters, like Larsen B in 2002. I built an open
> system that watches 12 ice shelves with three kinds of satellite — camera,
> radar, and laser — and measures not just where meltwater is, but whether
> it's *ponded* (the dangerous kind) and *how much water* it holds. It's
> blind-validated, it updates every season, it runs on free data, and it
> found something new: the dangerous ponded fraction collapses ~1,000×
> between warm and cool years even though the total wetted area barely
> changes. Ponding is nonlinear in melt — and that matters for sea-level
> projections."

If you memorize one paragraph, it's that one.

---

# DEPTH 2 — The whiteboard tour (eight layers)

Each layer: WHAT it does → WHY that way → the NUMBER it produces → its
HONEST WEAKNESS. You should be able to draw this as eight boxes with arrows.

## Layer 1 — The data (all free)

- **Sentinel-2** (ESA): optical camera, 10–30 m pixels, revisits every few
  days. Sees color — meltwater is dark blue, ice is bright. Blind to
  anything under cloud. Source: AWS Earth Search + Microsoft Planetary
  Computer (two independent backends; we proved they agree within 2%).
- **Sentinel-1** (ESA): radar. Sends microwave pulses, measures the echo
  ("backscatter"). Sees through clouds and darkness. Wet snow absorbs
  microwaves → echo drops. Cannot tell a puddle from damp slush.
- **ICESat-2** (NASA): laser altimeter. Fires green laser pulses; photons
  bounce back from surfaces. Over a pond, some photons reflect off the water
  surface and some off the pond BOTTOM → two echo populations → the gap
  between them is a direct depth measurement.
- **MEaSUReS**: the official shelf outline polygons, so "on-shelf" is
  defined by the community standard, not by us.

*Why three sensors?* Each one's weakness is another's strength: camera has
color but fears cloud; radar ignores cloud but lacks color; laser is
razor-accurate but only along thin lines. The system is designed so no
number depends on a single sensor's blind spot.

## Layer 2 — Optical detection (finding water in a photo)

**NDWI** (Normalized Difference Water Index) = (blue − red)/(blue + red) on
reflectance. Water absorbs red light and reflects blue → high NDWI. We use
the published Moussavi et al. (2020) recipe: NDWI > 0.19 = water seed, plus
a shadow test (green − red > 0.09) because cliff shadows also look "blue-ish"
but fail it. Then **hysteresis**: any pixel > 0.14 that TOUCHES a seed is
also water. Why: pond edges are shallow and less blue than centers; a single
hard threshold clips every pond's margins. Seeds must clear the strict bar;
margins may join only if connected to a real seed. (Same trick as Canny edge
detection in computer vision.)

- Produces: per-scene water masks; the per-shelf area record.
- Weakness: thin cloud/haze also reads water-blue → Layer 3 exists.

## Layer 3 — Clear-scene selection (the cloud war)

The single biggest error source in Antarctic optical melt mapping is not the
threshold — it's counting cloud as water. Metadata cloud % is tile-wide and
misses local haze sitting exactly over the shelf. Our test: the **halo
fraction** — how much of the shelf's surroundings our own cloud mask flags.
A season's area = the *peak-melt scene among genuinely clear scenes*
(halo < 0.08). Twist: in record melt years, melt itself darkens snow (lowers
NDSI) so no scene ever looks "clean" — the fallback branch accepts scenes
with metadata < 15% and halo < 0.30. Real haze fails one gate or the other.

We TRIED and REJECTED two popular alternatives, with receipts: persistence
compositing (haze recurs over the same spot → phantom 781 km²) and
clean-scene union (inflated a validated year 98 → 174). Selection windows
are evaluated PER MONTH for evidence series (decision 0005) because a clean
December was silently vetoing hazy record-January scenes.

- Produces: one defensible scene choice per tile per season.
- Weakness: one scene per season = a lower bound on the season's true peak.

## Layer 4 — Blind validation (how we know it's right)

For a validated tile: 80 points, stratified (40 in detected water, 26 near
the boundary, 14 far away), labeled against imagery WITHOUT seeing the
pipeline's answer — the answer key stays hidden until scoring, scored once.
This is design-based inference (the same logic as polling): precision = "of
what we call water, how much is water," recall = "of true water, how much we
caught," and **Horvitz-Thompson** weighting converts point results into an
unbiased area estimate with a confidence interval (each stratum's hit rate
is weighted by the area it represents).

Result: four validated REGIMES, each with its own correction factor —
George VI channel ponds (precision 0.625, correction ×1.0 — misses and
false alarms offset), Larsen C slush (×0.60), Amery frozen lake lids
(×0.50), Larsen B slush embayment (×0.35). Raw = "meltwater-affected area";
corrected = strict open water.

- Produces: error bars and per-regime corrections nobody else publishes.
- Weakness: ONE labeler (me→him). Leeson's 3-labeler protocol is the fix,
  on the roadmap.

## Layer 5 — Radar wet-snow record (the cloud-proof witness)

For each pixel and each satellite orbit track: build a WINTER BASELINE (the
median echo strength June–August, when everything is frozen solid). In melt
season, a pixel is WET when its echo drops ≥ 2 dB below its own winter self.
Per-orbit because look angle changes the echo; mixing tracks smears the
reference. Grid = the optical grid at 120 m (×4), same origin — so radar
and camera pixels align by construction, not by resampling luck.

Validated by four pre-registered gates: (1) 100.0% of blind-validated pond
cells read radar-wet (12,170 cells — alignment AND physics in one number);
(2) radar independently reproduces the optical season ranking (on
obs-normalized wet-exposure metrics; the naive max-extent metric SATURATES
because all of George VI wets every January — that failure is kept in the
gate record, decision 0002); (3) mid-winter false-wet 0.08% vs 2% allowed;
(4) melt onset Nov 2, end Feb 24 — climatologically sane.

- Produces: wet/dry state per 120 m cell per radar pass (~every 6 days).
- Weakness: wet ≠ ponded. Damp slush and a deep lake look identical. Hence:

## Layer 6 — Fusion (the two-witness rule)

Each radar pass classifies every cell: UNOBSERVED / DRY / WET / **PONDED**.
PONDED requires BOTH witnesses: radar says wet AND the camera saw a pond
there within ±6 days (evidence restricted to the Dec–Feb melt core —
outside it, optics sees frozen lake LIDS, which are stored water, not
ponding; decision 0004). By construction PONDED ⊆ WET — the code cannot
emit "ponded" where radar says dry. When optics says pond but radar says
dry, we COUNT it as a conflict instead of classifying it, and split it by
the cell's own radar season: in-season conflicts = real disagreement
(gated < 5%, measured 1.1–1.3%); off-season conflicts = the lid/refreeze
signal, reported as a product feature.

- Produces: the melt-state maps + the headline: wet:ponded = 10:1 in the
  record year, ~16,000:1 in the quiet year. **Ponding is nonlinear in melt.**
- Weakness: 120 m cells need ≥25% water fill → diffuse thin water is
  invisible to PONDED (deliberately conservative).

## Layer 7 — Depth & volume (the laser ruler)

Physics: light dims exponentially with depth (Beer–Lambert). For a pond
pixel: depth z = [ln(A_d − R_inf) − ln(R_w − R_inf)] / g, where R_w = what
the camera sees, A_d = the pond's floor brightness (estimated from its own
rim — the ice just outside the pond edge), R_inf = deep-water darkness
(0.05), and **g = the attenuation constant — the one unknown**. We never
guess g: ICESat-2 photon histograms over our ponds give two peaks (surface
& bottom); separation × 0.752 (light bends and slows in water — Parrish et
al. 2019) = measured depth. 111 such crossovers → fit g in closed form
(z = X/g is linear in 1/g), 3σ-trimmed.

The killer validation: George VI and Amery — 4,000 km apart — independently
fit g = 1.19 vs 1.15 (3% apart). The constant transfers across Antarctica.
Accuracy: RMSE 0.56 m per segment vs a 0.5 m target — a 12% near-miss,
reported as missed (decision 0003); tightening the time window below ±3
days doesn't help, so the floor is pairing geometry, not pond change.
Red band is only valid to ~3 m (deeper water goes black in red); green
extends deeper; the two bands' disagreement flags untrustworthy pixels.
Volume = Σ depth × pixel area, only where fused-PONDED (so lids and slush
never get fake depths), with Monte Carlo error bars (perturb g, rim albedo,
R_inf 500×).

- Produces: George VI 2019-20 = 0.129 km³ [0.108–0.166]; 2020-21 = 0.050
  [0.042–0.064]. Centre-deeper-than-rim sanity: 97–99% of ponds.
- Weakness: ±25% volumes; red-band saturation; segment-scale calibration.

## Layer 8 — Drainage catalogue & the honest nulls

Ponds that lose >80% of area between clear looks ≤14 days apart are
drainage/refreeze candidates, classified against radar (still wet after
loss = drained through the ice; frozen = benign refreeze). Coverage rule:
the pond must actually be OBSERVED on the second date — without it we got
784 phantom "drainages" that were just clouds and tile edges. With it, at
our ~2-week clear-scene cadence: ZERO detectable events, while literature
documents several. We report that as a cadence limitation, not a result.

- Produces: the method + the honest null.
- Weakness: the whole layer, and we say so. (Kingslake email asks the
  experts what cadence drainage detection truly needs.)

---

# DEPTH 3 — The hostile reviewer (and your answers)

**"Your precision is 0.35 on Amery. Why should I trust anything?"**
Because we MEASURED it — that's the point. 0.35 strict-open-water precision
on Amery means the raw number is meltwater-AFFECTED area (open water + lids
+ slush), and we publish the ×0.50 correction to strict open water with its
CI. Papers that don't do blind validation don't have better precision —
they have unknown precision. And lidded lakes ARE stored meltwater; for
hydrofracture, the affected area is arguably the relevant quantity.

**"Isn't one scene per season just a lower bound?"**
Yes, and we say so. The alternative — compositing many scenes — we tried,
and it inflated validated anchors (98→174, 781 phantom). A defensible lower
bound beats an indefensible bigger number. The radar record (every ~6 days,
cloud-proof) covers the sub-seasonal story.

**"Your RMSE missed its own gate."**
Correct: 0.56 m vs 0.5. We report it as missed, we widened the volume error
bars to carry it, and we located its source (segment-vs-pixel pairing, not
pond evolution — tightening time windows doesn't help). Same-day airborne
studies get ~0.3 m; for free-data automated harvesting at two sites, 0.56 m
with a 3% cross-site constant is honest state of the practice.

**"How much of this did AI build? Is this your work?"**
The code was AI-assisted throughout, documented in the README. Every
scientific decision — validation design, gate criteria, method rejections,
what got shipped vs killed — went through me and is recorded in five
decision documents, four of which exist because something FAILED and had to
be diagnosed. I can walk you through any layer on a whiteboard; ask me one.
(Then do it. This document is how.)

**"Why 2 dB and not 3?"**
Gate sweep: containment of blind-validated ponds was 1.000 / 0.999 / 0.994
at 2/3/4 dB, and mid-winter false-wet at 2 dB is 0.08% vs 2% allowed — the
permissive threshold costs nothing measurable, so we take the best
containment. Decision 0002.

**"Ranking ponded vs optical is apples and oranges — different units."**
Deliberately. Optical km² sums fractional water at 30 m; PONDED counts
120 m cells ≥25% full — a concentration measure. The claim is that the
ORDERING agrees (it does, exactly), not the magnitudes.

**"Why should anyone use this over Tuckett et al. 2025?"**
Complementary, not competing: they're continental and monthly but end in
2021 and are optical-only. We're 12 shelves but live through the current
season, blind-validated per regime, radar-cross-checked, with volumes.
(Stokes made this comparison for us — "perhaps not quite so efficiently as
your method.")

**"What's genuinely new here?"**
Three things: (1) the wet-vs-ponded decomposition at scale, and its
nonlinearity; (2) laser-calibrated ponded VOLUME time series with
uncertainty — essentially absent from the literature; (3) the observatory
model itself: living record, free compute, failures published.

---

# Where every headline number comes from

| Number | Origin file | One-line story |
|---|---|---|
| 12 shelves × 9 seasons | output/shelf/history*.json | region-agnostic driver over MEaSUReS boundaries |
| precision 0.625 [0.47–0.76] (GVI) | output/shelf_val/ + corrections.json | 80 blind stratified points, Wilson CI |
| 157±37 km² HT estimate | validation records | Horvitz-Thompson area correction vs 127 detected |
| 100.0% containment | output/sar/gates.json | validated ponds ∩ radar wet, ±6 d, 12,170 cells |
| 0.08% winter false-wet | gates.json gate3 | Jun–Jul 2021 median wet fraction at 2 dB |
| 796.6 / 186.2 / 0.5 km² ponded | output/sar/state_*_t2.json | fused PONDED peak extents, decisions 0004+0005 |
| conflicts 1.09% / 1.25% | state jsons | wet-season sensor disagreements / ponded instances |
| g = 1.62 (sites 1.19 vs 1.15) | output/depth/calibration.json | 111 ICESat-2 crossovers, closed-form fit |
| RMSE 0.56 m | calibration.json bands_tight | tight ±3 d pairs, red validity domain |
| 0.129 km³ [0.108–0.166] | output/depth/product_2019-20.json | Σ depth×area over PONDED, 500-draw Monte Carlo |
| 97–99% centre-deeper | product jsons | erosion-vs-rim depth comparison per pond |
| 784 → 0 drainage events | output/sar/drainage_*.json | coverage-blind vs coverage-aware detector |

---

# Glossary (terms you WILL be asked)

- **NDWI** — normalized difference water index; blue-vs-red contrast that
  makes water stand out.
- **NDSI** — snow index; melt depresses it, which is why record-melt scenes
  look "cloudy" to naive filters.
- **Hysteresis (detection)** — strict threshold to find pond cores, looser
  threshold to grow their connected margins.
- **Backscatter / gamma0 / dB** — radar echo strength; dB is log-scale, so
  "−2 dB" means the echo dropped ~37%.
- **Winter baseline** — a pixel's own frozen-season median echo; its
  personal "dry" reference.
- **Relative orbit** — a repeating satellite ground track; same track =
  same viewing geometry = comparable echoes.
- **STAC / COG** — the open catalog standard and cloud-optimized GeoTIFF
  format our products ship in; lets others stream our data without
  downloading everything.
- **Horvitz-Thompson** — survey-statistics estimator: weight each sampled
  point by the inverse of its inclusion probability → unbiased area
  estimates from stratified samples.
- **Wilson interval** — a confidence interval for proportions that behaves
  well at small n (better than the naive ±1.96√(pq/n)).
- **ATL03** — ICESat-2's geolocated photon product; every photon's lat/lon/
  elevation.
- **Strong beams** — ICESat-2 fires 3 beam pairs; the strong one of each
  pair has ~4× the photons — only those have enough bottom returns.
- **Refraction correction (×0.752)** — light bends and slows entering
  water; raw photon depth overestimates by ~33%.
- **Beer–Lambert** — exponential light attenuation with path length; the
  physics behind depth-from-darkness.
- **Monte Carlo uncertainty** — rerun the calculation hundreds of times
  with inputs perturbed within their plausible ranges; the spread of
  outputs IS the error bar.
- **EPSG:3031** — the polar stereographic map projection all our grids use;
  meters, not degrees, centered on the South Pole.
- **Halo fraction** — our cloud-mask coverage of the shelf's surroundings;
  the local-haze detector that tile-wide metadata misses.
- **Hydrofracture** — meltwater-filled crevasses wedging through the shelf;
  water is denser than ice, so a full crack keeps driving downward.

---

*Teaching sessions: work through Depth 2 one layer per session; close each
by answering that layer's hostile question from Depth 3 without the doc.
When all eight layers pass, run the full viva.*

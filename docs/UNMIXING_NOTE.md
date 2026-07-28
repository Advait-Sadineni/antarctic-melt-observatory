# Spectral unmixing — investigated, not adopted (yet)

Written during the accuracy-improvement session. This records why sub-pixel
spectral unmixing was investigated as the next accuracy step and why it was
**not** adopted, so the decision is not silently lost.

## The motivation

Blind reference-point validation showed the detector finds 0.89 of water
*features* (count-recall) but only 0.45 of water *area*, because it under-draws
pond **margins** — mixed water/ice pixels whose NDWI is diluted below
threshold. Linear spectral unmixing estimates each pixel's *fraction* of water
rather than a yes/no, which in principle recovers that sub-pixel margin area
and corrects the ~0.67x area under-report.

## What was tested

Two endmembers (water, ice) were derived from the blind-labelled points across
all three scenes, using the four 10 m bands (B2/B3/B4/B8) in reflectance:

```
water  (B2,B3,B4,B8) = [0.739, 0.618, 0.357, 0.177]   # dark, NIR-absorbing
ice    (B2,B3,B4,B8) = [0.967, 0.950, 0.905, 0.831]   # bright
```

Per-pixel water fraction was solved by constrained least squares. Findings:

- **Endmembers are physically correct** — water dark and falling toward the
  NIR, ice bright and flat.
- **Water fraction separates the labels well** — rank AUC 0.90 (water vs not).
- **Margins are only ~14% water.** Labelled water points in the margin stratum
  have median fraction 0.14, versus 0.84 for points in detected cores. The
  missed margin *area* is therefore much smaller than the missed *pixel count*
  suggests, which means the true area under-report is milder than the
  area-recall of 0.45 implies.

## Why it was not adopted

1. **It inherits the crevasse-shadow problem.** 10% of *not-water* points have
   water fraction > 0.72 — these are crevasse shadow, which unmixing confuses
   exactly as raw NDWI does. So unmixing would still need the shadow test; it
   is not an independent improvement to precision.
2. **Its one real benefit — fractional area — cannot be rigorously validated
   here.** Confirming a sub-pixel area estimate needs fractional ground truth,
   which point labels (binary, at pixel centres) do not provide. The cleanest
   proxy, cross-sensor fractional agreement, is not available either: only one
   of the three labelled scenes has a same-day Landsat pair, too few to derive
   robust Landsat TOA endmembers.

Adopting an area methodology that cannot be validated would break the rule
that has kept this project honest: no unvalidated number ships. So unmixing is
recorded as promising and deferred.

## When to revisit

- If fractional ground truth becomes available (hand-digitised pond outlines at
  sub-pixel precision, or higher-resolution imagery such as WorldView), the
  fractional area can be validated directly and adopted.
- If the pipeline moves to a learned model (roadmap item 3), a soft
  (probabilistic) output is the natural, trainable form of the same idea and
  can be validated the same way the current detector is.

The feasibility script is preserved in the session scratchpad; the endmembers
above are the reusable result.

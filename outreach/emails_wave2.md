# Outreach wave 2 — sent AFTER final numbers freeze

Placeholders: [WET_PONDED_RATIO], [VOL_2019], [VOL_2020], [TIGHT_N], [TIGHT_RMSE]
get filled from fusion_report.json / product_*.json before sending. Tone rules:
honest effort framing (built in ~2 weeks, CS undergrad), no inflation, one
specific question per email.

---

## 1. Amber Leeson (follow-up on her reply) — leeson thread, keep same subject

Hi Amber,

Thank you again for the labeling protocol — a quick follow-up with results,
since your note shaped the plan. The three-labeler weighted protocol is now
designed into the validation stack (recruiting the two additional labelers
this month), and in the meantime the pipeline finished its multi-sensor
phase:

- Fused Sentinel-1/-2 melt states for George VI: only ~[WET_PONDED_RATIO]%
  of radar-wetted area ever ponds — the hydrofracture-relevant fraction is
  far smaller than wet extent alone suggests.
- Pond depths calibrated against ICESat-2 photon-split crossovers at two
  sites 4,000 km apart (George VI and Amery) — attenuation coefficients
  agree to 3%, RMSE [TIGHT_RMSE] m over [TIGHT_N] tight crossovers, honestly
  reported against our 0.5 m target.

Everything is public (repo + DOI below). If any of it is useful to your
group's uncertainty work, I'd be glad to hear where it falls short.

Best,
Advait

## 2. Alison Banwell (CU Boulder / CIRES) — the George VI expert

Subject: ICESat-2-calibrated pond depths for George VI — student project,
would value your read

Dear Dr. Banwell,

I'm a computer science undergraduate at Penn State. Over the past two weeks
I've built an open pipeline that maintains a validated meltwater record for
12 Antarctic ice shelves from free satellite data, with George VI as the
blind-validated anchor — and because so much of the George VI literature is
yours, I wanted to share it and ask one question.

The pipeline fuses Sentinel-2 pond detection (blind stratified validation,
Horvitz-Thompson corrected) with Sentinel-1 wet-snow states, then retrieves
pond depths by dual-band attenuation calibrated against ICESat-2 photon
depths ([TIGHT_N] crossovers; George VI and Amery attenuation agree to 3%).
For 2019-20 George VI it gives [VOL_2019] km3 of ponded meltwater
([VOL_2020] km3 for 2020-21), with Monte Carlo uncertainty.

My question: is there any prospect of comparing against in-situ or
higher-fidelity depth data from your George VI field campaigns? Our laser
crossover RMSE is [TIGHT_RMSE] m at segment scale, and field-scale truth
would tell us what that number really means.

Repo: https://github.com/Advait-Sadineni/antarctic-melt-observatory
DOI: 10.5281/zenodo.21711608 · Live dashboard: [PAGES_URL]

Respectfully,
Advait Sadineni

## 3. Rajashree Tri Datta (CU Boulder) — ICESat-2 + Amery melt lakes

Subject: 82 ICESat-2 pond-depth crossovers on Amery/George VI from an open
pipeline — student project

Dear Dr. Datta,

Your ICESat-2 supraglacial lake work on Amery is the direct ancestor of
something I just built, so I wanted to share it. I'm a CS undergraduate at
Penn State; my open-source pipeline harvests ATL03 photon-split pond depths
automatically at optical crossovers (strong beams, Parrish refraction,
histogram peak splitting) and uses them to calibrate dual-band attenuation
depths — [TIGHT_N] tight crossovers so far across Amery and George VI, with
the two sites' attenuation coefficients agreeing to 3%.

One methodological question: at segment scale we plateau at ~[TIGHT_RMSE] m
RMSE against the laser, and tightening the crossover window below +/-3 days
stops helping — consistent with spatial pairing (photon-segment mean vs
150 m optical window) being the floor. Does that match your experience, and
is along-segment optical sampling the fix you'd recommend?

Repo: https://github.com/Advait-Sadineni/antarctic-melt-observatory
DOI: 10.5281/zenodo.21711608

Respectfully,
Advait Sadineni

## 4. CRYOLIST announcement (post-preprint ONLY)

Subject: [Data/Tool] Open Antarctic ice-shelf meltwater observatory —
12 shelves, 9 seasons, blind-validated, ICESat-2-calibrated depths

Dear all,

Announcing an open resource: a continuously updated meltwater record for 12
Antarctic ice shelves (2017-18 through the current season), built entirely
from free data and compute. Sentinel-2 detection with design-based blind
validation (per-regime precision/recall + Horvitz-Thompson corrections),
Sentinel-1 wet/dry cross-validation (4/4 gates), fused per-pixel melt states
(dry/wet/ponded with conflict accounting), and pond depths calibrated
against ICESat-2 photon crossovers at two sites.

Code, products, validation artifacts, and a browsable dashboard:
- Repo: https://github.com/Advait-Sadineni/antarctic-melt-observatory
- DOI: 10.5281/zenodo.21711608
- Dashboard: [PAGES_URL]
- Preprint: [PREPRINT_URL]

Feedback, use, and criticism all welcome — the validation design in
particular reports what fails as prominently as what passes.

Advait Sadineni (Penn State)

# Outreach drafts — Antarctic Meltwater Observatory

Six personalized cold emails. Rules applied: ≤ 6 sentences, one hook specific
to THEIR work, links (repo + DOI + one artifact), exactly one concrete
question, honest about student status and limitations. Placeholders:
`https://github.com/Advait-Sadineni/antarctic-melt-observatory` `https://doi.org/10.5281/zenodo.21711608` — filled at send time.

---

## 1. Alison Banwell (CIRES, CU Boulder) — the George VI hook

**Subject: Your 2019-20 George VI melt record — independently reproduced, and extended through 2025-26**

Dear Dr. Banwell,

I'm a computer science student who recently built an automated
meltwater-mapping pipeline on free Sentinel-2 data, and your George
VI work has been my ground truth throughout: my pipeline independently lands
on 19 January 2020 as the record melt date, and blind point-validation against
a hidden answer key puts its precision at 0.63 with a Horvitz-Thompson
bias-corrected area consistent with the detected 127 km² on the dominant tile.
The record now runs through the current season — and 2025-26 is George VI's
third-highest year in my series, with neighbouring Bach and Stange at their
all-time highs.

Code, per-season data, and the full validation trail (including what I got
wrong along the way): https://github.com/Advait-Sadineni/antarctic-melt-observatory · https://doi.org/10.5281/zenodo.21711608

One question, if you have 60 seconds: on Larsen C my blind validation shows
detections are ~40% slush rather than open water — does a ×0.6 open-water
correction for the eastern Peninsula match your field intuition, or is
counting slush the more defensible convention?

Thank you for the papers that made this buildable.
Advait Sadineni

---

## 2. Mariel Dirscherl (DLR) — the SAR fusion hook

**Subject: Sentinel-1 wet-snow detection independently confirming blind-validated Sentinel-2 ponds — 100% containment**

Dear Dr. Dirscherl,

Your Sentinel-1 supraglacial lake work convinced me the cloud problem in
optical melt mapping is solvable, so I built a per-pixel winter-baseline
(−3 dB, per relative orbit) wet-snow detector on Planetary Computer's S1-RTC
and gated it against my blind-validated Sentinel-2 pond record on George VI:
on first contact, radar confirmed 99.4-100% of the human-validated pond
pixels as wet across a 2-4 dB threshold sweep, and the retrieved melt season
(onset ~Nov 2, end ~Feb 24) matches the documented GVI calendar. Season-ranking
and winter-specificity gates are running now; no fusion is trusted until all
four pass.

Everything is open — code, gates, validation keys: https://github.com/Advait-Sadineni/antarctic-melt-observatory · https://doi.org/10.5281/zenodo.21711608

My one question: for HH-only IW scenes at ~71°S, is a fixed 3 dB threshold
below a per-orbit winter median defensible across melt regimes, or did your
work find the threshold needs to vary with incidence angle?

With thanks for lighting the path,
Advait Sadineni

---

## 3. Jennifer Arthur (Norwegian Polar Institute) — the East Antarctic hook

**Subject: 9-season Shackleton + Amery melt records (through 2025-26) — and they don't peak in the same years**

Dear Dr. Arthur,

Your Shackleton inventory and East Antarctic lakes review shaped how I built
my student project: an automated, blind-validated meltwater record that now
covers 10+ shelves × 9 seasons on fixed grids, through the current 2025-26
season. Two things I think might interest you: Shackleton and Amery both spike
in 2018-19 in my record (498 and ~2,957 km² meltwater-affected area
respectively), but Fimbul peaks in 2020-21 — and blind validation on Amery's
dominant tile shows the classic detection there is ~50% frozen lake lids and
slush by strict open-water standards, which I report as an explicit correction
rather than a single number.

Code, data, validation trail: https://github.com/Advait-Sadineni/antarctic-melt-observatory · https://doi.org/10.5281/zenodo.21711608

Question: is there a community convention you'd recommend for reporting
lidded-lake area vs open water — or is publishing both, with the blind
precision per regime, the right call?

Thank you for the review paper — it was my map into this field.
Advait Sadineni

---

## 4. Chris Stokes (Durham) — the record-currency hook

**Subject: East Antarctic shelves are melting hard right now — a live record through 2025-26**

Dear Professor Stokes,

Your circum-East-Antarctic lakes work is the reason my student project didn't
stop at one shelf: I've built an automated, blind-validated Sentinel-2 melt
record across 10+ shelves × 9 seasons on fixed comparable grids — and because
it runs on live archives, it extends through the *current* season. The 2025-26
signal is striking: Amery's second-largest year in my record (~1,000 km²),
with Bach and Stange (west Peninsula) at their all-time highs simultaneously.

Everything open, including the failure log and per-region blind validation:
https://github.com/Advait-Sadineni/antarctic-melt-observatory · https://doi.org/10.5281/zenodo.21711608

One question: for a student record like this to be genuinely useful to your
group, what's the one addition you'd prioritize — sub-seasonal time series,
more shelves, or independent labellers on the validation?

With thanks,
Advait Sadineni

---

## 5. Stewart Jamieson (Durham) — same team, complementary angle

**Subject: An open, region-agnostic pipeline for shelf meltwater — any Antarctic shelf from its MEaSUReS polygon**

Dear Professor Jamieson,

Following your and Chris Stokes' East Antarctic lakes work, I built a student
project that turns any Antarctic ice shelf's MEaSUReS boundary into a
9-season, blind-validated meltwater record automatically (Sentinel-2, free
archives, fixed EPSG:3031 grids) — currently 10+ shelves including Amery,
Shackleton and Fimbul, running through the present 2025-26 season. The
validation is the part I'm proudest of: stratified blind points with hidden
answer keys per region, which exposed three distinct melt regimes (channels /
slush / lidded lakes) needing different corrections.

Code and data: https://github.com/Advait-Sadineni/antarctic-melt-observatory · https://doi.org/10.5281/zenodo.21711608

Question: if you were pointing this at a scientific gap tomorrow, which
shelves or regions are most under-observed right now?

Best regards,
Advait Sadineni

---

## 6. Amber Leeson (Lancaster) — the data-science hook

**Subject: A data-science student's blind-validated Antarctic melt record — and its honest error budget**

Dear Dr. Leeson,

As someone at the lakes-and-data-science intersection, you may appreciate the
part of my student project I've worked hardest on: not the detection (published
Moussavi thresholds), but the error discipline around it — stratified blind
point validation with hidden keys per region, Horvitz-Thompson bias correction,
a documented catalogue of every way clouds fooled the pipeline before the
selection logic was hardened, and acceptance gates that have already correctly
rejected two "improvements" that looked good and weren't. The result: 10+
shelves × 9 seasons of comparable melt records through 2025-26, with per-region
precision (0.63 / 0.60 / 0.35) reported rather than hidden.

Code, data, and the full audit trail: https://github.com/Advait-Sadineni/antarctic-melt-observatory · https://doi.org/10.5281/zenodo.21711608

Question: my biggest known weakness is that the blind labeller and the author
are the same person — for a student project heading toward a preprint, what's
the lightest-weight path you'd suggest to independent labels?

Thank you,
Advait Sadineni

---

## Send checklist (do NOT send before all boxes tick)

- [ ] Repo public on GitHub, README rendering with hero figure
- [ ] Zenodo DOI minted (code + dataset v0.1)
- [ ] SAR gates 2-4 resolved (email #2 claims depend on them; soften if pending)
- [ ] Verify emails on staff pages: Banwell (CIRES), Leeson (Lancaster),
      Stokes/Jamieson (Durham)
- [x] Arthur jennifer.arthur@npolar.no confirmed
- [!] Dirscherl Mariel.Dirscherl@dlr.de BOUNCED 2026-07-30 (left DLR; ResearchGate
      stale). REPLACED with Katherine Deakin (Durham, katherine.a.deakin@durham.ac.uk),
      first author of the Frontiers 2025 GVI winter-SAR paper - our Gate-3 anchor;
      stronger hook. Backup co-author: Rebecca Dell (Cambridge SPRI).
- [ ] Advait reads each draft aloud once and edits anything that doesn't sound
      like him
- [ ] Send Tue-Thu, morning in recipient's timezone; no attachments, links only

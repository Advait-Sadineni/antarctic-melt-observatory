# 0005 — Evidence pools are evaluated per month

Date: 2026-08-03
Status: accepted

## Problem

Decision 0004's Dec-Feb evidence window broke 2019-20: the season-wide clean
pool preferred two dry, clean December scenes and excluded the hazy January
record-melt scenes (melt depresses NDSI, lifting halo into the 0.08-0.30
extreme-melt band that only branch 2 accepts - and branch 2 fires only when
branch 1 is empty). Production's Jan-Feb window had no clean scenes, so
branch 2 fired and captured the record; the widened evidence window silently
un-fired it. Result: fused PONDED evidence caught <=30 km2/date against
1,098 km2 of dense pond cells in the production product (97% missed), and
the season ranking inverted.

## Decision

`fusion.pond_series` evaluates `select_pool` PER MONTH (Dec, Jan, Feb; 2
scenes/month/tile), so each month's clean-vs-extreme-melt branching is
decided locally. Cache renamed ponds_<season>_corem.npz; the report's
construction check reads the same cache. All three pilot seasons re-walked
(third time; prior outputs archived pre0005/).

## Verification anchor

19DEA 2019-20: season-wide Dec-Feb pool = ['2019-12-21','2019-12-21'];
monthly January pool recovers ['2020-01-29','2020-01-19'] - production's
exact validated selections.

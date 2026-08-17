# A190502 a(57)–onward: values and timing for publication

See `data/b190502.txt` for the merged machine-readable table (published +
newly certified). Methodology notes:

## How each number was obtained

- **a(0)–a(56):** OEIS A190502 published data, taken from the official b-file
  (`data/b190502_known.txt`). Not recomputed here except as a validation check
  (see below).
- **a(57)–a(72): DONE.** Computed 2026-08-15 through 2026-08-17. Main anchors
  π(2ⁿ), π(2ⁿ⁻¹) independently cross-verified with two algorithms (Gourdon vs.
  Deleglise-Rivat) for every term. Computation was deliberately stopped after
  a(72) certified — this is the final term for now, a(73) onward is future
  work, not yet attempted (its prewarm header had just printed when the run
  was stopped; no computation happened for it and nothing was cached).
- **Reliability note:** three of these terms (a(69), and the two attempts
  toward a(72)) were interrupted mid-computation by unrelated machine
  reboots/crashes, unconnected to the pipeline itself. Every value below is
  still a clean, correct result — `pi_cache.json` persists every exact π(x)
  the moment it's computed, so a restart only ever re-did the one anchor call
  that was in flight, never corrupted or approximated anything. Timing for
  the affected terms should be read as "at least this long," not as clean
  single-pass numbers (see the Timing section below).

| n | a(n) | status |
|---|---|---|
| 57 | 1,838,131,803,685,114 | certified |
| 58 | 3,612,389,083,130,022 | certified |
| 59 | 7,101,392,627,457,670 | certified |
| 60 | 13,964,301,601,651,670 | certified |
| 61 | 27,467,389,548,130,805 | certified |
| 62 | 54,042,323,540,385,834 | certified |
| 63 | 106,356,804,267,969,409 | certified |
| 64 | 209,366,672,181,778,359 | certified |
| 65 | 412,246,861,431,389,466 | certified |
| 66 | 811,916,554,998,178,377 | certified |
| 67 | 1,599,434,686,587,771,626 | certified |
| 68 | 3,151,517,210,492,380,590 | certified |
| 69 | 6,211,093,021,215,415,074 | certified (crash-interrupted once, resumed cleanly) |
| 70 | 12,243,590,744,602,332,913 | certified |
| 71 | 24,140,116,236,766,772,570 | certified |
| 72 | 47,605,503,594,945,728,491 | certified (crash-interrupted once, resumed cleanly) |

## Timing

Two figures matter here: the "compute" time logged per term (anchors + sieve
walk + grid bracketing only) and the true wall-clock time including the
parallel checkpoint prewarm step, which is not included in the per-term
`seconds` field in `data/term_timings.csv` but is most of the real cost at
this magnitude:

| n | compute-only | true wall-clock (incl. prewarm) |
|---|---|---|
| 64 | 199s (3m19s) | ~22 min |
| 65 | 303s (5m03s) | ~36 min |
| 66 | 435s (7m15s) | ~48 min |
| 67 | 788s (13m08s) | ~90 min |
| 68 | 1091s (18m11s) | ~128 min |
| 69 | 1876s (31m16s) | not clean — crash-interrupted, excluded |
| 70 | 2419s (40m19s) | ~290 min (4h50m) |
| 71 | 3996s (66m36s) | ~469 min (7h49m) |
| 72 | 6008s (100m08s) | not clean — crash-interrupted once, excluded (the clean, post-restart prewarm alone was 576 min) |

Growth per term is noticeably faster than the naive $x^{2/3}$ estimate once
prewarm cost is counted — the checkpoint count itself is growing (14 → 19
grid points from n=64 to n=72), compounding with the per-`primecount`-call
cost increase. By n=71–72 a single term's wall-clock had grown to 8–11+ hours,
which is why the run was stopped at a(72).

## Validation against published data

Before trusting the pipeline on new terms, it was checked against OEIS's own
published values for **n = 20–52** (33 consecutive terms), freshly computed,
cache-bypassed for the anchors, all exact matches. 137 grid-bracketing
checkpoints were exercised across that range, so the bracketing lemma's loop
was genuinely tested, not skipped by short intervals. See
`data/term_timings.csv` for per-term wall times.

## Hardware

Laptop: 13th Gen Intel Core i5-13420H (8 cores / 12 threads), 13.7 GB RAM —
the same machine used for the A181671 (base-10) extension through a(20).

## Method

Certified via the bracketing lemma (see `README.md`): exact π(x) via
`primecount` at O(log) endpoints plus a small exact 128-bit sieve walk near Q,
backed by the published Johnston (2022) and Dusart (2010) analytic tail
bounds. Every new term additionally cross-verifies its record-defining
anchors with two independent `primecount` algorithms (Gourdon and
Deleglise-Rivat).

Base-2-specific: π(2ⁿ⁻¹) is exactly term (n−1)'s main anchor, so consecutive
runs get one of the two anchors for free from cache — but a cache hit is only
trusted as an anchor if it is marked cross-verified (`data/pi_verified.json`);
otherwise it is re-confirmed with the independent algorithm before use.

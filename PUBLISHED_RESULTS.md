# A190502 a(57)–onward: values and timing for publication

See `data/b190502.txt` for the merged machine-readable table (published +
newly certified). Methodology notes:

## How each number was obtained

- **a(0)–a(56):** OEIS A190502 published data, taken from the official b-file
  (`data/b190502_known.txt`). Not recomputed here except as a validation check
  (see below).
- **a(57)–a(67): DONE.** Freshly computed 2026-08-15, genuine single-pass
  timings, no interruptions, no cache warm start beyond what the validation
  sweep below legitimately produced. Main anchors π(2ⁿ), π(2ⁿ⁻¹) independently
  cross-verified with two algorithms (Gourdon vs. Deleglise-Rivat) for every
  term. The run was deliberately stopped after a(67) certified — a(68) onward
  is future work, not yet attempted.

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

Growth per term is noticeably faster than the naive $x^{2/3}$ estimate once
prewarm cost is counted — the checkpoint count itself is growing (14 → 16
grid points, 32 → 38 prewarm targets from n=64 to n=67), compounding with the
per-`primecount`-call cost increase.

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

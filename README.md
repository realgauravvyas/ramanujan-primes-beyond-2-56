<h1 align="center">Ramanujan Primes Beyond 2⁵⁶</h1>

<p align="center">
<b>Extending OEIS A190502 — computed entirely on a laptop, no institute, no cluster.</b>
</p>

<div align="center">

![Status](https://img.shields.io/badge/a(57)--a(67)-certified-6C63FF?style=flat-square)
![License](https://img.shields.io/badge/license-MIT-6C63FF?style=flat-square)
![Python](https://img.shields.io/badge/-Python-3776AB?style=flat-square&logo=python&logoColor=white)
![C++](https://img.shields.io/badge/-C%2B%2B-00599C?style=flat-square&logo=cplusplus&logoColor=white)

</div>

---

A Ramanujan prime $R_n$ is the least integer such that $\pi(x) - \pi(x/2) \geq n$
for all $x \geq R_n$ — Ramanujan's own 1919 strengthening of Bertrand's postulate.
[OEIS A190502](https://oeis.org/A190502) tabulates how many such primes lie below
each power of two, $A190502(n) = \pi_R(2^n)$, and publishes real data only through
$n=56$. This repo pushes that further — certified, cross-verified, and reproducible.

This is the base-2 sibling of
[`ramanujan-primes-beyond-1e19`](https://github.com/realgauravvyas/ramanujan-primes-beyond-1e19),
which extended the base-10 sequence [A181671](https://oeis.org/A181671) to $10^{23}$;
the certification method here is identical, only $Q=2^n$ instead of $Q=10^k$.

---

### 🔭 The result

| $n$ | $A190502(n)$ | status |
|---|---|---|
| 0–56 | — | OEIS A190502 published |
| 57 | 1,838,131,803,685,114 | **new**, cross-algorithm verified |
| 58 | 3,612,389,083,130,022 | **new**, cross-algorithm verified |
| 59 | 7,101,392,627,457,670 | **new**, cross-algorithm verified |
| 60 | 13,964,301,601,651,670 | **new**, cross-algorithm verified |
| 61 | 27,467,389,548,130,805 | **new**, cross-algorithm verified |
| 62 | 54,042,323,540,385,834 | **new**, cross-algorithm verified |
| 63 | 106,356,804,267,969,409 | **new**, cross-algorithm verified |
| 64 | 209,366,672,181,778,359 | **new**, cross-algorithm verified |
| 65 | 412,246,861,431,389,466 | **new**, cross-algorithm verified |
| 66 | 811,916,554,998,178,377 | **new**, cross-algorithm verified |
| 67 | 1,599,434,686,587,771,626 | **new**, cross-algorithm verified |

Full table, methodology, and per-term certificates: **[`certificates/`](certificates)**,
narrative in **[`PUBLISHED_RESULTS.md`](PUBLISHED_RESULTS.md)**.

Also independently checked: **33 consecutive published terms, $n=20$–$52$, all
exact matches** against the OEIS b-file, exercising 137 grid-bracketing
checkpoints — proof the method reproduces known data before it's trusted for
new data.

---

### 🧭 How it works

$\pi_R(x) = \min_{y \geq x} f(y)$ with $f(y) = \pi(y) - \pi(y/2)$, for *any* $x$
— nothing in that identity requires $x$ to be a power of ten, so the exact same
machinery that certified A181671 certifies A190502 with $Q=2^n$ instead. Direct
sieving is infeasible at this scale — the certification window spans far more
than the $2^{64}$ limit of general-purpose sieve libraries. Instead:

1. **Bracketing lemma** — exact $\pi(x)$ at $O(\log)$ grid endpoints
   (via [`primecount`](https://github.com/kimwalisch/primecount)) lower-bounds
   $f(x)$ across entire intervals, no sieving of the interior.
2. **Custom 128-bit segmented sieve** (`sieve128.cpp`) — handles the small region
   immediately above $Q$ exactly.
3. **Analytic tail bounds** — [Dusart (2010)](https://arxiv.org/abs/1002.0442) and
   [Johnston (2022)](https://arxiv.org/abs/2109.02249) certify everything beyond
   the sieved-and-bracketed region, for free (valid to $1.101\times10^{26}$, which
   caps this method at roughly $n\approx86$).
4. **Cross-algorithm verification** — every new term's defining anchors are
   computed independently by two `primecount` algorithms (Gourdon and
   Deleglise–Rivat) and must agree exactly.
5. **Free anchors in base 2** — $\pi(2^n/2) = \pi(2^{n-1})$ is *exactly* the
   previous term's main anchor, so consecutive terms reuse one of the two
   expensive anchors every time. A cached-but-not-yet-cross-verified anchor is
   still confirmed with the independent algorithm before it's trusted — see
   `pi_verified.json`.

---

### 📁 Repo layout

| path | contents |
|---|---|
| `src/` | `ramanujan2.py` (orchestration driver), `sieve128.cpp` (128-bit sieve), `build.sh` |
| `certificates/` | Machine-checkable JSON certificates, one per term |
| `data/` | `pi_cache.json` (every exact π(x) computed), `pi_verified.json` (which were cross-checked), b-file, timings |
| `PUBLISHED_RESULTS.md` | Term-by-term methodology and provenance notes |

---

### ⚙️ Reproducing

```bash
git clone --depth 1 https://github.com/kimwalisch/primecount
cmake -S primecount -B primecount/build -DCMAKE_BUILD_TYPE=Release
cmake --build primecount/build -j$(nproc)
pip install mpmath

g++ -O3 -march=native -fopenmp -std=c++17 src/sieve128.cpp -o sieve128 -lprimesieve

python3 src/ramanujan2.py validate --from 20 --to 52   # reproduce OEIS-published terms first
python3 src/ramanujan2.py run 57 --largest --threads 12 \
    --primecount primecount/build/primecount
```

`data/pi_cache.json` seeds every exact π(x) already computed, so reproducing any
term in `certificates/` should mostly be cache hits.

---

### ✅ Verified, not just computed

- Anchors cross-checked by two independent `primecount` algorithms — exact agreement, every term
- 33 consecutive published A190502 terms ($n=20$–$52$) independently reproduced, exact match
- The base-2 anchor-sharing optimization never bypasses cross-verification: cached
  values must be marked verified before they can anchor a new term (`pi_verified.json`)
- Dusart/Johnston bound formulas carried over unchanged from the already-verified
  A181671 pipeline, not re-derived from scratch

---

### 🤖 How this was built

The code (`ramanujan2.py`, adapted from the A181671 driver, and the orchestration
tooling) was developed with the assistance of Claude (Anthropic), directed and
reviewed throughout by the author, who verified the mathematical claims and
cross-checked results against the published OEIS b-file and independent
`primecount` algorithms.

---

### License

[MIT](LICENSE).

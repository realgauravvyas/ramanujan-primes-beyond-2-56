#!/usr/bin/env python3
"""
ramanujan2.py -- certified pi_R(2^n), i.e. OEIS A190502, for n far beyond 2^64.

    A190502(n) = number of Ramanujan primes <= 2^n
               = pi_R(2^n)
               = min_{y >= 2^n} f(y),      f(y) = pi(y) - pi(floor(y/2))

That last equality is the whole reason this is computable.  R_k (the k-th
Ramanujan prime) is by definition the least integer such that f(y) >= k for
every y >= R_k.  Hence R_k <= x  <=>  min_{y>=x} f(y) >= k, and therefore
pi_R(x) = #{k : R_k <= x} = min_{y>=x} f(y) exactly -- for ANY x, not just
powers of ten.  So the identical machinery that produced A181671 (pi_R(10^k))
produces A190502 with no mathematical change at all, only x = 2^n.

Three regions, three tools (unchanged from the A181671 pipeline):

  [Q, Q+D]      exact sieved walk         (sieve128, 128-bit)
  [Q+D, Y]      GRID BRACKETING           (exact pi at O(log) points only)
  [Y, inf)      Johnston / Dusart bounds  (analytic)

  BRACKETING LEMMA.  For a <= y <= b,
        f(y) = pi(y) - pi(y//2) >= pi(a) - pi(b//2)
  since pi is nondecreasing.  So exact pi values at the two ENDPOINTS
  lower-bound f on the whole interval -- no sieving inside at all.

WHAT IS DIFFERENT FOR BASE 2
----------------------------
1. Labelling.  The A181671 driver derived its term index as round(log10(Q)),
   which is correct only when Q is a power of ten.  Fed 2^57 it would have
   labelled the result "a(17)" and written "Q = 10^17" into the notes and
   certificate.  Here n is carried explicitly end to end and Q is always
   printed as an exact integer.

2. Free anchors.  For Q = 2^n the half-anchor is pi(Q/2) = pi(2^(n-1)),
   which is *exactly* the main anchor of term n-1.  Running terms
   consecutively therefore reuses one of the two expensive anchors every
   time -- a structural saving the powers-of-ten sequence cannot get
   (5*10^(k-1) is never a power of ten).

3. Because of (2), a cache hit can now legitimately supply a
   record-defining anchor.  The A181671 code returned early on any cache
   hit, which would silently skip the two-algorithm cross-check.  Cached
   values now carry a "was this cross-verified" flag (pi_verified.json) and
   an unverified hit is cross-checked before it is allowed to anchor a term.

Usage:
  python3 ramanujan2.py validate                 # reproduce published A190502
  python3 ramanujan2.py run 57                   # compute a(57), i.e. Q=2^57
  python3 ramanujan2.py run 57 --largest         # also find the record prime
  python3 ramanujan2.py explore 57               # 57, 58, 59, ... until stopped
  python3 ramanujan2.py bfile                    # emit merged b190502.txt
"""
import argparse, json, os, subprocess, sys, threading, time
from concurrent.futures import ThreadPoolExecutor
from mpmath import mp, mpf, li, log, sqrt, floor

mp.dps = 50
HERE = os.path.dirname(os.path.abspath(__file__))
SIEVE = os.path.join(HERE, "sieve128.exe" if os.name == "nt" else "sieve128")
CACHE = os.path.join(HERE, "pi_cache.json")
VERIFIED = os.path.join(HERE, "pi_verified.json")
TIMING_LOG = os.path.join(HERE, "term_timings.csv")
NOTES = os.path.join(HERE, "NOTES.md")
KNOWN_FILE = os.path.join(HERE, "b190502_known.txt")
JOHNSTON_LIMIT = mpf('1.101e26')

# The published terms of A190502, n -> a(n), loaded from the OEIS b-file that
# ships next to this script (b190502_known.txt, n = 0..56).  Anything computed
# beyond the largest n in that file is a new term.
def load_known():
    known = {}
    if os.path.exists(KNOWN_FILE):
        with open(KNOWN_FILE, encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                a, b = line.split()
                known[int(a)] = int(b)
    return known

KNOWN = load_known()
LAST_PUBLISHED = max(KNOWN) if KNOWN else -1


def _log_timing(label, Q, seconds, grid_points):
    new = not os.path.exists(TIMING_LOG)
    with open(TIMING_LOG, "a", encoding="utf-8") as f:
        if new:
            f.write("term,n,Q,seconds,minutes,grid_points,finished_at\n")
        n = label
        f.write(f"{label},{n},{Q},{seconds},{seconds/60:.2f},{grid_points},"
                f"{time.strftime('%Y-%m-%d %H:%M:%S')}\n")


def parse_exact_int(s):
    """Parse '1e23', '100000000000000000000000', '1.5e9' etc. EXACTLY.

    int(float(s)) silently corrupts large scientific notation: float('1e23')
    is 99999999999999991611392, not 10^23. Decimal parses the literal exactly
    at any size.
    """
    from decimal import Decimal
    v = Decimal(str(s).replace("_", ""))
    iv = int(v)
    if iv != v:
        raise ValueError(f"{s!r} is not an integer")
    return iv


# --------------------------------------------------------------- pi(x), cached
_cache = {}
_verified = set()
_cache_lock = threading.Lock()          # prewarming writes from several threads


def atomic_write_json(path, obj, **json_kwargs):
    """Write JSON so a kill mid-write can never leave a truncated/corrupt file."""
    d = os.path.dirname(os.path.abspath(path)) or "."
    fd, tmp = __import__("tempfile").mkstemp(dir=d, prefix=".tmp_", suffix=".json")
    try:
        with os.fdopen(fd, "w") as f:
            json.dump(obj, f, **json_kwargs)
        os.replace(tmp, path)          # atomic rename, same filesystem
    except Exception:
        try:
            os.remove(tmp)
        except OSError:
            pass
        raise


def _load():
    global _cache, _verified
    if os.path.exists(CACHE):
        _cache = json.load(open(CACHE))
    if os.path.exists(VERIFIED):
        _verified = set(json.load(open(VERIFIED)))


def _save():
    atomic_write_json(CACHE, _cache)


def _save_verified():
    atomic_write_json(VERIFIED, sorted(_verified))


def _run_primecount(x, primecount, threads, verbose, tag, extra_args=()):
    """One primecount subprocess call; returns the exact pi(x) value."""
    cmd = [primecount, str(x), "--status", *extra_args]
    if threads:
        cmd.append(f"--threads={threads}")
    t0 = time.time()
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            text=True, bufsize=1)
    last_num = None
    last_print = 0.0
    if verbose:
        print(f"      {tag}: counting exact primes (can take a while)...", end="", flush=True)
    for line in proc.stdout:
        s = line.strip()
        if s.isdigit():
            last_num = s
            continue
        now = time.time()
        if verbose and now - last_print >= 2.0:
            last_print = now
            print(f"\r      {tag}: still counting...  {now - t0:5.0f}s elapsed   ", end="", flush=True)
    proc.wait()
    if proc.returncode != 0 or last_num is None:
        if verbose:
            print()
        raise RuntimeError(f"primecount failed on {x} (exit {proc.returncode})")
    v = int(last_num)
    if verbose:
        print(f"\r      {tag}: done = {v:,}   ({time.time()-t0:.0f}s)" + " " * 25, flush=True)
    return v


def pi(x, primecount="primecount", threads=None, verbose=True, label=None,
      verify_independent=False):
    """Exact pi(x).  Cached, and cross-verifiable.

    If verify_independent, the value must be confirmed by TWO different
    primecount algorithms (Gourdon, the default, and Deleglise-Rivat).  A
    value already in the cache is NOT trusted for this purpose unless it is
    recorded in pi_verified.json -- otherwise a cheap prewarmed checkpoint
    could silently become a record-defining anchor without ever being
    cross-checked.  This matters much more in base 2 than it did in base 10,
    because pi(2^n / 2) = pi(2^(n-1)) really is a previous term's anchor and
    really will be a cache hit.
    """
    k = str(x)
    tag = label or f"pi({x:.3e})"
    cached = _cache.get(k)

    if cached is not None and (not verify_independent or k in _verified):
        if verbose:
            extra = " [previously cross-verified]" if k in _verified else ""
            print(f"      {tag}: already known (cached) = {int(cached):,}{extra}", flush=True)
        return int(cached)

    if cached is not None:
        # Cache hit, but this call needs a cross-verified value and this entry
        # has never been cross-checked.  Confirm it rather than recompute it.
        v = int(cached)
        if verbose:
            print(f"      {tag}: cached = {v:,}, but not yet cross-verified "
                  f"-- confirming with an independent algorithm", flush=True)
    else:
        v = _run_primecount(x, primecount, threads, verbose, tag)

    if verify_independent:
        v2 = _run_primecount(x, primecount, threads, verbose,
                             tag + " [cross-check, Deleglise-Rivat]",
                             extra_args=["--deleglise-rivat"])
        if v2 != v:
            raise RuntimeError(f"CROSS-CHECK MISMATCH for {tag}: "
                               f"Gourdon={v:,} vs Deleglise-Rivat={v2:,}")
        if verbose:
            print(f"      {tag}: cross-check PASSED -- Gourdon and Deleglise-Rivat agree",
                  flush=True)
        with _cache_lock:
            _verified.add(k)
            _save_verified()

    with _cache_lock:
        _cache[k] = str(v)
        _save()
    return v


# --------------------------------------------------- parallel checkpoint prewarm
# A single primecount process only saturates ~4 of 12 logical cores (its Sigma
# and Phi0 phases are largely serial), so most of the machine sits idle.  The
# grid checkpoints are deterministic -- offset_{n+1} = int(offset_n * growth) --
# and the two calls per checkpoint are independent, so they can be computed
# ahead of the main loop, several at a time, and handed over via _cache.
#
# Only CHECKPOINTS are prewarmed, never the anchors.  Prewarmed values are
# never added to _verified, so even if one did coincide with an anchor, pi()
# would still force the cross-check before using it.

def checkpoint_targets(Q, D, growth, n):
    """The exact pi() arguments certified_count's grid loop will ask for."""
    targets, y = [], Q + D
    for _ in range(n):
        step = int((y - Q) * growth)
        nxt = Q + step
        targets.append(nxt // 2)
        targets.append(nxt)
        y = nxt
    return targets


def predict_checkpoints(Q, D, growth, margin=2):
    """How many grid checkpoints the loop will need (li-estimate of the count).
    Margin covers estimate slop."""
    Q = mpf(Q)
    m = int(li(Q) - li(Q / 2))
    y = Q + D
    n = 0
    while int(floor(G(y))) < m and n < 90:
        n += 1
        y = Q + (y - Q) * mpf(growth)
    return n + margin


def _free_ram_gb():
    """Physically available RAM, or None if we cannot tell."""
    if os.name == "nt":
        try:
            import ctypes
            class MS(ctypes.Structure):
                _fields_ = [("dwLength", ctypes.c_ulong), ("dwMemoryLoad", ctypes.c_ulong),
                            ("ullTotalPhys", ctypes.c_ulonglong), ("ullAvailPhys", ctypes.c_ulonglong),
                            ("ullTotalPageFile", ctypes.c_ulonglong), ("ullAvailPageFile", ctypes.c_ulonglong),
                            ("ullTotalVirtual", ctypes.c_ulonglong), ("ullAvailVirtual", ctypes.c_ulonglong),
                            ("ullAvailExtendedVirtual", ctypes.c_ulonglong)]
            m = MS(); m.dwLength = ctypes.sizeof(MS)
            ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(m))
            return m.ullAvailPhys / (1024 ** 3)
        except Exception:
            return None
    try:
        info = {}
        with open("/proc/meminfo") as f:
            for line in f:
                k, _, rest = line.partition(":")
                info[k] = rest.strip()
        kb = info.get("MemAvailable") or info.get("MemFree")
        if kb:
            return int(kb.split()[0]) / (1024 ** 2)
    except Exception:
        pass
    try:
        import psutil
        return psutil.virtual_memory().available / (1024 ** 3)
    except Exception:
        return None


def safe_workers(Q, requested, verbose=True):
    """Cap concurrency so the batch cannot exhaust RAM.

    Footprint measured across two machines during the A181671 work: 0.36 GB at
    5e20, 0.95 GB at 5e21, 1.26 GB at 1e22 (laptop) and 5.9 GB at 1e23
    (desktop).  Exponent 0.67 anchored +15% margin at Q=1e22 matches both the
    1e22 and 1e23 measurements and stays conservative at smaller Q.  An OOM
    kill mid-call would waste hours.
    """
    free = _free_ram_gb()
    if free is None:
        return requested
    per = max(0.2, 1.45 * (float(Q) / 1e22) ** 0.67)
    fits = max(1, int(free * 0.95 / per))
    n = max(1, min(requested, fits))
    if verbose and n < requested:
        print(f"  [ram] {free:.1f} GB free, ~{per:.1f} GB per call -> "
              f"using {n} workers instead of {requested}", flush=True)
    return n


def prewarm(Q, D, growth, primecount, workers, threads, verbose=True):
    """Fill _cache with checkpoint values, `workers` primecount calls at once."""
    workers = safe_workers(Q, workers, verbose)
    n = predict_checkpoints(Q, D, growth)
    targets = [t for t in checkpoint_targets(Q, D, growth, n)
               if str(t) not in _cache]
    if not targets:
        return
    if verbose:
        print(f"  [0/3] prewarming ~{n} checkpoints "
              f"({len(targets)} values, {workers} at a time on {threads} threads each)",
              flush=True)
    t0 = time.time()
    done = {"n": 0}

    def one(x):
        r = subprocess.run([primecount, str(x), f"--threads={threads}"],
                           capture_output=True, text=True)
        v = r.stdout.strip()
        if r.returncode != 0 or not v.isdigit():
            return                      # leave it; the main loop will compute it
        with _cache_lock:
            _cache[str(x)] = v
            _save()                     # atomic (see atomic_write_json)
            done["n"] += 1
            k, tot = done["n"], len(targets)
            el = time.time() - t0
            if verbose:
                print(f"      prewarm {k}/{tot}   elapsed {el/60:.0f}m   "
                      f"ETA {el/k*(tot-k)/60:.0f}m", flush=True)

    with ThreadPoolExecutor(max_workers=workers) as ex:
        list(ex.map(one, targets))
    if verbose:
        print(f"  [0/3] prewarm done [{(time.time()-t0)/60:.0f}m]", flush=True)


# ------------------------------------------------------- analytic tail bounds
def G(y):
    """Johnston 2022 lower bound for f(y); valid 5314 <= y <= 1.101e26."""
    y = mpf(y)
    err = (sqrt(y) * log(y) + sqrt(y / 2) * log(y / 2)) / (8 * mp.pi)
    return li(y) - li(y / 2) - err


def dusart_ok(m):
    """f(y) >= m for all y >= 1.101e26, via Dusart 2010."""
    y = JOHNSTON_LIMIT; L = log(y); L2 = log(y / 2)
    lb = y / L * (1 + 1 / L + 2 / L**2)
    ub = (y / 2) / L2 * (1 + 1 / L2 + mpf('2.334') / L2**2)
    return (lb - ub) >= m


# ----------------------------------------------------------- exact local walk
def walk(Q, L, S, verbose=True):
    """Exact walk of f on (Q, Q+L]; returns (min_rel, min_off, end_rel)."""
    t0 = time.time()
    n_chunks = max(1, -(-L // S))
    if verbose:
        print(f"      scanning the {L:.2e}-wide region right after Q "
              f"(~{n_chunks} chunks)...", end="", flush=True)
    proc = subprocess.Popen([SIEVE, "walk", str(Q), str(L), str(S)],
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            text=True, bufsize=1)
    seg = 0
    for line in proc.stderr:
        if line.strip().startswith("seg"):
            seg += 1
            if verbose:
                print(f"\r      scanning near Q... chunk {seg}/{n_chunks}   "
                      f"({time.time()-t0:.0f}s)   ", end="", flush=True)
    out, _ = proc.communicate()
    if proc.returncode != 0:
        if verbose:
            print()
        raise RuntimeError(f"sieve128 walk failed (exit {proc.returncode})")
    d = dict(l.split() for l in out.strip().split("\n"))
    if verbose:
        print(f"\r      scanned near Q: done   ({time.time()-t0:.0f}s)" + " " * 25, flush=True)
    return int(d["min_rel"]), int(d["min_off"]), int(d["end_rel"])


# ------------------------------------------------------------ main certifier
def certified_count(n, D=10**9, S=10**8, growth=1.7, threads=None,
                    primecount="primecount", verbose=True,
                    verify_anchors=False):
    """Certified pi_R(2^n) = A190502(n).  Returns (count, info dict).

    verify_anchors=True cross-checks the two record-defining anchors
    (pi(Q), pi(Q/2)) with an independent algorithm -- use for genuinely new
    terms; not needed when replaying an already-published one.
    """
    Q = 2 ** n
    t0 = time.time()
    label = f"a({n})"
    if verbose:
        print()
        print("=" * 70)
        print(f"  NOW WORKING ON: A190502 {label}   (Q = 2^{n} = {Q:,})")
        if verify_anchors:
            print("  (main anchors will be independently cross-verified)")
        print("=" * 70)
        print("  step 1 of 3 -- exact prime counts near Q")
    # pi(Q/2) here is pi(2^(n-1)) -- literally term (n-1)'s main anchor, so on a
    # consecutive run this is a free, already-cross-verified cache hit.
    piQ = pi(Q, primecount, threads, verbose, label=f"pi(2^{n})",
             verify_independent=verify_anchors)
    piH = pi(Q // 2, primecount, threads, verbose, label=f"pi(2^{n-1}) [= Q/2]",
             verify_independent=verify_anchors)
    fQ = piQ - piH
    if verbose:
        print("  step 2 of 3 -- scanning the region right after Q")
    mn_rel, mn_off, end_rel = walk(Q, D, S, verbose)
    m = fQ + mn_rel                       # running minimum, exact on [Q, Q+D]
    y = Q + D
    f_y = fQ + end_rel                    # exact f at the right edge
    pi_half = {Q // 2: piH,
               y // 2: pi(y // 2, primecount, threads, verbose, label="pi((Q+D)/2)")}
    if verbose:
        print("  step 3 of 3 -- extending the safety net out to the analytic tail")
    grid = 0
    JLIM = 1101 * 10**23                  # 1.101e26, exact: G()'s validity limit
    while True:
        if y > JLIM:
            raise RuntimeError(
                "grid walked past the Johnston bound's validity limit "
                "(1.101e26); this n is beyond the method's ceiling")
        if int(floor(G(y))) >= m:
            break
        step = int((y - Q) * growth)
        nxt = Q + step
        # bracketing lemma: min over [y, nxt] >= f(y) - (pi(nxt//2) - pi(y//2))
        while True:
            if nxt // 2 not in pi_half:
                pi_half[nxt // 2] = pi(nxt // 2, primecount, threads, verbose,
                                       label=f"checkpoint {grid + 1} (pi half)")
            ph = pi_half[nxt // 2]
            bound = f_y - (ph - pi_half[y // 2])
            if bound >= m:
                break
            # bound too weak: refine by halving the increment BEYOND y (not the
            # offset from Q -- halving `step` itself would move nxt backward and
            # stall immediately).  Refinement lands the checkpoint somewhere
            # prewarm() did not predict; that only wastes a prewarmed value, it
            # cannot make the result wrong.
            inc = (nxt - y) // 2
            if inc < 10**6:
                raise RuntimeError("grid refinement stalled; raise D")
            if verbose:
                print(f"\n    [refine] bound {bound:,} < min {m:,} at Q+{nxt-Q:.3e}; "
                      f"halving increment", flush=True)
            nxt = y + inc
        pn = pi(nxt, primecount, threads, verbose, label=f"checkpoint {grid + 1}")
        f_n = pn - pi_half[nxt // 2]
        if f_n < m:                       # cannot happen if bound held, but check
            raise RuntimeError("grid endpoint below running minimum")
        grid += 1
        if verbose:
            print(f"    checkpoint {grid} verified safe out to Q+{nxt-Q:.2e}  [OK]", flush=True)
        y, f_y = nxt, f_n
    assert dusart_ok(m), "far-tail check failed"
    info = dict(sequence="A190502", n=n, Q=Q, Q_repr=f"2^{n}",
                pi_Q=piQ, pi_halfQ=piH, f_at_Q=fQ, count=m,
                min_offset=mn_off, walk_width=D, grid_points=grid,
                window_end=y, seconds=round(time.time() - t0),
                anchors_cross_verified=verify_anchors)
    _log_timing(label, Q, info['seconds'], grid)
    if verbose:
        mins, secs = divmod(info['seconds'], 60)
        print("-" * 70)
        print(f"  RESULT  A190502 {label} = {m:,}")
        print(f"  took {mins}m {secs}s, {grid} safety checkpoints")
        print("-" * 70)
    return m, info


def _write_notes(n, m, info):
    """Append a human-readable record of a newly-found term to NOTES.md."""
    is_new_file = not os.path.exists(NOTES)
    with open(NOTES, "a", encoding="utf-8") as f:
        if is_new_file:
            f.write("# A190502 -- Ramanujan primes <= 2^n -- new terms\n\n")
            f.write(
                f"OEIS A190502 publishes data through a({LAST_PUBLISHED}). "
                "Everything recorded below is a new, first-time result computed "
                "with `ramanujan2.py`, the base-2 driver for the certified "
                "pi_R pipeline previously used to extend A181671 (pi_R(10^k)) "
                "to k=23.\n\n"
                "Method: pi_R(x) = min_{y>=x} f(y) with f(y) = pi(y) - pi(y/2). "
                "A lower bound is proved with the bracketing lemma -- exact "
                "pi(x) (via `primecount`) at O(log) endpoints, plus an exact "
                "128-bit sieve walk (`sieve128`) in a small window just above "
                "Q -- backed by the published Johnston (2022) and Dusart (2010) "
                "analytic tail bounds.\n\n")
        f.write(f"## a({n}) = {m:,}\n\n")
        f.write(f"- Q = 2^{n} = {info['Q']:,}\n")
        f.write(f"- found: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        mins, secs = divmod(info['seconds'], 60)
        f.write(f"- wall time: {mins}m {secs}s ({info['seconds']}s)\n")
        f.write(f"- safety checkpoints (grid bracketing): {info['grid_points']}\n")
        f.write(f"- main anchors pi(2^{n}), pi(2^{n-1}) cross-verified with two "
                f"independent primecount algorithms (Gourdon vs Deleglise-Rivat): "
                f"{'yes, agreed' if info.get('anchors_cross_verified') else 'no'}\n")
        if info.get("largest_R"):
            f.write(f"- largest Ramanujan prime <= 2^{n}: {info['largest_R']:,}\n")
        f.write(f"- certificate file: {cert_name(n)}\n\n")


def cert_name(n):
    return f"piR_2pow{n}_certificate.json"


def largest_R(Q, count, L=10**9, S=10**8, verbose=True):
    """Largest Ramanujan prime <= Q (its index is `count`)."""
    piQ = int(_cache[str(Q)]); piH = int(_cache[str(Q // 2)])
    delta = count - 1 - (piQ - piH)
    r = subprocess.run([SIEVE, "back", str(Q), str(L), str(S), str(delta)],
                       stdout=subprocess.PIPE, stderr=None, text=True)
    if r.returncode != 0 or r.stdout.startswith("NOTFOUND"):
        return None
    x = int(r.stdout.split()[1])
    return x + 1


def finish_term(n, m, info, D, S, want_largest, primecount):
    """Shared post-processing for a completed term."""
    if want_largest:
        print("Searching for the actual largest Ramanujan prime <= Q "
              "(quick sieve scan, no more heavy counting)...")
        R = largest_R(2 ** n, m, L=D, S=S)
        if R:
            print(f"largest Ramanujan prime <= 2^{n}:  R_{m:,} = {R:,}")
            info["largest_R"] = R
    out = os.path.join(HERE, cert_name(n))
    atomic_write_json(out, info, indent=1, default=str)
    print("certificate ->", out)
    _write_notes(n, m, info)
    print(f"note appended -> {NOTES}\n")


def write_bfile(path=None):
    """Merge published terms with everything certified here into a b-file."""
    path = path or os.path.join(HERE, "b190502.txt")
    terms = dict(KNOWN)
    found = {}
    for fn in os.listdir(HERE):
        if fn.startswith("piR_2pow") and fn.endswith("_certificate.json"):
            info = json.load(open(os.path.join(HERE, fn)))
            found[int(info["n"])] = int(info["count"])
    for k, v in found.items():
        if k in terms and terms[k] != v:
            raise RuntimeError(f"certificate for a({k})={v} contradicts "
                               f"published value {terms[k]}")
        terms[k] = v
    with open(path, "w", encoding="utf-8") as f:
        for k in sorted(terms):
            f.write(f"{k} {terms[k]}\n")
    new = sorted(k for k in found if k > LAST_PUBLISHED)
    print(f"wrote {path}  ({len(terms)} terms, n = {min(terms)}..{max(terms)})")
    if new:
        print(f"new terms beyond published a({LAST_PUBLISHED}): "
              + ", ".join(f"a({k})" for k in new))
    else:
        print(f"no new terms yet (published through a({LAST_PUBLISHED}))")
    return path


# ------------------------------------------------------------------- CLI
def main():
    ap = argparse.ArgumentParser(
        description="Certified A190502(n) = number of Ramanujan primes <= 2^n")
    ap.add_argument("mode", choices=["validate", "run", "explore", "bfile"])
    ap.add_argument("n", nargs="?", type=int, default=None,
                    help="base-2 exponent (Q = 2^n)")
    ap.add_argument("--D", default="1e9", help="exact-walk width just above Q")
    ap.add_argument("--S", default="1e8", help="sieve segment size")
    ap.add_argument("--growth", type=float, default=1.7)
    ap.add_argument("--threads", type=int, default=None)
    ap.add_argument("--primecount",
                    default=os.path.join(HERE, "primecount.exe" if os.name == "nt"
                                         else "primecount"))
    ap.add_argument("--largest", action="store_true",
                    help="also locate the largest Ramanujan prime <= 2^n")
    ap.add_argument("--from", dest="lo", type=int, default=20,
                    help="validate: first n to check")
    ap.add_argument("--to", dest="hi", type=int, default=44,
                    help="validate: last n to check")
    ap.add_argument("--stop", type=int, default=None,
                    help="explore: last n to compute (default: run forever)")
    ap.add_argument("--workers", type=int, default=3,
                    help="concurrent primecount processes during prewarm")
    ap.add_argument("--pc-threads", type=int, default=4,
                    help="threads per primecount process during prewarm")
    ap.add_argument("--no-prewarm", action="store_true",
                    help="disable parallel checkpoint prewarming")
    a = ap.parse_args()
    _load()
    D, S = parse_exact_int(a.D), parse_exact_int(a.S)

    if a.mode == "bfile":
        write_bfile()
        return 0

    if a.mode == "validate":
        terms = [k for k in range(a.lo, a.hi + 1) if k in KNOWN]
        if not terms:
            print(f"no published A190502 values in range n={a.lo}..{a.hi} "
                  f"(b-file covers 0..{LAST_PUBLISHED})")
            return 1
        print(f"Checking A190502 a({terms[0]}) .. a({terms[-1]}) against the "
              f"published b-file ({len(terms)} terms).")
        print("These are real published OEIS values -- a mismatch means the")
        print("pipeline is wrong, not that a new term has been found.\n")
        ok = True
        results = []
        # Small terms: the walk width dominates and a huge D is pure waste.
        vD = min(D, 10**7)
        for i, k in enumerate(terms, 1):
            print(f"[[ term {i} of {len(terms)}: a({k}) ]]")
            m, info = certified_count(k, D=vD, S=min(S, vD),
                                      growth=a.growth, threads=a.threads,
                                      primecount=a.primecount)
            exp = KNOWN[k]
            good = (m == exp)
            ok &= good
            results.append((k, m, good))
            print(f"  >> a({k}) {'MATCH' if good else f'MISMATCH (expected {exp:,})'}\n")
        print("=" * 70)
        for k, m, good in results:
            print(f"  {'OK ' if good else 'XX '} a({k}) = {m:,}")
        print("=" * 70)
        print("VALIDATION", "PASSED - all terms matched" if ok
              else "FAILED - see MISMATCH above")
        return 0 if ok else 1

    if a.n is None:
        ap.error(f"{a.mode} needs an exponent n, e.g. `{a.mode} {LAST_PUBLISHED + 1}`")

    if a.mode == "explore":
        k = a.n
        print(f"Exploring A190502 from a({k}) upward"
              + (f" through a({a.stop})" if a.stop else ", continuing indefinitely")
              + " (each term ~2x the work of the last). Stop with Ctrl+C.")
        print(f"Published data ends at a({LAST_PUBLISHED}); "
              f"a({max(k, LAST_PUBLISHED + 1)}) onward is new.\n")
        while a.stop is None or k <= a.stop:
            if not a.no_prewarm:
                prewarm(2 ** k, D, a.growth, a.primecount, a.workers, a.pc_threads)
            m, info = certified_count(k, D=D, S=S, growth=a.growth,
                                      threads=a.threads, primecount=a.primecount,
                                      verify_anchors=True)
            if k in KNOWN:
                status = "MATCHES published value" if m == KNOWN[k] else \
                         f"*** DISAGREES with published {KNOWN[k]:,} ***"
                print(f"\n*** a({k}) = {m:,}   ({status}) ***")
                if m != KNOWN[k]:
                    return 1
            else:
                print(f"\n*** a({k}) = {m:,}   (NEW, previously unpublished term) ***")
            finish_term(k, m, info, D, S, True, a.primecount)
            k += 1
        write_bfile()
        return 0

    # single run
    n = a.n
    if n in KNOWN:
        print(f"Note: a({n}) is already published as {KNOWN[n]:,}; "
              f"this run will re-derive and check it.")
    else:
        print(f"Computing a NEW, previously unpublished term: "
              f"A190502 a({n}), Q = 2^{n} = {2**n:,}")
    if not a.no_prewarm:
        prewarm(2 ** n, D, a.growth, a.primecount, a.workers, a.pc_threads)
    m, info = certified_count(n, D=D, S=S, growth=a.growth, threads=a.threads,
                              primecount=a.primecount, verify_anchors=True)
    print(f"\n*** A190502 a({n}) = pi_R(2^{n}) = {m:,} ***")
    if n in KNOWN:
        print("MATCH against published value" if m == KNOWN[n]
              else f"*** MISMATCH: published value is {KNOWN[n]:,} ***")
    finish_term(n, m, info, D, S, a.largest, a.primecount)
    return 0


if __name__ == "__main__":
    sys.exit(main())

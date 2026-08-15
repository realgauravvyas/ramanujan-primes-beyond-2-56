// sieve128.cpp -- 128-bit segmented sieve and Ramanujan-prime walk engine.
//
// Works for arguments far above 2^64 (primesieve's hard limit), by using
// unsigned __int128 for positions while base primes (<= sqrt(x), always small
// enough) still come from libprimesieve.
//
// Modes:
//   primes  A B                 list primes in [A,B]                (testing)
//   count   A B                 count primes in [A,B]               (testing)
//   walk    Q L S               walk f(y)=pi(y)-pi(y/2) on (Q, Q+L]
//   back    Q L S DELTA         largest x in [Q-L,Q] with f(x)=f(Q)+DELTA
//
// Build:  g++ -O3 -march=native -fopenmp -std=c++17 sieve128.cpp \
//              -o sieve128 -lprimesieve
#include <primesieve.h>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif

typedef unsigned __int128 u128;
typedef uint64_t u64;

// ----------------------------------------------------------- 128-bit I/O
static u128 parse_u128(const char* s) {
    // accepts "123456", "1e20", "1.5e20", and underscores/commas
    std::string t;
    for (const char* p = s; *p; ++p)
        if (*p != '_' && *p != ',') t += *p;
    size_t e = t.find_first_of("eE");
    if (e == std::string::npos) {
        u128 v = 0;
        for (char c : t) if (c >= '0' && c <= '9') v = v * 10 + (u128)(c - '0');
        return v;
    }
    std::string mant = t.substr(0, e);
    int exp = atoi(t.c_str() + e + 1);
    size_t dot = mant.find('.');
    if (dot != std::string::npos) {
        int frac = (int)(mant.size() - dot - 1);
        mant.erase(dot, 1);
        exp -= frac;
    }
    u128 v = 0;
    for (char c : mant) if (c >= '0' && c <= '9') v = v * 10 + (u128)(c - '0');
    for (int i = 0; i < exp; i++) v *= 10;
    return v;
}
static std::string str_u128(u128 v) {
    if (v == 0) return "0";
    char b[64]; int i = 63; b[i--] = 0;
    while (v) { b[i--] = char('0' + (int)(v % 10)); v /= 10; }
    return std::string(b + i + 1);
}
static u64 isqrt_u128(u128 n) {
    if (n == 0) return 0;
    u64 x = (u64)sqrtl((long double)n);
    while (x > 0 && (u128)x * x > n) x--;
    while ((u128)(x + 1) * (x + 1) <= n) x++;
    return x;
}

// -------------------------------------------- odd-only segmented sieve
// Marks composites in [A,B]; returns primes as offsets from `base` (must be <= A).
static void sieve_range(u128 A, u128 B, u128 base, int nthreads,
                        std::vector<u64>& out) {
    out.clear();
    if (B < 2 || B < A) return;
    if (A < 2) A = 2;
    bool has2 = (A <= 2 && 2 <= B);
    u128 lo = (A & 1) ? A : A + 1;           // first odd >= A
    if (lo > B) { if (has2) out.push_back((u64)(2 - base)); return; }
    u64 n = (u64)((B - lo) / 2) + 1;         // count of odds in [lo,B]
    u64 words = (n + 63) / 64;
    std::vector<u64> bits(words, 0);
    u64 sq = isqrt_u128(B);
    const u64 CH = 1ull << 26;               // base-prime chunk

    #pragma omp parallel num_threads(nthreads)
    {
        std::vector<u64> loc(words, 0);
        #pragma omp for schedule(dynamic)
        for (long long c = 0; c <= (long long)(sq / CH); c++) {
            u64 clo = (u64)c * CH, chi = std::min(clo + CH - 1, sq);
            if (chi < 3) continue;
            if (clo < 3) clo = 3;
            size_t cnt = 0;
            u64* bp = (u64*)primesieve_generate_primes(clo, chi, &cnt, UINT64_PRIMES);
            for (size_t k = 0; k < cnt; k++) {
                u64 p = bp[k];
                u128 m0 = ((A + p - 1) / p) * p;      // first multiple >= A
                u128 pp = (u128)p * p;
                if (m0 < pp) m0 = pp;
                if ((m0 & 1) == 0) m0 += p;           // keep it odd
                if (m0 > B) continue;
                u64 idx = (u64)((m0 - lo) / 2);
                for (; idx < n; idx += p) loc[idx >> 6] |= (1ull << (idx & 63));
            }
            primesieve_free(bp);
        }
        #pragma omp critical
        for (u64 w = 0; w < words; w++) bits[w] |= loc[w];
    }
    if (lo == 1) bits[0] |= 1ull;            // 1 is not prime
    if (has2) out.push_back((u64)(2 - base));
    out.reserve(out.size() + n / 20);
    for (u64 i = 0; i < n; i++)
        if (!((bits[i >> 6] >> (i & 63)) & 1ull))
            out.push_back((u64)((lo + (u128)2 * i) - base));
}

// ------------------------------------------------------------ the walk
// Events on (a,b]:  +1 at each prime p;  -1 at each y=2q, q prime.
// Both event sets are produced relative to `base` so they merge as u64.
static void seg_events(u128 a, u128 b, int th,
                       std::vector<u64>& pos, std::vector<int8_t>& dlt) {
    std::vector<u64> P, Qh;
    sieve_range(a + 1, b, a, th, P);                       // +1 at p
    u128 a2 = a / 2, b2 = b / 2;
    sieve_range(a2 + 1, b2, a2, th, Qh);                   // -1 at 2q
    pos.clear(); dlt.clear();
    pos.reserve(P.size() + Qh.size()); dlt.reserve(P.size() + Qh.size());
    size_t i = 0, j = 0;
    while (i < P.size() || j < Qh.size()) {
        // positions relative to a:  p-a  and  2q-a = 2*(q-a2) + (2*a2-a)
        u64 vp = (i < P.size()) ? P[i] : ~0ull;
        u64 vq = (j < Qh.size()) ? (u64)(2 * (u128)Qh[j] + (2 * a2 - a)) : ~0ull;
        if (vp <= vq) { pos.push_back(vp); dlt.push_back(+1); i++; }
        else          { pos.push_back(vq); dlt.push_back(-1); j++; }
    }
}

int main(int argc, char** argv) {
    if (argc < 3) { fprintf(stderr, "usage: sieve128 MODE args...\n"); return 2; }
    std::string mode = argv[1];
    int th = 1;
#ifdef _OPENMP
    th = omp_get_max_threads();
#endif
    if (mode == "primes" || mode == "count") {
        u128 A = parse_u128(argv[2]), B = parse_u128(argv[3]);
        std::vector<u64> v; sieve_range(A, B, A, th, v);
        if (mode == "count") printf("%zu\n", v.size());
        else for (u64 o : v) printf("%s\n", str_u128(A + o).c_str());
        return 0;
    }
    if (mode == "walk") {
        u128 Q = parse_u128(argv[2]), L = parse_u128(argv[3]), S = parse_u128(argv[4]);
        long long f = 0, mn = 0; u128 mnpos = 0, np = 0, nm = 0;
        for (u128 off = 0; off < L; off += S) {
            u128 a = Q + off, b = std::min(Q + L, a + S);
            std::vector<u64> pos; std::vector<int8_t> dl;
            seg_events(a, b, th, pos, dl);
            for (size_t k = 0; k < pos.size(); k++) {
                f += dl[k];
                if (dl[k] > 0) np++; else nm++;
                if (f < mn) { mn = f; mnpos = a + pos[k] - Q; }
            }
            fprintf(stderr, "  seg %s .. %s  f=%lld min=%lld\n",
                    str_u128(a).c_str(), str_u128(b).c_str(), f, mn);
        }
        printf("min_rel %lld\nmin_off %s\nend_rel %lld\nn_plus %s\nn_minus %s\n",
               mn, str_u128(mnpos).c_str(), f, str_u128(np).c_str(),
               str_u128(nm).c_str());
        return 0;
    }
    if (mode == "back") {
        // largest x in [Q-L, Q] with f(x) = f(Q) + DELTA, walking downward
        u128 Q = parse_u128(argv[2]), L = parse_u128(argv[3]), S = parse_u128(argv[4]);
        long long target = atoll(argv[5]);
        long long v = 0;                       // f(Q) - f(Q) = 0
        u128 xhi = Q;
        u128 done = 0;
        while (done < L) {
            u128 b = Q - done, a = (b > S) ? b - S : 0;
            if (Q - a > L) a = Q - L;
            std::vector<u64> pos; std::vector<int8_t> dl;
            seg_events(a, b, th, pos, dl);
            for (size_t k = pos.size(); k-- > 0;) {
                u128 y = a + pos[k];
                if (v == target) { printf("x %s\n", str_u128(xhi).c_str()); return 0; }
                v -= dl[k]; xhi = y - 1;
            }
            done = Q - a;
        }
        if (v == target) { printf("x %s\n", str_u128(xhi).c_str()); return 0; }
        printf("NOTFOUND\n");
        return 1;
    }
    fprintf(stderr, "unknown mode %s\n", mode.c_str());
    return 2;
}

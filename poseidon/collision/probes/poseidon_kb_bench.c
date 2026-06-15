/*
 * Poseidon-16 KoalaBear permutation: optimized scalar C prototype + benchmark.
 *
 * Parameters (Poseidon Initiative 2026, q=4 partial-collision prize):
 *   Field:  KoalaBear p = 2^31 - 2^24 + 1 = 2130706433
 *   S-box:  x^3 (alpha = 3)
 *   Width:  t = 16
 *   Rounds: R_F = 8 (4 initial full + 4 terminal full), R_P = 20 partial
 *           (partial round applies S-box only to lane 0)
 *   MDS:    circulant, first row [1,1,51,1,11,17,2,1,101,63,15,2,67,22,13,3]
 *   RC:     Grain LFSR constants (448 = 28*16) matching the Python reference
 *
 * Arithmetic: Montgomery reduction, R = 2^32.
 *   p_inv'  = -p^{-1} mod 2^32 = 2130706431
 *   R2      = 2^64 mod p       = 402124772
 *
 * Build:
 *   clang -O3 -march=native -funroll-loops -o poseidon_kb_bench poseidon_kb_bench.c
 *
 * Verifies against the Python reference KAT, then benchmarks perms/sec/core.
 */

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>

#define P        2130706433u          /* KoalaBear prime */
#define P_INV    2130706431u          /* -p^{-1} mod 2^32 */
#define R2       402124772u           /* 2^64 mod p */
#define T        16
#define R_F      8
#define R_P      20
#define HALF_F   (R_F/2)
#define NROUNDS  (R_F + R_P)          /* 28 */

/* ---- Montgomery multiplication: returns a*b*R^{-1} mod p, R = 2^32 ---- */
static inline uint32_t mont_mul(uint32_t a, uint32_t b) {
    uint64_t x  = (uint64_t)a * b;          /* < p^2 < 2^62 */
    uint32_t m  = (uint32_t)x * P_INV;      /* low 32 bits */
    uint64_t t  = x + (uint64_t)m * P;      /* divisible by 2^32 */
    uint32_t r  = (uint32_t)(t >> 32);
    if (r >= P) r -= P;
    return r;
}

/* to/from Montgomery domain */
static inline uint32_t to_mont(uint32_t a)   { return mont_mul(a, R2); }
static inline uint32_t from_mont(uint32_t a) { return mont_mul(a, 1u); }

static inline uint32_t add_mod(uint32_t a, uint32_t b) {
    uint32_t s = a + b;
    if (s >= P) s -= P;     /* a,b < p so s < 2p < 2^32 */
    return s;
}

/* x^3 in Montgomery domain */
static inline uint32_t sbox(uint32_t x) {
    uint32_t x2 = mont_mul(x, x);
    return mont_mul(x2, x);
}

/* Circulant MDS first row, in normal (non-Montgomery) domain. */
static const uint32_t MDS_ROW[T] = {1,1,51,1,11,17,2,1,101,63,15,2,67,22,13,3};
/* In Montgomery domain (filled at init). */
static uint32_t MDS_ROW_M[T];

/* Round constants (normal domain) from Python Grain-LFSR reference. */
#include "rc_kb.h"
/* Round constants in Montgomery domain (filled at init). */
static uint32_t RC_M[NROUNDS][T];

/* Apply circulant MDS:  out[i] = sum_j MDS_ROW[(j - i) mod T] * state[j].
 * Equivalent to a cyclic convolution. We unroll the inner loop fully.
 * All values in Montgomery domain. */
static inline void apply_mds(uint32_t s[T]) {
    uint32_t out[T];
    for (int i = 0; i < T; i++) {
        uint32_t acc = 0;
        for (int j = 0; j < T; j++) {
            int idx = (j - i) & (T - 1);   /* T is a power of two */
            acc = add_mod(acc, mont_mul(MDS_ROW_M[idx], s[j]));
        }
        out[i] = acc;
    }
    memcpy(s, out, sizeof(out));
}

static inline void full_round(uint32_t s[T], const uint32_t rc[T]) {
    for (int i = 0; i < T; i++) s[i] = sbox(add_mod(s[i], rc[i]));
    apply_mds(s);
}

static inline void partial_round(uint32_t s[T], const uint32_t rc[T]) {
    for (int i = 0; i < T; i++) s[i] = add_mod(s[i], rc[i]);
    s[0] = sbox(s[0]);
    apply_mds(s);
}

/* Full permutation; input/output in Montgomery domain. */
static inline void permutation(uint32_t s[T]) {
    int r = 0;
    for (int k = 0; k < HALF_F; k++) full_round(s, RC_M[r++]);
    for (int k = 0; k < R_P;   k++) partial_round(s, RC_M[r++]);
    for (int k = 0; k < HALF_F; k++) full_round(s, RC_M[r++]);
}

static void init_montgomery(void) {
    for (int i = 0; i < T; i++) MDS_ROW_M[i] = to_mont(MDS_ROW[i] % P);
    for (int r = 0; r < NROUNDS; r++)
        for (int i = 0; i < T; i++) RC_M[r][i] = to_mont(RC[r*T + i] % P);
}

/* Expected outputs from the Python reference (circulant MDS + Grain LFSR). */
static const uint32_t KAT_ZERO[T] = {
    2067320972u, 506924172u, 1394794030u, 1814754695u, 169923386u, 1673494440u,
    1553037864u, 1678549726u, 238374927u, 1153010411u, 942253760u, 1034586261u,
    736793451u, 1744531091u, 1134330544u, 386040495u
};
static const uint32_t KAT_RANGE[T] = {
    1653397795u, 140286513u, 1428916292u, 1394470646u, 1355640352u, 2044647881u,
    774691762u, 743754183u, 970673709u, 1494356196u, 773929116u, 190527308u,
    19516788u, 1016605456u, 536974740u, 691034844u
};

static int verify(void) {
    uint32_t s[T];
    /* perm([0]*16) */
    for (int i = 0; i < T; i++) s[i] = to_mont(0);
    permutation(s);
    int ok = 1;
    for (int i = 0; i < T; i++) if (from_mont(s[i]) != KAT_ZERO[i]) ok = 0;
    /* perm(0..15) */
    for (int i = 0; i < T; i++) s[i] = to_mont((uint32_t)i);
    permutation(s);
    for (int i = 0; i < T; i++) if (from_mont(s[i]) != KAT_RANGE[i]) ok = 0;
    return ok;
}

int main(int argc, char** argv) {
    init_montgomery();

    if (!verify()) {
        fprintf(stderr, "KAT MISMATCH -- implementation is WRONG\n");
        /* print actual for debugging */
        uint32_t s[T];
        for (int i = 0; i < T; i++) s[i] = to_mont(0);
        permutation(s);
        fprintf(stderr, "perm([0]*16) =");
        for (int i = 0; i < T; i++) fprintf(stderr, " %u", from_mont(s[i]));
        fprintf(stderr, "\n");
        return 1;
    }
    printf("KAT OK: C output matches Python reference (circulant MDS + Grain LFSR).\n");

    /* Benchmark. Iterate the permutation in a tight dependency chain so the
     * compiler cannot hoist/eliminate it. We use the output to feed the next
     * input (a realistic VW iteration: state -> perm -> next state). */
    long long iters = (argc > 1) ? atoll(argv[1]) : 200000000LL;

    uint32_t s[T];
    for (int i = 0; i < T; i++) s[i] = to_mont((uint32_t)(i * 2654435761u % P));

    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);
    for (long long it = 0; it < iters; it++) {
        permutation(s);
        /* perturb lane to keep a data dependency / mimic VW step */
        s[0] = add_mod(s[0], to_mont(1));
    }
    clock_gettime(CLOCK_MONOTONIC, &t1);

    double secs = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) * 1e-9;
    double perms_per_sec = iters / secs;
    double ns_per_perm = secs * 1e9 / iters;

    /* sink to prevent dead-code elimination */
    volatile uint32_t sink = 0;
    for (int i = 0; i < T; i++) sink ^= from_mont(s[i]);

    printf("iters         : %lld\n", iters);
    printf("seconds       : %.4f\n", secs);
    printf("ns/perm       : %.3f\n", ns_per_perm);
    printf("perms/sec/core : %.6e\n", perms_per_sec);
    printf("log2(perms/sec/core) : %.3f\n", log2(perms_per_sec));
    printf("sink=%u\n", (unsigned)sink);
    return 0;
}

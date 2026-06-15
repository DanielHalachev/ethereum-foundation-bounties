/*
 * Poseidon-16 KoalaBear permutation: OPTIMIZED scalar C prototype + benchmark (v2).
 *
 * Same instance as v1 (KoalaBear p = 2^31-2^24+1, x^3, t=16, R_F=8, R_P=20,
 * circulant MDS [1,1,51,1,11,17,2,1,101,63,15,2,67,22,13,3], Grain-LFSR RC).
 *
 * Optimizations over v1:
 *  - State kept in NORMAL (non-Montgomery) domain.
 *  - MDS uses SMALL coefficients (max 101). out[i] = sum_j c_{(j-i) mod 16} * s[j].
 *    Each product c*s[j] <= 101*(p-1) < 2^38; the sum of 16 < 2^42 fits in u64,
 *    so we accumulate all 16 then do ONE modular reduction per output lane
 *    (lazy reduction) instead of 16 reductions.
 *  - x^3 S-box via Montgomery multiply (two mulmods).
 *
 * Build: clang -O3 -march=native -funroll-loops -o poseidon_kb_bench2 poseidon_kb_bench2.c -lm
 */

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>

#define P        2130706433u
#define P_INV    2130706431u   /* -p^{-1} mod 2^32 */
#define R2       402124772u    /* 2^64 mod p */
#define T        16
#define R_F      8
#define R_P      20
#define HALF_F   (R_F/2)
#define NROUNDS  (R_F + R_P)

/* reduce a 64-bit value < 2^42 into [0,p). Barrett-style via 128/64 not needed;
 * value < 2^42, p ~ 2^31, quotient < 2^11, so a couple of subtractions suffice
 * after one multiply-based estimate. We use a simple correct reduction. */
static inline uint32_t reduce64(uint64_t x) {
    return (uint32_t)(x % P);   /* one 64-bit div; correctness baseline */
}

/* Montgomery multiply for the cube (operands in Montgomery domain). */
static inline uint32_t mont_mul(uint32_t a, uint32_t b) {
    uint64_t x  = (uint64_t)a * b;
    uint32_t m  = (uint32_t)x * P_INV;
    uint64_t t  = x + (uint64_t)m * P;
    uint32_t r  = (uint32_t)(t >> 32);
    if (r >= P) r -= P;
    return r;
}
static inline uint32_t to_mont(uint32_t a)   { return mont_mul(a, R2); }
static inline uint32_t from_mont(uint32_t a) { return mont_mul(a, 1u); }

static inline uint32_t add_mod(uint32_t a, uint32_t b) {
    uint32_t s = a + b; if (s >= P) s -= P; return s;
}

/* x^3 in NORMAL domain: lift to Montgomery once is wasteful; instead compute
 * via Montgomery on the fly: x^3 = from_mont( mont(x)^3 * R^? ). Simpler and
 * fast: use a plain mulmod with reduce64 (x*x < 2^62, fits u64). */
static inline uint32_t mulmod(uint32_t a, uint32_t b) {
    return reduce64((uint64_t)a * b);
}
static inline uint32_t sbox(uint32_t x) {
    uint32_t x2 = mulmod(x, x);
    return mulmod(x2, x);
}

static const uint32_t MDS_ROW[T] = {1,1,51,1,11,17,2,1,101,63,15,2,67,22,13,3};

#include "rc_kb.h"  /* RC[448] in normal domain */
static uint32_t RC_R[NROUNDS][T];   /* reshaped */

/* MDS with small coefficients and lazy reduction. State in normal domain. */
static inline void apply_mds(uint32_t s[T]) {
    uint32_t out[T];
    for (int i = 0; i < T; i++) {
        uint64_t acc = 0;
        for (int j = 0; j < T; j++) {
            int idx = (j - i) & (T - 1);
            acc += (uint64_t)MDS_ROW[idx] * s[j];   /* <= 101*(p-1), 16 terms < 2^42 */
        }
        out[i] = reduce64(acc);
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
static inline void permutation(uint32_t s[T]) {
    int r = 0;
    for (int k = 0; k < HALF_F; k++) full_round(s, RC_R[r++]);
    for (int k = 0; k < R_P;   k++) partial_round(s, RC_R[r++]);
    for (int k = 0; k < HALF_F; k++) full_round(s, RC_R[r++]);
}

static const uint32_t KAT_ZERO[T] = {
    2067320972u,506924172u,1394794030u,1814754695u,169923386u,1673494440u,
    1553037864u,1678549726u,238374927u,1153010411u,942253760u,1034586261u,
    736793451u,1744531091u,1134330544u,386040495u };
static const uint32_t KAT_RANGE[T] = {
    1653397795u,140286513u,1428916292u,1394470646u,1355640352u,2044647881u,
    774691762u,743754183u,970673709u,1494356196u,773929116u,190527308u,
    19516788u,1016605456u,536974740u,691034844u };

static int verify(void) {
    uint32_t s[T];
    for (int i = 0; i < T; i++) s[i] = 0;
    permutation(s);
    int ok = 1;
    for (int i = 0; i < T; i++) if (s[i] != KAT_ZERO[i]) ok = 0;
    for (int i = 0; i < T; i++) s[i] = (uint32_t)i;
    permutation(s);
    for (int i = 0; i < T; i++) if (s[i] != KAT_RANGE[i]) ok = 0;
    return ok;
}

int main(int argc, char** argv) {
    for (int r = 0; r < NROUNDS; r++)
        for (int i = 0; i < T; i++) RC_R[r][i] = RC[r*T+i] % P;

    if (!verify()) {
        fprintf(stderr, "KAT MISMATCH\n");
        uint32_t s[T]; for (int i=0;i<T;i++) s[i]=0; permutation(s);
        for (int i=0;i<T;i++) fprintf(stderr," %u", s[i]); fprintf(stderr,"\n");
        return 1;
    }
    printf("KAT OK (v2 normal-domain + lazy-reduce MDS).\n");

    long long iters = (argc > 1) ? atoll(argv[1]) : 20000000LL;
    uint32_t s[T];
    for (int i = 0; i < T; i++) s[i] = (uint32_t)(i*2654435761u % P);

    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);
    for (long long it = 0; it < iters; it++) { permutation(s); s[0] = add_mod(s[0], 1); }
    clock_gettime(CLOCK_MONOTONIC, &t1);

    double secs = (t1.tv_sec-t0.tv_sec)+(t1.tv_nsec-t0.tv_nsec)*1e-9;
    double pps = iters/secs;
    volatile uint32_t sink=0; for (int i=0;i<T;i++) sink ^= s[i];
    printf("iters: %lld  seconds: %.4f  ns/perm: %.3f\n", iters, secs, secs*1e9/iters);
    printf("perms/sec/core: %.6e  log2: %.3f\n", pps, log2(pps));
    printf("sink=%u\n", (unsigned)sink);
    return 0;
}

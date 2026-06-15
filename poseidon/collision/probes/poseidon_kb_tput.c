/*
 * Poseidon-16 KoalaBear: THROUGHPUT benchmark (v4).
 *
 * Measures peak per-core perms/sec by running W INDEPENDENT permutation walks
 * interleaved, hiding per-instruction latency (this is how a real VW attack
 * saturates a core: many independent distinguished-point walks per core).
 *
 * Same instance & arithmetic as v3 (Barrett MDS, Montgomery available).
 * Build: clang -O3 -march=native -funroll-loops -o poseidon_kb_tput poseidon_kb_tput.c -lm
 */
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>

#define P 2130706433u
#define T 16
#define R_F 8
#define R_P 20
#define HALF_F (R_F/2)
#define NROUNDS (R_F+R_P)
#define W 8   /* independent walks interleaved per core */

static inline uint32_t reduce64(uint64_t x) {
    const uint64_t M = 8657571868ULL;          /* floor(2^64/p) */
    uint64_t q = (uint64_t)(((__uint128_t)x * M) >> 64);
    uint64_t r = x - q * P;
    if (r >= P) r -= P;
    return (uint32_t)r;
}
static inline uint32_t mulmod(uint32_t a, uint32_t b){ return reduce64((uint64_t)a*b); }
static inline uint32_t sbox(uint32_t x){ uint32_t x2=mulmod(x,x); return mulmod(x2,x); }
static inline uint32_t add_mod(uint32_t a,uint32_t b){ uint32_t s=a+b; if(s>=P)s-=P; return s; }

static const uint32_t MDS_ROW[T]={1,1,51,1,11,17,2,1,101,63,15,2,67,22,13,3};
#include "rc_kb.h"
static uint32_t RC_R[NROUNDS][T];

static inline void apply_mds(uint32_t s[T]){
    uint32_t out[T];
    for(int i=0;i<T;i++){ uint64_t a=0; for(int j=0;j<T;j++){ int idx=(j-i)&(T-1); a+=(uint64_t)MDS_ROW[idx]*s[j]; } out[i]=reduce64(a); }
    memcpy(s,out,sizeof(out));
}
static inline void full_round(uint32_t s[T],const uint32_t rc[T]){ for(int i=0;i<T;i++) s[i]=sbox(add_mod(s[i],rc[i])); apply_mds(s); }
static inline void partial_round(uint32_t s[T],const uint32_t rc[T]){ for(int i=0;i<T;i++) s[i]=add_mod(s[i],rc[i]); s[0]=sbox(s[0]); apply_mds(s); }
static inline void permutation(uint32_t s[T]){
    int r=0;
    for(int k=0;k<HALF_F;k++) full_round(s,RC_R[r++]);
    for(int k=0;k<R_P;k++)   partial_round(s,RC_R[r++]);
    for(int k=0;k<HALF_F;k++) full_round(s,RC_R[r++]);
}

static const uint32_t KAT_ZERO[T]={2067320972u,506924172u,1394794030u,1814754695u,169923386u,1673494440u,1553037864u,1678549726u,238374927u,1153010411u,942253760u,1034586261u,736793451u,1744531091u,1134330544u,386040495u};

int main(int argc,char**argv){
    for(int r=0;r<NROUNDS;r++) for(int i=0;i<T;i++) RC_R[r][i]=RC[r*T+i]%P;
    uint32_t v[T]; for(int i=0;i<T;i++) v[i]=0; permutation(v);
    int ok=1; for(int i=0;i<T;i++) if(v[i]!=KAT_ZERO[i]) ok=0;
    if(!ok){ fprintf(stderr,"KAT MISMATCH\n"); return 1; }
    printf("KAT OK (v4 throughput, W=%d interleaved walks).\n", W);

    long long iters = (argc>1)? atoll(argv[1]) : 4000000LL;  /* per walk */
    uint32_t st[W][T];
    for(int w=0;w<W;w++) for(int i=0;i<T;i++) st[w][i]=(uint32_t)((w*131+i)*2654435761u % P);

    struct timespec t0,t1; clock_gettime(CLOCK_MONOTONIC,&t0);
    for(long long it=0; it<iters; it++){
        for(int w=0; w<W; w++){ permutation(st[w]); st[w][0]=add_mod(st[w][0],1); }
    }
    clock_gettime(CLOCK_MONOTONIC,&t1);

    double secs=(t1.tv_sec-t0.tv_sec)+(t1.tv_nsec-t0.tv_nsec)*1e-9;
    long long total = iters*(long long)W;
    double pps = total/secs;
    volatile uint32_t sink=0; for(int w=0;w<W;w++) for(int i=0;i<T;i++) sink^=st[w][i];
    printf("total perms: %lld  seconds: %.4f  ns/perm: %.3f\n", total, secs, secs*1e9/total);
    printf("perms/sec/core: %.6e  log2: %.3f\n", pps, log2(pps));
    printf("sink=%u\n",(unsigned)sink);
    return 0;
}

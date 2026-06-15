"""
Probe bm_tmto_birthday -- PROBABILITY / TIME-MEMORY-TRADEOFF / COLLISION-SEARCH lens
on the Poseidon-1 CICO-2 bounty (RF=6, RP=10, p=2^31-2^24+1, t=16, x^3).

Map of interest
---------------
F : (x2,...,x15) in F_p^14  ->  (out0, out1) in F_p^2,   with x0=C1, x1=C2 fixed.
GOAL: a single preimage of (out0,out1) = (0,0).
"2-variable slice" F' : (x2,x3) in F_p^2 -> (out0,out1) in F_p^2, x4..x15 fixed.
Domain of F' ~ p^2 ~ 2^62.  Baseline brute = 2^62; best algebraic = ~2^55 (skip-1 resultant).

This probe answers, with cost formulas + tiny-instance numbers, whether ANY of:
  (1) Hellman / rainbow tables / distinguished points  (TMTO precompute->online)
  (2) Pollard rho / kangaroo                            (needs iterable self-map + group)
  (3) van Oorschot-Wiener parallel collision search     (golden collision / claw)
  (4) birthday meet-in-the-middle on V0={out0=0}, V1={out1=0}
beats 2^55 for a SINGLE preimage (no amortization).

We do NOT need the real prime to settle the STRUCTURE: random-function-ness,
absence of group structure, and the variety-intersection counting are field-size
robust, so we verify them on small primes and read off the asymptotic cost laws.
"""
import sys, os, random, math, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "reference"))
from poseidon.poseidon import Poseidon
from poseidon.mds_matrix import generate_circulant_mds_matrix

P_REAL = 2130706433            # 2^31 - 2^24 + 1
LOG2P  = math.log2(P_REAL)     # ~30.99
D = 3
ROW16 = [1, 1, 51, 1, 11, 17, 2, 1, 101, 63, 15, 2, 67, 22, 13, 3]
C1, C2 = 0xC09DE4, 0xEE6282


# ----------------------------------------------------------------------------
# small-instance Poseidon F' factory (scaled prime + t, real round structure)
# ----------------------------------------------------------------------------
def make_Fprime(p, t, rf, rp, c1, c2, fixed_tail, row=None, seed=0):
    """Return F'(a,b) -> (out0,out1) for the 2-var slice, real round structure."""
    rng = random.Random(seed)
    if row is None:
        # build an MDS-ish circulant first row for this t (reuse ROW16 prefix mod p, force coprime-ish)
        row = [(ROW16[i] % p) or 1 for i in range(t)]
    M = generate_circulant_mds_matrix(row[:t], p)
    pos = Poseidon(prime=p, alpha=D, t=t, r_f=rf, r_p=rp, mds=M)
    base_tail = list(fixed_tail)
    def Fp(a, b):
        st = [c1 % p, c2 % p, a % p, b % p] + [v % p for v in base_tail]
        out = pos.permutation_plus_linear(st)
        return (out[0], out[1])
    return Fp, p


# ============================================================================
# (A) Is F' a random function?  -> determines TMTO/collision cost laws.
#     For a SINGLE preimage with NO precomputation amortization, the relevant
#     quantity is: how many evaluations to hit the unique-ish target.
# ============================================================================
def test_random_function(p, t=4, rf=6, rp=4, trials_each=4):
    print(f"\n=== (A) F' random-function behaviour (p={p}, t={t}, RF={rf}, RP={rp}) ===")
    # A1: preimage-count distribution of F' over the WHOLE slice domain (enumerate p^2).
    #     For a random function F_p^2 -> F_p^2, #preimages of a fixed target ~ Poisson(1).
    rng = random.Random(7)
    tail = [rng.randrange(p) for _ in range(t - 4)]
    Fp, _ = make_Fprime(p, t, rf, rp, C1, C2, tail, seed=11)
    from collections import Counter
    img = Counter()
    for a in range(p):
        for b in range(p):
            img[Fp(a, b)] += 1
    # distribution of fibre sizes
    fib = Counter(img.values())               # size -> #targets having that many preimages
    total_targets = p * p
    hit_targets = sum(fib.values())
    miss = total_targets - hit_targets        # targets with 0 preimages
    mean_fib_over_image = sum(s * n for s, n in fib.items()) / hit_targets
    # Poisson(1) reference: P(0)=.3679 P(1)=.3679 P(2)=.1839 P(3)=.0613
    print(f"  domain p^2={p*p}, distinct images hit={hit_targets} ({hit_targets/total_targets:.3f} of codomain)")
    print(f"  fibre-size histogram (preimages per hit target): {dict(sorted(fib.items()))}")
    print(f"  P(0 preimages)={miss/total_targets:.3f}  [Poisson(1)=0.368]  -- random fn => ~37% targets unreachable")
    # specifically: is (0,0) reachable, and how many preimages?
    n00 = img.get((0, 0), 0)
    print(f"  #preimages of TARGET (0,0) in this slice = {n00}  (Poisson(1): expect ~0 or 1; many slices have 0)")
    return n00, miss / total_targets


# ============================================================================
# (B) Group-structure test: do rho / kangaroo / DLP-style walks apply?
#     They need F' (or some derived map) to be an ITERABLE self-map of a set
#     with a group operation the target relates to. Check: is F' even a self-map?
#     codomain F_p^2 == domain F_p^2 as sets, so iteration TYPE-CHECKS, but a
#     meaningful rho needs the cycle structure / a homomorphism to exploit.
#     We measure the rho cost to *collide* vs. to hit a *specific* point.
# ============================================================================
def test_group_structure(p, t=4, rf=6, rp=4):
    print(f"\n=== (B) group structure / rho applicability (p={p}, t={t}) ===")
    rng = random.Random(3)
    tail = [rng.randrange(p) for _ in range(t - 4)]
    Fp, _ = make_Fprime(p, t, rf, rp, C1, C2, tail, seed=21)
    # F' : F_p^2 -> F_p^2 is a self-map of the SET S=F_p^2 (size p^2).
    # Rho finds a COLLISION (two inputs same output) in ~sqrt(p^2)=p steps.
    # But we want a PREIMAGE of a SPECIFIC value (0,0), i.e. solve F'(z)=target.
    # Functional-graph iteration z->F'(z) reaches a fixed target only if target is
    # ON the path; for a random functional graph the rho-walk visits ~sqrt(p^2)=p
    # distinct nodes total, so probability the walk *passes through* a prescribed
    # node before cycling is ~ (#visited)/(p^2) = p / p^2 = 1/p. No leverage.
    # Demonstrate: run a rho walk, count distinct nodes till cycle (Floyd).
    def step(z):
        return Fp(z[0], z[1])
    z0 = (rng.randrange(p), rng.randrange(p))
    seen = {}
    z = z0; tail_len = None; rho_len = None
    for i in range(20 * p):          # cap
        if z in seen:
            tail_len = seen[z]; rho_len = i - tail_len; break
        seen[z] = i; z = step(z)
    visited = len(seen)
    print(f"  functional-graph rho walk: visited {visited} distinct nodes before cycle "
          f"(sqrt(p^2)=p={p}); tail={tail_len}, cycle={rho_len}")
    print(f"  P(a prescribed target lies on one random walk) ~ visited/p^2 = {visited/(p*p):.4f}  (~1/p)")
    print("  => no group homomorphism: target (0,0) is just one of p^2 nodes; iteration gives no")
    print("     shortcut to a SPECIFIC preimage. Rho/kangaroo solve DLP in a GROUP; here there is none.")
    return visited


# ============================================================================
# (C) Birthday meet-in-the-middle on the two output coordinates / varieties.
#     V0 = {x : out0(x)=0},  V1 = {x : out1(x)=0}.  Want a point in V0 ∩ V1.
#     Naive birthday idea: sample points of V0 and points of V1, hope two
#     samples coincide. Count how big V0∩V1 is and how big each Vi is.
# ============================================================================
def test_birthday_mitm(p, t=4, rf=6, rp=4):
    print(f"\n=== (C) birthday-MITM on varieties V0={{out0=0}}, V1={{out1=0}} (p={p}, t={t}) ===")
    rng = random.Random(5)
    tail = [rng.randrange(p) for _ in range(t - 4)]
    Fp, _ = make_Fprime(p, t, rf, rp, C1, C2, tail, seed=31)
    V0 = set(); V1 = set(); both = set()
    for a in range(p):
        for b in range(p):
            o0, o1 = Fp(a, b)
            if o0 == 0: V0.add((a, b))
            if o1 == 0: V1.add((a, b))
            if o0 == 0 and o1 == 0: both.add((a, b))
    print(f"  |slice domain|=p^2={p*p}")
    print(f"  |V0 (out0=0)|={len(V0)}  (expect ~ p^2/p = p = {p})   <- a curve, ~p points")
    print(f"  |V1 (out1=0)|={len(V1)}  (expect ~ p = {p})")
    print(f"  |V0 ∩ V1 (target)|={len(both)}  (expect ~ p^2/p^2 = 1)")
    # Key birthday accounting on the SLICE (2 free vars):
    # V0,V1 are subsets of a size-p^2 universe, each of size ~p. Two INDEPENDENT
    # random samples (one from V0, one from V1) collide w.p. ~ |V0∩V1| / (|V0|*|V1|)
    # if we just compare one-vs-one; the right birthday cost to find an element of
    # the (size m=|V0∩V1|) intersection by sampling both lists and hashing is:
    #   need to *enumerate* one full variety (cost ~p) and membership-test against the
    #   other.  That is the resultant/elimination cost, NOT a birthday speedup.
    # Demonstrate the FALSE hope: random one-from-each collision probability:
    if len(V0) and len(V1):
        p_one_each = len(both) / (len(V0) * len(V1))
        print(f"  P(one random V0-sample == one random V1-sample) = |∩|/(|V0||V1|) "
              f"= {p_one_each:.3e}  ~ 1/p^2  (NO birthday gain: V0,V1 nearly disjoint)")
    return len(V0), len(V1), len(both)


# ============================================================================
# (D) The CORRECT scaling: birthday on the FULL 14-var map vs the 2-var slice,
#     and why MITM across the partial-round middle has no low-degree meeting set.
#     We also state the exact cost formulas for the real instance.
# ============================================================================
def cost_formulas():
    print("\n=== (D) EXACT cost formulas at the real instance (p=2^31-2^24+1) ===")
    print(f"  log2 p = {LOG2P:.3f}")
    # brute / birthday baselines
    print("\n  -- Baselines --")
    print(f"   brute over 2-var slice            : p^2  ops  = 2^{2*LOG2P:.1f}")
    print(f"   best algebraic (skip-1 resultant) : D_I=3^20 -> ~2^55 time")
    # Hellman TMTO
    N = 2 ** (2 * LOG2P)
    print("\n  -- (1) Hellman/rainbow/DP TMTO (random-fn law: N = p^2 = 2^62) --")
    print("   coverage law: T * M^2 = N^2  (Hellman);  rainbow: T*M = N (roughly), with")
    print("   matrix-stopping-rule constraint m*t^2 <= N per table, ~t tables for full coverage.")
    print(f"   PRECOMPUTE cost = N = 2^{2*LOG2P:.1f}  (must evaluate ~the whole domain once).")
    print("   => For a SINGLE target with NO reuse, precompute (2^62) DOMINATES the 2^55 baseline.")
    print("      TMTO only wins when ONE precomputation is AMORTIZED over MANY inversions.")
    print("      CICO-2 bounty = ONE target (0,0) -> amortization factor = 1 -> TMTO is STRICTLY WORSE.")
    # rho / kangaroo
    print("\n  -- (2) Pollard rho / kangaroo --")
    print("   Require a GROUP + endomorphism (DLP). F' is F_p^2->F_p^2 with NO group law linking")
    print("   the target to iteration. Functional-graph rho finds a *collision* in ~sqrt(N)=2^31,")
    print("   but a collision is NOT a preimage of (0,0): P(walk hits prescribed node)~sqrt(N)/N=1/sqrt(N).")
    print("   To force a SPECIFIC preimage you'd need ~N/sqrt(N)=sqrt(N) restarts * sqrt(N) walk = N. No gain.")
    # VW parallel collision search
    print("\n  -- (3) van Oorschot-Wiener parallel collision (golden collision / claw) --")
    print("   VW finds a GOLDEN collision among W useful collisions in ~sqrt(N^3/W)/ (mem W) ... but")
    print("   that solves f(x)=g(y) CLAW problems where a *random* collision is the goal. Our goal is")
    print("   a preimage of a FIXED point (0,0), i.e. f(x)=const, which is a SEARCH not a collision.")
    print("   Recast as claw: define g0:a->out0(a,*) and g1:b->...; but out0,out1 share ALL 14 vars,")
    print("   they do not split into two independent halves whose outputs can be matched. No claw split.")
    # birthday-MITM
    print("\n  -- (4) birthday-MITM on V0={out0=0}, V1={out1=0} --")
    print("   On the 2-var slice: |V0|~p, |V1|~p, |V0∩V1|~1 (verified small-p above).")
    print("   To MITM you must ENUMERATE V0 (list ~p points) then test membership in V1.")
    print("   Enumerating V0 = solving the univariate out0(a,b)=0 for b given each a = the SAME")
    print("   elimination problem; listing V0 alone costs ~ p * (root-find of a deg-3^? poly) >> 2^55.")
    print("   Independent sampling cannot collide: P(rand V0 pt == rand V1 pt)=|∩|/(|V0||V1|)~1/p^2.")
    print("   The varieties are high-degree & 'in general position' (generic), so there is NO")
    print("   low-cost description of either to sample uniformly from.")
    # MITM across the middle
    print("\n  -- MITM across the partial-round middle (forward/backward) --")
    print("   Forward from input spans a coset reachable cheaply only up to the first S-box wall;")
    print("   backward from (0,0) output likewise. The meeting state is the FULL 16-word state in F_p,")
    print("   matching set size = p^16 = 2^496; you can only fix the 2 CICO words, leaving p^14 freedom")
    print("   on EACH side -> meet-set still p^14, birthday cost ~p^7=2^217. Catastrophically worse.")
    print("   No degree-collapse on either half (deg=3^(RF+RP), verified prior probes) to shrink it.")


def main():
    t0 = time.time()
    # Tiny exact-enumeration instances (p must be small: we enumerate p^2 fibres).
    # Use t=4 (minimal CICO-2 slice: words 0,1 fixed-ish, 2,3 vary), real round counts scaled.
    for p in (31, 61):
        test_random_function(p, t=4, rf=6, rp=4)
        test_group_structure(p, t=4, rf=6, rp=4)
        test_birthday_mitm(p, t=4, rf=6, rp=4)
    cost_formulas()
    print(f"\n[done in {time.time()-t0:.1f}s]")


if __name__ == "__main__":
    main()

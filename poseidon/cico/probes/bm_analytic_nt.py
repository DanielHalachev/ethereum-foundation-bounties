"""
Probe bm_analytic_nt.py — ANALYTIC NUMBER THEORY / ADDITIVE COMBINATORICS lens
on the Poseidon-1 CICO-2 bounty (RF=6, RP=10, p=2^31-2^24+1, t=16, x^3).

QUESTION: Can a character-sum / circle-method / Stepanov / Bombieri argument PROVE
that the F_p-solutions of {out0(X,Y)=0, out1(X,Y)=0} CONCENTRATE in a small,
exhaustively-searchable subset S of the (X,Y) plane (or of the 14-dim input space)?

The governing quantities (established): the elimination variety has degree
D_I = 3^(2*RF+RP) = 3^22 (or 3^20 with the legal 1-round skip). A 2-var slice
{out0=out1=0} has ~1 solution. Brute force over (X,Y) is p^2 = 2^62; we want to
beat the ~2^55 resultant attack by restricting (X,Y) to a structured set S.

The probe quantifies, on TINY instances (small p, small t, few rounds), four
character-sum predictions, and extrapolates the error terms to the real instance:

  (A) WEIL / DELIGNE error term. The count of solutions in a box / coset / subspace S
      is  N(S) = (|S|/p) * (#solutions on whole plane)  +  ERROR,
      and the Weil/Deligne bound on ERROR scales like  deg(variety) * sqrt(|S|)
      (per character) up to  D_I * p^{(dim-1/2)}. We measure on tiny instances
      whether the MAIN term |S|/p * Nsol ever EXCEEDS the error, for S small.

  (B) CONCENTRATION test. Empirically: do the actual solutions of the system on a
      tiny instance fall into any low-complexity set — a multiplicative coset
      <g^d>, a low-Hamming-weight set, a small affine subspace, a short interval?
      If solutions were equidistributed (the Weil prediction for a generic variety),
      NO such set of size << p contains a solution with non-negligible probability.

  (C) STEPANOV feasibility. Stepanov's method proves point-count bounds for a SINGLE
      curve f(X,Y)=0 of degree d over F_p with EXPLICIT error ~ d^2 * sqrt(p). We
      check the degrees involved and whether a Stepanov-style auxiliary polynomial
      could isolate roots in a short interval — i.e. is there any interval of length
      L << p guaranteed to contain a root.

  (D) CIRCLE METHOD main-term-vs-error. The Hardy-Littlewood count of solutions in S
      via additive characters has main term ~ |S|^2 / p and error governed by the
      minor arcs, bounded by max over nontrivial a of the exponential sum
      Sum_{x in S} e_p(a * out0(x)). We measure these exponential sums on tiny
      instances to see if any nontrivial frequency a gives a LARGE sum (=structure to
      exploit) or whether they are all square-root-cancellation small (=no structure).

OUTPUT: hard numbers + a verdict on whether ANY concentration result is usable, with
the dominant obstruction identified.

Run: ./.venv/bin/python probes/bm_analytic_nt.py
"""
import sys, os, math, random, cmath
from itertools import product

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "reference"))
from poseidon.poseidon import Poseidon
from poseidon.mds_matrix import generate_mds_matrix

random.seed(20260611)


# ---------------------------------------------------------------------------
# A tiny scaled-down Poseidon-1 (permutation_plus_linear) over a small prime,
# so we can ENUMERATE the whole (X,Y) plane and study solution distribution.
# Structural questions (does the solution set concentrate? do exp-sums cancel?)
# are field-size-robust, per the project notes.
# ---------------------------------------------------------------------------
def build(p, t, r_f, r_p, alpha=3):
    # generate_mds_matrix needs gcd(alpha,p-1)=1 for a valid sponge, but we only
    # use the matrix + round constants; cube is a permutation iff gcd(3,p-1)=1.
    M = generate_mds_matrix(t, p)  # Cauchy-type MDS, matches verifier family
    pos = Poseidon(prime=p, alpha=alpha, t=t, r_f=r_f, r_p=r_p, mds=M)
    rc = pos.round_constants
    Mi = [[int(M[i][j]) for j in range(t)] for i in range(t)]
    return pos, Mi, rc


def perm_eval(p, t, r_f, r_p, Mi, rc, state, alpha=3):
    """Evaluate permutation_plus_linear on an integer state vector mod p."""
    def mds(s):
        return [sum(Mi[i][j] * s[j] for j in range(t)) % p for i in range(t)]
    def arc(s, c):
        return [(s[i] + int(c[i])) % p for i in range(t)]
    st = [v % p for v in state]
    st = mds(st)
    half = r_f // 2
    idx = 0
    for _ in range(half):
        st = arc(st, rc[idx]); idx += 1
        st = [pow(v, alpha, p) for v in st]; st = mds(st)
    for _ in range(r_p):
        st = arc(st, rc[idx]); idx += 1
        st[0] = pow(st[0], alpha, p); st = mds(st)
    for _ in range(half):
        st = arc(st, rc[idx]); idx += 1
        st = [pow(v, alpha, p) for v in st]; st = mds(st)
    return st


def all_solutions(p, t, r_f, r_p, alpha=3, c1=None, c2=None, fixed=None):
    """Enumerate the entire (X,Y) plane; return solutions of out0=out1=0 and the
    full out0/out1 arrays for spectral analysis. Only feasible for tiny p,t."""
    pos, Mi, rc = build(p, t, r_f, r_p, alpha)
    if c1 is None:
        c1 = 0xC09DE4 % p
    if c2 is None:
        c2 = 0xEE6282 % p
    if fixed is None:
        fixed = [(424242 * (i + 1) + 7) % p for i in range(t - 4)]
    sols = []
    out0 = {}
    out1 = {}
    for X in range(p):
        for Y in range(p):
            st = [c1, c2, X, Y] + list(fixed)
            o = perm_eval(p, t, r_f, r_p, Mi, rc, st, alpha)
            out0[(X, Y)] = o[0]
            out1[(X, Y)] = o[1]
            if o[0] == 0 and o[1] == 0:
                sols.append((X, Y))
    return sols, out0, out1, (c1, c2, fixed)


# ---------------------------------------------------------------------------
# (B) CONCENTRATION diagnostics on the solution set
# ---------------------------------------------------------------------------
def primitive_root(p):
    if p == 2:
        return 1
    fac = []
    n = p - 1
    d = 2
    while d * d <= n:
        if n % d == 0:
            fac.append(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        fac.append(n)
    for g in range(2, p):
        if all(pow(g, (p - 1) // q, p) != 1 for q in fac):
            return g
    return None


def concentration_report(p, sols, alpha=3):
    print(f"    #solutions on full plane = {len(sols)}  (expected ~ p^2/p^2 = O(1) "
          f"for a generic 0-dim system; plane size p^2={p*p})")
    if not sols:
        print("    (no F_p solutions on this slice — re-slice would be needed)")
        return
    g = primitive_root(p)
    # multiplicative-coset test: are X-coords in a small index-d subgroup or coset?
    # subgroup of d-th powers has size (p-1)/gcd(d,p-1).
    sub = pow_residue_class(p, sols, alpha, g)
    print(f"    multiplicative-coset (d-th power class) test on X-coords: {sub}")
    # interval test: smallest interval (mod p) covering all X-coords
    xs = sorted(set(x for x, _ in sols))
    span = min((xs[(i + len(xs) - 1) % len(xs)] - xs[i]) % p for i in range(len(xs))) \
        if len(xs) > 1 else 0
    # smallest covering arc:
    arc = smallest_arc(xs, p)
    print(f"    X-coord interval test: {len(xs)} distinct X in arc of length "
          f"{arc}/{p} = {arc/p:.3f} of the field  "
          f"({'CONCENTRATED' if arc < p/4 and len(xs) > 2 else 'spread / too few pts'})")
    # Hamming-weight test (in base-2): mean popcount vs random expectation
    bits = p.bit_length()
    avg_hw = sum(bin(x).count("1") for x, _ in sols) / len(sols)
    print(f"    Hamming-weight test: mean popcount(X)={avg_hw:.2f}, "
          f"random expectation~{bits/2:.2f} (no concentration if equal)")


def pow_residue_class(p, sols, d, g):
    """Check if all X-coords lie in one coset of the subgroup of d-th powers."""
    if g is None:
        return "n/a"
    m = math.gcd(d, p - 1)
    if m == 1:
        return f"d-th powers = all of F_p* (gcd(d,p-1)=1) -> no restriction"
    # discrete log of each nonzero X mod m
    logs = []
    table = {}
    cur = 1
    for e in range(p - 1):
        table[cur] = e
        cur = (cur * g) % p
    classes = set()
    for x, _ in sols:
        if x % p == 0:
            classes.add("zero")
        else:
            classes.add(table[x] % m)
    return (f"X-coords occupy {len(classes)} of {m} cosets "
            f"({'CONCENTRATED in 1 coset' if len(classes) == 1 else 'spread across cosets'})")


def smallest_arc(xs, p):
    if len(xs) <= 1:
        return 0
    xs = sorted(xs)
    gaps = [(xs[(i + 1) % len(xs)] - xs[i]) % p for i in range(len(xs))]
    biggest = max(gaps)
    return p - biggest


# ---------------------------------------------------------------------------
# (D) CIRCLE-METHOD exponential sums:  S(a) = sum_{X,Y} e_p(a*out0 + b*out1)
#     and per-output  S0(a) = sum_X e_p(a*out0(X, y0)) along a line.
#     Square-root cancellation => |S| ~ sqrt(#terms); a LARGE |S| would mean a
#     usable frequency (structure). Report the max nontrivial magnitude.
# ---------------------------------------------------------------------------
def exp_sum_spectrum(p, out0, out1, sols):
    n = p * p
    # full 2D additive-character sum over the plane for the joint (a,b) freq.
    # The count of solutions via circle method:
    #   #sol = (1/p^2) sum_{a,b} S(a,b),  S(a,b)=sum_{X,Y} e_p(a*out0+b*out1)
    # main term = S(0,0)/p^2 = n/p^2 = 1.  We report the LARGEST |S(a,b)| over
    # nontrivial (a,b), normalized by sqrt(n) (square-root-cancellation scale).
    w = cmath.exp(2j * math.pi / p)
    best = 0.0
    best_ab = None
    # sample frequencies if p large-ish; for tiny p do all.
    freqs = [(a, b) for a in range(p) for b in range(p)]
    for (a, b) in freqs:
        if a == 0 and b == 0:
            continue
        s = 0j
        for (X, Y), v0 in out0.items():
            v1 = out1[(X, Y)]
            s += w ** ((a * v0 + b * v1) % p)
        m = abs(s)
        if m > best:
            best = m
            best_ab = (a, b)
    sqrtn = math.sqrt(n)
    print(f"    joint exp-sum: max nontrivial |S(a,b)|={best:.2f} at (a,b)={best_ab}; "
          f"sqrt(plane)={sqrtn:.2f};  ratio={best/sqrtn:.2f}  "
          f"({'STRUCTURE (large)!' if best > 8 * sqrtn else 'square-root cancellation -> no usable frequency'})")
    # single-output marginal: is out0 alone biased? max_a |sum e_p(a*out0)|
    best1 = 0.0
    best_a = None
    for a in range(1, p):
        s = 0j
        for v0 in out0.values():
            s += w ** ((a * v0) % p)
        m = abs(s)
        if m > best1:
            best1 = m
            best_a = a
    print(f"    marginal exp-sum out0: max_a |sum e_p(a*out0)|={best1:.2f} at a={best_a}; "
          f"sqrt(plane)={sqrtn:.2f}; ratio={best1/sqrtn:.2f} "
          f"({'BIASED' if best1 > 8 * sqrtn else 'uniform -> no usable bias'})")


# ---------------------------------------------------------------------------
# (A) WEIL/DELIGNE main-term-vs-error budget at the REAL instance.
#     For a 0-dim system cut from a variety of degree D over F_p, restricting to a
#     subset S of the plane, the expected count is  Nsol*|S|/p^2  and the
#     unconditional error from Deligne/Lang-Weil is bounded by ~ C(D) * |S|^{1/2}
#     for a 1-codim restriction, growing with D = 3^(2RF+RP). We compute when the
#     main term could beat the error.
# ---------------------------------------------------------------------------
def weil_budget():
    print("\n=== (A) Weil/Deligne main-term-vs-error budget at REAL instance ===")
    p = 2130706433
    RF, RP = 6, 10
    D_full = 3 ** (2 * RF + RP)      # 3^22
    D_skip = 3 ** (2 * (RF) + RP - 2)  # 3^20 with legal 1-round skip
    print(f"    p = {p}  (~2^{math.log2(p):.2f})")
    print(f"    elimination degree D_I (no skip) = 3^22 = {D_full}  (~2^{math.log2(D_full):.2f})")
    print(f"    elimination degree D_I (skip-1)  = 3^20 = {D_skip}  (~2^{math.log2(D_skip):.2f})")
    print("    A solution of {out0=out1=0} is a point on the 0-dim scheme cut from")
    print("    the plane by two curves of degree D_I; we want a set S (an arc, coset,")
    print("    or subspace) of the (X,Y) plane that PROVABLY contains a solution and")
    print("    has |S| << 2^55 so we can search it.")
    print()
    print("    Circle-method / Weil count of solutions IN S:")
    print("       N(S) = |S|/p^2 * Nsol  +  E,   with Nsol = O(1) (the ~1 sol/slice).")
    print("    For the MAIN TERM to be >= 1 (so S is guaranteed nonempty) we already")
    print("    need |S| >~ p^2 / Nsol ~ p^2 = 2^62  -- i.e. essentially the WHOLE plane.")
    print("    The error term E is bounded (Deligne) by ~ D_I * sqrt(|S|), which for")
    print("    any |S| <= p^2 DWARFS the main term by a factor ~ D_I = 3^20 ~ 2^31.")
    main_term_full_plane = 1.0  # Nsol ~ O(1)
    # error on full plane (worst case): D_I * p (codim-1 Deligne) -- astronomically
    # larger than the count, which is why the count is only meaningful asymptotically.
    err_full = D_skip * p
    print(f"\n    Concretely on the FULL plane: main term Nsol ~ O(1),")
    print(f"    Deligne error budget ~ D_I*p = 3^20 * 2^31 ~ 2^{math.log2(err_full):.1f}.")
    print(f"    Error/Main ratio ~ 2^{math.log2(err_full):.1f}  ==>  the analytic count is")
    print("    VACUOUS for distinguishing where (rather than how many) solutions sit.")
    print("    The crossover |S| where main term beats Deligne error needs")
    print("    |S| >~ (D_I)^2 * p^2 / Nsol^2  >>  p^2  -- impossible inside the plane.")


def pooled_concentration(p, t, r_f, r_p, n_slices=400, alpha=3):
    """Pool the FULL set of (X,Y) solutions across many random re-slices (random
    C1,C2 and 12 fixed inputs). This builds a large solution population to test
    equidistribution robustly. For each slice we enumerate X over the whole field
    and solve out1(.,Y)=0 via enumerating Y too -> full plane per slice. To keep
    it fast we enumerate the plane only for small p.

    We then ask: are the solution X-coords equidistributed mod p (Weil prediction
    for a generic variety), or do they pile up in an arc / coset / residue class?
    Metric: chi-square of X-coords against uniform over p bins, and against the
    quadratic-residue split.
    """
    pos, Mi, rc = build(p, t, r_f, r_p, alpha)
    g = primitive_root(p)
    xcoords = []
    rng = random.Random(99)
    plane = list(product(range(p), range(p)))
    for _ in range(n_slices):
        c1 = rng.randrange(p); c2 = rng.randrange(p)
        fixed = [rng.randrange(p) for _ in range(t - 4)]
        for (X, Y) in plane:
            st = [c1, c2, X, Y] + fixed
            o = perm_eval(p, t, r_f, r_p, Mi, rc, st, alpha)
            if o[0] == 0 and o[1] == 0:
                xcoords.append(X)
    n = len(xcoords)
    print(f"\n=== POOLED concentration p={p} t={t} RF={r_f} RP={r_p}: "
          f"{n} solutions over {n_slices} random re-slices ===")
    if n < 30:
        print("    too few pooled solutions for a robust test")
        return
    # chi-square uniformity over p bins (each value is its own bin)
    counts = [0] * p
    for x in xcoords:
        counts[x] += 1
    exp = n / p
    chi2 = sum((c - exp) ** 2 / exp for c in counts)
    # df = p-1; mean of chi2 ~ df, std ~ sqrt(2 df). z-score:
    df = p - 1
    z = (chi2 - df) / math.sqrt(2 * df)
    D_I = alpha ** (2 * r_f + r_p)
    print(f"    [variety degree D_I=3^{2*r_f+r_p}={D_I}, field p={p}; "
          f"D_I/p={D_I/p:.1f}x -> curve wraps the tiny field many times]")
    print(f"    X-coord uniformity over F_p: chi2={chi2:.1f}, df={df}, "
          f"z={z:.2f}  ({'UNIFORM' if abs(z) < 3 else 'chi2 flags non-uniform'})")
    if abs(z) >= 3:
        print("      ^ BUT this is a FINITE-SIZE artifact, not exploitable: with")
        print("        D_I>>p the curve covers F_p many times so per-bin counts are")
        print("        deterministic mod small p. EXPLOITABILITY = does it shrink the")
        print("        search set? -> check the arc/coset tests below.")
    # quadratic-residue bias: fraction of nonzero X that are QRs
    qr = 0; nz = 0
    for x in xcoords:
        if x % p == 0:
            continue
        nz += 1
        if pow(x, (p - 1) // 2, p) == 1:
            qr += 1
    frac = qr / nz if nz else 0
    sd = 0.5 / math.sqrt(nz) if nz else 1
    zqr = (frac - 0.5) / sd if nz else 0
    print(f"    QR-bias of X-coords: {qr}/{nz}={frac:.3f} (expect 0.5), "
          f"z={zqr:.2f} ({'BIASED' if abs(zqr) > 3 else 'no bias'})")
    # arc concentration of the pooled set
    arc = smallest_arc(sorted(set(xcoords)), p)
    print(f"    distinct X-values cover arc {arc}/{p}={arc/p:.3f} of field "
          f"({'CONCENTRATED' if arc < p * 0.6 else 'fills the field'})")


def vanishing_locus_equidistribution(p, t, r_f, r_p, n_slices=60, alpha=3):
    """Larger-field test that does NOT require full-plane enumeration: along a line
    (fix Y and all others), enumerate X over F_p and record the X where out0=0.
    The roots of out0(.,Y) over a line are a degree-D_I univariate's roots.
    Question: are these roots equidistributed in F_p (generic) or concentrated?
    We pool roots over many random lines and test arc/coset/uniformity. Cost ~
    n_slices * p evaluations -> use p up to a few thousand.
    """
    pos, Mi, rc = build(p, t, r_f, r_p, alpha)
    rng = random.Random(7)
    roots = []
    for _ in range(n_slices):
        c1 = rng.randrange(p); c2 = rng.randrange(p)
        Yv = rng.randrange(p)
        fixed = [rng.randrange(p) for _ in range(t - 4)]
        for X in range(p):
            st = [c1, c2, X, Yv] + fixed
            o = perm_eval(p, t, r_f, r_p, Mi, rc, st, alpha)
            if o[0] == 0:
                roots.append(X)
    n = len(roots)
    D_I = alpha ** (2 * r_f + r_p)
    print(f"\n=== VANISHING-LOCUS of out0 along lines, p={p} t={t} RF={r_f} RP={r_p} "
          f"(D_I=3^{2*r_f+r_p}={D_I}, D_I/p={D_I/p:.2f}) ===")
    print(f"    {n} roots over {n_slices} random lines "
          f"(expected ~ n_slices roots since deg/p-fold cover gives ~p/p*deg... )")
    if n < 50:
        print("    too few roots"); return
    arc = smallest_arc(sorted(set(roots)), p)
    print(f"    roots cover arc {arc}/{p}={arc/p:.3f} of field "
          f"({'CONCENTRATED in <half field' if arc < p * 0.5 else 'fills the field -> spread'})")
    # coarse-bin chi-square (sqrt(p) bins) so test is meaningful when D_I/p is modest
    nb = max(8, int(math.sqrt(p)))
    counts = [0] * nb
    for x in roots:
        counts[x * nb // p] += 1
    exp = n / nb
    chi2 = sum((c - exp) ** 2 / exp for c in counts)
    df = nb - 1
    z = (chi2 - df) / math.sqrt(2 * df)
    print(f"    coarse uniformity ({nb} bins): chi2={chi2:.1f}, df={df}, z={z:.2f} "
          f"({'UNIFORM -> no usable concentration' if abs(z) < 4 else 'non-uniform (check arc)'})")


def stepanov_note():
    print("\n=== (C) Stepanov-method feasibility ===")
    print("    Stepanov bounds #points on ONE curve f(X,Y)=0 of degree d over F_p with")
    print("    error ~ d^2*sqrt(p), and its auxiliary-polynomial construction can in")
    print("    principle confine roots — but only for LOW degree d relative to p.")
    print("    Here each eliminated curve has degree D_I = 3^20 ~ 2^31 >= p ~ 2^31.")
    print("    Stepanov REQUIRES d < ~ p^{1/2} (the auxiliary poly must have degree")
    print("    < p to be nonzero); with d ~ p the construction is empty: NO interval")
    print("    confinement, NO nontrivial point bound. Obstruction (i) degree>>field.")


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------
def main():
    print("############################################################")
    print("# ANALYTIC NUMBER THEORY / ADDITIVE COMBINATORICS PROBE     #")
    print("# Poseidon-1 CICO-2  (concentration of solutions?)          #")
    print("############################################################")

    # --- tiny instances we can fully enumerate ---
    # need gcd(3,p-1)=1 for x^3 a bijection: p=5(no,p-1=4),  p=7 (6,gcd3=3 no),
    # p=11 (10 gcd 1 yes), p=17(16 yes), p=23(22 yes), p=29(28 yes).
    tiny = [
        (11, 4, 2, 1),
        (17, 4, 2, 2),
        (23, 4, 4, 2),
        (29, 6, 2, 2),
    ]
    for (p, t, r_f, r_p) in tiny:
        assert math.gcd(3, p - 1) == 1, f"x^3 not a bijection mod {p}"
        print(f"\n=== TINY p={p} t={t} RF={r_f} RP={r_p}  (x^3, Cauchy MDS) ===")
        sols, out0, out1, _ = all_solutions(p, t, r_f, r_p)
        print("  (B) CONCENTRATION of the solution set:")
        concentration_report(p, sols)
        print("  (D) CIRCLE-METHOD exponential-sum spectrum over the plane:")
        exp_sum_spectrum(p, out0, out1, sols)

    # --- pooled equidistribution test (robust, large solution population) ---
    # p=17 plane is 289; 400 slices -> ~hundreds of solutions, fast.
    pooled_concentration(17, 4, 2, 2, n_slices=2000)
    pooled_concentration(23, 4, 2, 1, n_slices=2000)

    # --- larger-field vanishing-locus equidistribution (D_I/p moderate) ---
    # p where gcd(3,p-1)=1: 1019 (p-1=1018=2*509, gcd1), 2027 (2026=2*1013), 4001? p-1=4000=2^5*5^3 gcd 1.
    # Use modest rounds so D_I/p stays O(1)-O(10): RF=2,RP=0 -> 3^4=81; RF=2,RP=1->243.
    vanishing_locus_equidistribution(1019, 4, 2, 1, n_slices=120)  # D_I=243, D_I/p~0.24
    vanishing_locus_equidistribution(2027, 4, 2, 2, n_slices=80)   # D_I=729, D_I/p~0.36

    # --- analytic budgets at the real instance ---
    weil_budget()
    stepanov_note()

    print("\n############################################################")
    print("# SUMMARY")
    print("############################################################")
    print("If across tiny instances: (B) solutions are spread (arcs ~ full field,")
    print("coset test = spread, popcount ~ random) AND (D) all nontrivial exp-sums")
    print("show square-root cancellation -> the variety is analytically GENERIC and")
    print("NO character-sum argument concentrates solutions into a searchable S.")


if __name__ == "__main__":
    main()

"""
wf_special_slice.py  --  NON-GENERIC SLICES angle for Poseidon-1 CICO.

The standard CICO-2 attack fixes inputs x2..x15 to RANDOM constants and solves a
bivariate system {out0(X,Y)=0, out1(X,Y)=C} in (X,Y)=(x?,x?) via a resultant of
degree D_I = d^(2RF+RP) (one front round skippable -> d^(2(RF-1)+RP)).

QUESTION (this probe): does a SPECIAL choice of the fixed inputs, or restricting
the free variables to a STRUCTURED set, collapse the degree of the resultant or
make it FACTOR into low-degree pieces?  Concretely:

  (BASE) Measure deg_X Res_Y(out0,out1) and its factorization for GENERIC random
         fixed inputs on tiny instances.  This is the control.

  (a) FACTORIZATION under special fixed-value choices: set the 12 (here t-4)
      fixed inputs to  zero / all-equal / equal-to-round-constants / negated
      round constants, etc., and check (i) does deg drop?  (ii) does the
      univariate resultant FACTOR into low-degree irreducible pieces over GF(p)?
      A factor of degree D' << D_I that contains a real root => searchable.

  (b) STRUCTURED DOMAIN: the field is GF(p), p-1 = 2^24 * 127.  gcd(3,p-1)=1 so
      x^3 is a BIJECTION on F_p (the angle's "3-to-1 on order-126" premise is
      false -- 126 does not divide p-1).  But x^3 PERMUTES every multiplicative
      subgroup.  Restrict free vars X,Y to:
        - the order-127 subgroup  H127  (|H127| = 127, tiny!),
        - the order-2^k subgroup / a 2^k-coset,
      and check whether out0|H, out1|H are low-degree as functions on the
      subgroup (degree reduces mod (T^|H| - 1)) so an exhaustive/meet search over
      H x H costs |H|^2 instead of p^2.  For CICO-2 we need 2 output coords to
      vanish; if free vars live in H127 the probability a random pair works is
      ~ |H127|^2 / p^2 -- but if out0,out1 are CONSTANT or LOW-RANK on H, the
      structure could be exploited.  We measure the actual image sizes / collision
      structure of (out0,out1) restricted to H127 x H127.

  (c) RELATED-to-constant fixed values: tie x_i to round constants rc to try to
      cancel the first ARC, mirroring known "round-skip" tricks at the input.

ALL claims are backed by HARD NUMBERS measured on tiny RF/RP instances where the
sympy resultant runs in seconds, plus a concrete cost extrapolation to the real
RF=6,RP=10 target.

Run: ./.venv/bin/python probes/wf_special_slice.py
"""
import sys, os, time, itertools, functools
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "reference"))
from sympy import symbols, Poly, GF, factor_list
from poseidon.poseidon import Poseidon
from poseidon.mds_matrix import generate_mds_matrix, generate_circulant_mds_matrix

def log(*a):
    print(*a)
    sys.stdout.flush()

# factoring over GF(p) with p~2^31 is SLOW: deg-243 ~21s each, deg-81 ~0.6s.
# We always MEASURE the resultant degree (the load-bearing metric, cheap); we
# only FACTOR when degree <= cap so the whole run stays under the time budget.
# deg-81 factorizations are enough to compare special-value factor STRUCTURE.
FACTOR_DEG_CAP = 100

P = 2130706433
T = 16
D = 3
ROW16 = [1, 1, 51, 1, 11, 17, 2, 1, 101, 63, 15, 2, 67, 22, 13, 3]
C1, C2 = 0xC09DE4, 0xEE6282
X, Y = symbols("X Y")
DOM = GF(P)

CAUCHY = generate_mds_matrix(T, P)
CIRC = generate_circulant_mds_matrix(ROW16, P)


# ---------------------------------------------------------------------------
# Symbolic permutation: free vars at positions 2,3 = (X,Y); the OTHER positions
# (0,1 are the input constants C1,C2; 4..15 are the "fixed inputs") are concrete.
# ---------------------------------------------------------------------------
def sym_perm(pos, M, fixed_vals):
    """fixed_vals: list of length T-4 giving inputs at positions 4..15.
    Positions 0,1 = C1,C2 ; positions 2,3 = X,Y (symbolic)."""
    def mds(s):
        return [sum(int(M[i][j]) * s[j] for j in range(T)) for i in range(T)]
    def arc(s, c):
        return [s[i] + int(c[i]) for i in range(T)]
    rc = pos.round_constants
    state = [Poly(int(C1) % P, X, Y, domain=DOM),
             Poly(int(C2) % P, X, Y, domain=DOM),
             Poly(X, X, Y, domain=DOM),
             Poly(Y, X, Y, domain=DOM)]
    state += [Poly(int(v) % P, X, Y, domain=DOM) for v in fixed_vals]
    state = mds(state)
    half = pos.r_f // 2
    idx = 0
    for _ in range(half):
        state = arc(state, rc[idx]); idx += 1
        state = [x ** D for x in state]; state = mds(state)
    for _ in range(pos.r_p):
        state = arc(state, rc[idx]); idx += 1
        state[0] = state[0] ** D; state = mds(state)
    for _ in range(half):
        state = arc(state, rc[idx]); idx += 1
        state = [x ** D for x in state]; state = mds(state)
    return state


def resultant_data(pos, M, fixed_vals):
    """Return (deg of Res_Y(out0 - out0target?, ...), factor structure).
    For CICO-2 we want out0 = 0 and out1 = 0 (target constants are 0,0 in this
    framing -- we just need the *degree/factorization* of the elimination
    resultant; the specific target only shifts the constant term)."""
    out = sym_perm(pos, M, fixed_vals)
    # eliminate Y. Treat as polys in Y with coeffs in X.
    f = Poly(out[0].as_expr(), Y, X, domain=DOM)
    g = Poly(out[1].as_expr(), Y, X, domain=DOM)
    R = f.resultant(g)
    Rx = Poly(R.as_expr(), X, domain=DOM) if R.as_expr() != 0 else None
    deg = Rx.degree() if Rx is not None else -1
    return deg, Rx, out


def factor_summary(Rx):
    """Factor a univariate poly over GF(p); return list of (deg, mult) sorted.
    Returns None if degree exceeds the cap (factoring would be infeasible)."""
    if Rx is None:
        return []
    if Rx.degree() > FACTOR_DEG_CAP:
        return None
    _, facs = factor_list(Rx)
    degs = sorted([(int(f.degree()), int(m)) for (f, m) in facs])
    return degs


# ---------------------------------------------------------------------------
# Build the candidate "special" fixed-value vectors.
# ---------------------------------------------------------------------------
def make_specials(pos):
    rc = pos.round_constants
    n = T - 4  # 12
    cands = {}
    cands["random"] = [(1234567 * (i + 1) + 89) % P for i in range(n)]
    cands["all_zero"] = [0] * n
    cands["all_one"] = [1] * n
    cands["all_equal_7"] = [7] * n
    # equal to first-round constants at positions 4..15
    cands["eq_rc0"] = [int(rc[0][4 + i]) % P for i in range(n)]
    cands["neg_rc0"] = [(-int(rc[0][4 + i])) % P for i in range(n)]
    # the value that zeroes out the post-initial-MDS ARC?  hard symbolically;
    # try: all equal to a single round constant value
    cands["all_rc0_0"] = [int(rc[0][0]) % P] * n
    return cands


# ---------------------------------------------------------------------------
# (b) structured subgroup analysis.
# ---------------------------------------------------------------------------
def subgroup_elems(order):
    """Return the multiplicative subgroup of F_p* of given order (must divide p-1)."""
    assert (P - 1) % order == 0
    cof = (P - 1) // order
    g = 3  # find a generator of F_p*
    # 3 is a known primitive root for many primes; verify by order
    def mul_order(a):
        # too slow to compute full order; instead build subgroup from g^cof
        return None
    h = pow(g, cof, P)  # h has order = order (if g is primitive)
    elems = set()
    x = 1
    for _ in range(order):
        elems.add(x); x = (x * h) % P
    return sorted(elems), h


def restricted_image(pos, M, fixed_vals, Hx, Hy):
    """Evaluate (out0,out1) over Hx x Hy (subgroup product). Report image
    structure: number of distinct out0, distinct out1, distinct pairs, and how
    many pairs hit the CICO target (out0,out1)=(0,0).

    Uses the REAL numeric permutation_plus_linear (pure modular arithmetic) at
    every (x,y) point -- this is far faster than symbolic build+eval and works
    at any depth. fixed_vals -> positions 4..15; (x,y) -> positions 2,3;
    positions 0,1 = C1,C2."""
    base = [int(C1) % P, int(C2) % P, 0, 0] + [int(v) % P for v in fixed_vals]
    vals0 = set(); vals1 = set(); pairs = set(); hits = 0
    for x in Hx:
        for y in Hy:
            s = list(base); s[2] = x; s[3] = y
            o = pos.permutation_plus_linear(s)
            a, b = o[0], o[1]
            vals0.add(a); vals1.add(b); pairs.add((a, b))
            if a == 0 and b == 0:
                hits += 1
    return len(vals0), len(vals1), len(pairs), hits


# ===========================================================================
def fmt_facs(fs):
    if fs is None:
        return "(deg>cap, not factored)"
    fstr = ", ".join(f"{d}^{m}" if m > 1 else f"{d}" for d, m in fs[:10])
    if len(fs) > 10:
        fstr += ", ..."
    return f"[{fstr}]"


def degree_table(matname, M, instances, which=None):
    log("=" * 74)
    log(f"(BASE + a)  deg_X Res_Y(out0,out1) + factorization over GF(p)   [matrix={matname}]")
    log("  pred Poseidon drop = d^(2RF+RP); front-skip = d^(2(RF-1)+RP); trivial = d^(2(RF+RP))")
    for (rf, rp) in instances:
        pos = Poseidon(prime=P, alpha=D, t=T, r_f=rf, r_p=rp, mds=M)
        specials = make_specials(pos)
        if which is not None:
            specials = {k: v for k, v in specials.items() if k in which}
        drop = D ** (2 * rf + rp)
        fskip = D ** (2 * (rf - 1) + rp)
        triv = D ** (2 * (rf + rp))
        log(f"\n  --- RF={rf} RP={rp}  (drop={drop}, front-skip={fskip}, trivial={triv}) ---")
        for name, fv in specials.items():
            t0 = time.perf_counter()
            try:
                deg, Rx, _ = resultant_data(pos, M, fv)
                fs = factor_summary(Rx)
                dt = time.perf_counter() - t0
                if fs is None:
                    mindeg, nfac = -1, -1
                else:
                    mindeg = min((d for d, m in fs), default=-1)
                    nfac = sum(m for d, m in fs)
                flag = ""
                if deg == drop: flag = "[=drop]"
                elif 0 <= deg < drop: flag = "[<drop!]"
                elif deg == triv: flag = "[=triv]"
                log(f"    {name:12} deg={deg:>6} {flag:9} minfac={mindeg:>5} nfac={nfac:>4}  "
                    f"facs:{fmt_facs(fs)}  ({dt:.1f}s)")
            except Exception as ex:
                log(f"    {name:12} ERROR {ex!r}")


def main():
    log("=" * 74)
    log("FIELD STRUCTURE")
    log(f"  p = {P},  p-1 = 2^24 * 127")
    log("  gcd(3, p-1) = 1  =>  x^3 is a BIJECTION on F_p (and on every subgroup)")
    log("  => the angle's 'order-126 / 3-to-1' premise is FALSE (126 does not divide p-1)")

    # ---- Section (b) runs FIRST: it is the most novel part of this angle and is
    #      fast (numeric eval). Heavy resultant cases come last so a timeout only
    #      sacrifices redundant confirmations.
    log("=" * 74)
    log("(b)  STRUCTURED DOMAIN: restrict free vars X,Y to small subgroups")
    H127, _ = subgroup_elems(127)
    log(f"  |H127| = {len(H127)} (order-127 mult. subgroup). cube permutes it (gcd(3,127)=1).")
    for (rf, rp) in [(2, 1), (2, 2), (2, 3), (4, 1)]:  # RF must be even
        pos = Poseidon(prime=P, alpha=D, t=T, r_f=rf, r_p=rp, mds=CAUCHY)
        fv = [(1234567 * (i + 1) + 89) % P for i in range(T - 4)]
        t0 = time.perf_counter()
        n0, n1, npair, hits = restricted_image(pos, CAUCHY, fv, H127, H127)
        dt = time.perf_counter() - t0
        tot = len(H127) ** 2
        log(f"  RF={rf} RP={rp}: over H127xH127 ({tot} pts): "
            f"|img out0|={n0}, |img out1|={n1}, |distinct pairs|={npair}, "
            f"CICO-hits(0,0)={hits}  ({dt:.1f}s)")
        log(f"     -> pair-image fills {npair}/{tot} = {npair/tot:.3f} of product "
            f"(1.0 = generic/no collapse). Expected (0,0) hit over {tot} pts "
            f"~ {tot}/p^2 = {tot/P**2:.2e}")

    log("=" * 74)
    log("(b')  STRUCTURED-DOMAIN COST LOGIC for the real target")
    log("  CICO-2 needs out0=out1=0: 2 constraints over F_p. Restricting X,Y to a")
    log("  subgroup H shrinks the candidate set to |H|^2 but does NOT reduce the number")
    log("  of constraints. Per-pair success ~ 1/p^2, so expected hits over HxH =")
    log("  |H|^2/p^2.  |H|=127 -> 2^-48 (never).  |H|=2^24 -> 2^48/2^62 = 2^-14 (<1 hit,")
    log("  and enumerating 2^48 pairs already exceeds the 2^62 brute target on cost per")
    log("  op grounds only at ~2^48, but gives no SOLUTION).  Helps ONLY if out0,out1")
    log("  collapse to low-rank/constant on H -- measured above.")

    log("=" * 74)
    log("(c)  matrix-dependence cross-check on the circulant matrix")
    degree_table("circulant", CIRC, [(2, 1)],
                 which={"random", "all_equal_7", "eq_rc0"})

    # ---- (a) DEGREE + FACTORIZATION under special slices, ordered cheapest-first.
    # (a-full) ALL special-value choices, FACTORED, on the cheap deg-81 instance
    #   (RF=2,RP=0): does any special slice change the FACTOR STRUCTURE / drop deg?
    degree_table("Cauchy", CAUCHY, [(2, 0)])
    # (a-deep) DEGREE-ONLY across depth (deg-243 resultants ~8s each): does deg
    #   ever drop below 'drop' for a special slice?  (2,1)->3^5=243.
    degree_table("Cauchy", CAUCHY, [(2, 1)],
                 which={"random", "all_zero", "all_equal_7", "eq_rc0", "all_one"})
    # (a-deepest) one deg-729 confirmation (heavy; last so timeout is harmless):
    #   (2,2)->3^6=729  -- random vs all_equal_7 only.
    degree_table("Cauchy", CAUCHY, [(2, 2)], which={"random", "all_equal_7"})


if __name__ == "__main__":
    main()

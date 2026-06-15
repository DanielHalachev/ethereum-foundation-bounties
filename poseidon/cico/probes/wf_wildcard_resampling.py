"""
Wildcard angle: SMARTER-THAN-BRUTE search via cheap univariate root-finding.

Observation under test:
  Fix 13 of the 14 free CICO inputs.  Then out0(X) (output word 0 as a function
  of the single remaining free input X) is a UNIVARIATE polynomial over F_p of
  degree 3^(RF+RP).  For the full instance 3^16 ~= 2^25.4, which is "small" as a
  univariate.  Root-finding out0(X)=0 over F_p costs ~ O(D log D log p) field ops
  (build via Frobenius x^p mod out0, gcd, then factor).  This finds X with out0=0.
  But CICO-2 needs out0=0 AND out1=0 simultaneously.

Design candidates costed here:
  (A) "1D resample": repeatedly randomize the 13 fixed inputs; root-find out0(X)=0;
      among the F_p-roots, check whether any also gives out1=0.  Per fixing, P(root
      gives out1=0) ~ (#roots)/p.  This is the cost we measure.
  (B) "2D via two univariates": keep X and Y free; root-find out0 in X for the
      special structure -- but out0 depends on both, so this reduces to the
      resultant (no free lunch).  We argue why.

This probe:
  1. MEASURES the true univariate degree of out0(X) and out1(X) on tiny instances,
     confirming deg = 3^(RF+RP) in one free input (matrix-independent).
  2. MEASURES the expected number of F_p-roots of out0(X)=0 (should be ~1 by the
     "random degree-D poly over F_p has ~1 root" heuristic, Poisson(1)).
  3. END-TO-END demo on a SMALL instance (RF=2,RP=2,t=16): runs strategy (A),
     i.e. resample 13 inputs, root-find out0=0 numerically (eval+root over F_p via
     small-degree GF factorization), test out1, until a real CICO-2 solution is
     found, then VERIFIES it with verify_cico_solution.
  4. COSTS strategy (A) for the full instance and compares to 2^62 brute and 2^55
     resultant.
"""
import sys, os, time, random
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "reference"))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "reference", "bounties"))

from sympy import symbols, Poly, GF
from poseidon.poseidon import Poseidon
from poseidon.mds_matrix import generate_mds_matrix, generate_circulant_mds_matrix

P = 2130706433
T = 16
D = 3
ROW16 = [1, 1, 51, 1, 11, 17, 2, 1, 101, 63, 15, 2, 67, 22, 13, 3]
C1, C2 = 0xC09DE4, 0xEE6282
X = symbols("X")
DOM = GF(P)
random.seed(12345)


# ---------------------------------------------------------------------------
# Symbolic single-variable permutation (only x2 = X is symbolic; rest constant)
# ---------------------------------------------------------------------------
def sym_perm_1var(pos, x2_is_X, fixed_free):
    """permutation_plus_linear with input [C1,C2,X,f3,...,f15] (X = free word #2).
       fixed_free is list of 13 ints for words 3..15."""
    p, t, mds, rc = pos.prime, pos.t, pos.mds, pos.round_constants

    def mds_mul(s):
        return [sum(int(mds[i][j]) * s[j] for j in range(t)) for i in range(t)]

    def add_rc(s, c):
        return [s[i] + int(c[i]) for i in range(t)]

    def sb(x):
        return x ** D

    state = [Poly(int(C1) % P, X, domain=DOM),
             Poly(int(C2) % P, X, domain=DOM),
             Poly(X, X, domain=DOM) if x2_is_X else Poly(int(fixed_free[-1]) % P, X, domain=DOM)]
    state += [Poly(int(v) % P, X, domain=DOM) for v in fixed_free]
    state = state[:t]
    assert len(state) == t

    state = mds_mul(state)
    half = pos.r_f // 2
    idx = 0
    for _ in range(half):
        state = add_rc(state, rc[idx]); idx += 1
        state = [sb(x) for x in state]
        state = mds_mul(state)
    for _ in range(pos.r_p):
        state = add_rc(state, rc[idx]); idx += 1
        state[0] = sb(state[0])
        state = mds_mul(state)
    for _ in range(half):
        state = add_rc(state, rc[idx]); idx += 1
        state = [sb(x) for x in state]
        state = mds_mul(state)
    return state


def measure_univariate_degree():
    print("=== (1) Univariate degree of out0(X), out1(X) in ONE free input ===")
    print(f"{'matrix':12} {'RF':>2} {'RP':>2} {'deg out0':>9} {'deg out1':>9} {'3^(RF+RP)':>10} {'build_s':>8}")
    mats = {
        "circulant": generate_circulant_mds_matrix(ROW16, P),
        "cauchy": generate_mds_matrix(T, P),
    }
    rows = []
    for (rf, rp) in [(2, 0), (2, 1), (2, 2), (2, 3), (4, 2)]:
        for name, M in mats.items():
            pos = Poseidon(prime=P, alpha=D, t=T, r_f=rf, r_p=rp, mds=M)
            fixed = [(987654 * (i + 1) + 13) % P for i in range(13)]
            t0 = time.perf_counter()
            out = sym_perm_1var(pos, True, fixed)
            d0 = out[0].degree(); d1 = out[1].degree()
            bt = time.perf_counter() - t0
            pred = D ** (rf + rp)
            print(f"{name:12} {rf:>2} {rp:>2} {d0:>9} {d1:>9} {pred:>10} {bt:>8.2f}")
            rows.append((name, rf, rp, d0, d1, pred))
    ok = all(d0 == pred and d1 == pred for (_, _, _, d0, d1, pred) in rows)
    print(f"  -> univariate degree == 3^(RF+RP) for all: {ok}")
    print(f"  -> for FULL instance RF=6,RP=10: deg out0 = 3^16 = {3**16} = 2^{(3**16).bit_length()-1:.0f}.. bits")
    return ok


# ---------------------------------------------------------------------------
# Root counting: how many F_p roots does a random out0(X)=0 have?
# Use the GF(p)-factor / number-of-roots = deg gcd(X^p - X, f).
# Computing X^p mod f costs ~ log2(p) squarings of a deg-D poly: ~31 * M(D).
# We just COUNT roots on small instances (small deg) to validate Poisson(1).
# ---------------------------------------------------------------------------
def count_fp_roots(poly_in_X):
    f = Poly(poly_in_X, X, domain=DOM)
    if f.degree() <= 0:
        return 0
    # number of distinct F_p-roots = deg gcd(f, X^p - X)
    # compute X^p mod f by repeated squaring
    xp = Poly(X, X, domain=DOM)
    base = Poly(X, X, domain=DOM)
    e = P
    acc = Poly(1, X, domain=DOM)
    # acc = X^P mod f
    b = base
    while e > 0:
        if e & 1:
            acc = (acc * b) % f
        e >>= 1
        if e:
            b = (b * b) % f
    g = (acc - Poly(X, X, domain=DOM)) % f
    from sympy import gcd as sgcd
    gg = sgcd(g, f)
    return gg.degree() if gg.degree() > 0 else 0


def measure_root_count():
    print("\n=== (2) #F_p-roots of out0(X)=0 (target: ~Poisson(1), mean ~1) ===")
    pos = Poseidon(prime=P, alpha=D, t=T, r_f=2, r_p=2,
                   mds=generate_circulant_mds_matrix(ROW16, P))
    counts = []
    N = 30
    t0 = time.perf_counter()
    for trial in range(N):
        fixed = [random.randrange(P) for _ in range(13)]
        out = sym_perm_1var(pos, True, fixed)
        # solve out0(X) = 0  (out0 already includes no target; CICO wants out0==0)
        nr = count_fp_roots(out[0].as_expr())
        counts.append(nr)
    dt = time.perf_counter() - t0
    mean = sum(counts) / len(counts)
    frac_pos = sum(1 for c in counts if c > 0) / len(counts)
    print(f"  RF=2,RP=2 (deg {3**4}): {N} trials, mean #roots = {mean:.2f}, "
          f"P(>=1 root) = {frac_pos:.2f}, time {dt:.1f}s")
    print(f"  counts: {counts}")
    return mean, frac_pos


# ---------------------------------------------------------------------------
# (3) END-TO-END strategy (A) on a small instance: find a real CICO-2 solution.
#     For each randomization of the 13 fixed words: form out0(X), out1(X) as
#     univariate polys; find common F_p-roots of {out0=0, out1=0} = roots of
#     gcd over F_p; if any X-root exists, we have a CICO-2 solution (x2=X, rest
#     = the 13 fixed words).  But that's just the resultant in 1 var pre-fixed.
#     The PURE wildcard variant: root-find out0=0 only, then TEST out1 at those
#     roots.  We implement the pure-wildcard variant and verify.
# ---------------------------------------------------------------------------
def fp_roots_list(poly_in_X):
    """Return the list of distinct F_p roots of poly_in_X (small degree only)."""
    f = Poly(poly_in_X, X, domain=DOM)
    if f.degree() <= 0:
        return []
    xp = Poly(X, X, domain=DOM)
    b = Poly(X, X, domain=DOM)
    acc = Poly(1, X, domain=DOM)
    e = P
    while e > 0:
        if e & 1:
            acc = (acc * b) % f
        e >>= 1
        if e:
            b = (b * b) % f
    from sympy import gcd as sgcd
    g = (acc - Poly(X, X, domain=DOM)) % f
    gg = sgcd(g, f)
    if gg.degree() <= 0:
        return []
    # split gg (squarefree, all roots in F_p) into linear factors
    fac = gg.factor_list()[1]
    roots = []
    for (fct, mult) in fac:
        if fct.degree() == 1:
            a = int(fct.nth(1)); b0 = int(fct.nth(0))
            r = (-b0 * pow(a, -1, P)) % P
            roots.append(r)
    return roots


def end_to_end_small():
    print("\n=== (3) END-TO-END pure-wildcard CICO-2 on small instance RF=2,RP=2 ===")
    from cico_verifier import verify_cico_solution
    rf, rp = 2, 2
    M = generate_circulant_mds_matrix(ROW16, P)
    pos = Poseidon(prime=P, alpha=D, t=T, r_f=rf, r_p=rp, mds=M)
    constants = [C1, C2, 0, 0]  # need out0==0, out1==0
    tries = 0
    t0 = time.perf_counter()
    found = None
    while tries < 200000 and time.perf_counter() - t0 < 120:
        tries += 1
        fixed = [random.randrange(P) for _ in range(13)]  # words 3..15
        out = sym_perm_1var(pos, True, fixed)
        roots = fp_roots_list(out[0].as_expr())  # X s.t. out0 = 0
        for r in roots:
            # test out1 at X=r by plugging into out1 poly (cheap eval)
            v1 = int(Poly(out[1].as_expr(), X, domain=DOM).eval(r)) % P
            if v1 == 0:
                found = (r, fixed)
                break
        if found:
            break
    dt = time.perf_counter() - t0
    if not found:
        print(f"  no solution in {tries} tries / {dt:.1f}s (expected since "
              f"P(out1=0 | out0=0) ~ #roots/p ~ 2^-31)")
        return False, tries, dt
    r, fixed = found
    free_inputs = [r] + fixed  # words 2..15  (length 14)
    ok = verify_cico_solution(free_inputs, r_f=rf, r_p=rp, mds=M, constants=constants)
    print(f"  FOUND after {tries} tries / {dt:.1f}s ; verify_cico_solution = {ok}")
    print(f"  x2={r}, fixed words 3..15 = {fixed[:3]}...")
    return ok, tries, dt


# ---------------------------------------------------------------------------
# (3b) Control: cross-check that the SAME small instance is solvable cheaply by
#      the bivariate (resultant-in-one-var) variant, to confirm the instance is
#      not degenerate, and to MEASURE how the wildcard compares.
#      Here: keep X free, also let word 3 = second symbolic var? No -- to stay 1D
#      and demonstrate a *genuine* shortcut would need out1's root set to overlap.
# ---------------------------------------------------------------------------
def end_to_end_resultant_small():
    print("\n=== (3b) Control: bivariate resultant CICO-2 (1 free word root-find both) ===")
    # Strategy: fix 12 words, keep words 2 (X) and 3 (Z) free; build out0(X,Z),
    # out1(X,Z); resultant over Z gives R0(X) of degree d^(2RF+RP); root-find;
    # back-substitute Z.  This is the KNOWN attack; we run it small to confirm a
    # solution exists & measure D_I, for honest comparison to wildcard.
    from sympy import symbols as syms, resultant
    from cico_verifier import verify_cico_solution
    rf, rp = 2, 2
    M = generate_circulant_mds_matrix(ROW16, P)
    pos = Poseidon(prime=P, alpha=D, t=T, r_f=rf, r_p=rp, mds=M)
    Xs, Zs = syms("Xs Zs")
    DOM2 = GF(P)
    mds, rc, t = pos.mds, pos.round_constants, pos.t

    def mds_mul(s): return [sum(int(mds[i][j]) * s[j] for j in range(t)) for i in range(t)]
    def add_rc(s, c): return [s[i] + int(c[i]) for i in range(t)]
    def sb(x): return x ** D

    fixed12 = [random.randrange(P) for _ in range(12)]  # words 4..15
    state = [Poly(int(C1) % P, Xs, Zs, domain=DOM2),
             Poly(int(C2) % P, Xs, Zs, domain=DOM2),
             Poly(Xs, Xs, Zs, domain=DOM2),
             Poly(Zs, Xs, Zs, domain=DOM2)]
    state += [Poly(int(v) % P, Xs, Zs, domain=DOM2) for v in fixed12]
    state = mds_mul(state)
    half = rf // 2; idx = 0
    for _ in range(half):
        state = add_rc(state, rc[idx]); idx += 1
        state = [sb(x) for x in state]; state = mds_mul(state)
    for _ in range(rp):
        state = add_rc(state, rc[idx]); idx += 1
        state[0] = sb(state[0]); state = mds_mul(state)
    for _ in range(half):
        state = add_rc(state, rc[idx]); idx += 1
        state = [sb(x) for x in state]; state = mds_mul(state)
    P0 = Poly(state[0].as_expr(), Xs, Zs, domain=DOM2)  # out0
    Q0 = Poly(state[1].as_expr(), Xs, Zs, domain=DOM2)  # out1
    t0 = time.perf_counter()
    f = Poly(P0.as_expr(), Zs, Xs, domain=DOM2)
    g = Poly(Q0.as_expr(), Zs, Xs, domain=DOM2)
    R = f.resultant(g)
    DI = R.degree()
    # root-find R(X)=0
    rootsX = fp_roots_list(R.as_expr().subs(syms("Zs"), 0) if False else R.as_expr())
    # R is in Xs; adapt fp_roots_list to Xs
    Rx = Poly(R.as_expr(), Xs, domain=DOM2)
    # reuse generic root finder
    def roots_in(var, poly):
        f2 = Poly(poly, var, domain=DOM2)
        if f2.degree() <= 0: return []
        b = Poly(var, var, domain=DOM2); acc = Poly(1, var, domain=DOM2); e = P
        while e > 0:
            if e & 1: acc = (acc * b) % f2
            e >>= 1
            if e: b = (b * b) % f2
        from sympy import gcd as sgcd
        gg = sgcd((acc - Poly(var, var, domain=DOM2)) % f2, f2)
        if gg.degree() <= 0: return []
        rs = []
        for (fct, m) in gg.factor_list()[1]:
            if fct.degree() == 1:
                a = int(fct.nth(1)); b0 = int(fct.nth(0))
                rs.append((-b0 * pow(a, -1, P)) % P)
        return rs
    Xroots = roots_in(Xs, R.as_expr())
    sol = None
    for xr in Xroots:
        # solve out0(xr, Z)=0 for Z, intersect with out1(xr,Z)=0
        p0z = Poly(P0.as_expr(), Xs, Zs, domain=DOM2).eval({Xs: xr})
        q0z = Poly(Q0.as_expr(), Xs, Zs, domain=DOM2).eval({Xs: xr})
        z0 = set(roots_in(Zs, Poly(p0z, Zs, domain=DOM2).as_expr()))
        z1 = set(roots_in(Zs, Poly(q0z, Zs, domain=DOM2).as_expr()))
        common = z0 & z1
        if common:
            sol = (xr, next(iter(common)))
            break
    dt = time.perf_counter() - t0
    if sol is None:
        print(f"  resultant D_I={DI}, no common root this fixing ({dt:.1f}s).")
        return False, DI
    xr, zr = sol
    free_inputs = [xr, zr] + fixed12
    ok = verify_cico_solution(free_inputs, r_f=rf, r_p=rp, mds=M, constants=[C1, C2, 0, 0])
    print(f"  resultant D_I={DI} (pred d^(2RF+RP)={D**(2*rf+rp)}); SOLVED, "
          f"verify={ok}, time {dt:.1f}s")
    return ok, DI


# ---------------------------------------------------------------------------
# (4) Cost model for the FULL instance under strategy (A) pure-wildcard.
# ---------------------------------------------------------------------------
def cost_model(mean_roots, frac_pos):
    print("\n=== (4) Cost model: pure-wildcard strategy (A) for FULL RF=6,RP=10 ===")
    Dfull = D ** (6 + 10)          # univariate degree of out0
    Dbi = D ** (2 * 6 + 10)        # resultant degree (1 front skip -> 2*(6-1)+10)
    Dbi_skip = D ** (2 * (6 - 1) + 10)
    log2 = lambda x: (x).bit_length() - 1
    print(f"  univariate deg out0 = 3^16 = {Dfull} ~ 2^{log2(Dfull)}")
    # cost of ONE root-find of out0(X)=0 over F_p:
    #   build poly out0(X): need its 3^16 coefficients. THIS is the hidden cost.
    #   Either (a) symbolic build: O(3^16) terms -> ~2^25 monomials but building
    #   them through the rounds is ~ (#rounds)*M(3^16) ~ 16 * 2^25*25 ~ 2^33 ops
    #   AND 2^25 * 4 bytes = ~256 MB to STORE the dense coeff vector. Feasible-ish.
    #   (b) root-find = Frobenius X^p mod out0: log2(p)=31 squarings of deg-3^16
    #   poly. Each squaring via NTT ~ M(D)=D log D ~ 2^25 * 25 = 2^30. Times 31 =
    #   ~2^35. Plus gcd ~ M(D) log D ~ 2^30. So per root-find ~ 2^35 ops.
    per_rootfind = 31 * (Dfull * log2(Dfull))   # ~ Frobenius cost in field ops
    print(f"  per out0=0 root-find (Frobenius) ~ 31 * D*log2 D = 2^{log2(per_rootfind)} field ops")
    print(f"  storage for dense out0 coeffs ~ 2^25 * 8B = {Dfull*8/1e9:.2f} GB")
    # Each root-find yields ~mean_roots roots. For each, P(out1 root)=1/p.
    # So per FIXING, P(CICO hit) = mean_roots / p  (each independent root ~ unif).
    # Expected #fixings = p / mean_roots ~ 2^31.
    print(f"  measured mean #out0-roots per fixing = {mean_roots:.2f} (Poisson(1)),")
    print(f"  P(a given out0-root also gives out1=0) = 1/p = 2^-31")
    print(f"  => P(CICO hit per fixing) = mean_roots/p ~ 2^{log2(int(P/max(mean_roots,1e-9)))*-1 if False else ''}")
    exp_fixings = P / max(mean_roots, 1e-9)
    print(f"  => expected #fixings ~ p/mean = 2^{log2(int(exp_fixings))}")
    total = exp_fixings * per_rootfind
    print(f"  => TOTAL strategy-(A) cost ~ 2^{log2(int(exp_fixings))} * 2^{log2(per_rootfind)} "
          f"= 2^{log2(int(total))} field ops")
    print(f"  COMPARE: brute = 2^62 ; resultant(skip-1) D_I=2^{log2(Dbi_skip)} solve ~2^55")
    # The honest punchline:
    print("\n  HONEST ANALYSIS:")
    print("  * Per-fixing we resample ALL 13 fixed words -> out0 changes -> must")
    print("    REBUILD out0 (2^25 coeffs) AND re-run a 2^35 root-find. The root-find")
    print("    only enforces out0=0; out1=0 is still a fresh 1/p gamble per root.")
    print("  * mean_roots ~ 1, so each fixing gives ~1 candidate X with out0=0,")
    print("    and that candidate passes out1=0 with prob 1/p. => ~p ~ 2^31 fixings,")
    print("    each costing >= 2^35 (root-find) + 2^25 (rebuild). Total >= 2^66.")
    print("  * Pure brute force = evaluate perm directly: 2^31 fixings * 2^31 X-values")
    print("    but per (fixing) the X-sweep is the SAME 2^31 evals; total 2^62, and")
    print("    each eval is ONE permutation (~16*16 ~ 2^8 ops) not 2^35.")
    print("  * So replacing a brute X-sweep (2^31 cheap perm evals) by a root-find")
    print("    (2^35 field ops) is STRICTLY WORSE per fixing, AND yields only the")
    print("    out0=0 slice (~1 point) which we'd get from the sweep anyway.")
    return total


def main():
    ok_deg = measure_univariate_degree()
    mean, frac = measure_root_count()
    e2e_ok, tries, dt = end_to_end_small()
    res_ok, DI = end_to_end_resultant_small()
    cost_model(mean, frac)
    print("\n=== SUMMARY ===")
    print(f"univariate-degree law confirmed: {ok_deg}")
    print(f"mean out0-roots/fixing: {mean:.2f}; P(>=1): {frac:.2f}")
    print(f"pure-wildcard small e2e found+verified: {e2e_ok} (tries={tries}, {dt:.1f}s)")
    print(f"resultant control solved+verified: {res_ok} (D_I={DI})")


if __name__ == "__main__":
    main()

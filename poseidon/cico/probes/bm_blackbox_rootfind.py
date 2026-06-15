"""
Probe: black-box univariate root-finding for the CICO resultant R(X).

CONTEXT
-------
After pinning x0=C1, x1=C2 and fixing x4..x15, the two output words become
bivariate polynomials out0(X,Y), out1(X,Y) over F_p with X=x2, Y=x3, each of
total degree delta = 3^(RF+RP) in (X,Y). Eliminating Y gives a univariate
    R(X) = Res_Y(out0, out1)   of degree  D_I = delta^2 = 3^(2RF+RP).
For the real instance (skip-1) D_I = 3^20 ~ 2^31.7. R has ~1 root in F_p.

Two known cost regimes:
  (A) INTERPOLATE+FACTOR: get the D_I+1 coefficients of R by evaluating it at
      D_I+1 points (each eval = a univariate resultant in Y, cost ~Otilde(delta)
      via half-gcd ~ delta*polylog), then Cantor-Zassenhaus. Cost ~ Otilde(D_I*delta).
  (B) SCAN: try x in F_p one at a time; for each, gcd(out0(x,.),out1(x,.)) and
      read off the shared root. Cost ~ p*Otilde(delta).

QUESTIONS PROBED
----------------
Q1. Is there a black-box root-finder for a degree-D poly with ~1 F_p-root that
    beats min( O(D) evals to interpolate, O(p) evals to scan )?  -> theory + test.
Q2. Does evaluating R at one point really cost ~Otilde(delta) (half-gcd resultant
    of two univariate-in-Y polys of degree delta) -- i.e. << D_I?  -> measure.
Q3. The 10 partial rounds are rank-1 updates of the state map (only word 0 is
    cubed). Does this give R a low DISPLACEMENT RANK / quasi-Toeplitz structure
    that lets us form R in o(D_I*delta)?  -> measure displacement rank of the
    Sylvester/Bezout structure on tiny instances.
Q4. Lower bound: can an adversary force Omega(D) black-box queries to pin a root?

We instantiate at small primes / small (RF,RP) so everything is computable, and
read off the structural quantities that are field-size robust.
"""
import sys, os, time, random
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "reference"))
from sympy import symbols, Poly, GF, ZZ, Matrix
from poseidon.poseidon import Poseidon
from poseidon.mds_matrix import generate_mds_matrix

X, Y = symbols("X Y")


def build_outs(p, t, rf, rp, fixed_tail, C1, C2, seed=0):
    """Symbolic out0,out1 in X=word2, Y=word3 for permutation_plus_linear."""
    M = generate_mds_matrix(t, p)
    pos = Poseidon(prime=p, alpha=3, t=t, r_f=rf, r_p=rp, mds=M)
    DOM = GF(p)
    rc = pos.round_constants
    def mds(s): return [sum(int(M[i][j]) * s[j] for j in range(t)) for i in range(t)]
    def arc(s, c): return [s[i] + int(c[i]) for i in range(t)]
    state = [Poly(C1 % p, X, Y, domain=DOM), Poly(C2 % p, X, Y, domain=DOM),
             Poly(X, X, Y, domain=DOM), Poly(Y, X, Y, domain=DOM)]
    state += [Poly(int(v) % p, X, Y, domain=DOM) for v in fixed_tail]
    state = mds(state)
    def cube(x): return x * x * x
    half = rf // 2; idx = 0
    for _ in range(half):
        state = arc(state, rc[idx]); idx += 1
        state = [cube(x) for x in state]; state = mds(state)
    for _ in range(rp):
        state = arc(state, rc[idx]); idx += 1
        state[0] = cube(state[0]); state = mds(state)
    for _ in range(half):
        state = arc(state, rc[idx]); idx += 1
        state = [cube(x) for x in state]; state = mds(state)
    return state[0].as_expr(), state[1].as_expr(), p


def total_deg(expr):
    pol = Poly(expr, X, Y)
    return pol.total_degree()


# ---------- Q2: cost of a SINGLE black-box evaluation of R(x*) -------------
def eval_R_at(o0, o1, p, x0):
    """R(x0) = Res_Y( o0(x0,Y), o1(x0,Y) ). Univariate resultant -> half-gcd cost."""
    DOM = GF(p)
    f = Poly(o0.subs(X, x0), Y, domain=DOM)
    g = Poly(o1.subs(X, x0), Y, domain=DOM)
    if f.is_zero or g.is_zero:
        return 0
    return int(f.resultant(g))


def full_resultant(o0, o1, p):
    DOM = GF(p)
    f = Poly(o0, Y, X, domain=DOM)
    g = Poly(o1, Y, X, domain=DOM)
    return f.resultant(g)  # Poly in X over GF(p)


# ---------- Q3: displacement rank of the Sylvester resultant structure -----
def sylvester_in_Y(o0, o1, p):
    """Return Sylvester matrix (entries are polys in X) coefficients to inspect
       structure. We return the two coeff-vectors (in Y) whose coeffs are polys in X."""
    f = Poly(o0, Y)  # coeffs are polys in X
    g = Poly(o1, Y)
    return f, g


def main():
    random.seed(1)
    print("=" * 78)
    print("Q2/Q3 STRUCTURE & COST: tiny instances, degrees / displacement rank")
    print("=" * 78)
    # small prime, small (rf, rp); keep degrees tiny so resultants are instant
    # r_f must be even (reference enforces). Use rf=2 (1 front+1 back) and rf=4.
    cases = [
        (1009, 4, 2, 0),  # delta = 3^2 = 9,   D_I = 3^4 = 81
        (1009, 4, 2, 1),  # delta = 3^3 = 27,  D_I = 3^5 = 243
        (1009, 6, 2, 2),  # delta = 3^4 = 81,  D_I = 3^6 = 729
    ]
    C1, C2 = 0xC09DE4, 0xEE6282
    for (p, t, rf, rp) in cases:
        tail = [random.randrange(p) for _ in range(t - 4)]
        t0 = time.perf_counter()
        o0, o1, _ = build_outs(p, t, rf, rp, tail, C1, C2)
        d0, d1 = total_deg(o0), total_deg(o1)
        delta_pred = 3 ** (rf + rp)
        # full resultant
        R = full_resultant(o0, o1, p)
        degR = R.degree()
        DI_pred = 3 ** (2 * rf + rp)
        # how many F_p roots does R have?
        try:
            roots = R.ground_roots()
            nroot = sum(roots.values())
        except Exception:
            nroot = "?"
        # degree of R in X vs delta^2
        # measure: is R full-degree?  is delta achieved?
        # displacement / Sylvester structure: degrees of the coeff polys (in X) of o_i(Y)
        f = Poly(o0, Y); g = Poly(o1, Y)
        degsf = [Poly(c, X).degree() if not Poly(c, X).is_zero else -1 for c in f.all_coeffs()]
        degsg = [Poly(c, X).degree() if not Poly(c, X).is_zero else -1 for c in g.all_coeffs()]
        dt = time.perf_counter() - t0
        print(f"\np={p} t={t} RF={rf} RP={rp}  ({dt:.2f}s)")
        print(f"  total_deg(o0,o1) = ({d0},{d1})   delta_pred=3^{rf+rp}={delta_pred}")
        print(f"  deg_X R = {degR}   D_I_pred=3^{2*rf+rp}={DI_pred}   full-degree={degR==DI_pred}")
        print(f"  #F_p-roots(R) = {nroot}")
        print(f"  deg_X of o0 coeffs-in-Y: {degsf}")
        print(f"  deg_X of o1 coeffs-in-Y: {degsg}")
        # Q2 sanity: eval R at a point and confirm it matches R(x*)
        xt = random.randrange(p)
        val_direct = eval_R_at(o0, o1, p, xt)
        val_fromR = int(R.eval(xt)) % p
        print(f"  Q2 check: R({xt}) via 1-point resultant = {val_direct}, via full R = {val_fromR}, "
              f"match={val_direct==val_fromR}", flush=True)

        # ---- Q3: displacement rank of the Sylvester matrix S(x) (entries in X) ----
        # If S(x) had small displacement rank as a structured matrix, R could be
        # formed in o(D_I) per the structured-resultant theory. Build S at a numeric
        # x and measure displacement rank D = rank(S - Z S Z^T).
        dr = displacement_rank_sylvester(f, g, p, xt)
        n = len(f.all_coeffs()) - 1 + len(g.all_coeffs()) - 1  # Sylvester dim ~ degf+degg
        print(f"  Q3 displacement rank of Sylvester(in Y) at x={xt}: {dr}  (matrix dim {n})", flush=True)


def displacement_rank_sylvester(fY, gY, p, xval):
    """Sylvester matrix of fY,gY (polys in Y, coeffs in X) evaluated at X=xval.
       Return rank over GF(p) of S - Z S Z^T  (Z = down-shift), the displacement rank.
       Small (O(1)) displacement rank => Toeplitz-like => fast structured solve."""
    fc = [int(Poly(c, X).eval(xval)) % p for c in fY.all_coeffs()]  # high..low in Y
    gc = [int(Poly(c, X).eval(xval)) % p for c in gY.all_coeffs()]
    m, nn = len(fc) - 1, len(gc) - 1
    N = m + nn
    if N == 0:
        return 0
    S = [[0] * N for _ in range(N)]
    # standard Sylvester layout
    for i in range(nn):
        for j, c in enumerate(fc):
            S[i][i + j] = c % p
    for i in range(m):
        for j, c in enumerate(gc):
            S[nn + i][i + j] = c % p
    # Z down-shift NxN
    def shift(Mrows):
        out = [[0] * N for _ in range(N)]
        for i in range(1, N):
            for j in range(1, N):
                out[i][j] = Mrows[i - 1][j - 1]
        return out
    ZSZ = shift(S)
    Dmat = [[(S[i][j] - ZSZ[i][j]) % p for j in range(N)] for i in range(N)]
    return rank_gfp(Dmat, p)


def rank_gfp(rows, p):
    rows = [r[:] for r in rows]
    nr = len(rows); nc = len(rows[0]) if nr else 0
    rank = 0; pr = 0
    for pc in range(nc):
        piv = -1
        for r in range(pr, nr):
            if rows[r][pc] % p != 0:
                piv = r; break
        if piv < 0:
            continue
        rows[pr], rows[piv] = rows[piv], rows[pr]
        inv = pow(rows[pr][pc], -1, p)
        rows[pr] = [(x * inv) % p for x in rows[pr]]
        for r in range(nr):
            if r != pr and rows[r][pc] % p:
                f = rows[r][pc]
                rows[r] = [(rows[r][c] - f * rows[pr][c]) % p for c in range(nc)]
        pr += 1; rank += 1
        if pr == nr:
            break
    return rank


if __name__ == "__main__":
    main()


# ============================================================================
# Q3b: Is displacement-rank-2 generic (true of EVERY Sylvester matrix)?
#      And does it give leverage beyond the half-gcd we already use?
# Q4 : Lower bound -- forcing Omega(D) queries.
# Run: ./.venv/bin/python -c "import probes.bm_blackbox_rootfind as m; m.q3b_q4()"
# ============================================================================
def q3b_q4():
    import random as _r
    print("=" * 78)
    print("Q3b: displacement rank of Sylvester for RANDOM (unstructured) poly pairs")
    print("=" * 78)
    p = 1009
    for (df, dg) in [(9, 9), (27, 27), (50, 50), (81, 81)]:
        fc = [_r.randrange(1, p)] + [_r.randrange(p) for _ in range(df)]
        gc = [_r.randrange(1, p)] + [_r.randrange(p) for _ in range(dg)]
        dr = _disp_rank_from_coeffs(fc, gc, p)
        print(f"  random pair deg ({df},{dg}): Sylvester displacement rank = {dr}")
    print("  -> If this is also ~2, displacement-rank-2 is GENERIC (no special")
    print("     structure from the partial rounds): it's just the universal fact")
    print("     'Sylvester = 2 stacked Toeplitz bands', already used by half-gcd.")


def _disp_rank_from_coeffs(fc, gc, p):
    m, nn = len(fc) - 1, len(gc) - 1
    N = m + nn
    S = [[0] * N for _ in range(N)]
    for i in range(nn):
        for j, c in enumerate(fc):
            S[i][i + j] = c % p
    for i in range(m):
        for j, c in enumerate(gc):
            S[nn + i][i + j] = c % p
    Z = [[0] * N for _ in range(N)]
    for i in range(1, N):
        Z[i][i - 1] = 1
    def mm(A, B):
        return [[sum(A[i][k] * B[k][j] for k in range(N)) % p for j in range(N)] for i in range(N)]
    ZS = mm(Z, S)
    ZT = [[Z[j][i] for j in range(N)] for i in range(N)]
    ZSZ = mm(ZS, ZT)
    D = [[(S[i][j] - ZSZ[i][j]) % p for j in range(N)] for i in range(N)]
    return rank_gfp(D, p)

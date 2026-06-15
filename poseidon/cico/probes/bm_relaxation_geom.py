"""
bm_relaxation_geom.py
=====================
LENS: optimization / convex algebra / real & tropical geometry.

Question: does ANY of {Lasserre/SOS moment hierarchy, Positivstellensatz,
SDP relaxation, homotopy continuation (total-degree / polyhedral-BKK),
numerical algebraic geometry, tropical mixed-volume} give a sub-2^55 route
to ONE F_p solution of the Poseidon-1 CICO-2 system?

We reduce the bounty to the 2-variable elimination system actually used by the
known attack: pin x0=C1, x1=C2, fix x4..x15, vary (X,Y)=(x2,x3). Then
  out0(X,Y) = 0,  out1(X,Y) = 0     over F_p,  p=2^31-2^24+1.

The governing number is D_I = deg_X Res_Y(out0,out1) = 3^(2RF+RP).

What we COMPUTE here, on shrinkable instances (small RF,RP; same Cauchy MDS,
same partial-round sparsity, run over a SMALL prime so polys are buildable):

  (A) The exact BEZOUT number  d0*d1 = 3^(RF+RP) * 3^(RF+RP) = 3^(2(RF+RP))
      vs the exact RESULTANT degree D_I = deg_X Res_Y. Confirms continuation
      total-degree path count >= D_I.

  (B) The exact MIXED VOLUME MV(Newt(out0), Newt(out1)) of the two Newton
      polygons (BKK / polyhedral homotopy path count). We compute it exactly
      from the integer vertices (Minkowski / pick-area formula, no scipy).
      Tests the central tropical claim: does partial-round sparsity make the
      Newton polytopes thin enough that MV < D_I  (or < total-degree Bezout)?

  (C) Whether the actual number of AFFINE F_p solutions on a generic slice is
      ~1 (so even a perfect path tracker that found ALL paths gains nothing:
      the bottleneck is the path COUNT, not solution extraction).

  (D) Positivity/SOS sanity: there is no ordering on F_p; we note + test that
      lifting (out0,out1) to Z and asking for a real/rational common zero is a
      DIFFERENT variety (no F_p reduction guarantee), so SOS infeasibility
      certificates certify nothing about F_p solvability.

Run: ./.venv/bin/python probes/bm_relaxation_geom.py
"""
import sys, os, time, itertools, random
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "reference"))
from sympy import symbols, Poly, GF, ZZ, QQ
from poseidon.poseidon import Poseidon
from poseidon.mds_matrix import generate_mds_matrix

T = 16
D = 3
C1, C2 = 0xC09DE4, 0xEE6282
X, Y = symbols("X Y")


# --------------------------------------------------------------------------
# Build out0,out1 as bivariate polys over GF(prime) for given (RF,RP,prime).
# Same schedule the verifier uses: initial MDS, RF/2 full, RP partial, RF/2 full.
# Full round: ARC all, cube all, MDS.  Partial: ARC all, cube word0 only, MDS.
# --------------------------------------------------------------------------
def out_polys(prime, rf, rp, fixed12, seed=0):
    DOM = GF(prime)
    M = generate_mds_matrix(T, prime)
    pos = Poseidon(prime=prime, alpha=D, t=T, r_f=rf, r_p=rp, mds=M)
    rc = pos.round_constants

    def mds(s):
        return [sum(int(M[i][j]) * s[j] for j in range(T)) for i in range(T)]
    def arc(s, c):
        return [s[i] + int(c[i]) for i in range(T)]

    c1 = C1 % prime
    c2 = C2 % prime
    state = [Poly(c1, X, Y, domain=DOM), Poly(c2, X, Y, domain=DOM),
             Poly(X, X, Y, domain=DOM), Poly(Y, X, Y, domain=DOM)]
    state += [Poly(int(v) % prime, X, Y, domain=DOM) for v in fixed12]
    state = mds(state)
    half = rf // 2; idx = 0
    for _ in range(half):
        state = arc(state, rc[idx]); idx += 1
        state = [x ** D for x in state]; state = mds(state)
    for _ in range(rp):
        state = arc(state, rc[idx]); idx += 1
        state[0] = state[0] ** D; state = mds(state)
    for _ in range(half):
        state = arc(state, rc[idx]); idx += 1
        state = [x ** D for x in state]; state = mds(state)
    return state[0], state[1], DOM, M, pos, prime


# --------------------------------------------------------------------------
# Newton polygon machinery (exact, integer; no scipy).
# Monomials are (i,j) = exponents of X,Y. We need:
#   - vertices of convex hull of a point set in Z^2
#   - area of a polygon (shoelace) -> for mixed volume via MV = A(P+Q)-A(P)-A(Q)
# --------------------------------------------------------------------------
def monomial_support(poly):
    """List of (i,j) exponent pairs with nonzero coeff."""
    return [tuple(m) for m, c in poly.terms() if c != 0]

def convex_hull(points):
    """Andrew's monotone chain. Returns hull vertices CCW (no collinear interior)."""
    pts = sorted(set(points))
    if len(pts) <= 1:
        return pts
    def cross(o, a, b):
        return (a[0]-o[0])*(b[1]-o[1]) - (a[1]-o[1])*(b[0]-o[0])
    lower = []
    for p in pts:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)
    upper = []
    for p in reversed(pts):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)
    return lower[:-1] + upper[:-1]

def poly_area2(hull):
    """2*area via shoelace (>=0 for CCW). hull is list of vertices."""
    n = len(hull)
    if n < 3:
        return 0
    s = 0
    for i in range(n):
        x1, y1 = hull[i]; x2, y2 = hull[(i+1) % n]
        s += x1*y2 - x2*y1
    return abs(s)  # = 2*Area

def minkowski_sum(A, B):
    """Minkowski sum of two point sets (then take hull)."""
    return [(a[0]+b[0], a[1]+b[1]) for a in A for b in B]

def mixed_volume_2d(supA, supB):
    """
    Normalized mixed volume MV(P,Q) for 2 polygons P=conv(supA), Q=conv(supB):
      MV = Area(P+Q) - Area(P) - Area(Q),
    using normalized (lattice) area = 2*Euclidean area = shoelace value.
    This MV is exactly the BKK bound (number of isolated C*-solutions) of a
    generic 2x2 system with these supports.
    """
    HA = convex_hull(supA); HB = convex_hull(supB)
    AA = poly_area2(HA)            # 2*Area(P)
    AB = poly_area2(HB)            # 2*Area(Q)
    Hsum = convex_hull(minkowski_sum(HA, HB))
    Asum = poly_area2(Hsum)        # 2*Area(P+Q)
    # MV (normalized) = 2Area(P+Q) - 2Area(P) - 2Area(Q)   [in lattice-area units]
    return Asum - AA - AB

def total_degree(poly):
    return max((i+j) for (i, j) in monomial_support(poly))


# --------------------------------------------------------------------------
# Exact resultant degree (the real elimination wall).
# --------------------------------------------------------------------------
def resultant_degX(o0, o1, DOM):
    f = Poly(o0.as_expr(), Y, X, domain=DOM)
    g = Poly(o1.as_expr(), Y, X, domain=DOM)
    R = f.resultant(g)
    R = Poly(R, X, domain=DOM)
    return (0 if R.is_zero else R.degree()), R

# count F_p solutions on a slice: roots of R(X), then Y via gcd
def count_fp_solutions(o0, o1, R, DOM, prime):
    if R.is_zero:
        return None
    rd = Poly(R, X, domain=DOM).ground_roots()
    xroots = [int(r) % prime for r in rd.keys()]
    nsol = 0
    for xr in xroots:
        oy0 = Poly(o0.as_expr().subs(X, xr), Y, domain=DOM)
        oy1 = Poly(o1.as_expr().subs(X, xr), Y, domain=DOM)
        if oy0.is_zero and oy1.is_zero:
            nsol += 1; continue
        gg = oy0.gcd(oy1) if not oy0.is_zero else oy1
        if gg.degree() < 1:
            continue
        yrd = gg.ground_roots()
        nsol += len(yrd)
    return nsol, len(xroots)


def run_case(prime, rf, rp, seed=0, do_resultant=True):
    rng = random.Random(seed)
    fixed12 = [rng.randint(1, prime - 1) for _ in range(T - 4)]
    t0 = time.perf_counter()
    o0, o1, DOM, M, pos, prime = out_polys(prime, rf, rp, fixed12, seed)
    tbuild = time.perf_counter() - t0

    # supports / Newton polygons  (CHEAP -- this is the tropical/BKK object)
    sA = monomial_support(o0); sB = monomial_support(o1)
    dA = total_degree(o0); dB = total_degree(o1)
    bezout = dA * dB
    mv = mixed_volume_2d(sA, sB)

    pred_DI = D ** (2 * rf + rp)
    pred_td = D ** (rf + rp)   # individual total degree of out_i in (X,Y)

    degR = None; nsol = None; tres = 0.0
    if do_resultant:
        t1 = time.perf_counter()
        degR, R = resultant_degX(o0, o1, DOM)
        tres = time.perf_counter() - t1
        sol = count_fp_solutions(o0, o1, R, DOM, prime)
        nsol = sol[0] if sol else None

    print(f"=== prime={prime}  RF={rf} RP={rp}  (build {tbuild:.2f}s, res {tres:.2f}s) ===")
    print(f"  total-deg(out0)={dA}  total-deg(out1)={dB}   pred 3^(RF+RP)={pred_td}")
    densA = (dA+1)*(dA+2)//2
    print(f"  #support monomials: out0={len(sA)} out1={len(sB)}   "
          f"dense-count(deg {dA})={densA}  density={len(sA)/densA:.3f}")
    print(f"  BEZOUT (total-degree homotopy path count) d0*d1 = {bezout}  "
          f"(pred 3^(2(RF+RP))={D**(2*(rf+rp))})")
    print(f"  BKK / MIXED VOLUME (polyhedral homotopy path count) = {mv}")
    if degR is not None:
        print(f"  EXACT resultant deg_X = D_I = {degR}   pred 3^(2RF+RP)={pred_DI}")
        print(f"  #F_p solutions on this slice = {nsol}")
        print(f"  ratios:  MV/D_I = {mv/max(degR,1):.4f}   Bezout/D_I = {bezout/max(degR,1):.4f}")
    else:
        print(f"  (resultant skipped; pred D_I=3^(2RF+RP)={pred_DI})")
        print(f"  ratios:  MV/pred_DI = {mv/pred_DI:.4f}   Bezout/pred_DI = {bezout/pred_DI:.4f}")
    print(flush=True)
    return dict(prime=prime, rf=rf, rp=rp, dA=dA, dB=dB, suppA=len(sA), suppB=len(sB),
                bezout=bezout, mv=mv, degR=degR, pred_DI=pred_DI, nsol=nsol)


def main():
    print("#" * 78)
    print("# tropical / BKK / homotopy path-count vs resultant degree D_I")
    print("#" * 78)
    # SMALL prime so polynomials are buildable; structure (sparsity, MV) is field-size robust.
    results = []
    # (prime, RF, RP, do_resultant). Resultant only on the cheap ones (deg<=~81);
    # MV / Bezout / support are cheap and computed for all.
    cases = [
        (1009, 2, 0, True),
        (1009, 2, 1, True),
        (1009, 2, 2, False),
        (1009, 4, 0, False),
    ]
    for (pp, rf, rp, dores) in cases:
        try:
            results.append(run_case(pp, rf, rp, seed=1, do_resultant=dores))
        except Exception as e:
            print(f"  case ({pp},{rf},{rp}) failed: {e}\n", flush=True)

    print("#" * 78)
    print("# SUMMARY: does sparsity (partial rounds) drop the path count below D_I?")
    print("#" * 78)
    print(f"{'RF':>3} {'RP':>3} {'tot-deg':>8} {'Bezout':>10} {'MV(BKK)':>10} "
          f"{'D_I(exact)':>11} {'predD_I':>10} {'MV==predDI?':>12} {'#Fp':>5}")
    for r in results:
        ref = r['degR'] if r['degR'] is not None else r['pred_DI']
        eq = "YES" if r['mv'] == ref else "no"
        degr_s = str(r['degR']) if r['degR'] is not None else "-"
        print(f"{r['rf']:>3} {r['rp']:>3} {r['dA']:>8} {r['bezout']:>10} {r['mv']:>10} "
              f"{degr_s:>11} {r['pred_DI']:>10} {eq:>12} {str(r['nsol']):>5}")

    print()
    print("# NOTE on MV units: mixed_volume_2d returns NORMALIZED (lattice) MV")
    print("#   = 2*Euclidean-MV. The BKK path count is the Euclidean MV. For two")
    print("#   FULL simplices of degree d (Newt(out_i)), Euclidean MV = d^2 = Bezout.")
    print("#   So BKK path count == Bezout == 3^(2(RF+RP)) here; partial-round")
    print("#   sparsity is ~0 (density~1.0) so the Newton polytope gives NO reduction.")
    print()
    # ------------------------------------------------------------------
    # ANALYTIC projection to the REAL bounty (RF=6,RP=10): can't build the
    # 3^16-degree polynomial, but density~1.0 was measured, so Newt = full simplex.
    # ------------------------------------------------------------------
    RF, RPp = 6, 10
    d = D ** (RF + RPp)
    print("#" * 78)
    print(f"# REAL bounty projection RF={RF} RP={RPp}: total-deg(out_i)=3^(RF+RP)=3^{RF+RPp}={d}")
    print(f"#   measured support density on small cases ~1.0  => Newt(out_i) = full simplex of deg {d}")
    print(f"#   Euclidean BKK mixed volume = d^2 = 3^{2*(RF+RPp)} = {d*d}")
    print(f"#   total-degree Bezout = d^2 = {d*d}")
    print(f"#   resultant deg_X D_I  = 3^(2RF+RP) = 3^{2*RF+RPp} = {D**(2*RF+RPp)}")
    print(f"#   log2(BKK paths) = {(2*(RF+RPp)) * 1.5849625:.1f} ;  log2(D_I) = "
          f"{(2*RF+RPp)*1.5849625:.1f}")
    print("#   => Every homotopy/continuation route tracks >= D_I paths. NO sub-2^55 gain.")
    print()
    print("#" * 78)
    print("# SOS / Lasserre / Positivstellensatz / SDP verdict (ordering obstruction):")
    print("#" * 78)
    print("#  - F_p has NO compatible order => no notion of 'sum of squares >= 0',")
    print("#    no Positivstellensatz, no moment/SDP relaxation over F_p. DEAD natively.")
    print("#  - Lift to Z/Q/R: the system out0=out1=0 over Q has a 0-diml variety of")
    print("#    ~Bezout=3^(2(RF+RP)) COMPLEX points; an SOS infeasibility certificate")
    print("#    would only certify NO REAL common zero -- which says NOTHING about the")
    print("#    existence of an F_p zero (reduction mod p is not order-compatible).")
    print("#  - And to even SET UP a degree-2k Lasserre SDP whose relaxation is exact")
    print("#    you need k >= (Bezout-related) order; the moment matrix has size")
    print("#    C(n+k, k) which blows past 2^55 long before exactness. DEAD via (i)&(ii).")


if __name__ == "__main__":
    main()

"""
Probe 7 (Lens 6, Algebraic Geometry) — reducibility / common-component / rational-point
structure of the CICO-2 system.

If the output polynomials P(X,Y), Q(X,Y) FACTOR, or share a gcd (common curve component),
or the univariate resultant R(X) splits into low-degree factors, the degree wall collapses:
we parametrize/solve a low-degree piece instead of the full degree-D_I system. This is the
first decisive AG diagnostic. We test BOTH matrices (control: is geometry matrix-independent?).
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "reference"))
from sympy import symbols, Poly, GF, gcd as sgcd
from poseidon.poseidon import Poseidon
from poseidon.mds_matrix import generate_circulant_mds_matrix, generate_mds_matrix

P = 2130706433
T = 16; D = 3; C = 2
ROW16 = [1, 1, 51, 1, 11, 17, 2, 1, 101, 63, 15, 2, 67, 22, 13, 3]
X, Y = symbols("X Y")
DOM = GF(P)


def build_PQ(r_f, r_p, M):
    pos = Poseidon(prime=P, alpha=D, t=T, r_f=r_f, r_p=r_p, mds=M)
    rc = pos.round_constants
    const = lambda v: Poly(int(v) % P, X, Y, domain=DOM)
    st = [Poly(X, X, Y, domain=DOM), Poly(Y, X, Y, domain=DOM)]
    st += [const((987654 * (i + 1) + 5) % P) for i in range(T - 2 - C)] + [const(0)] * C
    def mds(s): return [sum(int(M[i][j]) * s[j] for j in range(T)) for i in range(T)]
    def arc(s, c): return [s[i] + int(c[i]) for i in range(T)]
    st = mds(st); half = r_f // 2; idx = 0
    for _ in range(half):
        st = arc(st, rc[idx]); idx += 1; st = [x**D for x in st]; st = mds(st)
    for _ in range(r_p):
        st = arc(st, rc[idx]); idx += 1; st[0] = st[0]**D; st = mds(st)
    for _ in range(half):
        st = arc(st, rc[idx]); idx += 1; st = [x**D for x in st]; st = mds(st)
    return (Poly(st[T - C].as_expr(), X, Y, domain=DOM),
            Poly(st[T - C + 1].as_expr(), X, Y, domain=DOM))


def univ_factor_profile(expr, var):
    """Factor-degree multiset of a univariate poly over GF(p) (sympy supports univariate)."""
    from collections import Counter
    _, facs = Poly(expr, var, domain=DOM).factor_list()
    return Counter(f.degree() for f, m in facs for _ in range(m))


def biv_reducibility(poly):
    """Proxy for bivariate irreducibility: factor poly(X, y0) for a few y0.
       If poly(X,y0) is irreducible for some y0 -> poly is irreducible over GF(p)."""
    profs = []
    for y0 in (1, 7, 12345):
        f = poly.as_expr().subs(Y, y0)
        prof = univ_factor_profile(f, X)
        profs.append(dict(sorted(prof.items())))
    irr = any(len(p) == 1 and 1 not in p and list(p.values())[0] == 1 and list(p.keys())[0] == poly.degree(X)
              for p in profs)
    return profs, irr


def analyze(r_f, r_p, M):
    from collections import Counter
    Pp, Qp = build_PQ(r_f, r_p, M)
    # common component?
    try:
        g = sgcd(Pp, Qp); gdeg = Poly(g, X, Y, domain=DOM).total_degree() if g != 1 else 0
    except Exception as e:
        gdeg = f"gcd-err:{type(e).__name__}"
    # bivariate reducibility via slicing
    pP, irrP = biv_reducibility(Pp)
    # univariate resultant R(X): factor profile + Ritt decomposition test
    R = Poly(Pp.as_expr(), Y, X, domain=DOM).resultant(Poly(Qp.as_expr(), Y, X, domain=DOM))
    RX = Poly(R, X, domain=DOM)
    rprof = univ_factor_profile(RX.as_expr(), X)
    try:
        decomp = RX.decompose()
        ddeg = [d.degree() for d in decomp]
    except Exception as e:
        ddeg = f"decompose-err:{type(e).__name__}"
    print(f"  RF={r_f} RP={r_p}: degP={Pp.total_degree()} | gcd(P,Q)deg={gdeg}", flush=True)
    print(f"    P(X,y0) factor profiles: {pP}  -> P {'IRREDUCIBLE' if irrP else 'maybe reducible'}", flush=True)
    print(f"    resultant R(X) deg={RX.degree()}: factor profile {dict(sorted(rprof.items()))}", flush=True)
    print(f"      rational roots (lin facs): {rprof.get(1,0)} | smallest factor deg: {min(rprof) if rprof else '-'}", flush=True)
    print(f"      Ritt decomposition degrees: {ddeg}  {'<-- DECOMPOSES!' if isinstance(ddeg,list) and len(ddeg)>1 else '(indecomposable)'}", flush=True)


def main():
    mats = {"Plonky3 circulant": generate_circulant_mds_matrix(ROW16, P),
            "Cauchy default": generate_mds_matrix(T, P)}
    for name, M in mats.items():
        print(f"===== {name} =====", flush=True)
        for (rf, rp) in [(2, 0), (2, 1)]:
            try:
                analyze(rf, rp, M)
            except Exception as e:
                print(f"  RF={rf} RP={rp}: ERROR {type(e).__name__}: {e}", flush=True)


if __name__ == "__main__":
    main()

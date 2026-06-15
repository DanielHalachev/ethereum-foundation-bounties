"""
Probe 6 (Lens 2) — structure of the CICO-2 multiplication matrix T_X.

The surviving constructive lead: a Krylov/Wiedemann change-of-order uses O(D_I) memory
(breaks the resultant's delta^2 / PB wall) but O(D_I^2) time UNLESS T_X has low
*displacement rank* in a cheaply-available basis -> then a superfast structured solve
gives ~Õ(D_I) time = personal-hardware feasible.

We build T_X (multiply-by-X on the quotient ring R/I) for tiny CICO-2 ideals and measure:
  - sparsity of T_X (nonzeros) -> Wiedemann matvec cost / memory
  - displacement rank  rank(Z*T_X - T_X*Z)  for the down-shift Z in the basis ordering,
    and for the degree-grading diagonal -> is T_X Toeplitz/Hankel-like (low) or generic (~D_I)?
Low displacement rank in this cheap basis = the superfast dream is alive. Full = it's dead
(structure, if any, lives only in shape position, which is itself the expensive step).
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "reference"))
from sympy import symbols, Poly, GF, groebner, ZZ
from poseidon.poseidon import Poseidon
from poseidon.mds_matrix import generate_circulant_mds_matrix, generate_mds_matrix

P = 2130706433
T = 16
D = 3
C = 2
ROW16 = [1, 1, 51, 1, 11, 17, 2, 1, 101, 63, 15, 2, 67, 22, 13, 3]
X, Y = symbols("X Y")
DOM = GF(P)


def build_PQ(r_f, r_p, M):
    """CICO-2 bivariate system (no skip): input (X,Y,consts...,0,0), out last 2 coords."""
    pos = Poseidon(prime=P, alpha=D, t=T, r_f=r_f, r_p=r_p, mds=M)
    rc = pos.round_constants
    const = lambda v: Poly(int(v) % P, X, Y, domain=DOM)
    st = [Poly(X, X, Y, domain=DOM), Poly(Y, X, Y, domain=DOM)]
    st += [const((987654 * (i + 1) + 5) % P) for i in range(T - 2 - C)] + [const(0)] * C

    def mds(s): return [sum(int(M[i][j]) * s[j] for j in range(T)) for i in range(T)]
    def arc(s, c): return [s[i] + int(c[i]) for i in range(T)]
    st = mds(st)
    half = r_f // 2; idx = 0
    for _ in range(half):
        st = arc(st, rc[idx]); idx += 1; st = [x**D for x in st]; st = mds(st)
    for _ in range(r_p):
        st = arc(st, rc[idx]); idx += 1; st[0] = st[0]**D; st = mds(st)
    for _ in range(half):
        st = arc(st, rc[idx]); idx += 1; st = [x**D for x in st]; st = mds(st)
    return st[T - C].as_expr(), st[T - C + 1].as_expr()


def rank_modp(M, p):
    M = [[x % p for x in row] for row in M]; rows = len(M); cols = len(M[0]) if rows else 0
    rk = 0; pr = 0
    for c in range(cols):
        piv = next((r for r in range(pr, rows) if M[r][c] % p), None)
        if piv is None: continue
        M[pr], M[piv] = M[piv], M[pr]
        inv = pow(M[pr][c], -1, p); M[pr] = [(x * inv) % p for x in M[pr]]
        for r in range(rows):
            if r != pr and M[r][c] % p:
                f = M[r][c]; M[r] = [(M[r][k] - f * M[pr][k]) % p for k in range(cols)]
        pr += 1; rk += 1
    return rk


def analyze(r_f, r_p, M):
    Pe, Qe = build_PQ(r_f, r_p, M)
    G = groebner([Pe, Qe], X, Y, order="grevlex", domain=DOM)
    # leading exponent pair (a,b) of each basis poly under grevlex
    lead = [g.monoms(order="grevlex")[0] for g in G.polys]

    def in_LTideal(a, b):
        return any(a >= la and b >= lb for (la, lb) in lead)
    amax = max(la for la, _ in lead) + max(lb for _, lb in lead) + 2
    basis = [(a, b) for a in range(amax) for b in range(amax) if not in_LTideal(a, b)]
    Dq = len(basis)
    index = {m: i for i, m in enumerate(basis)}

    def reduce_to_col(expr):
        res = G.reduce(expr)
        r = res[1] if isinstance(res, (list, tuple)) else res
        rp = Poly(r, X, Y, domain=DOM)
        col = [0] * Dq
        for monom, coeff in rp.terms():
            col[index[tuple(monom)]] = int(coeff) % P
        return col

    TX = [[0] * Dq for _ in range(Dq)]
    for (a, b), j in index.items():
        col = reduce_to_col((X ** (a + 1)) * (Y ** b))
        for i in range(Dq):
            TX[i][j] = col[i]
    nz = sum(1 for i in range(Dq) for j in range(Dq) if TX[i][j] % P)

    # displacement w.r.t. down-shift Z (subdiagonal ones) in the basis ordering
    Z = [[1 if i == j + 1 else 0 for j in range(Dq)] for i in range(Dq)]
    def matmul(A, B): return [[sum(A[i][k] * B[k][j] for k in range(Dq)) % P for j in range(Dq)] for i in range(Dq)]
    ZT = matmul(Z, TX); TZ = matmul(TX, Z)
    disp = [[(ZT[i][j] - TZ[i][j]) % P for j in range(Dq)] for i in range(Dq)]
    dr = rank_modp(disp, P)
    print(f"RF={r_f} RP={r_p}: D_I={Dq}  nnz(T_X)={nz} ({100*nz/Dq/Dq:.0f}% dense)  "
          f"displacement_rank(down-shift)={dr}  (Toeplitz-like<=2, generic~D_I={Dq})", flush=True)


def main():
    mats = {
        "Plonky3 circulant": generate_circulant_mds_matrix(ROW16, P),
        "Cauchy default":    generate_mds_matrix(T, P),
    }
    for name, M in mats.items():
        print(f"===== {name} =====", flush=True)
        for (rf, rp) in [(2, 0), (2, 1)]:
            try:
                analyze(rf, rp, M)
            except Exception as e:
                print(f"RF={rf} RP={rp}: ERROR {type(e).__name__}: {e}", flush=True)


if __name__ == "__main__":
    main()

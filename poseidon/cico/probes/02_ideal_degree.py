"""
Probe 2 — measure the CICO-2 ideal degree D_I empirically.

Builds the bivariate system {P(X,Y)=0, Q(X,Y)=0} for permutation_plus_linear over
GF(p) symbolically (sympy), then D_I := deg_X Res_Y(P,Q). Checks:
  - the Poseidon-specific degree drop: D_I should equal d^(2*R_F + R_P),
    NOT the trivial Bezout d^(2*(R_F+R_P)).  (This factor-d^R_P drop is what makes
    the resultant attack cheap; 2026/150 sec 4.3.)
  - matrix independence: same D_I for the Plonky3 circulant and the Cauchy default.

Run with the project venv:  ./.venv/bin/python probes/02_ideal_degree.py
"""
import sys, os, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "reference"))

from sympy import symbols, Poly, GF, resultant
from poseidon.poseidon import Poseidon
from poseidon.mds_matrix import generate_circulant_mds_matrix, generate_mds_matrix

P = 2130706433
T = 16
D = 3
ROW16 = [1, 1, 51, 1, 11, 17, 2, 1, 101, 63, 15, 2, 67, 22, 13, 3]
X, Y = symbols("X Y")
DOM = GF(P)   # coefficients reduced mod p


def sym_perm_plus_linear(pos, state):
    """permutation_plus_linear evaluated on a symbolic state (list of sympy Polys)."""
    p, t, mds, rc = pos.prime, pos.t, pos.mds, pos.round_constants

    def mds_mul(s):
        return [sum(int(mds[i][j]) * s[j] for j in range(t)) for i in range(t)]

    def add_rc(s, c):
        return [s[i] + int(c[i]) for i in range(t)]

    def sbox(x):
        return x ** D

    state = mds_mul(state)                      # initial linear layer ("+")
    half = pos.r_f // 2
    idx = 0
    for _ in range(half):                       # front full rounds
        state = add_rc(state, rc[idx]); idx += 1
        state = [sbox(x) for x in state]
        state = mds_mul(state)
    for _ in range(pos.r_p):                    # partial rounds
        state = add_rc(state, rc[idx]); idx += 1
        state[0] = sbox(state[0])
        state = mds_mul(state)
    for _ in range(half):                       # back full rounds
        state = add_rc(state, rc[idx]); idx += 1
        state = [sbox(x) for x in state]
        state = mds_mul(state)
    return state


def measure_DI(r_f, r_p, matrix):
    pos = Poseidon(prime=P, alpha=D, t=T, r_f=r_f, r_p=r_p, mds=matrix)
    # CICO-2 input: words 0,1 pinned (C1,C2); two free words become X,Y; rest fixed.
    C1, C2 = 0xC09DE4, 0xEE6282
    zero = Poly(0, X, Y, domain=DOM)
    def const(v): return Poly(int(v) % P, X, Y, domain=DOM)
    state = [const(C1), const(C2),
             Poly(X, X, Y, domain=DOM), Poly(Y, X, Y, domain=DOM)]
    # remaining 12 free words fixed to fixed pseudo-random constants (deterministic)
    fillers = [ (1234567 * (i + 1) + 89) % P for i in range(T - 4) ]
    state += [const(v) for v in fillers]

    tb = time.perf_counter()
    out = sym_perm_plus_linear(pos, state)
    Pp = out[0].as_expr()
    Qp = out[1].as_expr()
    degP = Poly(Pp, X, Y, domain=DOM).total_degree()
    build = time.perf_counter() - tb
    t0 = time.perf_counter()
    # Y as main generator -> Poly.resultant eliminates Y, returns a Poly in X over GF(p)
    f = Poly(Pp, Y, X, domain=DOM)
    g = Poly(Qp, Y, X, domain=DOM)
    R = f.resultant(g)
    DI = R.degree()        # degree in X
    dt = time.perf_counter() - t0
    return degP, DI, build, dt


def main():
    import sys as _sys
    cases = [(2, 0), (2, 1), (2, 2)]
    mats = {
        "Plonky3 circulant": generate_circulant_mds_matrix(ROW16, P),
        "Cauchy default":    generate_mds_matrix(T, P),
    }
    print(f"{'matrix':18} {'RF':>2} {'RP':>2} {'degP':>5} {'measured D_I':>13} "
          f"{'d^(2RF+RP)':>11} {'trivial':>9} {'build_s':>8} {'res_s':>7}  flag", flush=True)
    for (rf, rp) in cases:
        for name, M in mats.items():
            degP, DI, build, dt = measure_DI(rf, rp, M)
            pred = D ** (2 * rf + rp)
            triv = D ** (2 * (rf + rp))
            flag = "DROP-OK" if DI == pred else ("trivial" if DI == triv else "??")
            print(f"{name:18} {rf:>2} {rp:>2} {degP:>5} {DI:>13} {pred:>11} {triv:>9} "
                  f"{build:>8.1f} {dt:>7.1f}  {flag}", flush=True)
            _sys.stdout.flush()


if __name__ == "__main__":
    main()

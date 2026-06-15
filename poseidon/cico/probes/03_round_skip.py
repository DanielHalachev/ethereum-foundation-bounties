"""
Probe 3 — implement the BBLP22 / 2026-150 sec4.2 first-round skip and measure its
effect on the CICO-2 ideal degree D_I.

Convention (matching the paper): capacity = LAST c coords. CICO-2: c=k=2.
Permutation = permutation_plus_linear: M_init, then [r_f/2 full][r_p partial][r_f/2 full].
The skip absorbs M_init + round-1 (ARC,S-box,M) into an affine-in-(X,Y) state.

Construction (sec 4.2, Fig.2), c=2, d=3, t=16, Minv = M^{-1}:
  a_k (k=1,2): kernel of the c x (c+1) submatrix  Minv[rows t-c..t-1, cols (k-1)(c+1) .. (k-1)(c+1)+c].
  b (len c):   solve  Minv[t-c.., t-c..] . b = (Minv . C^(1))[t-c..]
  S-box input s2 = (a_{1,*} X1^{1/d}, a_{2,*} X2^{1/d}, 0...0 (extra), b_1, b_2)
  => post-S-box (cube) = (a_{1,*}^3 X1, a_{2,*}^3 X2, 0...0, b_1^3, b_2^3)   [affine in X1,X2]
  => state before round 2 = M . post-S-box.

Checks:
  (1) round-trip: recovered input = Minv.(s2 - C^(1)) must have its last c coords = 0
      for random X1,X2  -> construction is a valid CICO input.
  (2) D_I(skip) should equal d^(2(R_F-1)+R_P), vs D_I(no skip) = d^(2 R_F + R_P).

Run:  ./.venv/bin/python probes/03_round_skip.py
"""
import sys, os, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "reference"))

from sympy import symbols, Poly, GF, Matrix
from poseidon.poseidon import Poseidon
from poseidon.mds_matrix import generate_circulant_mds_matrix, _right_null_space

P = 2130706433
T = 16
D = 3
C = 2                       # capacity / number of CICO vars (CICO-2)
ROW16 = [1, 1, 51, 1, 11, 17, 2, 1, 101, 63, 15, 2, 67, 22, 13, 3]
X, Y = symbols("X Y")
DOM = GF(P)


def matinv(M):
    return [[int(x) % P for x in row] for row in Matrix(M).inv_mod(P).tolist()]


def matvec(M, v):
    return [sum(int(M[i][j]) * int(v[j]) for j in range(len(v))) % P for i in range(len(M))]


def build_skip_params(M, rc1):
    """Return (a_list, b) for the first-round skip. a_list[k] is length c+1; b is length c."""
    Minv = matinv(M)
    # a_k = kernel of Minv[rows t-c..t-1, cols (k-1)(c+1).. (k-1)(c+1)+c]   (c x (c+1))
    a_list = []
    for k in range(C):
        cols = [(k) * (C + 1) + j for j in range(C + 1)]          # block of c+1 columns
        sub = [[Minv[T - C + r][cols[j]] for j in range(C + 1)] for r in range(C)]
        ker = _right_null_space([row[:] for row in sub], C + 1, P)
        assert ker, f"no kernel for variable {k} (matrix not MDS here?)"
        a_list.append([x % P for x in ker[0]])
    # b: solve  Minv[t-c.., t-c..] . b = (Minv . C^(1))[t-c..]
    rhs_full = matvec(Minv, [int(x) % P for x in rc1])
    A = Matrix([[Minv[T - C + r][T - C + j] for j in range(C)] for r in range(C)])
    rhsv = Matrix([rhs_full[T - C + r] for r in range(C)])
    b = list((A.inv_mod(P) * rhsv) % P)
    b = [int(x) % P for x in b]
    return a_list, b


def roundtrip_ok(M, rc1, a_list, b):
    """Pick random-ish X1,X2, build s2, recover input, check last c coords are 0."""
    Minv = matinv(M)
    e = pow(D, -1, P - 1)                 # X^{1/d} = X^e
    ok = True
    for (x1, x2) in [(7, 11), (123456, 7654321), (1, 2)]:
        x1e, x2e = pow(x1, e, P), pow(x2, e, P)
        s2 = [a_list[0][0] * x1e % P, a_list[0][1] * x1e % P, a_list[0][2] * x1e % P,
              a_list[1][0] * x2e % P, a_list[1][1] * x2e % P, a_list[1][2] * x2e % P]
        s2 += [0] * (T - 2 * (C + 1) - C)         # extra branches = 0
        s2 += [b[0], b[1]]                         # capacity branches
        s1 = [(s2[i] - int(rc1[i])) % P for i in range(T)]
        inp = matvec(Minv, s1)
        if any(inp[T - C + r] != 0 for r in range(C)):
            ok = False
    return ok


# ---------- symbolic evaluation ----------
def sbox(x): return x ** D

def mds_mul_sym(s, M):
    return [sum(int(M[i][j]) * s[j] for j in range(T)) for i in range(T)]

def add_rc_sym(s, c):
    return [s[i] + int(c[i]) for i in range(T)]

def full_round(s, M, c):
    s = add_rc_sym(s, c); s = [sbox(x) for x in s]; return mds_mul_sym(s, M)

def partial_round(s, M, c):
    s = add_rc_sym(s, c); s[0] = sbox(s[0]); return mds_mul_sym(s, M)


def DI_noskip(pos, M):
    """Baseline: input=(X,Y,consts...,0,0); run full permutation_plus_linear; D_I from out[-2],out[-1]."""
    rc = pos.round_constants
    const = lambda v: Poly(int(v) % P, X, Y, domain=DOM)
    state = [Poly(X, X, Y, domain=DOM), Poly(Y, X, Y, domain=DOM)]
    state += [const((987654 * (i + 1) + 5) % P) for i in range(T - 2 - C)]
    state += [const(0)] * C                          # capacity input = 0 (last c)
    state = mds_mul_sym(state, M)                    # initial linear layer
    half = pos.r_f // 2; idx = 0
    for _ in range(half): state = full_round(state, M, rc[idx]); idx += 1
    for _ in range(pos.r_p): state = partial_round(state, M, rc[idx]); idx += 1
    for _ in range(half): state = full_round(state, M, rc[idx]); idx += 1
    return resultant_deg(state[T - C].as_expr(), state[T - C + 1].as_expr())


def DI_skip(pos, M, a_list, b):
    """Skipped: post-S-box state affine in X,Y; close round1 MDS; run remaining R-1 rounds."""
    rc = pos.round_constants
    A = [[int(pow(a_list[k][j], D, P)) for j in range(C + 1)] for k in range(C)]  # a_{k,j}^d
    Bc = [int(pow(bk, D, P)) for bk in b]                                          # b_k^d
    post = [A[0][0] * Poly(X, X, Y, domain=DOM), A[0][1] * Poly(X, X, Y, domain=DOM), A[0][2] * Poly(X, X, Y, domain=DOM),
            A[1][0] * Poly(Y, X, Y, domain=DOM), A[1][1] * Poly(Y, X, Y, domain=DOM), A[1][2] * Poly(Y, X, Y, domain=DOM)]
    post += [Poly(0, X, Y, domain=DOM)] * (T - 2 * (C + 1) - C)
    post += [Poly(int(Bc[0]), X, Y, domain=DOM), Poly(int(Bc[1]), X, Y, domain=DOM)]
    state = mds_mul_sym(post, M)                     # close round-1 MDS -> state before round 2
    half = pos.r_f // 2; idx = 1                     # rc[0] absorbed into b
    for _ in range(half - 1): state = full_round(state, M, rc[idx]); idx += 1
    for _ in range(pos.r_p): state = partial_round(state, M, rc[idx]); idx += 1
    for _ in range(half): state = full_round(state, M, rc[idx]); idx += 1
    return resultant_deg(state[T - C].as_expr(), state[T - C + 1].as_expr())


def resultant_deg(Pe, Qe):
    f = Poly(Pe, Y, X, domain=DOM); g = Poly(Qe, Y, X, domain=DOM)
    return f.resultant(g).degree()


def main():
    M = generate_circulant_mds_matrix(ROW16, P)
    print("Plonky3 circulant, CICO-2, d=3, t=16\n")
    for (rf, rp) in [(2, 0), (2, 1), (4, 0)]:
        pos = Poseidon(prime=P, alpha=D, t=T, r_f=rf, r_p=rp, mds=M)
        a_list, b = build_skip_params(M, pos.round_constants[0])
        rt = roundtrip_ok(M, pos.round_constants[0], a_list, b)
        t0 = time.perf_counter(); di_ns = DI_noskip(pos, M); t1 = time.perf_counter()
        di_sk = DI_skip(pos, M, a_list, b); t2 = time.perf_counter()
        pred_ns = D ** (2 * rf + rp)
        pred_sk = D ** (2 * (rf - 1) + rp)
        print(f"RF={rf} RP={rp} | roundtrip(input cap==0): {rt}")
        print(f"   no-skip D_I = {di_ns:>6}  (pred d^(2RF+RP)     = {pred_ns:>6})  "
              f"{'OK' if di_ns==pred_ns else 'MISMATCH'}  [{t1-t0:.1f}s]")
        print(f"   skip-1  D_I = {di_sk:>6}  (pred d^(2(RF-1)+RP) = {pred_sk:>6})  "
              f"{'OK' if di_sk==pred_sk else 'MISMATCH'}  [{t2-t1:.1f}s]")
        print(f"   --> skip divides D_I by {di_ns/di_sk:.1f} (expect d^2={D**2})\n", flush=True)


if __name__ == "__main__":
    main()

"""
Probe 10 — INDEPENDENT ground-truth re-verification (do not trust prior probes).

Checks, all computed fresh here:
  (A) Which matrix does the *verifier* actually use? (cico_verifier default = Cauchy.)
  (B) D_I = deg_X Res_Y(P,Q) on tiny instances, for BOTH Cauchy and circulant.
      Compare to d^(2RF+RP) (Poseidon drop) vs d^(2(RF+RP)) (trivial Bezout).
  (C) Round-constant anomaly scan for the REAL RF=6,RP=10 instance:
      zeros, duplicates, small values, partial-round-constant structure.
  (D) Cauchy-matrix structure: Krylov dim of e0,e1 (subspace-trail capacity);
      any invariant subspace aligned with output coords {0,1}.
  (E) Inverse permutation sanity (P+ is a bijection; we can invert it).

Run: ./.venv/bin/python probes/10_ground_truth.py
"""
import sys, os, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "reference"))
from sympy import symbols, Poly, GF, Matrix
from poseidon.poseidon import Poseidon
from poseidon.mds_matrix import generate_mds_matrix, generate_circulant_mds_matrix

P = 2130706433
T = 16
D = 3
ROW16 = [1, 1, 51, 1, 11, 17, 2, 1, 101, 63, 15, 2, 67, 22, 13, 3]
X, Y = symbols("X Y")
DOM = GF(P)


def sym_perm(pos, M, state):
    def mds(s): return [sum(int(M[i][j]) * s[j] for j in range(T)) for i in range(T)]
    def arc(s, c): return [s[i] + int(c[i]) for i in range(T)]
    rc = pos.round_constants
    state = mds(state)
    half = pos.r_f // 2; idx = 0
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


def measure_DI(rf, rp, M):
    pos = Poseidon(prime=P, alpha=D, t=T, r_f=rf, r_p=rp, mds=M)
    C1, C2 = 0xC09DE4, 0xEE6282
    def c(v): return Poly(int(v) % P, X, Y, domain=DOM)
    state = [c(C1), c(C2), Poly(X, X, Y, domain=DOM), Poly(Y, X, Y, domain=DOM)]
    state += [c((1234567 * (i + 1) + 89) % P) for i in range(T - 4)]
    out = sym_perm(pos, M, state)
    f = Poly(out[0].as_expr(), Y, X, domain=DOM)
    g = Poly(out[1].as_expr(), Y, X, domain=DOM)
    R = f.resultant(g)
    return R.degree()


def part_B():
    print("=" * 70)
    print("(B) D_I on tiny instances  [pred drop = d^(2RF+RP), trivial = d^(2(RF+RP))]")
    mats = {"Cauchy(verifier)": generate_mds_matrix(T, P),
            "circulant":        generate_circulant_mds_matrix(ROW16, P)}
    for (rf, rp) in [(2, 0), (2, 1), (2, 2), (4, 0)]:
        for name, M in mats.items():
            t0 = time.perf_counter()
            di = measure_DI(rf, rp, M)
            dt = time.perf_counter() - t0
            drop = D ** (2 * rf + rp); triv = D ** (2 * (rf + rp))
            flag = "DROP" if di == drop else ("trivial" if di == triv else "??")
            print(f"  {name:18} RF={rf} RP={rp}: D_I={di:>6}  drop={drop:>6} triv={triv:>7}  {flag}  [{dt:.1f}s]")


def part_C():
    print("=" * 70)
    print("(C) Round-constant anomaly scan, REAL instance RF=6, RP=10")
    pos = Poseidon(prime=P, alpha=D, t=T, r_f=6, r_p=10)
    rc = pos.round_constants
    flat = [x for row in rc for x in row]
    print(f"  total constants = {len(flat)} (= (6+10)*16 = 256)")
    print(f"  zeros: {sum(1 for x in flat if x == 0)}   duplicates: {len(flat) - len(set(flat))}")
    print(f"  min={min(flat)}  max={max(flat)}  (<2^16: {sum(1 for x in flat if x < (1<<16))})")
    # partial-round word-0 constants (the only S-box that fires in partial rounds)
    part0 = [rc[3 + i][0] for i in range(10)]
    print(f"  partial-round word-0 constants (10): {part0}")
    print(f"    any zero among partial word-0 ARC: {any(v == 0 for v in part0)}")


def part_D():
    print("=" * 70)
    print("(D) Cauchy matrix structure: Krylov dim of e0,e1; invariant alignment")
    M = Matrix(generate_mds_matrix(T, P))
    for name, e in [("e0", [1] + [0]*15), ("e1", [0,1] + [0]*14)]:
        v = Matrix(e); basis = [v]
        cur = v
        for _ in range(T + 2):
            cur = (M * cur).applyfunc(lambda z: z % P)
            stacked = Matrix.hstack(*basis, cur)
            if stacked.rank(iszerofunc=lambda z: z % P == 0) > len(basis):
                basis.append(cur)
            else:
                break
        print(f"  Krylov dim of {name} under Cauchy M: {len(basis)} (full = {T})")
    # smallest M-invariant subspace containing both e0,e1
    basis = [Matrix([1]+[0]*15), Matrix([0,1]+[0]*14)]
    def rank(bs): return Matrix.hstack(*bs).rank(iszerofunc=lambda z: z % P == 0)
    changed = True
    while changed and len(basis) < T:
        changed = False
        for b in list(basis):
            nb = (M * b).applyfunc(lambda z: z % P)
            if rank(basis + [nb]) > len(basis):
                basis.append(nb); changed = True
    print(f"  smallest M-invariant subspace containing e0,e1: dim {len(basis)} (full={T})")


def part_E():
    print("=" * 70)
    print("(E) Inverse-permutation sanity (P+ bijective, invertible)")
    M = generate_mds_matrix(T, P)
    pos = Poseidon(prime=P, alpha=D, t=T, r_f=6, r_p=10, mds=M)
    Minv = [[int(x) % P for x in r] for r in Matrix(M).inv_mod(P).tolist()]
    einv = pow(D, -1, P - 1)
    def mv(Mx, v): return [sum(Mx[i][j]*v[j] for j in range(T)) % P for i in range(T)]
    def inv_full(s, c):
        s = mv(Minv, s); s = [pow(x, einv, P) for x in s]
        return [(s[i]-int(c[i])) % P for i in range(T)]
    def inv_part(s, c):
        s = mv(Minv, s); s[0] = pow(s[0], einv, P)
        return [(s[i]-int(c[i])) % P for i in range(T)]
    s_in = [0xC09DE4, 0xEE6282] + [(7*i+3) % P for i in range(14)]
    out = pos.permutation_plus_linear(s_in)
    rc = pos.round_constants; s = list(out)
    for i in range(2, -1, -1): s = inv_full(s, rc[13+i]) if False else s  # placeholder
    # invert: 3 back full, 10 partial, 3 front full, then initial M
    idx = 15
    for _ in range(3): s = inv_full(s, rc[idx]); idx -= 1
    for _ in range(10): s = inv_part(s, rc[idx]); idx -= 1
    for _ in range(3): s = inv_full(s, rc[idx]); idx -= 1
    s = mv(Minv, s)
    print(f"  inverse round-trip exact: {s == [x % P for x in s_in]}")


if __name__ == "__main__":
    part_C(); part_D(); part_E(); part_B()

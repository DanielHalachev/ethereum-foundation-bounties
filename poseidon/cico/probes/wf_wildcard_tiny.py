"""
TINY end-to-end demonstration of the pure-wildcard CICO-2 route, on a small
field where the ~p search actually terminates, plus a head-to-head cost count
against brute force on the SAME instance.

Instance: p=193, t=8, alpha=3, RF=2, RP=2, Cauchy MDS, Grain RCs.
CICO-2: input [C1,C2,x2..x7], require out0==0 and out1==0.

We compare three solvers on the identical instance and COUNT work:
  (A) pure-wildcard: for each randomization of x3..x7, build out0(X) (X=x2),
      root-find out0(X)=0 over F_p, test out1 at each root.  WORK = #perm-equiv
      evaluations (we count root-finds and out1-tests).
  (B) brute: sweep all (x2,...,x7) until out0==out1==0.  WORK = #perm evals.
  (C) bivariate resultant: fix x4..x7, keep x2=X,x3=Z, Res_Z(out0,out1)=R(X),
      root-find, back-substitute -- the known cheap attack.

All three call verify_cico_solution at the end to confirm a real solution.
"""
import sys, os, time, random
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "reference"))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "reference", "bounties"))

from sympy import symbols, Poly, GF, gcd as sgcd, resultant
from poseidon.poseidon import Poseidon
from poseidon.mds_matrix import generate_mds_matrix

random.seed(2026)
p = 193
T = 8
D = 3
C1, C2 = 5, 7        # tiny "pinned" inputs
TARGET = (0, 0)      # out0, out1 must be 0
X = symbols("X")
DOM = GF(p)
M = generate_mds_matrix(T, p)


def make_pos(rf, rp):
    return Poseidon(prime=p, alpha=D, t=T, r_f=rf, r_p=rp, mds=M)


def perm(pos, state):
    return pos.permutation_plus_linear(state)


# ---- symbolic out0(X),out1(X) with x2=X, words3..7 fixed ----
def sym_out01(pos, fixed5):
    mds, rc, t = pos.mds, pos.round_constants, pos.t
    def mm(s): return [sum(int(mds[i][j]) * s[j] for j in range(t)) for i in range(t)]
    def arc(s, c): return [s[i] + int(c[i]) for i in range(t)]
    def sb(x): return x ** D
    state = [Poly(C1, X, domain=DOM), Poly(C2, X, domain=DOM), Poly(X, X, domain=DOM)]
    state += [Poly(int(v) % p, X, domain=DOM) for v in fixed5]
    state = state[:t]
    state = mm(state)
    half = pos.r_f // 2; idx = 0
    for _ in range(half):
        state = arc(state, rc[idx]); idx += 1; state = [sb(x) for x in state]; state = mm(state)
    for _ in range(pos.r_p):
        state = arc(state, rc[idx]); idx += 1; state[0] = sb(state[0]); state = mm(state)
    for _ in range(half):
        state = arc(state, rc[idx]); idx += 1; state = [sb(x) for x in state]; state = mm(state)
    return state[0].as_expr(), state[1].as_expr()


def fp_roots(expr, var=X):
    f = Poly(expr, var, domain=DOM)
    if f.degree() <= 0:
        return []
    b = Poly(var, var, domain=DOM); acc = Poly(1, var, domain=DOM); e = p
    while e > 0:
        if e & 1: acc = (acc * b) % f
        e >>= 1
        if e: b = (b * b) % f
    g = sgcd((acc - Poly(var, var, domain=DOM)) % f, f)
    if g.degree() <= 0:
        return []
    rs = []
    for (fct, m) in g.factor_list()[1]:
        if fct.degree() == 1:
            a = int(fct.nth(1)); b0 = int(fct.nth(0))
            rs.append((-b0 * pow(a, -1, p)) % p)
    return rs


def solver_wildcard(pos):
    from cico_verifier import verify_cico_solution
    work_rootfind = 0; work_test = 0; fixings = 0
    t0 = time.perf_counter()
    while time.perf_counter() - t0 < 60:
        fixings += 1
        fixed5 = [random.randrange(p) for _ in range(5)]  # x3..x7
        o0, o1 = sym_out01(pos, fixed5)
        work_rootfind += 1
        roots = fp_roots(o0)
        for r in roots:
            work_test += 1
            v1 = int(Poly(o1, X, domain=DOM).eval(r)) % p
            if v1 == TARGET[1]:
                free = [r] + fixed5
                ok = verify_cico_solution(free, prime=p, k=2, t=T, r_f=pos.r_f,
                                          r_p=pos.r_p, mds=M, constants=[C1, C2, 0, 0])
                return dict(found=True, ok=ok, fixings=fixings,
                            rootfinds=work_rootfind, tests=work_test, sol=(r, fixed5),
                            time=time.perf_counter() - t0)
    return dict(found=False, fixings=fixings, rootfinds=work_rootfind,
                tests=work_test, time=time.perf_counter() - t0)


def solver_brute(pos):
    from cico_verifier import verify_cico_solution
    evals = 0; t0 = time.perf_counter()
    # brute sweep over (x2..x7): full space p^6 too big; but we only need to match
    # 2 outputs -> expected p^2 evals.  Random search counting evals.
    while time.perf_counter() - t0 < 60:
        free = [random.randrange(p) for _ in range(T - 2)]  # x2..x7
        evals += 1
        st = [C1, C2] + free
        out = perm(pos, st)
        if out[0] % p == 0 and out[1] % p == 0:
            ok = verify_cico_solution(free, prime=p, k=2, t=T, r_f=pos.r_f,
                                      r_p=pos.r_p, mds=M, constants=[C1, C2, 0, 0])
            return dict(found=True, ok=ok, evals=evals, time=time.perf_counter() - t0)
    return dict(found=False, evals=evals, time=time.perf_counter() - t0)


def solver_resultant(pos):
    from cico_verifier import verify_cico_solution
    Xs, Zs = symbols("Xs Zs"); DOM2 = GF(p)
    mds, rc, t = pos.mds, pos.round_constants, pos.t
    def mm(s): return [sum(int(mds[i][j]) * s[j] for j in range(t)) for i in range(t)]
    def arc(s, c): return [s[i] + int(c[i]) for i in range(t)]
    def sb(x): return x ** D
    fixings = 0; t0 = time.perf_counter()
    def roots_in(var, expr):
        f = Poly(expr, var, domain=DOM2)
        if f.degree() <= 0: return []
        b = Poly(var, var, domain=DOM2); acc = Poly(1, var, domain=DOM2); e = p
        while e > 0:
            if e & 1: acc = (acc * b) % f
            e >>= 1
            if e: b = (b * b) % f
        g = sgcd((acc - Poly(var, var, domain=DOM2)) % f, f)
        rs = []
        if g.degree() > 0:
            for (fct, m) in g.factor_list()[1]:
                if fct.degree() == 1:
                    a = int(fct.nth(1)); b0 = int(fct.nth(0))
                    rs.append((-b0 * pow(a, -1, p)) % p)
        return rs
    DI = None
    while time.perf_counter() - t0 < 60:
        fixings += 1
        fixed4 = [random.randrange(p) for _ in range(4)]  # x4..x7
        state = [Poly(C1, Xs, Zs, domain=DOM2), Poly(C2, Xs, Zs, domain=DOM2),
                 Poly(Xs, Xs, Zs, domain=DOM2), Poly(Zs, Xs, Zs, domain=DOM2)]
        state += [Poly(int(v) % p, Xs, Zs, domain=DOM2) for v in fixed4]
        state = mm(state)
        half = pos.r_f // 2; idx = 0
        for _ in range(half):
            state = arc(state, rc[idx]); idx += 1; state = [sb(x) for x in state]; state = mm(state)
        for _ in range(pos.r_p):
            state = arc(state, rc[idx]); idx += 1; state[0] = sb(state[0]); state = mm(state)
        for _ in range(half):
            state = arc(state, rc[idx]); idx += 1; state = [sb(x) for x in state]; state = mm(state)
        P0 = state[0].as_expr(); Q0 = state[1].as_expr()
        f = Poly(P0, Zs, Xs, domain=DOM2); g = Poly(Q0, Zs, Xs, domain=DOM2)
        R = f.resultant(g); DI = Poly(R.as_expr(), Xs, domain=DOM2).degree()
        for xr in roots_in(Xs, R.as_expr()):
            p0z = Poly(P0, Xs, Zs, domain=DOM2).eval({Xs: xr})
            q0z = Poly(Q0, Xs, Zs, domain=DOM2).eval({Xs: xr})
            z0 = set(roots_in(Zs, Poly(p0z, Zs, domain=DOM2).as_expr()))
            z1 = set(roots_in(Zs, Poly(q0z, Zs, domain=DOM2).as_expr()))
            common = z0 & z1
            if common:
                zr = next(iter(common)); free = [xr, zr] + fixed4
                ok = verify_cico_solution(free, prime=p, k=2, t=T, r_f=pos.r_f,
                                          r_p=pos.r_p, mds=M, constants=[C1, C2, 0, 0])
                return dict(found=True, ok=ok, fixings=fixings, DI=DI,
                            time=time.perf_counter() - t0)
    return dict(found=False, fixings=fixings, DI=DI, time=time.perf_counter() - t0)


def main():
    print(f"TINY instance: p={p}, t={T}, alpha={D}, CICO-2 [C1={C1},C2={C2}] -> out0=out1=0")
    for (rf, rp) in [(2, 2)]:
        pos = make_pos(rf, rp)
        print(f"\n--- RF={rf}, RP={rp} ---")
        rW = solver_wildcard(pos)
        print("WILDCARD :", rW)
        rB = solver_brute(pos)
        print("BRUTE    :", rB)
        rR = solver_resultant(pos)
        print("RESULTANT:", rR)
        print("\nINTERPRETATION (this instance):")
        if rW.get("found"):
            print(f"  wildcard needed ~{rW['fixings']} fixings & {rW['rootfinds']} root-finds "
                  f"(each root-find ~ a degree-{D**(rf+rp)} univariate solve = many perm-evals)")
        if rB.get("found"):
            print(f"  brute needed ~{rB['evals']} perm evals (theory: p^2={p*p})")
        if rR.get("found"):
            print(f"  resultant: 1 fixing usually suffices, D_I={rR.get('DI')} "
                  f"(pred d^(2RF+RP)={D**(2*rf+rp)})")


if __name__ == "__main__":
    main()

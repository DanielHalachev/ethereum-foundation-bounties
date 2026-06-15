"""
Probe 11 — a CORRECT, end-to-end resultant-based CICO-2 solver, validated against
the REAL verifier (reference/bounties/cico_verifier.py :: verify_cico_solution).

Goal: prove the attack pipeline is correct and scalable. We solve the EXACT bounty
problem for given (RF, RP) using the verifier's actual Cauchy matrix and the actual
fixed constants C1=0xC09DE4, C2=0xEE6282. We demonstrate on small (RF,RP) that fit
sympy in-session; the same pipeline scales to RF=6,RP=10 modulo compute (D_I=3^22).

Pipeline (no skip; straightforward and obviously correct):
  1. Pin x0=C1, x1=C2. Leave x2=X, x3=Y free. Fix x4..x15 to deterministic constants.
  2. Symbolically evaluate permutation_plus_linear over GF(p)[X,Y] -> out0(X,Y), out1(X,Y).
  3. R(X) = Res_Y(out0, out1). Find F_p-roots of R (sympy ground_roots / factor over GF(p)).
  4. For each root x*, recover Y via gcd(out0(x*,Y), out1(x*,Y)); take its F_p root y*.
  5. Build free_inputs = [x*, y*, fixed...]; check with verify_cico_solution.
  If no F_p solution on this slice, re-randomize the 12 fixed inputs and retry.

Run: ./.venv/bin/python probes/11_resultant_solver.py
"""
import sys, os, time, random
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "reference"))
from sympy import symbols, Poly, GF
from poseidon.poseidon import Poseidon
from poseidon.mds_matrix import generate_mds_matrix
from bounties.cico_verifier import verify_cico_solution

P = 2130706433
T = 16
D = 3
C1, C2 = 0xC09DE4, 0xEE6282
X, Y = symbols("X Y")
DOM = GF(P)  # symmetric=False would matter for printing only


def out_polys(pos, M, fixed12):
    """Return (out0, out1) as sympy expressions in X,Y for permutation_plus_linear."""
    def mds(s): return [sum(int(M[i][j]) * s[j] for j in range(T)) for i in range(T)]
    def arc(s, c): return [s[i] + int(c[i]) for i in range(T)]
    rc = pos.round_constants
    state = [Poly(C1, X, Y, domain=DOM), Poly(C2, X, Y, domain=DOM),
             Poly(X, X, Y, domain=DOM), Poly(Y, X, Y, domain=DOM)]
    state += [Poly(int(v) % P, X, Y, domain=DOM) for v in fixed12]
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
    return state[0].as_expr(), state[1].as_expr()


def fp_roots_univariate(poly_expr, var):
    """All F_p roots (with the X^p-X trick implicit via sympy ground_roots over GF(p))."""
    pol = Poly(poly_expr, var, domain=DOM)
    if pol.is_zero:
        return None  # identically zero -> infinitely many (shouldn't happen generically)
    # ground_roots returns roots in the domain GF(p)
    rd = pol.ground_roots()
    return [int(r) % P for r in rd.keys()]


def solve(rf, rp, max_slices=40, seed=0, verbose=True):
    M = generate_mds_matrix(T, P)
    pos = Poseidon(prime=P, alpha=D, t=T, r_f=rf, r_p=rp, mds=M)
    rng = random.Random(seed)
    DI_pred = D ** (2 * rf + rp)
    if verbose:
        print(f"--- solve RF={rf} RP={rp}  (predicted D_I = 3^{2*rf+rp} = {DI_pred}) ---")
    for slc in range(max_slices):
        fixed12 = [rng.randint(0, P - 1) for _ in range(T - 4)]
        t0 = time.perf_counter()
        o0, o1 = out_polys(pos, M, fixed12)
        f = Poly(o0, Y, X, domain=DOM)
        g = Poly(o1, Y, X, domain=DOM)
        R = f.resultant(g)               # univariate in X over GF(p)
        tb = time.perf_counter() - t0
        if R.is_zero:
            continue
        xroots = fp_roots_univariate(R.as_expr(), X)
        for xr in xroots:
            # recover Y: common root of o0(xr,Y), o1(xr,Y)
            oy0 = Poly(o0.subs(X, xr), Y, domain=DOM)
            oy1 = Poly(o1.subs(X, xr), Y, domain=DOM)
            if oy0.is_zero and oy1.is_zero:
                yroots = [0]
            else:
                gg = oy0.gcd(oy1) if not oy0.is_zero else oy1
                if gg.degree() < 1:
                    continue
                yroots = fp_roots_univariate(gg.as_expr(), Y)
            for yr in yroots:
                free = [xr, yr] + [int(v) % P for v in fixed12]
                if verify_cico_solution(free, r_f=rf, r_p=rp):
                    if verbose:
                        print(f"  [slice {slc}] SOLVED in {tb:.2f}s resultant + rootfind.")
                        print(f"  deg R = {R.degree()} (=D_I), #Fp-roots(R) = {len(xroots)}")
                        print(f"  free_inputs = {free}")
                        print(f"  verify_cico_solution -> True  ✓")
                    return free
        if verbose and slc < 3:
            print(f"  [slice {slc}] deg R={R.degree()}, Fp-roots(R)={len(xroots)}, no full Fp sol; retry")
    if verbose:
        print(f"  no F_p solution found in {max_slices} slices")
    return None


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(); ap.add_argument('--cases', default='2,0;2,1')
    cases = [tuple(int(x) for x in c.split(',')) for c in ap.parse_args().cases.split(';')]
    for (rf, rp) in cases:
        t0 = time.perf_counter()
        sol = solve(rf, rp)
        print(f"  total {time.perf_counter()-t0:.1f}s\n")

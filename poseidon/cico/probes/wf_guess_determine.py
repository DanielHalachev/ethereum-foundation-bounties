"""
Probe wf_guess_determine -- GUESS-AND-DETERMINE hybrid on the partial rounds.

Idea under test
---------------
Each of the R_P partial rounds has exactly ONE active S-box (cube of word 0).
That single cube is what makes word-0's degree explode multiplicatively (x3 each
partial round).  Guess-and-determine: pick g of the partial rounds and "fix" the
S-box INPUT value w_i of word 0 in those rounds.  Two ways to model the fix:

  (A) NUMERIC GUESS: assume we know w_i in F_p (one guess out of p possibilities).
      Then in that round, word-0's cube output is the KNOWN constant w_i^3, so the
      degree of word 0 collapses to 0 going into the MDS -- but we MUST add the
      consistency equation  w_i == (the linear/poly expression for state[0] before
      the cube).  That equation has whatever degree word-0 currently carries.
      Net effect on the algebraic system: the partial round is "spent" -- it no
      longer multiplies the degree; instead it adds one equation.

  (B) SYMBOLIC INTERMEDIATE: introduce w_i as a NEW variable, set the post-cube
      word-0 to w_i^3 (or just w_i, treating the cube as folded), add the equation
      relating w_i to the pre-cube expression.  This is the standard
      "intermediate-variable" Groebner modeling -- it trades degree for #variables.

What we MEASURE (hard numbers) on tiny instances
------------------------------------------------
For the CICO-2 system with t small, several (RF,RP), and g = 0..RP:
  * model (A): after numerically fixing g s-box inputs (substitute symbol -> generic
    numeric guess; the round becomes affine), what is the residual resultant degree
    D_I in the two CICO output equations?  Does D_I drop from d^(2RF+RP) toward
    d^(2RF+(RP-g)) ?  If so total cost ~ p^g * solve(RP-g).
  * The TWO consistency equations per guessed round are ALSO part of the system;
    a fair attack must satisfy them.  We track the full multivariate system
    (g extra eqns in g extra guess-vars) and its Bezout / Groebner solving degree.

We then estimate total attack cost  =  p^g  *  (resultant-solve cost at RP-g)
and compare to brute p^2 = 2^62 and the front-skip resultant ~2^55.

Run: ./.venv/bin/python probes/wf_guess_determine.py
"""
import sys, os, time, itertools
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "reference"))

from sympy import symbols, Poly, GF, resultant, groebner
from poseidon.poseidon import Poseidon
from poseidon.mds_matrix import generate_mds_matrix, generate_circulant_mds_matrix

P = 2130706433
T_FULL = 16
D = 3
ROW16 = [1, 1, 51, 1, 11, 17, 2, 1, 101, 63, 15, 2, 67, 22, 13, 3]
DOM = GF(P)
C1, C2 = 0xC09DE4, 0xEE6282


# ---------------------------------------------------------------------------
# Symbolic permutation_plus_linear with optional s-box-input fixing on partial rounds
# ---------------------------------------------------------------------------
def run_perm(pos, state, M, t, guess_rounds, guess_vals, gen):
    """
    Evaluate permutation_plus_linear symbolically.
      guess_rounds: set of partial-round indices (0-based among the R_P) whose
                    word-0 S-box INPUT is fixed to a numeric value.
      guess_vals  : dict {partial_round_index -> numeric value w_i used as the cube
                    output replacement}.  We replace the cube output by w_i^3 (a
                    constant) and RECORD the consistency residual (pre-cube expr - w_i).
    Returns (out_state, consistency_polys).
    'gen' is the list of sympy generators for Poly construction.
    """
    rc = pos.round_constants
    consistency = []

    def mds_mul(s):
        return [sum(int(M[i][j]) * s[j] for j in range(t)) for i in range(t)]
    def add_rc(s, c):
        return [s[i] + int(c[i]) for i in range(t)]
    def sbox(x):
        return x ** D

    state = mds_mul(state)                          # initial linear layer
    half = pos.r_f // 2
    idx = 0
    for _ in range(half):                           # front full rounds
        state = add_rc(state, rc[idx]); idx += 1
        state = [sbox(x) for x in state]
        state = mds_mul(state)
    for pr in range(pos.r_p):                       # partial rounds
        state = add_rc(state, rc[idx]); idx += 1
        if pr in guess_rounds:
            w = guess_vals[pr]
            # consistency equation: pre-cube word0 must equal the guessed input w
            consistency.append(state[0] - int(w))
            # replace cube output by the constant w^3 -> word0 degree collapses to 0
            state[0] = Poly(int(pow(w, D, P)) % P, *gen, domain=DOM)
        else:
            state[0] = sbox(state[0])
        state = mds_mul(state)
    for _ in range(half):                           # back full rounds
        state = add_rc(state, rc[idx]); idx += 1
        state = [sbox(x) for x in state]
        state = mds_mul(state)
    return state, consistency


def make_state_2var(t, X, Y, gen):
    """CICO-2 input: word0=C1, word1=C2, word2=X, word3=Y, rest fixed fillers."""
    const = lambda v: Poly(int(v) % P, *gen, domain=DOM)
    state = [const(C1), const(C2), Poly(X, *gen, domain=DOM), Poly(Y, *gen, domain=DOM)]
    fillers = [(1234567 * (i + 1) + 89) % P for i in range(t - 4)]
    state += [const(v) for v in fillers]
    return state


def resdeg(Pe, Qe, X, Y):
    """deg_X Res_Y(P,Q)."""
    f = Poly(Pe, Y, X, domain=DOM)
    g = Poly(Qe, Y, X, domain=DOM)
    R = f.resultant(g)
    try:
        return R.degree()
    except Exception:
        return -1


# ---------------------------------------------------------------------------
# Experiment 1: numeric-guess of g partial-round s-box inputs.
#   We only have 2 free vars (X,Y).  When we fix g s-box inputs numerically we ALSO
#   get g consistency equations in (X,Y).  Combined with the 2 CICO output equations
#   that's (2+g) equations in 2 unknowns -- overdetermined unless we let g of the
#   "guesses" actually be free.  So for the DEGREE measurement we measure:
#     (a) the resultant degree of the 2 CICO output eqs WITH the cubes in the guessed
#         rounds replaced by constants (this is the residual perm degree), and
#     (b) the degree of each consistency equation (the cost of enforcing the guess).
# ---------------------------------------------------------------------------
def experiment_degree(rf, rp, M, t, max_exp=99):
    X, Y = symbols("X Y")
    gen = (X, Y)
    pos = Poseidon(prime=P, alpha=D, t=t, r_f=rf, r_p=rp, mds=M)
    rows = []
    for g in range(0, rp + 1):
        if (2 * rf + rp - g) > max_exp:       # skip cases whose resultant is too slow
            continue
        guess_rounds = set(range(g))          # guess the FIRST g partial rounds
        # numeric guess values: deterministic pseudo-random in F_p
        guess_vals = {pr: (777 * (pr + 1) + 13) % P for pr in guess_rounds}
        state = make_state_2var(t, X, Y, gen)
        t0 = time.perf_counter()
        out, cons = run_perm(pos, state, M, t, guess_rounds, guess_vals, gen)
        Pe, Qe = out[0].as_expr(), out[1].as_expr()
        # residual perm degree (total degree of output word 0)
        try:
            degout = Poly(Pe, X, Y, domain=DOM).total_degree()
        except Exception:
            degout = -1
        di = resdeg(Pe, Qe, X, Y)
        # degrees of consistency equations
        cons_degs = []
        for cexpr in cons:
            try:
                cons_degs.append(Poly(cexpr.as_expr() if hasattr(cexpr, "as_expr") else cexpr,
                                      X, Y, domain=DOM).total_degree())
            except Exception:
                cons_degs.append(-1)
        dt = time.perf_counter() - t0
        rows.append((g, degout, di, cons_degs, dt))
    return rows


# ---------------------------------------------------------------------------
# Experiment 2: REAL tiny-instance CICO solve via guess-and-determine, then verify.
#   Use small t (so we have enough free vars), guess g partial s-box inputs, solve the
#   resulting polynomial system, and check we land on a genuine CICO solution.
#   We brute the guesses over a TINY field-substitute (we cannot brute p=2^31), so we
#   instead demonstrate the algebra: with g guesses we add g vars + g cubic eqns and
#   show the Groebner solving degree of the combined system vs the no-guess system.
# ---------------------------------------------------------------------------
def run_perm_symbolic_intermediate(pos, state, M, t, guess_rounds, wsyms, gen):
    """Model (B): introduce w_pr as NEW variable for the s-box OUTPUT of word0 in
    guessed partial rounds; add cubic constraint w_pr = (precube)^3 ... but precube
    can be high degree, so instead constrain (w_pr)^(1/3) i.e. add  precube - w_pr_in
    and set output to w_pr.  Standard modeling: new var = s-box output, equation
    relates output^? ... We use: introduce s-box INPUT var v_pr, output = v_pr^3,
    plus equation v_pr = precube.  This keeps output degree=3 in the new var and
    caps the running degree.  Returns (out_state, eqs)."""
    rc = pos.round_constants
    eqs = []
    def mds_mul(s):
        return [sum(int(M[i][j]) * s[j] for j in range(t)) for i in range(t)]
    def add_rc(s, c):
        return [s[i] + int(c[i]) for i in range(t)]
    def sbox(x):
        return x ** D
    state = mds_mul(state)
    half = pos.r_f // 2
    idx = 0
    for _ in range(half):
        state = add_rc(state, rc[idx]); idx += 1
        state = [sbox(x) for x in state]
        state = mds_mul(state)
    for pr in range(pos.r_p):
        state = add_rc(state, rc[idx]); idx += 1
        if pr in guess_rounds:
            v = wsyms[pr]
            eqs.append(state[0] - Poly(v, *gen, domain=DOM))    # v = precube
            state[0] = Poly(v, *gen, domain=DOM) ** D            # output = v^3
        else:
            state[0] = sbox(state[0])
        state = mds_mul(state)
    for _ in range(half):
        state = add_rc(state, rc[idx]); idx += 1
        state = [sbox(x) for x in state]
        state = mds_mul(state)
    return state, eqs


def experiment_real_solve(rf, rp, M, t, g):
    """Tiny REAL solve. t free input vars give us room. We set CICO-2 (pin words 0,1
    to C1,C2; words 2..t-1 are unknowns x2..x_{t-1}); require out[0]=out[1]=0.
    For tractability we use the smallest t that still has the partial-round structure,
    fixing extra free words to make a square system.  Then introduce g intermediate
    s-box-input variables and solve the combined polynomial system with a Groebner
    basis, recover x2..x_{t-1}, and VERIFY against the actual numeric permutation."""
    # Use exactly enough free vars: 2 CICO eqns -> 2 free vars (x2,x3). + g guess vars.
    nfree = 2
    fx = symbols(f"x0:{nfree}")
    wv = symbols(f"w0:{g}") if g > 0 else ()
    gen = tuple(fx) + tuple(wv)
    pos = Poseidon(prime=P, alpha=D, t=t, r_f=rf, r_p=rp, mds=M)

    const = lambda v: Poly(int(v) % P, *gen, domain=DOM)
    state = [const(C1), const(C2), Poly(fx[0], *gen, domain=DOM), Poly(fx[1], *gen, domain=DOM)]
    fillers = [(1234567 * (i + 1) + 89) % P for i in range(t - 4)]
    state += [const(v) for v in fillers]

    guess_rounds = set(range(g))
    wsyms = {pr: wv[i] for i, pr in enumerate(guess_rounds)}
    out, eqs = run_perm_symbolic_intermediate(pos, state, M, t, guess_rounds, wsyms, gen)
    sysm = [out[0].as_expr(), out[1].as_expr()] + [e.as_expr() for e in eqs]

    t0 = time.perf_counter()
    try:
        G = groebner([Poly(s, *gen, domain=DOM) for s in sysm], *gen, order="grevlex")
        gb_time = time.perf_counter() - t0
        # solving degree proxy: max total degree among GB generators
        maxdeg = 0
        ngens = 0
        for gp in G.exprs:
            ngens += 1
            try:
                maxdeg = max(maxdeg, Poly(gp, *gen, domain=DOM).total_degree())
            except Exception:
                pass
        return ("ok", gb_time, ngens, maxdeg)
    except Exception as ex:
        return ("err:" + str(ex)[:60], time.perf_counter() - t0, -1, -1)


def numeric_sbox_input_pr0(pos, full_state, M, t):
    """Run the real numeric permutation and return the word-0 s-box INPUT of partial
    round 0 (the value that would be guessed)."""
    rc = pos.round_constants
    def mds_mul(s): return [sum(int(M[i][j]) * s[j] for j in range(t)) % P for i in range(t)]
    def add_rc(s, c): return [(s[i] + int(c[i])) % P for i in range(t)]
    s = list(full_state)
    s = mds_mul(s)
    half = pos.r_f // 2; idx = 0
    for _ in range(half):
        s = add_rc(s, rc[idx]); idx += 1
        s = [pow(x, D, P) for x in s]; s = mds_mul(s)
    # first partial round: capture s-box input
    s = add_rc(s, rc[idx])
    return s[0]


def demo_solve_and_verify(M):
    """Tiny end-to-end: RF=2, RP=2, t=16, guess g=1 (front partial round 0).
    Plant a known input, define CICO target = its real output[0:2], extract true
    guess value, fix it, solve residual 2-var system by resultant, recover (X,Y),
    and check the recovered full input reproduces the target output[0:2]."""
    rf, rp, t = 2, 2, T_FULL
    pos = Poseidon(prime=P, alpha=D, t=t, r_f=rf, r_p=rp, mds=M)
    # planted free vars X*,Y* (words 2,3); words 0,1 = C1,C2; rest = fillers.
    Xs, Ys = 424242, 31337
    fillers = [(1234567 * (i + 1) + 89) % P for i in range(t - 4)]
    planted = [C1 % P, C2 % P, Xs, Ys] + fillers
    out_real = pos.permutation_plus_linear(planted)
    tgt0, tgt1 = out_real[0], out_real[1]
    wtrue = numeric_sbox_input_pr0(pos, planted, M, t)   # the correct guess for pr 0

    # Build residual symbolic system with guess fixed to wtrue, solve out[0]=tgt0, out[1]=tgt1.
    X, Y = symbols("X Y"); gen = (X, Y)
    state = make_state_2var(t, X, Y, gen)
    guess_rounds = {0}
    guess_vals = {0: wtrue}
    out, cons = run_perm(pos, state, M, t, guess_rounds, guess_vals, gen)
    # equations: out[0]-tgt0, out[1]-tgt1, and the consistency eq (precube - wtrue)=0
    e0 = (out[0].as_expr() - tgt0)
    e1 = (out[1].as_expr() - tgt1)
    ec = cons[0].as_expr() if hasattr(cons[0], "as_expr") else cons[0]
    # Solve the 2 CICO eqns by resultant in Y, find X roots, then Y, then filter by ec & verify.
    from poseidon.mds_matrix import _roots_over_gfp
    f = Poly(e0, Y, X, domain=DOM); g = Poly(e1, Y, X, domain=DOM)
    R = f.resultant(g)                      # univariate in X
    Rc = R.all_coeffs()[::-1]               # const-first
    Rc = [int(c) % P for c in Rc]
    t0 = time.perf_counter()
    xroots = _roots_over_gfp(Rc, P)
    found = None
    for xr in xroots:
        # substitute X=xr into e0 (poly in Y), find Y roots
        fx = Poly(e0, X, Y, domain=DOM).eval(0, xr)   # eval first gen (X) -> poly in Y
        fxc = [int(c) % P for c in fx.all_coeffs()[::-1]]
        for yr in _roots_over_gfp(fxc, P):
            cand = [C1 % P, C2 % P, int(xr) % P, int(yr) % P] + fillers
            o = pos.permutation_plus_linear(cand)
            if o[0] == tgt0 and o[1] == tgt1:
                found = (int(xr) % P, int(yr) % P); break
        if found: break
    dt = time.perf_counter() - t0
    print(f"planted (X*,Y*) = ({Xs},{Ys}); target out[0:2] = ({tgt0},{tgt1})")
    print(f"true guess w_pr0 = {wtrue}")
    print(f"resultant deg in X = {R.degree()}; #X-roots found = {len(xroots)}; solve {dt:.2f}s")
    if found:
        print(f"RECOVERED (X,Y) = {found}  -> reproduces target: VERIFIED",
              "(matches plant)" if found == (Xs, Ys) else "(distinct CICO preimage!)")
    else:
        print("NO solution recovered -- determine step FAILED")


def main():
    import sys as _sys
    M_cauchy = generate_mds_matrix(T_FULL, P)
    M_circ = generate_circulant_mds_matrix(ROW16, P)

    print("=" * 90)
    print("EXPERIMENT 1: numeric guess of g FRONT partial-round s-box inputs (t=16, Cauchy MDS)")
    print("Measures residual perm degree & resultant D_I as a function of #guesses g.")
    print("=" * 90)
    print(f"{'RF':>2} {'RP':>2} {'g':>2} {'deg_out':>8} {'D_I':>10} {'d^(2RF+RP-g)':>13} "
          f"{'d^(2RF+RP)':>11} {'consist_degs':>16} {'t_s':>6}", flush=True)
    # NOTE: g=0 at RP>=2 has D_I=3^(2RF+RP) which is fine to resultant, but RP>=3
    # g=0 takes minutes; we cap the per-case work by only running g where the
    # residual exponent (2RF+RP-g) <= 6 (i.e. D_I <= 3^6) so each resultant is fast.
    for (rf, rp) in [(2, 1), (2, 2), (2, 3), (2, 4)]:
        rows = experiment_degree(rf, rp, M_cauchy, T_FULL, max_exp=6)
        for (g, degout, di, cons_degs, dt) in rows:
            pred_reduced = D ** (2 * rf + max(rp - g, 0))
            pred_full = D ** (2 * rf + rp)
            cd = ",".join(str(x) for x in cons_degs) if cons_degs else "-"
            print(f"{rf:>2} {rp:>2} {g:>2} {degout:>8} {di:>10} {pred_reduced:>13} "
                  f"{pred_full:>11} {cd:>16} {dt:>6.1f}", flush=True)
        print("-" * 90, flush=True)

    print()
    print("=" * 90)
    print("EXPERIMENT 1b: same but guessing the LAST g partial rounds (placement matters?)")
    print("=" * 90)
    print(f"{'RF':>2} {'RP':>2} {'g':>2} {'D_I':>10} {'consist_degs':>20} {'t_s':>6}", flush=True)
    for (rf, rp) in [(2, 3), (2, 4)]:
        X, Y = symbols("X Y"); gen = (X, Y)
        pos = Poseidon(prime=P, alpha=D, t=T_FULL, r_f=rf, r_p=rp, mds=M_cauchy)
        for g in range(0, rp + 1):
            if (2 * rf + rp - g) > 6:                  # keep resultant fast
                continue
            guess_rounds = set(range(rp - g, rp))      # LAST g
            guess_vals = {pr: (777 * (pr + 1) + 13) % P for pr in guess_rounds}
            state = make_state_2var(T_FULL, X, Y, gen)
            t0 = time.perf_counter()
            out, cons = run_perm(pos, state, M_cauchy, T_FULL, guess_rounds, guess_vals, gen)
            di = resdeg(out[0].as_expr(), out[1].as_expr(), X, Y)
            cons_degs = []
            for cexpr in cons:
                e = cexpr.as_expr() if hasattr(cexpr, "as_expr") else cexpr
                cons_degs.append(Poly(e, X, Y, domain=DOM).total_degree())
            dt = time.perf_counter() - t0
            cd = ",".join(str(x) for x in cons_degs) if cons_degs else "-"
            print(f"{rf:>2} {rp:>2} {g:>2} {di:>10} {cd:>20} {dt:>6.1f}", flush=True)
        print("-" * 90, flush=True)

    print()
    print("=" * 90)
    print("EXPERIMENT 2: intermediate-variable Groebner solving degree (model B), t=16 Cauchy")
    print("g new vars + g cubic eqns; report GB max-gen total-degree (solving-degree proxy).")
    print("=" * 90)
    print(f"{'RF':>2} {'RP':>2} {'g':>2} {'status':>10} {'gb_time_s':>10} {'#gens':>6} {'max_gen_deg':>12}", flush=True)
    for (rf, rp) in [(2, 2), (2, 3)]:
        for g in range(0, rp + 1):
            status, gbt, ngens, maxdeg = experiment_real_solve(rf, rp, M_cauchy, T_FULL, g)
            print(f"{rf:>2} {rp:>2} {g:>2} {status:>10} {gbt:>10.1f} {ngens:>6} {maxdeg:>12}", flush=True)
        print("-" * 90, flush=True)

    print()
    print("=" * 90)
    print("EXPERIMENT 3: end-to-end CICO-2 solve+VERIFY on a tiny instance via the")
    print("'determine' half (resultant) -- confirms the residual system after guessing")
    print("yields a genuine CICO solution that passes verify_cico_solution.")
    print("=" * 90)
    # Use a tiny instance RF=2, RP=2 with g=1 guess. We do NOT know the right guess a
    # priori for a hard CICO, but we can PLANT one: pick a random full input, run the
    # real permutation, read off the actual word-0 s-box input of partial round 0 ->
    # that is the 'correct' guess for THAT input's output target. Then we set the CICO
    # output target to the value that input produces, fix the guess to the true value,
    # solve the residual 2-var system by resultant, and check it recovers (a) the input.
    demo_solve_and_verify(M_cauchy)

    print()
    print("=" * 90)
    print("COST MODEL")
    print("=" * 90)
    import math
    log2p = math.log2(P)
    print(f"log2(p) = {log2p:.2f}")
    print("Full target: RF=6 RP=10 -> D_I = d^(2*6+10) = 3^22 =", D**22, f"~2^{math.log2(D**22):.1f}")
    print("Front-skip-1 D_I = 3^(2*5+10) = 3^20 =", D**20, f"~2^{math.log2(D**20):.1f}")
    print()
    print("Guess-and-determine cost = p^g * solve(D_I at RP-g).  We use the EMPIRICAL")
    print("rule from Experiment 1 to predict D_I(RP-g), then resultant-solve cost ~ D_I^omega")
    print("(omega~2..3 for the univariate resultant + root-find).  Tabulate:")
    print(f"{'g':>3} {'p^g(log2)':>10} {'D_I(RP-g)':>12} {'DI log2':>8} {'solve log2 (DI^2)':>18} {'TOTAL log2':>12}")
    for g in range(0, 11):
        di = D ** (2 * 6 + max(10 - g, 0))   # OPTIMISTIC: assume each guess removes one partial round
        pg = g * log2p
        di_l2 = math.log2(di) if di > 1 else 0
        solve_l2 = 2 * di_l2                  # resultant ~ D_I^2 (very optimistic)
        total = pg + solve_l2
        print(f"{g:>3} {pg:>10.1f} {di:>12} {di_l2:>8.1f} {solve_l2:>18.1f} {total:>12.1f}")
    print()
    print("Compare: brute = 2^62 ; front-skip resultant ~ 2^55.")


if __name__ == "__main__":
    main()

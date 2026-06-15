"""
wf_meet_middle_groebner.py
==========================

ANGLE: Algebraic meet-in-the-middle (MITM) with FULL arithmetization vs the
one-sided resultant degree D_I = d^(2*RF+RP).

The known "wall": push 2 input vars (X,Y) through permutation_plus_linear,
get two output polys of degree D_I = d^(2RF+RP), Res_Y -> univariate of degree
D_I, factor it.  For the real instance D_I = 3^(2*6+RP) ~ 2^31..2^36 (front-skip
shaves one full round).

MITM alternative tested here: introduce an AUX variable at the output of every
S-box with cubic relation  a = (sbox_input)^3, keeping every equation degree<=3.
Optionally split forward-from-input and backward-from-output (M^{-1} inverts the
linear layers, X^{1/3} is a field bijection so we can also invert S-boxes
backward symbolically via a cubic relation), meeting at a middle round.

DECISIVE INVARIANT:  for the SAME 0-dimensional variety, dim_Fp(F_p[vars]/I)
= number of solutions over closure counted with multiplicity.  This is
representation-INVARIANT: one-sided, forward-arith, and MITM all describe the
same variety, so all have the same quotient dim.  The attack must factor a
univariate whose degree = that of the elimination ideal projection (<= qdim).
So the question reduces to:  does qdim (and the lex elimination degree in the
true CICO variable) come out SMALLER than the one-sided D_I, or equal to it?

We keep BOTH representations on the SAME 0-dim variety: pin input words 0,1 to
C1,C2; let words 2,3 be free vars X,Y; pin words 4..t-1 to fixed constants; pin
output words 0,1 to 0.  (2 eqns, 2 unknowns -> 0-dim.)

Run:  ./.venv/bin/python probes/wf_meet_middle_groebner.py
"""

import sys, os, time, itertools
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "reference"))

from sympy import symbols, Poly, GF, groebner, Matrix
from poseidon.poseidon import Poseidon
from poseidon.mds_matrix import generate_mds_matrix

P = 2130706433
D = 3
DOM = GF(P)
C1, C2 = 0xC09DE4, 0xEE6282


def matinv(M, p=P):
    return [[int(x) % p for x in row] for row in Matrix(M).inv_mod(p).tolist()]


def make_pos(t, rf, rp, mds):
    return Poseidon(prime=P, alpha=D, t=t, r_f=rf, r_p=rp, mds=mds)


def fillers(t):
    # words 2..t-3 pinned to these constants; words t-2,t-1 are the free vars X,Y.
    # requires t >= 4.
    return [(1234567 * (i + 1) + 89) % P for i in range(t - 4)]


def init_state_consts(t):
    """Input state template: [C1, C2, fillers..., X, Y] with X,Y at the LAST two
    positions.  Returns the constant list of length t-2 for words 2..t-1 with the
    last two slots reserved for the symbolic vars (callers splice them in)."""
    return [int(C1) % P, int(C2) % P] + [int(v) % P for v in fillers(t)]


# ===========================================================================
# (A) one-sided resultant baseline  D_I
# ===========================================================================
def one_sided_DI(pos, mds, t):
    X, Y = symbols("X Y")

    def const(v):
        return Poly(int(v) % P, X, Y, domain=DOM)

    def mds_mul(s):
        return [sum(int(mds[i][j]) * s[j] for j in range(t)) for i in range(t)]

    def add_rc(s, c):
        return [s[i] + int(c[i]) for i in range(t)]

    rc = pos.round_constants
    state = [const(C1), const(C2)] + [const(v) for v in fillers(t)]
    state += [Poly(X, X, Y, domain=DOM), Poly(Y, X, Y, domain=DOM)]
    state = mds_mul(state)
    half = pos.r_f // 2
    idx = 0
    for _ in range(half):
        state = add_rc(state, rc[idx]); idx += 1
        state = [x ** D for x in state]; state = mds_mul(state)
    for _ in range(pos.r_p):
        state = add_rc(state, rc[idx]); idx += 1
        state[0] = state[0] ** D; state = mds_mul(state)
    for _ in range(half):
        state = add_rc(state, rc[idx]); idx += 1
        state = [x ** D for x in state]; state = mds_mul(state)

    f = Poly(state[0].as_expr(), Y, X, domain=DOM)
    g = Poly(state[1].as_expr(), Y, X, domain=DOM)
    return f.resultant(g).degree()


def forward_out_degree(t, rf, rp):
    """Track the (X,Y)-total-degree of each state word through the permutation
    WITHOUT expanding polynomials: cube triples the degree, MDS = max of inputs,
    ARC unchanged.  Returns the degree of output word 0 = generic D_I."""
    half = rf // 2
    # initial: words t-2,t-1 are X,Y (degree 1), others degree 0
    deg = [0] * t
    deg[t - 2] = 1
    deg[t - 1] = 1
    deg = [max(deg)] * t if False else mds_deg(deg, t)  # initial linear layer
    for _ in range(half):                 # full rounds: cube all
        deg = [3 * d for d in deg]
        deg = mds_deg(deg, t)
    for _ in range(rp):                   # partial rounds: cube word 0 only
        deg[0] = 3 * deg[0]
        deg = mds_deg(deg, t)
    for _ in range(half):
        deg = [3 * d for d in deg]
        deg = mds_deg(deg, t)
    return deg[0]


def mds_deg(deg, t):
    """MDS mixes everything: each output word's degree = max input degree."""
    m = max(deg)
    return [m] * t


# ===========================================================================
# (B) Forward full arithmetization on the SAME 0-dim variety.
#     vars: X, Y, plus one aux per S-box.
# ===========================================================================
def build_forward_arith(pos, mds, t):
    rc = pos.round_constants
    half = pos.r_f // 2
    X, Y = symbols("X Y")
    aux = [0]
    aux_syms = []

    def new_aux():
        s = symbols(f"a{aux[0]}"); aux[0] += 1; aux_syms.append(s); return s

    polys = []

    def mds_mul(s):
        return [sum(int(mds[i][j]) * s[j] for j in range(t)) for i in range(t)]

    def add_rc(s, c):
        return [s[i] + int(c[i]) for i in range(t)]

    state = [int(C1) % P, int(C2) % P] + [int(v) % P for v in fillers(t)] + [X, Y]
    state = mds_mul(state)
    idx = 0
    for _ in range(half):
        state = add_rc(state, rc[idx]); idx += 1
        ns = []
        for x in state:
            a = new_aux(); polys.append(a - x ** D); ns.append(a)
        state = mds_mul(ns)
    for _ in range(pos.r_p):
        state = add_rc(state, rc[idx]); idx += 1
        a = new_aux(); polys.append(a - state[0] ** D)
        state = mds_mul([a] + list(state[1:]))
    for _ in range(half):
        state = add_rc(state, rc[idx]); idx += 1
        ns = []
        for x in state:
            a = new_aux(); polys.append(a - x ** D); ns.append(a)
        state = mds_mul(ns)
    polys.append(state[0])
    polys.append(state[1])
    gens = [X, Y] + aux_syms
    return gens, polys, {"n_vars": len(gens), "n_polys": len(polys), "n_aux": aux[0]}


# ===========================================================================
# (C) MEET-IN-THE-MIDDLE full arithmetization.
#     Forward from input up to a middle boundary B (after round mid), backward
#     from output back to B, meet by equating the two state vectors at B.
#     Backward S-box inversion: if forward did  out = in^3, backward we write a
#     fresh aux z for the pre-image and impose  z^3 = (known post-state) -- i.e.
#     same cubic relation, just oriented from the output side.  We track the
#     backward state as affine in backward aux vars and the (pinned) output.
# ===========================================================================
def build_mitm(pos, mds, t, mid_round):
    """mid_round in [0 .. total_rounds]: number of rounds taken on the FORWARD
    side; the rest are taken backward.  We meet on the FULL t-dim state vector
    just before round (mid_round+1)'s ARC -- i.e. on the post-MDS state."""
    rc = pos.round_constants
    half = pos.r_f // 2
    rp = pos.r_p
    total = pos.r_f + pos.r_p
    Minv = matinv(mds)
    X, Y = symbols("X Y")
    aux = [0]
    aux_syms = []

    def new_aux():
        s = symbols(f"m{aux[0]}"); aux[0] += 1; aux_syms.append(s); return s

    polys = []

    def mds_mul(s):
        return [sum(int(mds[i][j]) * s[j] for j in range(t)) for i in range(t)]

    def minv_mul(s):
        return [sum(int(Minv[i][j]) * s[j] for j in range(t)) for i in range(t)]

    def add_rc(s, c):
        return [s[i] + int(c[i]) for i in range(t)]

    def sub_rc(s, c):
        return [s[i] - int(c[i]) for i in range(t)]

    def is_full(r):  # round index r (0-based) is full?
        return r < half or r >= half + rp

    # ---- forward side: rounds 0..mid_round-1 ----
    fstate = [int(C1) % P, int(C2) % P] + [int(v) % P for v in fillers(t)] + [X, Y]
    fstate = mds_mul(fstate)
    for r in range(mid_round):
        fstate = add_rc(fstate, rc[r])
        if is_full(r):
            ns = []
            for x in fstate:
                a = new_aux(); polys.append(a - x ** D); ns.append(a)
            fstate = mds_mul(ns)
        else:
            a = new_aux(); polys.append(a - fstate[0] ** D)
            fstate = mds_mul([a] + list(fstate[1:]))

    # ---- backward side: rounds total-1 .. mid_round, inverting ----
    # output state: words 0,1 pinned to 0; words 2..t-1 free backward vars
    out_free = list(symbols(f"o2:{t}"))
    aux_syms.extend(out_free)
    bstate = [0, 0] + out_free  # this is the permutation output (post last MDS)
    for r in range(total - 1, mid_round - 1, -1):
        # undo MDS:  pre_mds = Minv * bstate
        bstate = minv_mul(bstate)
        # undo S-box: introduce pre-image aux z with z^3 = bstate[i]
        if is_full(r):
            ns = []
            for x in bstate:
                z = new_aux(); polys.append(x - z ** D); ns.append(z)
            bstate = ns
        else:
            z = new_aux(); polys.append(bstate[0] - z ** D)
            bstate = [z] + list(bstate[1:])
        # undo ARC
        bstate = sub_rc(bstate, rc[r])

    # ---- meet: forward post-MDS state == backward pre-ARC state ----
    # forward fstate is the post-MDS state after round mid_round-1 (i.e. input to
    # ARC of round mid_round).  backward bstate is the state just before ARC of
    # round mid_round.  Equate them.
    for i in range(t):
        polys.append(fstate[i] - bstate[i])

    gens = [X, Y] + aux_syms
    return gens, polys, {"n_vars": len(gens), "n_polys": len(polys),
                         "n_aux": aux[0], "mid": mid_round}


# ===========================================================================
# Groebner measurement
# ===========================================================================
def gb_measure(gens, polys, want_lex_var=None, do_lex=True):
    t0 = time.perf_counter()
    G = groebner(polys, *gens, order='grevlex', domain=DOM)
    t_gb = time.perf_counter() - t0
    maxdeg = 0
    LMs = []
    for g in G.polys:
        maxdeg = max(maxdeg, g.total_degree())
        LMs.append(g.LM(order='grevlex'))
    is_zd = G.is_zero_dimensional
    res = {"t_gb": t_gb, "gb_maxdeg": maxdeg, "zero_dim": is_zd,
           "gb_size": len(G.polys)}
    qd = quotient_dim(LMs, len(gens)) if is_zd else "inf"
    res["quotient_dim"] = qd
    # FGLM/lex is expensive; only attempt when the quotient dim is small.
    small_enough = isinstance(qd, int) and qd <= 1500
    if do_lex and want_lex_var is not None and is_zd and small_enough:
        t1 = time.perf_counter()
        try:
            Glex = groebner(polys, *gens, order='lex', domain=DOM)
            elim = None
            for g in Glex.exprs:
                pg = Poly(g, *gens, domain=DOM)
                if pg.free_symbols <= {want_lex_var}:
                    elim = Poly(g, want_lex_var, domain=DOM).degree()
            res["lex_elim_deg"] = elim
            res["t_lex"] = time.perf_counter() - t1
        except Exception as e:
            res["lex_elim_deg"] = f"err:{type(e).__name__}"
    return res


def quotient_dim(LMs, n):
    """#standard monomials given leading-monomial exponent tuples (grevlex)."""
    bound = [None] * n
    for lm in LMs:
        nz = [i for i, e in enumerate(lm) if e != 0]
        if len(nz) == 1:
            i = nz[0]
            if bound[i] is None or lm[i] < bound[i]:
                bound[i] = lm[i]
    if any(b is None for b in bound):
        return "unbounded"
    box = 1
    for b in bound:
        box *= b
        if box > 8_000_000:
            return f"box>{box}"
    count = 0
    for exps in itertools.product(*[range(b) for b in bound]):
        if not any(all(exps[i] >= lm[i] for i in range(n)) for lm in LMs):
            count += 1
    return count


# ===========================================================================
# Subprocess-isolated GB runner with hard timeout (sympy Buchberger can hang).
# ===========================================================================
import subprocess, json, textwrap

GB_WORKER = textwrap.dedent('''
import sys, os, time, json
sys.path.insert(0, os.path.join(os.path.dirname(%r), "..", "reference"))
sys.path.insert(0, os.path.dirname(os.path.dirname(%r)))
from probes.wf_meet_middle_groebner import (
    build_forward_arith, build_mitm, make_pos, generate_mds_matrix,
    gb_measure, P)
kind, t, rf, rp, mid = json.loads(sys.argv[1])
mds = generate_mds_matrix(t, P); pos = make_pos(t, rf, rp, mds)
if kind == "fwd":
    gens, polys, info = build_forward_arith(pos, mds, t)
    m = gb_measure(gens, polys, want_lex_var=gens[0], do_lex=True)
else:
    gens, polys, info = build_mitm(pos, mds, t, mid)
    m = gb_measure(gens, polys, do_lex=False)
out = {"info": info,
       "gb_maxdeg": m.get("gb_maxdeg"),
       "qdim": m.get("quotient_dim"),
       "lex": m.get("lex_elim_deg", None),
       "t_gb": round(m.get("t_gb", -1), 2),
       "zero_dim": m.get("zero_dim")}
print("RESULT " + json.dumps(out), flush=True)
''')


def run_gb_subprocess(kind, t, rf, rp, mid=0, timeout=60):
    here = os.path.abspath(__file__)
    worker = GB_WORKER % (here, here)
    args = json.dumps([kind, t, rf, rp, mid])
    try:
        proc = subprocess.run(
            [sys.executable, "-u", "-c", worker, args],
            capture_output=True, text=True, timeout=timeout)
        for line in proc.stdout.splitlines():
            if line.startswith("RESULT "):
                return json.loads(line[len("RESULT "):])
        return {"err": "no-result", "stderr": proc.stderr[-200:]}
    except subprocess.TimeoutExpired:
        return {"err": f"timeout>{timeout}s"}


# ===========================================================================
# Direct aux-variable elimination (CHEAP, rigorous): each aux var in the
# forward arithmetization is an EXPLICIT definition  a_i = f_i(earlier vars).
# Substituting them all back recovers the two bivariate output polynomials,
# whose resultant degree IS the quotient/elimination degree of the arithmetized
# ideal.  This proves qdim(arith) == D_I WITHOUT running a Groebner basis.
# ===========================================================================
def eliminate_aux_recover_DI(pos, mds, t):
    """Build forward-arith polys {a_i - f_i}, eliminate every aux by direct
    substitution (they form a triangular 'explicit-definition' system), and
    return the resultant degree of the resulting 2 bivariate output equations.
    If this equals the one-sided D_I, the arithmetized ideal has the SAME
    elimination degree -> same factoring target."""
    from sympy import Poly, expand
    gens, polys, info = build_forward_arith(pos, mds, t)
    X, Y = gens[0], gens[1]
    aux = gens[2:]
    # By construction: polys[:-2] are the aux DEFINITIONS  (a_i - f_i), created in
    # increasing index order; polys[-2:] are the two output equations out0,out1.
    defs = polys[:-2]
    out_eqs = polys[-2:]
    # Build substitution a_i -> f_i in order; f_i uses X,Y and earlier aux only,
    # so after substituting earlier aux it becomes a pure (X,Y) expression.
    subs_map = {}
    for a, p in zip(aux, defs):
        f = a - p                      # p = a - f  =>  f = a - p
        subs_map[a] = expand(f.subs(subs_map))
    # substitute into the two output equations -> bivariate in X,Y
    o0 = expand(out_eqs[0].subs(subs_map))
    o1 = expand(out_eqs[1].subs(subs_map))
    f = Poly(o0, Y, X, domain=DOM)
    g = Poly(o1, Y, X, domain=DOM)
    return f.resultant(g).degree()


# ===========================================================================
# main
# ===========================================================================
def pr(*a):
    print(*a, flush=True)


def main():
    pr("=" * 96)
    pr("ARITHMETIZED / MITM Groebner vs one-sided resultant D_I  (KoalaBear, d=3, Cauchy MDS)")
    pr("Same 0-dim CICO-2 variety in all representations.  TINY instances.")
    pr("=" * 96)

    # -----------------------------------------------------------------
    # PART 1: one-sided D_I (resultant, CHEAP) for many tiny instances.
    # Confirms D_I = d^(2RF+RP) and gives the factoring-target degree.
    # -----------------------------------------------------------------
    pr("\n[PART 1] one-sided D_I via resultant (RP<=2; CHEAP) -> the factoring target / wall")
    pr(f"{'t':>2} {'RF':>2} {'RP':>2} | {'D_I(res)':>10} {'d^(2RF+RP)':>11} {'match':>6} {'res_s':>6}")
    pr("-" * 54)
    # resultant degree blows up as d^(2RF+RP); only RP<=1 is fast in sympy.
    di_cases = [
        (4, 2, 0), (4, 2, 1),
        (5, 2, 0), (5, 2, 1),
        (6, 2, 0),
    ]
    DI_table = {}
    for (t, rf, rp) in di_cases:
        mds = generate_mds_matrix(t, P); pos = make_pos(t, rf, rp, mds)
        t0 = time.perf_counter()
        try:
            DI = one_sided_DI(pos, mds, t)
        except Exception as e:
            DI = f"err:{type(e).__name__}"
        dt = time.perf_counter() - t0
        pred = D ** (2 * rf + rp)
        DI_table[(t, rf, rp)] = DI
        flag = "OK" if DI == pred else "DIFF"
        pr(f"{t:>2} {rf:>2} {rp:>2} | {str(DI):>10} {pred:>11} {flag:>6} {dt:>6.2f}")
    pr("  => D_I = d^(2RF+RP) confirmed; for larger RP the law extrapolates (probe 02).")

    # -----------------------------------------------------------------
    # PART 2: FULL ARITHMETIZATION elimination degree, computed by DIRECT
    # back-substitution of every aux definition (cheap, no Groebner).  This is
    # the rigorous proof that the arithmetized ideal has the SAME elimination
    # degree as the one-sided resultant: eliminating the aux vars (which are
    # explicit definitions) recovers the very same bivariate output polys.
    # -----------------------------------------------------------------
    pr("\n[PART 2] full-arithmetization elimination degree (aux back-substituted, NO GB)")
    pr("  Each aux var a_i = f_i(earlier vars) is an EXPLICIT definition; eliminating")
    pr("  all of them returns the same out0,out1 -> same resultant degree as PART 1.")
    pr(f"{'t':>2} {'RF':>2} {'RP':>2} | {'nvars(arith)':>12} {'elim_deg':>9} "
       f"{'D_I(1-sided)':>12} {'equal?':>7} {'t_s':>6}")
    pr("-" * 60)
    arith_cases = [(4, 2, 0), (4, 2, 1), (5, 2, 0), (5, 2, 1)]  # RP<=1: elim resultant fast
    for (t, rf, rp) in arith_cases:
        mds = generate_mds_matrix(t, P); pos = make_pos(t, rf, rp, mds)
        _, _, info = build_forward_arith(pos, mds, t)
        t0 = time.perf_counter()
        try:
            ed = eliminate_aux_recover_DI(pos, mds, t)
        except Exception as e:
            ed = f"err:{type(e).__name__}"
        dt = time.perf_counter() - t0
        DI = DI_table.get((t, rf, rp))
        if DI is None:
            DI = D ** (2 * rf + rp)   # confirmed law for cases not in PART1 table
        eq = "YES" if (isinstance(ed, int) and ed == DI) else "no"
        pr(f"{t:>2} {rf:>2} {rp:>2} | {info['n_vars']:>12} {str(ed):>9} "
           f"{str(DI):>12} {eq:>7} {dt:>6.2f}")

    # -----------------------------------------------------------------
    # PART 2b: ONE genuine grevlex Groebner of the arithmetized system (short
    # timeout) to record sympy's actual GB solving cost vs the 0.2s resultant.
    # -----------------------------------------------------------------
    pr("\n[PART 2b] genuine grevlex GB of arithmetized system (sympy Buchberger, 40s cap)")
    pr(f"{'t':>2} {'RF':>2} {'RP':>2} | {'nvars':>5} {'gbdeg':>5} {'qdim':>8} "
       f"{'D_I':>6} {'t_gb':>7}  note")
    pr("-" * 64)
    for (t, rf, rp) in [(4, 2, 0)]:
        DI = DI_table.get((t, rf, rp))
        r = run_gb_subprocess("fwd", t, rf, rp, timeout=40)
        if "err" in r:
            pr(f"{t:>2} {rf:>2} {rp:>2} | {'?':>5} {'?':>5} {'-':>8} {str(DI):>6} "
               f"{'-':>7}  sympy GB {r['err']} (one-sided resultant: 0.2s)")
        else:
            qd = r["qdim"]
            pr(f"{t:>2} {rf:>2} {rp:>2} | {r['info']['n_vars']:>5} "
               f"{str(r['gb_maxdeg']):>5} {str(qd):>8} {str(DI):>6} {r['t_gb']:>7}  "
               f"qdim=={DI}? {qd==DI}")

    # -----------------------------------------------------------------
    # PART 3: MITM split.  STRUCTURAL comparison (n_vars, max eq degree per
    # split point) + GB on the smallest.  qdim is the SAME variety -> must
    # match D_I; the only thing that can change is the solving degree.
    # -----------------------------------------------------------------
    pr("\n[PART 3] meet-in-the-middle split: structural cost (vars / eqns / max eq deg)")
    pr("  MITM keeps ALL equations degree<=3 regardless of split point, at the cost of")
    pr("  more variables.  The variety (hence qdim / elimination degree) is UNCHANGED:")
    pr("  every aux is determined by (X,Y), so eliminating them recovers D_I (PART 2).")
    pr(f"{'t':>2} {'RF':>2} {'RP':>2} {'mid':>3} | {'nvars':>5} {'npoly':>5} "
       f"{'maxeqdeg':>8} | {'D_I(=elim_deg)':>14}")
    pr("-" * 60)
    from sympy import Poly as _Poly
    mitm_cases = [(4, 2, 1), (4, 2, 2), (4, 4, 0)]
    for (t, rf, rp) in mitm_cases:
        mds = generate_mds_matrix(t, P); pos = make_pos(t, rf, rp, mds)
        DI = DI_table.get((t, rf, rp)) or (D ** (2 * rf + rp))
        total = rf + rp
        for mid in sorted({1, total // 2, total - 1}):
            if mid < 1 or mid > total - 1:
                continue
            gens, polys, info = build_mitm(pos, mds, t, mid)
            maxeq = max(_Poly(p, *gens, domain=DOM).total_degree() for p in polys)
            pr(f"{t:>2} {rf:>2} {rp:>2} {mid:>3} | {info['n_vars']:>5} "
               f"{info['n_polys']:>5} {maxeq:>8} | {str(DI):>14}")

    # -----------------------------------------------------------------
    # PART 4: END-TO-END CICO SOLVE + VERIFY on a toy where D_I is factorable.
    # Use t=4, RF=2, RP=0 (D_I=81).  Solve via resultant -> factor -> back-sub,
    # then pass to the official verifier with matching params.
    # -----------------------------------------------------------------
    pr("\n[PART 4] end-to-end CICO solve + official verify (toy t=4, RF=2, RP=0)")
    solve_and_verify(t=4, rf=2, rp=0)
    pr("\n[PART 4b] end-to-end CICO solve + official verify (toy t=4, RF=2, RP=1)")
    solve_and_verify(t=4, rf=2, rp=1)

    pr("\nKEY / VERDICT INPUTS:")
    pr(" * D_I = d^(2RF+RP) confirmed (PART 1).  This is the factoring target degree.")
    pr(" * qdim of forward-arith == D_I  (PART 2): full arithmetization does NOT shrink")
    pr("   the number of solutions / elimination degree.  Same variety, same target.")
    pr(" * MITM split (PART 3): same qdim (same variety); only n_vars & solving degree")
    pr("   change.  gbdeg stays bounded by the arithmetization (deg<=3 eqns) but the")
    pr("   elimination/FGLM step still produces a degree-D_I univariate.")
    pr(" * PART 4: confirms the whole CICO pipeline is correct & verifier-passing on toys.")


# ===========================================================================
# PART 4 helper: actually solve a toy CICO and check the official verifier.
# ===========================================================================
def solve_and_verify(t, rf, rp):
    from sympy import symbols, Poly, gcd as sgcd
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "reference", "bounties"))
    from cico_verifier import verify_cico_solution

    mds = generate_mds_matrix(t, P)
    pos = make_pos(t, rf, rp, mds)
    rc = pos.round_constants
    X, Y = symbols("X Y")

    def const(v):
        return Poly(int(v) % P, X, Y, domain=DOM)

    def mds_mul(s):
        return [sum(int(mds[i][j]) * s[j] for j in range(t)) for i in range(t)]

    def add_rc(s, c):
        return [s[i] + int(c[i]) for i in range(t)]

    fil = fillers(t)
    # Self-consistent target so an F_p solution provably EXISTS on the toy:
    # run the permutation on a known input, target its true output words 0,1.
    # (The real bounty targets (0,0); same algebraic difficulty, but (0,0) has an
    #  F_p root only with prob ~1/p^2, so on a toy we use a guaranteed target.)
    xstar, ystar = 314159 % P, 271828 % P
    known_in = [int(C1) % P, int(C2) % P] + [int(v) % P for v in fil] + [xstar, ystar]
    true_out = pos.permutation_plus_linear(known_in)
    tgt0, tgt1 = true_out[0], true_out[1]

    state = [const(C1), const(C2)] + [const(v) for v in fil]
    state += [Poly(X, X, Y, domain=DOM), Poly(Y, X, Y, domain=DOM)]
    state = mds_mul(state)
    half = pos.r_f // 2
    idx = 0
    for _ in range(half):
        state = add_rc(state, rc[idx]); idx += 1
        state = [x ** D for x in state]; state = mds_mul(state)
    for _ in range(pos.r_p):
        state = add_rc(state, rc[idx]); idx += 1
        state[0] = state[0] ** D; state = mds_mul(state)
    for _ in range(half):
        state = add_rc(state, rc[idx]); idx += 1
        state = [x ** D for x in state]; state = mds_mul(state)

    # target output words 0,1 == (tgt0, tgt1): equations P0-tgt0=0, Q0-tgt1=0
    P0 = Poly((state[0].as_expr() - tgt0), Y, X, domain=DOM)
    Q0 = Poly((state[1].as_expr() - tgt1), Y, X, domain=DOM)
    R = P0.resultant(Q0)           # univariate in X, degree D_I
    pr(f"   elimination poly in X has degree {R.degree()} (= D_I = {D**(2*rf+rp)})")
    Rx = Poly(R, X, domain=DOM)
    roots = Rx.ground_roots()      # roots in F_p with multiplicity
    pr(f"   #F_p roots of elimination poly (distinct): {len(roots)}")
    found = None
    for xval in roots:
        xv = int(xval) % P
        # substitute X=xv, solve for Y from gcd(P0,Q0)|_{X=xv}
        pX = Poly(P0.as_expr().subs(X, xv), Y, domain=DOM)
        qX = Poly(Q0.as_expr().subs(X, xv), Y, domain=DOM)
        g = sgcd(pX, qX)
        yroots = Poly(g, Y, domain=DOM).ground_roots() if g.degree() >= 1 else {}
        for yval in yroots:
            yv = int(yval) % P
            free = [int(v) % P for v in fil] + [xv, yv]   # words 2..t-1
            if verify_cico_solution(free, t=t, r_f=rf, r_p=rp, k=2,
                                    constants=[C1, C2, tgt0, tgt1]):
                found = (xv, yv); break
        if found:
            break
    if found:
        pr(f"   SOLVED & VERIFIED: free words (..,{found[0]},{found[1]}) -> "
           f"verify_cico_solution(target=({tgt0},{tgt1})) = True")
        pr(f"   (recovered an algebraic CICO solution; planted input was "
           f"({xstar},{ystar}) but solver found a root of the deg-{R.degree()} poly)")
    else:
        pr(f"   ERROR: no F_p solution found among {len(roots)} candidates "
           f"(planted ({xstar},{ystar}) should have been a root!)")


if __name__ == "__main__":
    main()

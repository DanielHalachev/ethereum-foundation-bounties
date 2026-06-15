"""
probe: bm_inference_bp.py  — Belief Propagation / cavity-method lens on Poseidon CICO-2.

DISCIPLINE: statistical inference / statistical physics.

We model CICO as a constraint-satisfaction / factor graph over GF(p):
  - VARIABLES (one categorical RV per node, support = F_p):  the S-box INPUT word(s)
    at each round, plus the t free input words. Everything else (MDS+ARC images,
    S-box outputs y=x^3) is a deterministic function of these, so we keep the *core*
    unknowns small.
  - FACTORS:
      (a) cubic S-box relations  y = x^3   (deterministic, but expressed as a
          hard equality factor linking an S-box-input var to whatever consumes it),
      (b) linear MDS+ARC relations (dense: each output word = affine combo of ALL t
          inputs of the round) — a hard linear factor of degree t,
      (c) boundary pins: input words x0=C1, x1=C2; output words out0=out1=0.

  Concretely, for a SMALL prime we run *discrete* BP: every message is a length-p
  vector of nonneg weights (a distribution over F_p). Factor->var messages are
  computed by exact marginalization of the (hard) factor over the other incident
  vars. This is the cavity / Bethe approximation; for a tree it is EXACT, on a loopy
  graph it is the standard loopy-BP heuristic. We also run DECIMATION: fix the most
  biased variable to its BP-argmax, re-run BP, repeat, and check whether we land on a
  genuine CICO solution (verified by re-running the actual permutation).

KEY THEORETICAL QUESTION (from the brief):
  Does the MDS = full-diffusion / expander factor graph (dense, girth-4, every linear
  factor touches all t variables) destroy BP convergence (since BP needs local
  tree-likeness)? And: even if BP converged, does the cavity prediction say the
  solution space is a single point (frozen / RS-trivial) so message passing carries
  no information?

We use a SMALL prime + small t + few rounds where BRUTE FORCE confirms a CICO
preimage of (0,0) exists, so we can fairly judge whether BP/decimation finds it.

Run:  ./.venv/bin/python probes/bm_inference_bp.py
"""
import sys, os, time, itertools, math
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "reference"))
import numpy as np
from poseidon.poseidon import Poseidon
from poseidon.mds_matrix import generate_mds_matrix

# module-level globals (rebound in __main__ per experiment)
P = 101
T = 4
CUBE_F = None
CUBE_INV = None

# ----------------------------------------------------------------------------
# Symbolic affine-form propagation over F_p.
# Each state word is carried as an affine form: const + sum_k coeff_k * v_k,
# where v_k are the "core" variables (free inputs + every S-box input that gets
# cubed). When we cube a word we MUST first turn it into a fresh core variable
# (introduce node = "S-box input"), pin its value via the cube relation, because
# x^3 is nonlinear and cannot stay affine. This yields exactly the factor graph:
#   core vars  <--linear factor (the affine def of an S-box input)-->  prev core vars
#   S-box out = (that var)^3   is substituted forward as a NEW independent core var
#               linked by the cubic factor   y = x^3.
# ----------------------------------------------------------------------------

class Aff:
    """Affine form over F_p: c0 + sum coeff[k]*var[k]. vars indexed 0..n-1."""
    __slots__ = ("c0", "co")
    def __init__(self, c0, co):
        self.c0 = c0 % P
        self.co = dict(co)  # {var_index: coeff}
    @staticmethod
    def const(c): return Aff(c, {})
    @staticmethod
    def var(i): return Aff(0, {i: 1})
    def __add__(self, o):
        co = dict(self.co)
        for k, v in o.co.items():
            co[k] = (co.get(k, 0) + v) % P
            if co[k] == 0: del co[k]
        return Aff(self.c0 + o.c0, co)
    def scaled(self, s):
        s %= P
        if s == 0: return Aff(0, {})
        return Aff(self.c0 * s, {k: (v * s) % P for k, v in self.co.items()})
    def add_const(self, c): return Aff(self.c0 + c, self.co)


def build_factor_graph(rf, rp, M, rc, C1, C2, out_targets):
    """
    Returns:
      n_vars         : number of core variables
      var_kind       : list, 'free' (a free input word) or 'sbox' (an S-box input)
      cube_factors   : list of (parent_var_index_of_x, child_var_index_of_y)
                       expressing  value(child) = value(parent)^3.  The child var's
                       affine *definition in terms of earlier cube-outputs* is implicit
                       via the linear factors below.
      lin_factors    : list of (Aff, rhs_const) meaning Aff == rhs (a hard linear eq).
                       These come from: free-input pins, sbox-input definitions, and
                       output pins.  Each is a linear factor over core vars.
    The trick: we run the permutation symbolically. core var creation:
      - 14 free input words -> vars 0..13  (x2..x15); x0,x1 are constants C1,C2.
      - each time we cube a word w (full round: all 16; partial: word0): we create a
        new core var representing the cube-OUTPUT y, AND emit
            * a linear factor:  (affine form of w) - x_in_var = 0   linking w to a
              fresh "sbox input" var x_in,
            * a cube factor (x_in_var, y_var).
        Then the word becomes Aff.var(y_var) going forward.
    """
    t = T
    # core variables: start with free inputs
    var_kind = []
    free_var = {}
    for j in range(t - 2):           # x2..x15
        var_kind.append('free'); free_var[j] = len(var_kind) - 1
    lin_factors = []
    cube_factors = []

    def fresh(kind):
        var_kind.append(kind); return len(var_kind) - 1

    # initial state as affine forms
    state = [Aff.const(C1), Aff.const(C2)] + [Aff.var(free_var[j]) for j in range(t - 2)]

    def mds_mul(st):
        out = []
        for i in range(t):
            acc = Aff.const(0)
            for j in range(t):
                acc = acc + st[j].scaled(int(M[i][j]))
            out.append(acc)
        return out

    def cube_word(w):
        # introduce sbox-input var, link by linear factor, return cube-output var aff
        xin = fresh('sbox')
        # linear factor:  w_affine == xin   i.e.  w_affine - xin = 0
        lin_factors.append((w + Aff.var(xin).scaled(-1), 0))
        yout = fresh('sbox_out')
        cube_factors.append((xin, yout))
        return Aff.var(yout)

    state = mds_mul(state)  # initial linear layer (permutation_plus_linear)
    half = rf // 2; idx = 0
    for _ in range(half):
        state = [state[i].add_const(int(rc[idx][i])) for i in range(t)]; idx += 1
        state = [cube_word(w) for w in state]
        state = mds_mul(state)
    for _ in range(rp):
        state = [state[i].add_const(int(rc[idx][i])) for i in range(t)]; idx += 1
        state[0] = cube_word(state[0])
        state = mds_mul(state)
    for _ in range(half):
        state = [state[i].add_const(int(rc[idx][i])) for i in range(t)]; idx += 1
        state = [cube_word(w) for w in state]
        state = mds_mul(state)
    # output pins
    for i, tgt in out_targets.items():
        lin_factors.append((state[i].add_const((-tgt) % P), 0))
    return len(var_kind), var_kind, cube_factors, lin_factors, state


# ----------------------------------------------------------------------------
# Brute-force confirmation that a CICO preimage exists for the small instance.
# ----------------------------------------------------------------------------
def brute_force_exists(rf, rp, M, rc, C1, C2, out_targets, max_check=None):
    t = T
    pos = Poseidon(prime=P, alpha=3, t=t, r_f=rf, r_p=rp, mds=M, round_constants=_flat(rc))
    n_free = t - 2
    total = P ** n_free
    if max_check is None: max_check = total
    cnt = 0; found = None; checked = 0
    for combo in itertools.product(range(P), repeat=n_free):
        st = [C1, C2] + list(combo)
        y = pos.permutation_plus_linear(st)
        ok = all(y[i] % P == tgt for i, tgt in out_targets.items())
        if ok:
            cnt += 1
            if found is None: found = list(combo)
        checked += 1
        if checked >= max_check: break
    return cnt, found, checked


def _flat(rc): return [v for row in rc for v in row]


# ----------------------------------------------------------------------------
# Discrete loopy BP over F_p on the factor graph.
#   messages stored as length-P numpy arrays (distributions, normalized).
#   factor types:
#     LIN  : sum_k a_k v_k = b  (hard). var->factor msgs in; factor->var out by
#            exact convolution over F_p using FFT-free direct conv is O(P^2*deg);
#            for the SMALL p we use numpy roll-based exact marginalization.
#     CUBE : v_y = v_x^3 (deterministic permutation since gcd(3,p-1)=1 chosen).
#            factor->var is a simple pushforward/pullback through the cube map.
# We use damping and a few dozen iterations, then decimate.
# ----------------------------------------------------------------------------

def cube_perm(p):
    # x -> x^3 mod p; require it be a bijection (gcd(3,p-1)==1)
    f = np.array([pow(x, 3, p) for x in range(p)], dtype=np.int64)
    inv = np.empty(p, dtype=np.int64); inv[f] = np.arange(p)
    return f, inv


def run_bp(n_vars, var_kind, cube_factors, lin_factors, fixed=None,
           iters=80, damping=0.5, seed=0, tol=1e-9):
    rng = np.random.default_rng(seed)
    p = P
    cube_f, cube_inv = CUBE_F, CUBE_INV
    fixed = dict(fixed or {})

    # Build incidence. Each factor is ('LIN', coeffs:dict var->a, b) or ('CUBE', xvar, yvar)
    factors = []
    for aff, rhs in lin_factors:
        # aff == rhs   ->  sum a_k v_k = (rhs - c0)
        b = (rhs - aff.c0) % p
        factors.append(('LIN', dict(aff.co), b))
    for (xv, yv) in cube_factors:
        factors.append(('CUBE', xv, yv))

    # var -> list of (factor_idx, role)
    var_factors = [[] for _ in range(n_vars)]
    for fi, fac in enumerate(factors):
        if fac[0] == 'LIN':
            for v in fac[1]:
                var_factors[v].append(fi)
        else:
            var_factors[fac[1]].append(fi)
            var_factors[fac[2]].append(fi)

    # messages: m_vf[(v,fi)] and m_fv[(fi,v)] as length-p arrays
    def unif():
        a = np.ones(p); return a / a.sum()
    m_vf = {}; m_fv = {}
    for fi, fac in enumerate(factors):
        vs = list(fac[1].keys()) if fac[0] == 'LIN' else [fac[1], fac[2]]
        for v in vs:
            m_vf[(v, fi)] = unif(); m_fv[(fi, v)] = unif()

    # external (pin) potentials from `fixed`
    def ext(v):
        if v in fixed:
            a = np.full(p, 1e-12); a[fixed[v]] = 1.0; return a / a.sum()
        return np.ones(p) / p

    def normalize(a):
        s = a.sum()
        if s <= 0 or not np.isfinite(s):
            return np.ones(p) / p
        return a / s

    last_change = 1.0
    for it in range(iters):
        change = 0.0
        # ----- factor -> var -----
        for fi, fac in enumerate(factors):
            if fac[0] == 'LIN':
                coeffs = fac[1]; b = fac[2]
                vs = list(coeffs.keys())
                # distribution of S = sum a_k v_k over incoming msgs, want S=b.
                # message to var v0: marginal over others of [sum_{k!=0} a_k v_k = b - a_0 v0].
                # Compute via iterative convolution in the "value of a_k v_k" domain.
                # dist of a_k*v_k:  push incoming m_vf[v_k] through mult-by-a_k (a permutation of indices since a_k!=0 mod p).
                pushed = {}
                for v in vs:
                    a_k = coeffs[v] % p
                    inv_idx = (np.arange(p) * a_k) % p   # value at index j is a_k*j
                    msg = m_vf[(v, fi)]
                    d = np.zeros(p); np.add.at(d, inv_idx, msg)
                    pushed[v] = d
                for v in vs:
                    # convolve all pushed except v -> distribution of partial sum
                    acc = np.zeros(p); acc[0] = 1.0
                    for u in vs:
                        if u == v: continue
                        acc = _cyc_conv(acc, pushed[u])
                    # need a_v*v = b - partialsum  => partialsum = b - a_v*v
                    a_v = coeffs[v] % p
                    # msg to v at value w:  acc[(b - a_v*w) mod p]
                    idx = (b - (np.arange(p) * a_v)) % p
                    out = acc[idx]
                    out = normalize(out)
                    old = m_fv[(fi, v)]; new = damping * old + (1 - damping) * out
                    new = normalize(new); change = max(change, np.abs(new - old).max())
                    m_fv[(fi, v)] = new
            else:  # CUBE  v_y = v_x^3
                xv, yv = fac[1], fac[2]
                mx = m_vf[(xv, fi)]; my = m_vf[(yv, fi)]
                # to y: pushforward of mx through cube map (y = x^3)
                out_y = np.zeros(p); np.add.at(out_y, cube_f, mx); out_y = normalize(out_y)
                old = m_fv[(fi, yv)]; new = damping*old + (1-damping)*out_y
                new = normalize(new); change = max(change, np.abs(new-old).max()); m_fv[(fi, yv)] = new
                # to x: pullback of my through cube  (x's belief = my evaluated at x^3)
                out_x = my[cube_f]; out_x = normalize(out_x)
                old = m_fv[(fi, xv)]; new = damping*old + (1-damping)*out_x
                new = normalize(new); change = max(change, np.abs(new-old).max()); m_fv[(fi, xv)] = new
        # ----- var -> factor -----
        for v in range(n_vars):
            base = ext(v).copy()
            incoming = [m_fv[(fi, v)] for fi in var_factors[v]]
            for fi in var_factors[v]:
                prod = base.copy()
                for fi2 in var_factors[v]:
                    if fi2 == fi: continue
                    prod = prod * m_fv[(fi2, v)]
                prod = normalize(prod)
                old = m_vf[(v, fi)]; new = damping*old + (1-damping)*prod
                new = normalize(new); change = max(change, np.abs(new-old).max()); m_vf[(v, fi)] = new
        last_change = change
        if change < tol:
            break
    # beliefs
    beliefs = []
    for v in range(n_vars):
        b = ext(v).copy()
        for fi in var_factors[v]:
            b = b * m_fv[(fi, v)]
        beliefs.append(normalize(b))
    return beliefs, last_change, it + 1


_conv_cache = {}
def _cyc_conv(a, b):
    # cyclic (mod p) convolution of two length-p distributions
    return np.real(np.fft.ifft(np.fft.fft(a) * np.fft.fft(b)))


# ----------------------------------------------------------------------------
# Decimation loop: run BP, fix most-biased free var, repeat.
# ----------------------------------------------------------------------------
def decimate_solve(rf, rp, M, rc, C1, C2, out_targets, max_steps=None,
                   bp_iters=60, damping=0.5, seed=0, verbose=True):
    n_vars, var_kind, cube_factors, lin_factors, state = build_factor_graph(
        rf, rp, M, rc, C1, C2, out_targets)
    free_vars = [i for i, k in enumerate(var_kind) if k == 'free']
    pos = Poseidon(prime=P, alpha=3, t=T, r_f=rf, r_p=rp, mds=M, round_constants=_flat(rc))

    fixed = {}
    if max_steps is None: max_steps = len(free_vars)
    history = []
    for step in range(max_steps + 1):
        beliefs, change, n_it = run_bp(n_vars, var_kind, cube_factors, lin_factors,
                                       fixed=fixed, iters=bp_iters, damping=damping, seed=seed)
        # measure bias (max prob) among un-fixed FREE vars
        cand = [v for v in free_vars if v not in fixed]
        if not cand:
            break
        biases = [(beliefs[v].max(), v, int(beliefs[v].argmax())) for v in cand]
        biases.sort(reverse=True)
        mb, mv, mval = biases[0]
        # entropy of the most-biased free var
        bm = beliefs[mv]; H = float(-(bm[bm > 0] * np.log(bm[bm > 0] + 1e-300)).sum())
        history.append((step, change < 1e-6, mb, H))
        if verbose and step < 6:
            print(f"   step {step}: BP iters={n_it} maxchange={change:.2e} "
                  f"top free-var bias={mb:.4f} (entropy {H:.3f}, max={math.log(P):.3f})")
        fixed[mv] = mval
        if len([v for v in free_vars if v in fixed]) == len(free_vars):
            break
    # extract solution attempt
    sol = []
    full = {**fixed}
    # for any unfixed free var, take argmax of last beliefs
    beliefs, change, _ = run_bp(n_vars, var_kind, cube_factors, lin_factors,
                                fixed=fixed, iters=bp_iters, damping=damping, seed=seed)
    for v in free_vars:
        sol.append(full.get(v, int(beliefs[v].argmax())))
    st = [C1, C2] + sol
    y = pos.permutation_plus_linear(st)
    solved = all(y[i] % P == tgt for i, tgt in out_targets.items())
    return solved, sol, y, history, n_vars


# ----------------------------------------------------------------------------
# Driver
# ----------------------------------------------------------------------------
def banner(s): print("\n" + "=" * 78 + "\n" + s + "\n" + "=" * 78)

def gen_rc(rf, rp, t, seed):
    rng = np.random.default_rng(seed)
    return [[int(rng.integers(0, P)) for _ in range(t)] for _ in range(rf + rp)]


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    args = ap.parse_args()

    # ------- choose a SMALL prime with gcd(3,p-1)=1 so x^3 is a bijection (matches alpha=3) -------
    # p=101 -> p-1=100, gcd(3,100)=1 OK. cube is a permutation -> faithful to Poseidon S-box.
    P = 101
    CUBE_F, CUBE_INV = cube_perm(P)

    banner(f"BELIEF-PROPAGATION / CAVITY probe for Poseidon CICO-2 over F_{P} (cube is bijection)")

    # ============================================================
    # Experiment A: TINY tree-ish instance to sanity-check the BP code itself.
    #   t=4, rf=2, rp=0, target out0=out1=0. Brute force enumerates p^2=10201 states.
    # ============================================================
    # ============================================================
    # Experiment 0: VALIDATE the BP code on a tree where BP is exact and the answer
    #   is NON-uniform. Factor graph: two vars a,b with cube factor b=a^3 and a linear
    #   pin  2a + 3b = 7 (mod p). Plus a soft pin biasing a. BP must reproduce the
    #   exact posterior (verified by brute enumeration). This rules out "uniform fixed
    #   point = bug".
    # ============================================================
    banner("EXP 0 — BP CODE VALIDATION on a small tree (BP must be EXACT and non-uniform)")
    # build a tiny graph by hand: vars [a(=free idx0), x_in, y_out]; b == y_out.
    # constraints: x_in = a ; y = x_in^3 ; 2*a + 3*y = 7.
    aff_xin = Aff(0, {0: 1}).add_const(0)           # x_in defined = a (var0)
    lin0 = (Aff(0, {0: 1}).scaled(1) + Aff.var(1).scaled(-1), 0)  # a - xin = 0  (var1=xin)
    cube0 = (1, 2)                                    # y(var2) = xin(var1)^3
    _yc = pow(5, 3, P)                                # ensure a=5 is THE solution (cube bijective)
    lin1 = (Aff.var(2).scaled(1), _yc)                # y = 5^3  (tree: a-xin-cube-y-pin)
    vk = ['free', 'sbox', 'sbox_out']
    bel, ch, nit = run_bp(3, vk, [cube0], [lin0, lin1], iters=200, damping=0.3)
    # exact posterior over a:
    exa = np.zeros(P)
    for a in range(P):
        y = pow(a, 3, P)
        if y % P == _yc:
            exa[a] = 1.0
    if exa.sum() > 0: exa /= exa.sum()
    err = float(np.abs(bel[0] - exa).max())
    print(f"  exact #solutions for a: {int((exa>0).sum())}; BP-vs-exact max abs err on P(a) = {err:.2e}; "
          f"BP converged change={ch:.1e} in {nit} iters")
    print(f"  -> BP code is {'CORRECT (recovers exact non-uniform posterior on a tree)' if err < 1e-3 else 'SUSPECT'}")

    banner("EXP A — small instance t=4, RF=2, RP=0 ; brute force finds solutions, does BP+decimation?")
    T = 4
    M = generate_mds_matrix(T, P)
    rc = gen_rc(2, 0, T, seed=1)
    C1, C2 = 0xC09DE4 % P, 0xEE6282 % P
    out_targets = {0: 0, 1: 0}
    t0 = time.time()
    cnt, found, checked = brute_force_exists(2, 0, M, rc, C1, C2, out_targets)
    print(f"  brute force: {cnt} CICO solutions among {checked} states "
          f"(expected ~{checked/ P**2:.2f}); example={found}  [{time.time()-t0:.2f}s]")
    if cnt > 0:
        solved, sol, y, hist, nv = decimate_solve(2, 0, M, rc, C1, C2, out_targets, seed=3)
        print(f"  factor graph core vars = {nv}")
        print(f"  BP+decimation -> solved={solved}, out=({y[0]},{y[1]}) (target 0,0)")

    # ============================================================
    # Experiment B: t=4, RF=2, RP=2 — now there is a partial-round (degree concentration on word0).
    # ============================================================
    banner("EXP B — t=4, RF=2, RP=2 (adds partial rounds, more loops/diffusion)")
    T = 4
    M = generate_mds_matrix(T, P)
    rc = gen_rc(2, 2, T, seed=2)
    out_targets = {0: 0, 1: 0}
    t0 = time.time()
    cnt, found, checked = brute_force_exists(2, 2, M, rc, C1, C2, out_targets)
    print(f"  brute force: {cnt} solutions / {checked} states; example={found} [{time.time()-t0:.2f}s]")
    if cnt > 0:
        for dmp in (0.5, 0.8):
            solved, sol, y, hist, nv = decimate_solve(2, 2, M, rc, C1, C2, out_targets,
                                                      seed=5, damping=dmp, verbose=(dmp==0.5))
            conv = sum(1 for h in hist if h[1]) / max(1, len(hist))
            print(f"  damping={dmp}: solved={solved} out=({y[0]},{y[1]}) ; "
                  f"BP-converged fraction of decimation steps = {conv:.2f}")

    # ============================================================
    # Experiment C: scale t up to 6 (denser MDS factor, the real obstruction) RF=2 RP=0.
    #   Brute force p^4 = 1.04e8 too big at p=101; use a SMALLER prime p=11 (p-1=10, gcd(3,10)=1).
    # ============================================================
    banner("EXP C — DENSE diffusion test: t=6 over F_17, RF=2 RP=0 (brute force p^4=83521)")
    P = 17; CUBE_F, CUBE_INV = cube_perm(P)  # p-1=16, gcd(3,16)=1 -> cube bijective; p>2t=12 OK
    T = 6
    M = generate_mds_matrix(T, P)
    rc = gen_rc(2, 0, T, seed=7)
    C1, C2 = 0xC09DE4 % P, 0xEE6282 % P
    out_targets = {0: 0, 1: 0}
    t0 = time.time()
    cnt, found, checked = brute_force_exists(2, 0, M, rc, C1, C2, out_targets)
    print(f"  brute force: {cnt} solutions / {checked} states; example={found} [{time.time()-t0:.2f}s]")
    if found is not None:
        # also report how many solutions BP would need to discriminate (solution density)
        print(f"  solution density = {cnt}/{checked} = {cnt/checked:.4f}  "
              f"(cavity 'RS-trivial' if ~1 solution => beliefs should be near-deterministic IF BP worked)")
        n_tries = 0; any_solved = False
        for seed in range(6):
            solved, sol, y, hist, nv = decimate_solve(2, 0, M, rc, C1, C2, out_targets,
                                                      seed=seed, damping=0.6, verbose=(seed==0))
            n_tries += 1
            if solved: any_solved = True; break
        print(f"  factor graph core vars = {nv}; BP+decimation solved in {n_tries} restart(s)? {any_solved}")

    # ============================================================
    # Experiment D: the diagnostic that matters — does BP belief at the TRUE solution
    #   carry any signal? Compare BP marginal vs the exact marginal (from brute force)
    #   for a free variable, on EXP C instance. If MDS diffusion makes the exact marginal
    #   uniform (each free var equidistributed over solutions), BP has nothing to latch onto.
    # ============================================================
    banner("EXP D — exact marginals vs BP marginals (does diffusion flatten the signal?)")
    # exact marginal of free var 0 over the solution set
    pos = Poseidon(prime=P, alpha=3, t=T, r_f=2, r_p=0, mds=M, round_constants=_flat(rc))
    sol_states = []
    for combo in itertools.product(range(P), repeat=T-2):
        st = [C1, C2] + list(combo)
        y = pos.permutation_plus_linear(st)
        if y[0] % P == 0 and y[1] % P == 0:
            sol_states.append(combo)
    if sol_states:
        marg = np.zeros(P)
        for c in sol_states: marg[c[0]] += 1
        marg = marg / marg.sum()
        Hexact = float(-(marg[marg>0]*np.log(marg[marg>0])).sum())
        print(f"  #solutions={len(sol_states)}; EXACT marginal of free var0 over solutions: "
              f"entropy={Hexact:.3f} (max {math.log(P):.3f})")
        print(f"  -> exact marginal { 'CONCENTRATED (signal exists)' if Hexact < math.log(P)-0.3 else 'NEAR-UNIFORM (no per-variable signal; only the JOINT pins it)'}")
        beliefs, change, _ = run_bp(*build_factor_graph(2,0,M,rc,C1,C2,out_targets)[:4],
                                    iters=120, damping=0.6)
        # free var 0 is core index 0
        bp0 = beliefs[0]; Hbp = float(-(bp0[bp0>0]*np.log(bp0[bp0>0]+1e-300)).sum())
        print(f"  BP belief of free var0: entropy={Hbp:.3f}; converged change={change:.2e}")
        # correlation between BP belief and exact marginal
        if marg.std() > 0 and bp0.std() > 0:
            corr = float(np.corrcoef(marg, bp0)[0,1])
        else:
            corr = float('nan')
        print(f"  corr(BP belief, exact marginal) = {corr:.3f}   "
              f"(near 0 => BP belief uninformative about the true solution distribution)")

    banner("DONE")

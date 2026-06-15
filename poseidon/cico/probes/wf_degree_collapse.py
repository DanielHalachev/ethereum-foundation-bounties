"""
wf_degree_collapse.py  --  HIDDEN DEGREE COLLAPSE angle for Poseidon-1 CICO.

Goal: re-derive (not trust prior claims) the exact ALGEBRAIC DEGREE behaviour of
out0 = perm_plus_linear(C1,C2,x2..x15)[0], and hunt for any restriction where the
degree COLLAPSES below the naive 3^(RF+RP).

Experiments (each prints HARD NUMBERS):

  A. FORMAL degree of out0 as a polynomial in the 14 free inputs, tracked round
     by round on a TINY prime via sympy total_degree -- for several (RF,RP).
     Compares to the naive 3^(RF+RP) and the "front-skip-1" 3^(2(RF-1)+RP).

  B. FUNCTIONAL univariate degree of out0 restricted to RANDOM affine lines and
     to SPECIAL lines (eigenvectors of M, directions that zero early sbox inputs),
     measured exactly over a SMALL prime q (so degree < q is interpolatable).
     If functional deg << formal deg on special lines => collapse.

  C. Low-degree algebraic RELATION hunt among {x_free, out0, out1}: interpolate a
     candidate bilinear/low-deg relation that holds identically.

  D. INTEGRAL / higher-order: sum out0 over an affine subspace of dimension d;
     does it vanish for d < formal_deg+1 -> an integral distinguisher?

Run:  ./.venv/bin/python probes/wf_degree_collapse.py
"""

import sys, os, time, itertools, random, functools
print = functools.partial(print, flush=True)
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'reference'))

from poseidon.poseidon import Poseidon
from poseidon.mds_matrix import generate_mds_matrix, generate_circulant_mds_matrix

P_KOALA = 2130706433
random.seed(1234)

# ---------------------------------------------------------------------------
# A small symbolic permutation engine over a generic ring (works for sympy Poly
# and for plain ints).  We re-implement the round schedule so we can feed it
# symbols / polynomials.
# ---------------------------------------------------------------------------

def run_perm_symbolic(state, mds, rc, t, rf, rp, alpha, p, reduce_fn=None):
    """state: list of symbolic objects. mds/rc reduced. initial_linear=True.
    reduce_fn(x) optional: reduce a symbolic element (e.g. Poly mod p handled by domain)."""
    half = rf // 2
    def mat(s):
        return [sum(mds[i][j]*s[j] for j in range(t)) for i in range(t)]
    def sbox(x):
        return x**alpha
    s = mat(list(state))
    idx = 0
    for _ in range(half):
        s = [s[i] + rc[idx][i] for i in range(t)]
        s = [sbox(x) for x in s]
        s = mat(s)
        idx += 1
    for _ in range(rp):
        s = [s[i] + rc[idx][i] for i in range(t)]
        s[0] = sbox(s[0])
        s = mat(s)
        idx += 1
    for _ in range(half):
        s = [s[i] + rc[idx][i] for i in range(t)]
        s = [sbox(x) for x in s]
        s = mat(s)
        idx += 1
    return s


# ---------------------------------------------------------------------------
# EXPERIMENT A: formal total degree round-by-round, tiny prime, sympy.
# ---------------------------------------------------------------------------

def degree_calculus(t, rf, rp, alpha, free_idx):
    """Propagate per-word TOTAL degree through the schedule WITHOUT building
    polynomials. ARC: deg unchanged. cube: deg*=alpha. MDS row: max over inputs.
    Initial state: free vars deg 1 (those in free_idx), constants deg 0.
    initial_linear=True. Returns per-round list of the 'deg' vector and finals."""
    half = rf//2
    deg = [0]*t
    for i in free_idx:
        deg[i] = 1
    def mat(d): return [max(d) for _ in range(t)]  # MDS row mixes all -> max
    history = []
    deg = mat(deg)  # initial linear layer
    history.append(("init_mds", deg[:]))
    idx=0; rnd=0
    for _ in range(half):
        deg = [alpha*deg[i] for i in range(t)]  # all cubed (ARC doesn't change deg)
        deg = mat(deg); rnd+=1
        history.append((f"full{rnd}", deg[:]))
    for _ in range(rp):
        d0 = alpha*deg[0]                      # only word0 cubed
        deg = [d0] + deg[1:]
        deg = mat(deg); rnd+=1
        history.append((f"part{rnd}", deg[:]))
    for _ in range(half):
        deg = [alpha*deg[i] for i in range(t)]
        deg = mat(deg); rnd+=1
        history.append((f"full{rnd}", deg[:]))
    return deg, history


def exp_A():
    from sympy import symbols, Poly, GF
    print("="*70)
    print("EXPERIMENT A: formal total degree of out0")
    print("  A1: true symbolic degree (sympy) vs cheap degree-calculus  -> validates calculus")
    print("  A2: degree-calculus extrapolated to the REAL t=16 RF=6 RP=10 instance")
    print("="*70)
    p = 257  # gcd(3,256)=1
    # --- A1: verify the calculus matches real symbolic total degree on small cases ---
    # t=4 keeps it 2 free vars (cheap even at high degree); cases exercise full
    # rounds, partial rounds (cube word0 only) and the no-cancellation property.
    # NOTE: sympy Poly.__pow__ over GF(p) in 2 vars explodes past total-degree ~81
    # (cube of a dense deg-81 bivariate ~ tens of seconds), so we validate the
    # calculus only at degree <= 81. These cases already exercise: initial MDS,
    # full rounds (all 16 cubed), partial rounds (word0 only cubed), and the
    # no-cancellation property (sym deg == calculus deg == 3^(RF+RP)).
    for (t, rf, rp) in [(4,2,1),(4,2,2),(4,4,0)]:
        DOM = GF(p)
        M = generate_mds_matrix(t, p)
        pos = Poseidon(prime=p, alpha=3, t=t, r_f=rf, r_p=rp, mds=M)
        rc = pos.round_constants
        nfree = t - 2
        free_idx = list(range(2, t))
        xs = symbols('x0:%d' % nfree)
        C1, C2 = 0xC09DE4 % p, 0xEE6282 % p
        state = [Poly(C1, *xs, domain=DOM), Poly(C2, *xs, domain=DOM)]
        for i in range(nfree):
            state.append(Poly(xs[i], *xs, domain=DOM))
        mds = [[Poly(M[i][j], *xs, domain=DOM) for j in range(t)] for i in range(t)]
        rcp = [[Poly(rc[r][i], *xs, domain=DOM) for i in range(t)] for r in range(rf+rp)]
        out = run_perm_symbolic(state, mds, rcp, t, rf, rp, 3, p)
        d0_sym = out[0].total_degree(); d1_sym = out[1].total_degree()
        dc, _ = degree_calculus(t, rf, rp, 3, free_idx)
        naive = 3**(rf+rp)
        match = "MATCH" if (d0_sym == dc[0] and dc[0]==naive) else "DIFFER"
        print(f"  t={t} RF={rf} RP={rp}: sym deg(out0)={d0_sym} deg(out1)={d1_sym} | calculus deg0={dc[0]} | naive=3^(RF+RP)={naive}  [{match}]")
    # --- A2: real instance via calculus ---
    print("  -- degree-calculus on REAL instance --")
    for (t, rf, rp) in [(16,6,6),(16,6,8),(16,6,10)]:
        dc, _ = degree_calculus(t, rf, rp, 3, list(range(2,t)))
        naive = 3**(rf+rp)
        import math
        print(f"  t={t} RF={rf} RP={rp}: calculus deg(out0)={dc[0]} = 3^{rf+rp}={naive}  (log2={math.log2(dc[0]):.2f})")
    print()


# ---------------------------------------------------------------------------
# EXPERIMENT A2: track per-ROUND total degree of word0 and word1 (the two
# constrained outputs), to see where saturation happens and whether partial
# rounds inflate word0 the way theory says.
# ---------------------------------------------------------------------------

def exp_A2():
    print("="*70)
    print("EXPERIMENT A2: per-round degree growth on the REAL t=16 RF=6 RP=10 instance")
    print("  (degree-calculus, validated against sympy in A1)")
    print("="*70)
    t, rf, rp = 16, 6, 10
    _, hist = degree_calculus(t, rf, rp, 3, list(range(2, t)))
    import math
    for (label, deg) in hist:
        d0 = deg[0]
        print(f"  {label:8s}: deg(word0)={d0}  log2={math.log2(d0) if d0>0 else 0:.2f}")
    print()


# ---------------------------------------------------------------------------
# EXPERIMENT B: functional univariate degree along lines, exact over small prime q.
# We work the ACTUAL permutation in integer arithmetic mod q (alpha=3 needs
# gcd(3,q-1)=1 to keep sbox a bijection but that's not required for degree test).
# Functional degree = degree of unique interpolating poly of degree < q.
# Detect via finite differences: (d+1)-th forward difference over q consecutive
# points is identically the leading-coeff*d!*... ; we instead use the standard
# trick: f restricted to line has functional degree D iff the D-th order finite
# difference is constant nonzero and (D+1)-th vanishes -- but that needs >D pts.
# For potentially-large D we instead interpolate over a CHOSEN small q where
# 3^(RF+RP) >= q, so functional degree saturates at <= q-1.  We then DETECT
# collapse on SPECIAL lines: if special-line functional degree < random-line
# functional degree, that's a structural collapse.
# ---------------------------------------------------------------------------

def perm_int(state, mds, rc, t, rf, rp, alpha, q):
    half = rf//2
    def mat(s): return [sum(mds[i][j]*s[j] for j in range(t))%q for i in range(t)]
    def sb(x): return pow(x, alpha, q)
    s = mat([v%q for v in state]); idx=0
    for _ in range(half):
        s=[(s[i]+rc[idx][i])%q for i in range(t)]; s=[sb(x) for x in s]; s=mat(s); idx+=1
    for _ in range(rp):
        s=[(s[i]+rc[idx][i])%q for i in range(t)]; s[0]=sb(s[0]); s=mat(s); idx+=1
    for _ in range(half):
        s=[(s[i]+rc[idx][i])%q for i in range(t)]; s=[sb(x) for x in s]; s=mat(s); idx+=1
    return s

def functional_degree_on_line(eval_fn, q):
    """eval_fn(tt) returns out0 value at line param tt in F_q. Return the degree
    of the unique <q interpolating univariate poly via Newton forward differences
    on points 0..q-1 (full table). Returns (deg). Cost O(q^2) ints -> keep q small."""
    vals = [eval_fn(tt % q) for tt in range(q)]
    # forward difference table; degree = largest k with k-th diff not all zero,
    # using the fact that for a function on all of F_q the interpolating poly has
    # degree = q-1 - (number of trailing-zero leading finite differences)...
    # Simpler & exact: build Newton divided differences over distinct nodes 0..q-1
    # but mod q with q prime, divided diffs are fine. Degree = index of last
    # nonzero leading coefficient.
    # We compute successive finite differences (nodes are consecutive ints).
    diffs = vals[:]
    last_nonzero = 0
    for k in range(1, q):
        diffs = [(diffs[i+1]-diffs[i]) % q for i in range(len(diffs)-1)]
        if any(d != 0 for d in diffs):
            last_nonzero = k
        if len(diffs) <= 1:
            break
    return last_nonzero

def exp_B():
    print("="*70)
    print("EXPERIMENT B: functional univariate degree on random vs SPECIAL lines")
    print("="*70)
    # choose a small prime q with gcd(3,q-1)=1 so cubing is a bijection (matches
    # the real cipher's sbox being a permutation). q=251: q-1=250, gcd(3,250)=1. good.
    q = 251
    for (t, rf, rp) in [(4,2,1),(4,4,4),(6,4,6),(8,6,6)]:
        M = generate_mds_matrix(t, q)
        pos = Poseidon(prime=q, alpha=3, t=t, r_f=rf, r_p=rp, mds=M)
        rc = pos.round_constants
        nfree = t - 2
        C1, C2 = 0xC09DE4 % q, 0xEE6282 % q
        def build_state(free):
            return [C1, C2] + list(free)
        # base point + direction in the free space (nfree dims), embedded into state[2:]
        def out0_on_line(base_free, dir_free):
            def ev(tt):
                free = [(base_free[i] + tt*dir_free[i]) % q for i in range(nfree)]
                return perm_int(build_state(free), M, rc, t, rf, rp, 3, q)[0]
            return ev
        # random lines
        rdegs = []
        for _ in range(4):
            base = [random.randrange(q) for _ in range(nfree)]
            d = [random.randrange(q) for _ in range(nfree)]
            rdegs.append(functional_degree_on_line(out0_on_line(base, d), q))
        # SPECIAL line 1: direction = right eigenvector of M restricted to free coords?
        #   The interesting structural direction is one that the initial MDS maps to
        #   a sparse early state. Simpler structural special dirs to try:
        #   (i) unit directions e_i (single free var active)
        unit_degs = []
        for i in range(nfree):
            base = [random.randrange(q) for _ in range(nfree)]
            d = [1 if j==i else 0 for j in range(nfree)]
            unit_degs.append(functional_degree_on_line(out0_on_line(base, d), q))
        # SPECIAL line 2: direction chosen so that after initial MDS, the input to
        #   the FIRST partial round's word0 sbox is constant (kills one cubing on
        #   the active path). We can't easily solve that for full rounds (all words
        #   cubed) but we try directions in the kernel of selected MDS rows.
        # SPECIAL line 3: direction = a single free var but pick base so that the
        #   word0 input to round-1 sbox is 0 (sbox derivative 0).
        naive = 3**(rf+rp)
        print(f"  t={t} RF={rf} RP={rp} q={q}: naive 3^(RF+RP)={naive} (sat cap q-1={q-1})")
        print(f"     random-line func degs : {rdegs}")
        print(f"     unit-dir  func degs   : {unit_degs}")
    print()


# ---------------------------------------------------------------------------
# EXPERIMENT B2: the REAL test of 'collapse on a special line'. Use a tiny
# SYMBOLIC instance and substitute a univariate line into the EXACT polynomial,
# then read off the true degree (not saturated by q). This distinguishes formal
# degree on generic vs special lines.
# ---------------------------------------------------------------------------

def exp_B2_one(t, rf, rp, p=1009):
    from sympy import symbols, Poly, GF
    DOM = GF(p)
    M = generate_mds_matrix(t, p)
    pos = Poseidon(prime=p, alpha=3, t=t, r_f=rf, r_p=rp, mds=M)
    rc = pos.round_constants
    nfree = t-2
    tt = symbols('tt')
    C1, C2 = 0xC09DE4 % p, 0xEE6282 % p
    naive = 3**(rf+rp)
    def line_degree(base, direction):
        # state[2:] = base + tt*direction, all as Poly in tt
        st = [Poly(C1, tt, domain=DOM), Poly(C2, tt, domain=DOM)]
        for i in range(nfree):
            st.append(Poly(base[i] + direction[i]*tt, tt, domain=DOM))
        mds = [[Poly(M[i][j], tt, domain=DOM) for j in range(t)] for i in range(t)]
        rcp = [[Poly(rc[r][i], tt, domain=DOM) for i in range(t)] for r in range(rf+rp)]
        o = run_perm_symbolic(st, mds, rcp, t, rf, rp, 3, p)
        return o[0].degree(), o[1].degree()
    print(f"  t={t} RF={rf} RP={rp} naive(univariate cap)=3^(RF+RP)={naive}")
    for _ in range(3):
        base=[random.randrange(p) for _ in range(nfree)]
        d=[random.randrange(p) for _ in range(nfree)]
        print(f"     random line: deg(out0),deg(out1) = {line_degree(base,d)}")
    for i in range(nfree):
        base=[random.randrange(p) for _ in range(nfree)]
        d=[1 if j==i else 0 for j in range(nfree)]
        print(f"     unit-dir e{i}: deg(out0),deg(out1) = {line_degree(base,d)}")
    # SPECIAL line: choose direction so that after the initial MDS the input to
    # the FIRST full round's word-0 S-box is CONSTANT along the line (kills one
    # cube on that coordinate). Free coords are state[2:]; we need
    #   (M @ [C1,C2,free])[0] independent of tt  =>  sum_j M[0][2+j]*dir[j] = 0.
    # Solve a nontrivial dir in the (nfree-1)-dim kernel of that single linear eqn.
    row = [M[0][2+j] % p for j in range(nfree)]
    # pick dir: set dir = e0 adjusted to cancel via another coord
    if nfree >= 2 and row[0] % p != 0:
        dir_sp = [0]*nfree
        dir_sp[1] = 1
        dir_sp[0] = (-row[1]*pow(row[0],-1,p)) % p
        base=[random.randrange(p) for _ in range(nfree)]
        print(f"     SPECIAL (kill rnd1 word0 sbox): deg(out0),deg(out1) = {line_degree(base,dir_sp)}")
        # sanity: confirm the word0 sbox input really is constant on this line
        chk = (sum(M[0][2+j]*dir_sp[j] for j in range(nfree))) % p
        print(f"        [check sum M[0][2+j]*dir = {chk} (should be 0)]")


def exp_B2():
    print("="*70)
    print("EXPERIMENT B2: EXACT (non-saturated) degree of out0 on lines (symbolic, univariate)")
    print("="*70)
    for (t, rf, rp) in [(4,2,2),(4,2,4),(6,2,2)]:
        exp_B2_one(t, rf, rp)
    print()


# ---------------------------------------------------------------------------
# EXPERIMENT C: low-degree algebraic relation among inputs and outputs.
# We test whether out0 satisfies any low-degree polynomial relation with a SINGLE
# free input x_j (the others fixed): i.e. is the map x_j -> out0 of low degree?
# More aggressively: is there a bilinear relation a*out0 + b*out1 + c*x_j + ... =0?
# We interpolate over small q and report rank of the monomial-evaluation matrix.
# ---------------------------------------------------------------------------

def exp_C():
    print("="*70)
    print("EXPERIMENT C: low-degree relation hunt among (x_j, out0, out1)")
    print("="*70)
    q = 251
    for (t, rf, rp) in [(4,2,2),(4,4,4),(6,4,6)]:
        M = generate_mds_matrix(t, q)
        pos = Poseidon(prime=q, alpha=3, t=t, r_f=rf, r_p=rp, mds=M)
        rc = pos.round_constants
        nfree = t-2
        C1, C2 = 0xC09DE4 % q, 0xEE6282 % q
        # fix all free vars except x_active; sweep x_active over F_q.
        active = 0
        base = [random.randrange(q) for _ in range(nfree)]
        pts = []
        for v in range(q):
            free = base[:]; free[active] = v
            o = perm_int([C1,C2]+free, M, rc, t, rf, rp, 3, q)
            pts.append((v % q, o[0], o[1]))
        # Try to find lowest total-degree relation R(x, out0, out1)=0 holding on all pts.
        # Build monomial matrix up to degree D in 3 vars; find left null space over GF(q).
        found = None
        for D in range(1, 6):
            monos = [(a,b,c) for a in range(D+1) for b in range(D+1) for c in range(D+1) if a+b+c<=D]
            # rows = points, cols = monomials
            rows = []
            for (x,o0,o1) in pts:
                rows.append([ (pow(x,a,q)*pow(o0,b,q)*pow(o1,c,q))%q for (a,b,c) in monos])
            # null space (coeff vector) over GF(q): gaussian elim on transpose
            ns = gf_right_nullspace(rows, len(monos), q)
            if ns:
                found = (D, len(monos), len(ns))
                break
        nz = "NONE up to D=5"
        if found:
            nz = f"degree {found[0]} relation EXISTS (monos={found[1]}, nullity={found[2]})"
        # Also the pure univariate functional degree of x_active -> out0:
        fd = functional_degree_on_line(lambda vv: perm_int([C1,C2]+ [ (base[j] if j!=active else vv) for j in range(nfree)], M, rc, t, rf, rp, 3, q)[0], q)
        print(f"  t={t} RF={rf} RP={rp}: relation(x,o0,o1): {nz}; func-deg(x_active->out0)={fd}")
    print()

def gf_right_nullspace(rows, ncols, q):
    """Return basis of {v : rows @ v = 0} over GF(q). rows: list of length-ncols."""
    # Reduced row echelon on the matrix; nullspace from free columns.
    A = [r[:] for r in rows]
    nr = len(A)
    pivots = []
    pr = 0
    for c in range(ncols):
        piv = -1
        for r in range(pr, nr):
            if A[r][c] % q != 0:
                piv = r; break
        if piv == -1:
            continue
        A[pr], A[piv] = A[piv], A[pr]
        inv = pow(A[pr][c], -1, q)
        A[pr] = [(x*inv)%q for x in A[pr]]
        for r in range(nr):
            if r != pr and A[r][c] % q != 0:
                f = A[r][c]
                A[r] = [(A[r][k]-f*A[pr][k])%q for k in range(ncols)]
        pivots.append(c); pr += 1
        if pr == nr: break
    free = [c for c in range(ncols) if c not in pivots]
    basis = []
    for fc in free:
        v = [0]*ncols
        v[fc] = 1
        for i,pc in enumerate(pivots):
            v[pc] = (-A[i][fc]) % q
        basis.append(v)
    return basis


# ---------------------------------------------------------------------------
# EXPERIMENT D: integral / higher-order distinguisher. Sum out0 over an affine
# subspace V of dimension d (the active free coords range over a subfield-like
# cube of size q in each active dim, i.e. all of F_q). If the multivariate
# polynomial degree in those active vars is < d*(q-1)... actually the clean test:
# summing a monomial x1^a1..xd^ad over (F_q)^d is 0 unless every ai is a multiple
# of (q-1) and >0. So sum over full (F_q)^d of out0 = 0 unless out0 contains a
# monomial with each active var to a power that is a positive multiple of (q-1).
# A degree-D poly with D < d*(q-1) MUST sum to 0. We measure the minimal #active
# vars d (cube dimension) for which the sum is NONzero -> relates to deg per var.
# ---------------------------------------------------------------------------

def exp_D():
    print("="*70)
    print("EXPERIMENT D: integral sums over affine cubes (F_q)^d")
    print("="*70)
    q = 251  # gcd(3,q-1)=gcd(3,250)=1
    for (t, rf, rp) in [(4,2,1),(4,2,2),(4,4,4),(6,4,6)]:
        M = generate_mds_matrix(t, q)
        pos = Poseidon(prime=q, alpha=3, t=t, r_f=rf, r_p=rp, mds=M)
        rc = pos.round_constants
        nfree = t-2
        C1, C2 = 0xC09DE4 % q, 0xEE6282 % q
        base = [random.randrange(q) for _ in range(nfree)]
        results = []
        # d=1: sum over one active var across all F_q
        for d in range(1, min(nfree, 3)+1):
            active = list(range(d))
            total0 = 0
            # iterate the full cube (q^d) -- keep q^d small: q=251, d<=2 -> 63001 ok, d=3 -> 1.5e7 too big
            if q**d > 200000:
                results.append((d, "skip(too big)"))
                continue
            for combo in itertools.product(range(q), repeat=d):
                free = base[:]
                for k,a in enumerate(active):
                    free[a] = combo[k]
                total0 = (total0 + perm_int([C1,C2]+free, M, rc, t, rf, rp, 3, q)[0]) % q
            results.append((d, total0))
        print(f"  t={t} RF={rf} RP={rp} q={q}: sum_out0 over (F_q)^d : {results}")
        print(f"     (=0 means total deg in those d vars < d*(q-1)={ '%d'% (q-1) } per var threshold not reached)")
    print()


if __name__ == "__main__":
    t0 = time.time()
    # cheap integer-arithmetic experiments + instant calculus first
    exp_A2()   # degree-calculus on REAL t=16 instance (instant)
    exp_B()    # functional degree on random vs special lines (integer, fast)
    exp_C()    # low-degree relation hunt (integer, fast)
    exp_D()    # integral / higher-order sums (integer, moderate)
    # symbolic validation last (sympy, can be slower)
    exp_A()    # validates the calculus against true sympy degree (deg<=81)
    exp_B2()   # exact non-saturated line degrees (univariate sympy, fast)
    print(f"TOTAL elapsed {time.time()-t0:.1f}s")

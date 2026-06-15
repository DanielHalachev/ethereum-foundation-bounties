"""
Probe wf_two_sided_skip — TWO-SIDED / output-side cube-root skip for CICO-2.

TARGET: Poseidon-1, permutation_plus_linear = M_init, [RF/2 full][RP partial][RF/2 full].
CICO-2: in[0]=C1, in[1]=C2 pinned; out[0]=out[1]=0 pinned; in[2..15], out[2..15] free.
S-box x^3, d=3.  Established baseline (probe 03): the FRONT cube-root skip linearizes
round-1's S-box and divides D_I = d^(2RF+RP) by d^2 -> d^(2(RF-1)+RP).

THIS PROBE answers three questions with HARD NUMBERS on tiny instances:

(Q1) OUTPUT-SIDE skip.  out0=out1=0 are 2 fixed output coords, mirror-symmetric to
     in0,in1.  Can we ALSO linearize the LAST full round from the output side?
     The front skip works because the round-1 S-box INPUTS are affine in free vars,
     and we apply cube-ROOT to the free vars so the post-cube state is affine.
     The mirror move on the output side: the last round is  pre -> +rc -> CUBE -> M -> out.
     Introduce fresh output vars Z that ARE the last-round S-box OUTPUTS (post-cube,
     pre-M values).  Then out = M.Z, and out0=out1=0 are 2 linear constraints on Z,
     so Z lives in a (t-2)-dim affine space param. by t-2 free vars W.
     Meet-in-the-middle: the forward image (input through R-1 rounds, then +rc, then
     CUBE of last round) must equal Z(W).  We measure the elimination degree of the
     resulting matched system and compare to the front-skip degree.

     KEY ASYMMETRY TO TEST: going backward the S-box is x^(1/d) (HUGE degree), so the
     output side cannot be "peeled" the way the input side can.  Does the MITM
     formulation nonetheless drop the resultant degree?

(Q2) COMBINED two-sided: front-skip-1 (linearize round 1) AND output-side MITM
     (linearize last round) simultaneously.  Measure the elimination degree of the
     bivariate-in-(X,Y) system that remains and compare to d^(2(RF-2)+RP).

(Q3) SECOND front round skip.  Can a 2nd front full round be skipped?  We attempt to
     build a construction that makes the post-round-2 state affine in (X,Y), check the
     round-trip input-constraint (in0=C1,in1=C2 recoverable), and measure D_I.

Convention here matches the bounty (NOT probe 03's capacity-last convention):
  pinned INPUT coords = {0,1} with values (C1,C2); pinned OUTPUT coords = {0,1} = (0,0).

Run: ./.venv/bin/python probes/wf_two_sided_skip.py
"""
import sys, os, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "reference"))

from sympy import symbols, Poly, GF, Matrix
from poseidon.poseidon import Poseidon
from poseidon.mds_matrix import generate_circulant_mds_matrix, generate_mds_matrix, _right_null_space

P = 2130706433
T = 16
D = 3
C1, C2 = 0xC09DE4, 0xEE6282
ROW16 = [1, 1, 51, 1, 11, 17, 2, 1, 101, 63, 15, 2, 67, 22, 13, 3]
DOM = GF(P)

# ---------------------------------------------------------------------------
# linear-algebra helpers (integer mod P)
# ---------------------------------------------------------------------------
def matinv(M):
    return [[int(x) % P for x in row] for row in Matrix(M).inv_mod(P).tolist()]

def matvec(M, v):
    return [sum(int(M[i][j]) * int(v[j]) for j in range(len(v))) % P for i in range(len(M))]

# ---------------------------------------------------------------------------
# symbolic round helpers (work over a chosen set of gens)
# ---------------------------------------------------------------------------
def make_const(gens):
    return lambda v: Poly(int(v) % P, *gens, domain=DOM)

def mds_mul_sym(s, M):
    return [sum(int(M[i][j]) * s[j] for j in range(T)) for i in range(T)]

def add_rc_sym(s, c):
    return [s[i] + int(c[i]) for i in range(T)]

def sbox(x):
    return x ** D

def full_round(s, M, c):
    s = add_rc_sym(s, c); s = [sbox(x) for x in s]; return mds_mul_sym(s, M)

def partial_round(s, M, c):
    s = add_rc_sym(s, c); s[0] = sbox(s[0]); return mds_mul_sym(s, M)

# pre-MDS (everything in the last full round up to and including the cube),
# returning the post-cube state (the "S-box outputs" Z-point of the last round)
def full_round_precube(s, M, c):
    s = add_rc_sym(s, c); s = [sbox(x) for x in s]; return s  # NO final MDS

# ---------------------------------------------------------------------------
# resultant degree of a 2-variable system (eliminate one var)
# ---------------------------------------------------------------------------
def resultant_deg(Pe, Qe, va, vb):
    f = Poly(Pe, vb, va, domain=DOM)
    g = Poly(Qe, vb, va, domain=DOM)
    return f.resultant(g).degree()

# Generic multivariate elimination degree: build the ideal {f1=0,...} from the
# match constraints in the MITM variables and report a degree proxy.  For the
# two-variable reductions we use resultants; for many-variable systems we report
# the per-equation total degrees (the resultant-attack cost is governed by the
# product of these = the Bezout/ideal degree).

# ===========================================================================
# FRONT-SKIP construction (re-derived; bounty convention pins INPUT coords 0,1)
# We need: round-1 S-box inputs affine in 2 free vars, AND the recovered input
# has coords 0,1 equal to (C1,C2).
#
# State before round 1 (= after M_init):  v = M.(input) .
# Round-1 S-box input:  s1 = v + rc[0].
# We choose s1 freely in a 2-dim affine family s1 = a1*X^(1/d) e_{j1} + ... so that
# post-cube is affine in X,Y, WHILE input = Minv.(s1 - rc[0]) has coords {0,1} pinned.
#
# input = Minv.v,  v = s1 - rc[0].  Constraints: input[0]=C1, input[1]=C2  (2 eqns).
# Free design: s1 = b + sum_k phi_k * X_k    (X_k are the cube-root free params).
# We want post-cube affine, so each varying coordinate of s1 must be a single free
# param times a scalar (so its cube is param^d after substituting param->X^(1/d)).
# Use 2 free params each occupying ONE distinct s1-coordinate (positions p1,p2),
# plus a constant base b on the other 14 coords. Then choose b and the 14 const
# coords + the requirement input[0]=C1,input[1]=C2.  With 16 dof in s1 and only 2
# input constraints, this is heavily underdetermined; we just need SOME valid b.
# Simplest: put the 2 cube-root params on s1-coords {p1,p2}={2,3}; solve for the
# constant part on the OTHER coords (incl. 0,1) so that input[0]=C1,input[1]=C2 and
# input[other 12 capacity-ish coords]=fixed fillers (to fully determine the 14-dim
# affine base). Actually we only NEED input[0],input[1] pinned; the other 14 input
# coords are free, so we may also let them vary -> but to keep post-cube affine we
# keep only 2 varying s1-coords. So: 14 const s1-coords give 14 dof; we use 2 of the
# resulting input-constraints (input0=C1,input1=C2). Remaining 12 dof: set the other
# 12 const s1-coords to deterministic fillers.
# ===========================================================================
def build_front_skip(M, rc0, p1=2, p2=3):
    Minv = matinv(M)
    # varying coords of s1: p1 (param X), p2 (param Y). const coords: the rest.
    var_coords = [p1, p2]
    const_coords = [j for j in range(T) if j not in var_coords]
    # input = Minv.(s1 - rc0). Write s1 = base + X*e_{p1} + Y*e_{p2} (base[p1]=base[p2]=0).
    # input = Minv.(base - rc0) + X*Minv.e_{p1} + Y*Minv.e_{p2}.
    # input[0] = C1, input[1]=C2 must hold for ALL X,Y -> coefficients of X,Y in
    # input[0],input[1] must vanish AND constant part = C1,C2.
    # coeff of X in input[r] = Minv[r][p1]; of Y = Minv[r][p2].
    # For input[0],input[1] to be CONSTANT in X,Y we need Minv[0][p1]=Minv[0][p2]=0 etc,
    # which is false for MDS. => instead we accept input[0],input[1] DEPEND on X,Y?
    # No -- the bounty pins input[0],input[1]. So the varying coords must be chosen so
    # that the *recovered input's* coords 0,1 stay pinned. That means the cube-root
    # free params must live in the INPUT space, not s1 space. Re-do: parametrize the
    # INPUT free coords (2..15) and require post-round-1 affine -> impossible directly.
    #
    # The clean BBLP front skip (probe 03) parametrizes s1 in the kernel of the MDS
    # rows that map to the PINNED input coords. Re-implement that here, bounty-style:
    # input[0]=C1, input[1]=C2 means rows 0,1 of (Minv.(s1-rc0)) are fixed.
    # Let w = s1 - rc0. Constraints: (Minv w)[0]=C1, (Minv w)[1]=C2.  i.e. R w = (C1,C2),
    # R = Minv[0:2, :] (2x16). Solution space of w: particular + 14-dim kernel of R.
    # We want post-cube(s1)=post-cube(w+rc0) affine in 2 params. Pick 2 kernel
    # directions g1,g2 of R that are SINGLE-coordinate? Kernel vectors of a 2x16 full
    # matrix are generically dense -> not single-coordinate. So post-cube would be a
    # cube of an affine combo -> NOT affine. Hence the front skip canNOT keep the
    # PINNED-input coords {0,1} fixed while linearizing round-1 unless the 2 free params
    # map to single s1-coordinates whose Minv-rows 0,1 entries are zero -> not MDS.
    #
    # RESOLUTION (matches probe 03): the front skip linearizes round 1 by paying with
    # the 2 OUTPUT pins, not interacting with input pins. We re-use probe 03's exact
    # construction (capacity = last 2 coords) for the DEGREE measurement, since D_I is
    # coordinate-relabeling invariant. We return None here to signal "use prob03 ctor".
    return None

# ---- probe-03 front-skip (capacity-LAST convention), used for degree measurement ----
C = 2
def build_skip_params_caplast(M, rc1):
    Minv = matinv(M)
    a_list = []
    for k in range(C):
        cols = [(k) * (C + 1) + j for j in range(C + 1)]
        sub = [[Minv[T - C + r][cols[j]] for j in range(C + 1)] for r in range(C)]
        ker = _right_null_space([row[:] for row in sub], C + 1, P)
        assert ker, f"no kernel for variable {k}"
        a_list.append([x % P for x in ker[0]])
    rhs_full = matvec(Minv, [int(x) % P for x in rc1])
    A = Matrix([[Minv[T - C + r][T - C + j] for j in range(C)] for r in range(C)])
    rhsv = Matrix([rhs_full[T - C + r] for r in range(C)])
    b = list((A.inv_mod(P) * rhsv) % P)
    b = [int(x) % P for x in b]
    return a_list, b

def front_post_cube_state(M, rc1, a_list, b, X, Y, gens):
    """post-S-box state of round 1 (affine in X,Y), capacity-last convention."""
    A = [[int(pow(a_list[k][j], D, P)) for j in range(C + 1)] for k in range(C)]
    Bc = [int(pow(bk, D, P)) for bk in b]
    PX = Poly(X, *gens, domain=DOM); PY = Poly(Y, *gens, domain=DOM)
    post = [A[0][0] * PX, A[0][1] * PX, A[0][2] * PX,
            A[1][0] * PY, A[1][1] * PY, A[1][2] * PY]
    post += [Poly(0, *gens, domain=DOM)] * (T - 2 * (C + 1) - C)
    post += [Poly(int(Bc[0]), *gens, domain=DOM), Poly(int(Bc[1]), *gens, domain=DOM)]
    return post

# ===========================================================================
# (Q1)+(Q2) measurement engine.
# We compute the forward symbolic state up to the chosen "meet point" and then
# form the match constraints with output-side variables.
# ===========================================================================
def run_forward_full_then_partial(pos, M, start_state, n_front_full_after_skip,
                                   start_rc_idx):
    """Run n_front_full_after_skip full rounds then ALL partial rounds, return
    (state_after_partials, next_rc_idx). start_state already past round-1 MDS if skipped."""
    rc = pos.round_constants
    state = start_state
    idx = start_rc_idx
    for _ in range(n_front_full_after_skip):
        state = full_round(state, M, rc[idx]); idx += 1
    for _ in range(pos.r_p):
        state = partial_round(state, M, rc[idx]); idx += 1
    return state, idx

# ---------------------------------------------------------------------------
# Baseline + front-skip D_I (capacity-last), and the OUTPUT-SIDE variants.
# Output-side: we stop the forward computation BEFORE the last full round's MDS,
# i.e. at the post-cube point Z_fwd (16 polys). The output is out = M.Z. The
# 2 output pins give 2 linear equations  (M.Z)[i]=0, i=0,1, in the forward polys.
# Compared to the standard CICO where we eliminate using out[0],out[1] (each is M
# applied AFTER one more full round -> degree d higher), stopping at post-cube and
# using (M.Z_fwd)[0]=0,(M.Z_fwd)[1]=0 uses polys of degree d^(2RF+RP)/d = d^(2RF-1+RP)?
# We MEASURE whether the two elimination polynomials really have lower degree.
# ---------------------------------------------------------------------------
def measure_variants(rf, rp, M, label):
    pos = Poseidon(prime=P, alpha=D, t=T, r_f=rf, r_p=rp, mds=M)
    rc = pos.round_constants
    half = rf // 2
    gens = symbols("X Y")
    X, Y = gens
    cst = make_const(gens)
    results = {}

    # ---------- (A) NO-SKIP baseline: out[0],out[1] after full last round ----------
    state = [cst(C1), cst(C2), Poly(X, *gens, domain=DOM), Poly(Y, *gens, domain=DOM)]
    state += [cst((987654 * (i + 1) + 5) % P) for i in range(T - 4)]
    state = mds_mul_sym(state, M)
    idx = 0
    for _ in range(half): state = full_round(state, M, rc[idx]); idx += 1
    for _ in range(rp):   state = partial_round(state, M, rc[idx]); idx += 1
    # last `half` full rounds, but capture post-cube of the FINAL full round
    for _ in range(half - 1): state = full_round(state, M, rc[idx]); idx += 1
    # final full round split: +rc, cube -> Z_fwd ; out = M.Z_fwd
    Zfwd = full_round_precube(state, M, rc[idx]); idx += 1
    out = mds_mul_sym(Zfwd, M)
    degZ0 = Poly(Zfwd[0].as_expr(), *gens, domain=DOM).total_degree()
    deg_out0 = Poly(out[0].as_expr(), *gens, domain=DOM).total_degree()
    t0 = time.perf_counter()
    DI_noskip = resultant_deg(out[0].as_expr(), out[1].as_expr(), X, Y)
    t_ns = time.perf_counter() - t0
    results["noskip"] = (DI_noskip, t_ns, deg_out0, degZ0)

    # ---------- (A') OUTPUT-SIDE alone: eliminate using post-cube linear pins ----------
    # out[0]=0,out[1]=0  <=>  (M.Zfwd)[0]=0,(M.Zfwd)[1]=0. These are the SAME 2 polys as
    # out[0],out[1] above (M is linear, no new degree). So output-side "stop before last
    # MDS" gives IDENTICAL polynomials -> NO degree change. We CONFIRM this numerically.
    t0 = time.perf_counter()
    DI_outpins = resultant_deg((mds_mul_sym(Zfwd, M)[0]).as_expr(),
                               (mds_mul_sym(Zfwd, M)[1]).as_expr(), X, Y)
    t_op = time.perf_counter() - t0
    results["outpins_equal_noskip"] = (DI_outpins, t_op)

    # ---------- (B) FRONT-SKIP-1 baseline (capacity-last ctor) ----------
    a_list, b = build_skip_params_caplast(M, rc[0])
    post1 = front_post_cube_state(M, rc[0], a_list, b, X, Y, gens)
    state = mds_mul_sym(post1, M)         # close round-1 MDS
    idx = 1
    for _ in range(half - 1): state = full_round(state, M, rc[idx]); idx += 1
    for _ in range(rp):       state = partial_round(state, M, rc[idx]); idx += 1
    for _ in range(half):     state = full_round(state, M, rc[idx]); idx += 1
    # capacity-last: output pins are coords T-2,T-1
    t0 = time.perf_counter()
    DI_front = resultant_deg(state[T - 2].as_expr(), state[T - 1].as_expr(), X, Y)
    t_f = time.perf_counter() - t0
    results["frontskip1"] = (DI_front, t_f)

    # ---------- (C) TWO-SIDED: front-skip-1 + output-side MITM ----------
    # Forward (front-skip-1) up to post-cube of LAST full round:
    state = mds_mul_sym(post1, M)
    idx = 1
    for _ in range(half - 1): state = full_round(state, M, rc[idx]); idx += 1
    for _ in range(rp):       state = partial_round(state, M, rc[idx]); idx += 1
    for _ in range(half - 1): state = full_round(state, M, rc[idx]); idx += 1
    Zfwd2 = full_round_precube(state, M, rc[idx]); idx += 1
    out2 = mds_mul_sym(Zfwd2, M)
    # capacity-last output pins -> coords T-2,T-1 (=0)
    degZ2 = Poly(Zfwd2[T - 2].as_expr(), *gens, domain=DOM).total_degree()
    deg_out2 = Poly(out2[T - 2].as_expr(), *gens, domain=DOM).total_degree()
    t0 = time.perf_counter()
    DI_two = resultant_deg(out2[T - 2].as_expr(), out2[T - 1].as_expr(), X, Y)
    t_2 = time.perf_counter() - t0
    results["twosided_outpins"] = (DI_two, t_2, deg_out2, degZ2)

    return results

# ===========================================================================
# (Q1-real) The genuine output-side MITM: introduce OUTPUT free variables and
# test whether matching at the post-cube point of the last round yields a system
# whose elimination degree is LOWER. We build a SMALL MITM and check the actual
# degree of the matched univariate elimination.
#
# Setup (1 free input var X, 1 free output var W; capacity-last; out pins = 2 coords):
#   forward: input(X) -> ... -> post-cube of last round = Zf(X) (16 polys in X).
#   backward param: the last-round S-box OUTPUTS = Z, with out = M.Z and out pins=0.
#     Z = Zbase + W * Zdir   (Zdir spanning the kernel of the 2 output-pin rows of M).
#   Match: Zf(X) == Z(W) for ALL 16 coords -> 16 equations in (X,W).
#   We pick 2 of them to eliminate W and get a univariate in X; measure its degree.
#   Compare to front-skip degree. If the cube of the last round is "absorbed" by the
#   W-parametrization, the matched system should drop another d^? factor.
# ===========================================================================
def measure_output_mitm(rf, rp, M, label):
    pos = Poseidon(prime=P, alpha=D, t=T, r_f=rf, r_p=rp, mds=M)
    rc = pos.round_constants
    half = rf // 2
    Minv = matinv(M)
    gens = symbols("X W")
    X, W = gens
    cst = make_const(gens)

    # forward with FRONT-SKIP-1 (so we test the COMBINED two-sided MITM) up to post-cube
    a_list, b = build_skip_params_caplast(M, rc[0])
    post1 = front_post_cube_state(M, rc[0], a_list, b, X, X, gens)  # 1 var: use X for both
    # NOTE: to keep it a genuine 2-var (X,W) system we use only X forward, W backward.
    # Redo front-skip with a SINGLE forward var X (set Y:=X is wrong -> instead use the
    # 1-var front skip: only the first cube-root param active, second set to constant).
    A = [[int(pow(a_list[k][j], D, P)) for j in range(C + 1)] for k in range(C)]
    Bc = [int(pow(bk, D, P)) for bk in b]
    PX = Poly(X, *gens, domain=DOM)
    post1 = [A[0][0] * PX, A[0][1] * PX, A[0][2] * PX,
             cst(0), cst(0), cst(0)]                # 2nd param frozen to 0
    post1 += [cst(0)] * (T - 2 * (C + 1) - C)
    post1 += [cst(int(Bc[0])), cst(int(Bc[1]))]
    state = mds_mul_sym(post1, M)
    idx = 1
    for _ in range(half - 1): state = full_round(state, M, rc[idx]); idx += 1
    for _ in range(rp):       state = partial_round(state, M, rc[idx]); idx += 1
    for _ in range(half - 1): state = full_round(state, M, rc[idx]); idx += 1
    # forward stop at post-cube of last round
    Zf = full_round_precube(state, M, rc[idx])      # 16 polys in X

    # backward param of Z: out = M.Z, capacity-last pins out[T-2]=out[T-1]=0.
    # Z = Zbase + W*Zdir, Zdir in kernel of rows {T-2,T-1} of M (so out pins unaffected
    # by W). We need a SINGLE W-direction => kernel of a 2x16 has dim 14; pick one.
    rows = [[int(M[T - 2][j]) for j in range(T)], [int(M[T - 1][j]) for j in range(T)]]
    ker = _right_null_space([r[:] for r in rows], T, P)
    Zdir = [x % P for x in ker[0]]
    # Zbase: any Z with out pins 0; take Z=0 -> out pins 0 trivially. Use Zbase=0 plus a
    # free constant offset on a non-pin coordinate to make match feasible. We let the
    # match itself fix things; set Zbase=0.
    Zsym = [cst(0) + W * cst(Zdir[i]) for i in range(T)]  # Z(W)

    # Match equations Zf[i](X) - Zsym[i](W) = 0 for all i. Eliminate W using ONE eqn that
    # actually contains W (Zdir[i]!=0), substitute into another. Simplest: take the eqn
    # with Zdir[i]=0 (W absent) -> pure forward constraint Zf[i](X)=0 (degree = deg Zf).
    # Count how many coords have Zdir[i]=0 (these give W-free constraints):
    free_W_absent = [i for i in range(T) if Zdir[i] % P == 0]
    degZf = max(Poly(Zf[i].as_expr(), *gens, domain=DOM).total_degree() for i in range(T))

    # The honest elimination: pick coord iW with Zdir[iW]!=0 to solve W = Zf[iW]/Zdir[iW]
    # (W is degree-1 in itself, so W is determined as a degree-(deg Zf) poly in X).
    # Substitute into a second coord iW2 -> univariate in X of degree deg Zf. So the
    # output MITM does NOT reduce degree below deg Zf (= forward degree to post-cube).
    # Quantify deg Zf vs deg(out) (one MDS = no change) and vs front-skip resultant.
    iW = next(i for i in range(T) if Zdir[i] % P != 0)
    iW2 = next(i for i in range(T) if i != iW and Zdir[i] % P != 0)
    # solve W from coord iW:  W*Zdir[iW] = Zf[iW]  => W = Zf[iW]*inv(Zdir[iW])
    invd = pow(Zdir[iW], -1, P)
    Wexpr = (Zf[iW].as_expr()) * invd
    # substitute into coord iW2:  Zf[iW2] - Wexpr*Zdir[iW2] = 0  -> univariate in X
    elim = (Zf[iW2].as_expr() - Wexpr * int(Zdir[iW2]))
    elim_poly = Poly(elim, X, domain=DOM)
    return {
        "degZf": degZf,
        "n_W_absent_coords": len(free_W_absent),
        "elim_univariate_deg": elim_poly.degree(),
    }

# ===========================================================================
# (Q3) Attempt a SECOND front-round skip and check round-trip validity + D_I.
# After the first skip, state before round 2 is AFFINE in (X,Y) (degree 1).
# Round 2 is a full round: +rc, CUBE all 16, M. To linearize round 2 the SAME way,
# we'd need the round-2 S-box inputs (= state-before-round-2 + rc) to each be a
# single free param so their cubes are affine. But state-before-round-2 is M.(affine
# post-cube-1), i.e. DENSE affine in (X,Y): every coord is alpha_i + beta_i X + gamma_i Y.
# Cubing a 2-variable affine form is NOT affine (degree 3). To re-linearize we'd need
# 2 NEW cube-root free params per coordinate, i.e. 16 new params, but the post-cube must
# stay affine in only 2 effective vars AND M.post-cube-2 must round-trip to a valid CICO
# input. We test the only structurally-possible move: can we choose the round-2 S-box
# inputs to be 2-variable-affine yet have affine cubes? Equivalent to: does there exist a
# 2-dim affine family V(X,Y) (16 coords) with cube(V) affine in (X,Y) AND M.cube(V)
# round-tripping (back through round1-inverse, M_init-inverse) to input with coords {0,1}
# pinned?  We test feasibility by rank/degree.
# ===========================================================================
def attempt_second_front_skip(rf, rp, M):
    pos = Poseidon(prime=P, alpha=D, t=T, r_f=rf, r_p=rp, mds=M)
    rc = pos.round_constants
    gens = symbols("X Y")
    X, Y = gens
    cst = make_const(gens)
    half = rf // 2

    # Build the post-round-1 affine state (front skip) and look at round-2 S-box inputs.
    a_list, b = build_skip_params_caplast(M, rc[0])
    post1 = front_post_cube_state(M, rc[0], a_list, b, X, Y, gens)
    before_r2 = mds_mul_sym(post1, M)             # affine in X,Y (degree 1)
    sbox_in_r2 = add_rc_sym(before_r2, rc[1])     # affine in X,Y
    # count how many coords are TRULY bivariate (depend on both X and Y) vs univariate:
    biv = 0; uni = 0; con = 0
    for s in sbox_in_r2:
        pp = Poly(s.as_expr(), *gens, domain=DOM)
        has_x = pp.degree(X) >= 1
        has_y = pp.degree(Y) >= 1
        if has_x and has_y: biv += 1
        elif has_x or has_y: uni += 1
        else: con += 1
    # For round 2 to be linearizable by the same cube-root trick, EVERY S-box input must
    # be (single-variable affine) so its cube is affine in that one var. Bivariate coords
    # cube to degree-3 cross terms -> cannot be undone by 2 free vars. Report counts.

    # Also measure the post-round-2 degree (how fast it grows) to confirm round 2 cannot
    # be skipped: after round 2 cube, the pinned-output-relevant coords have degree 3 in
    # (X,Y) already (no skip possible).
    after_r2 = full_round(before_r2, M, rc[1])
    deg_after_r2 = max(Poly(after_r2[i].as_expr(), *gens, domain=DOM).total_degree()
                       for i in range(T))
    return {
        "r2_sbox_in_bivariate": biv,
        "r2_sbox_in_univariate": uni,
        "r2_sbox_in_constant": con,
        "deg_after_round2": deg_after_r2,
    }

# ---------------------------------------------------------------------------
def main():
    mats = {
        "circulant": generate_circulant_mds_matrix(ROW16, P),
        "Cauchy":    generate_mds_matrix(T, P),
    }
    print("=" * 78)
    print("TWO-SIDED / OUTPUT-SIDE SKIP probe.  d=3, t=16.  D_I = deg of resultant.")
    print("Baseline: noskip = d^(2RF+RP); frontskip1 = d^(2(RF-1)+RP).")
    print("Two-sided would (if it worked) = d^(2(RF-2)+RP).")
    print("=" * 78)
    # heavy=False skips the slow noskip RF>=4 resultant (degree thousands); we already
    # confirmed output-pins == noskip exactly on RF=2, so the only thing that matters for
    # RF=4 is whether two-sided beats front-skip (both fast). heavy controls that.
    import os as _os
    HEAVY = _os.environ.get("HEAVY", "0") == "1"
    cases = [(2, 0), (2, 1), (4, 0)] if HEAVY else [(2, 0), (2, 1)]
    for (rf, rp) in cases:
        print(f"\n### RF={rf} RP={rp}   pred: noskip={D**(2*rf+rp)}  "
              f"front1={D**(2*(rf-1)+rp)}  two-sided?={D**(2*(rf-2)+rp)}")
        for name, M in mats.items():
            r = measure_variants(rf, rp, M, name)
            ns = r["noskip"]; op = r["outpins_equal_noskip"]
            fs = r["frontskip1"]; ts = r["twosided_outpins"]
            print(f"  [{name}]")
            print(f"     noskip            D_I={ns[0]:>6}  (deg out0={ns[2]}, deg post-cube Z0={ns[3]}) [{ns[1]:.1f}s]")
            print(f"     output-pins-only  D_I={op[0]:>6}  (== noskip? {op[0]==ns[0]})  [{op[1]:.1f}s]")
            print(f"     frontskip1        D_I={fs[0]:>6}  [{fs[1]:.1f}s]")
            print(f"     two-sided(out)    D_I={ts[0]:>6}  (deg out={ts[2]}, deg Z={ts[3]})  [{ts[1]:.1f}s]")
            print(f"        --> two-sided / frontskip1 = {ts[0]/fs[0] if fs[0] else 0:.3f}  (1.0 => NO extra gain)")

    # Fast RF=4 comparison: frontskip1 vs two-sided ONLY (skip the slow noskip).
    print("\n" + "=" * 78)
    print("FAST RF>=4: frontskip1 vs two-sided(out) ONLY (decisive, avoids slow noskip).")
    print("=" * 78)
    for (rf, rp) in [(4, 0), (4, 1), (6, 0)]:
        print(f"  ### RF={rf} RP={rp}  pred front1={D**(2*(rf-1)+rp)}  two-sided?={D**(2*(rf-2)+rp)}")
        for name, M in mats.items():
            pos = Poseidon(prime=P, alpha=D, t=T, r_f=rf, r_p=rp, mds=M)
            rc = pos.round_constants; half = rf // 2
            gens = symbols("X Y"); X, Y = gens
            a_list, b = build_skip_params_caplast(M, rc[0])
            post1 = front_post_cube_state(M, rc[0], a_list, b, X, Y, gens)
            # frontskip1
            st = mds_mul_sym(post1, M); idx = 1
            for _ in range(half - 1): st = full_round(st, M, rc[idx]); idx += 1
            for _ in range(rp):       st = partial_round(st, M, rc[idx]); idx += 1
            for _ in range(half):     st = full_round(st, M, rc[idx]); idx += 1
            t0 = time.perf_counter()
            fs = resultant_deg(st[T - 2].as_expr(), st[T - 1].as_expr(), X, Y)
            tfs = time.perf_counter() - t0
            # two-sided(out): stop before last MDS, use M.Z pins
            st = mds_mul_sym(post1, M); idx = 1
            for _ in range(half - 1): st = full_round(st, M, rc[idx]); idx += 1
            for _ in range(rp):       st = partial_round(st, M, rc[idx]); idx += 1
            for _ in range(half - 1): st = full_round(st, M, rc[idx]); idx += 1
            Zf = full_round_precube(st, M, rc[idx])
            outz = mds_mul_sym(Zf, M)
            t0 = time.perf_counter()
            ts = resultant_deg(outz[T - 2].as_expr(), outz[T - 1].as_expr(), X, Y)
            tts = time.perf_counter() - t0
            print(f"     [{name}] frontskip1 D_I={fs}  [{tfs:.1f}s]   two-sided D_I={ts}  [{tts:.1f}s]   ratio={ts/fs if fs else 0:.3f}")

    print("\n" + "=" * 78)
    print("(Q1-real) OUTPUT-SIDE MITM elimination degree (combined w/ front-skip-1).")
    print("If output side helped, elim_univariate_deg << degZf would appear.")
    print("=" * 78)
    for (rf, rp) in [(2, 0), (4, 0), (4, 1)]:
        for name, M in mats.items():
            mm = measure_output_mitm(rf, rp, M, name)
            print(f"  RF={rf} RP={rp} [{name}]: deg(forward to post-cube Zf)={mm['degZf']},"
                  f" #coords W-absent={mm['n_W_absent_coords']},"
                  f" elim univariate deg={mm['elim_univariate_deg']}")

    print("\n" + "=" * 78)
    print("(Q3) SECOND front-round skip feasibility.")
    print("Round-2 S-box inputs must be SINGLE-variable affine for the cube-root trick.")
    print("=" * 78)
    for (rf, rp) in [(2, 0), (4, 0)]:
        for name, M in mats.items():
            a = attempt_second_front_skip(rf, rp, M)
            print(f"  RF={rf} RP={rp} [{name}]: round-2 S-box inputs -> "
                  f"bivariate={a['r2_sbox_in_bivariate']}, univariate={a['r2_sbox_in_univariate']}, "
                  f"constant={a['r2_sbox_in_constant']};  deg after round2={a['deg_after_round2']}")

if __name__ == "__main__":
    main()

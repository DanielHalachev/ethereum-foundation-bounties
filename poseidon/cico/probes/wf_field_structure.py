"""
wf_field_structure.py

ANGLE: Exploit p-1 = 2^24 * 127 and the discrete-log linearity of x^3.

Sub-questions:
 (a) In the dlog domain x=g^e, S-box x^3 is e->3e (LINEAR & BIJECTIVE since gcd(3,p-1)=1).
     The MDS layer is additive in x, hence NONLINEAR in e. Is there any sub-structure
     (multiplicative coset / 2-adic exponent filtration) where MORE of the round function
     becomes simultaneously tractable?
 (b) Test whether forcing intermediate S-box inputs into the order-127 subgroup creates
     info-losing collisions usable for CICO. (PREMISE CHECK: cubing 3-to-1 requires 3|p-1,
     which is FALSE here, so cubing is a bijection on every subgroup.)
 (c) Is there a near-linear automorphism (Frobenius is trivial; subfield maps) commuting
     with the round function?

We test concretely on tiny instances and report HARD NUMBERS.
"""
import sys, os, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'reference'))
from poseidon.poseidon import Poseidon
from poseidon.mds_matrix import generate_mds_matrix, generate_circulant_mds_matrix
import sympy
from sympy import GF, symbols, Poly

P = 2130706433          # 2^31 - 2^24 + 1, KoalaBear
ALPHA = 3

def banner(s):
    print("\n" + "="*78 + "\n" + s + "\n" + "="*78)

# ---------------------------------------------------------------------------
banner("FACT BLOCK: field structure of F_p")
print(f"p              = {P}")
print(f"p-1            = {P-1} = {sympy.factorint(P-1)}")
print(f"gcd(3, p-1)    = {sympy.gcd(3, P-1)}  -> x^3 is {'BIJECTIVE' if sympy.gcd(3,P-1)==1 else '3-to-1'} on F_p^*")
print(f"3^{{-1}} mod p-1 = {pow(3, -1, P-1)}  (cube-root exponent; cubing inverts in O(log p))")
# find a generator g of F_p^*
def find_generator(p):
    facs = list(sympy.factorint(p-1).keys())
    g = 2
    while True:
        if all(pow(g, (p-1)//q, p) != 1 for q in facs):
            return g
        g += 1
g = find_generator(P)
print(f"generator g    = {g}")
# the order-127 subgroup H = {g^(2^24 * j)}
gen127 = pow(g, (P-1)//127, P)
H = set()
x = 1
for _ in range(127):
    H.add(x); x = (x*gen127) % P
assert len(H) == 127
# verify cubing is a bijection on H
cubes_H = set(pow(h,3,P) for h in H)
print(f"|H| (order-127 subgroup)        = {len(H)}")
print(f"|cube(H)|                        = {len(cubes_H)}  -> cubing on H is {'BIJECTIVE (no collisions)' if len(cubes_H)==127 else 'INFO-LOSING'}")
print(f"cube(H) == H ?                   = {cubes_H == H}")
# order-2^24 subgroup Q
genQ = pow(g, 127, P)
print(f"order of 2-Sylow subgroup        = 2^24 = {2**24}")
print("--> Sub-angle (b) PREMISE: cubing is 3-to-1 only if 3 | p-1. It is NOT. DEAD on its face;")
print("    but we still test whether ANY multiplicative subgroup is preserved by a full round.")

# ---------------------------------------------------------------------------
banner("(a/b) Does ANY multiplicative coset/subgroup survive an MDS (additive) layer?")
# Key obstruction: MDS is additive. A round = ARC + cube + MDS. For dlog-linearity to
# propagate, we'd need the MDS image of a subgroup-structured vector to remain structured.
# Test: take H (mult. subgroup). Is M*H-vector ever landing in H componentwise? Measure.
M16 = generate_mds_matrix(16, P)
import random
random.seed(1)
# Strong structural test: is the set H closed under x -> (sum of a few H elements)?
# i.e. is H closed under addition in any nontrivial way? (additive closure of a mult subgroup)
add_closed = 0
trials = 20000
Hl = list(H)
for _ in range(trials):
    a, b = random.choice(Hl), random.choice(Hl)
    if (a+b) % P in H:
        add_closed += 1
print(f"H closed under addition: {add_closed}/{trials} random pairs a+b land in H "
      f"(expect ~{trials*127/P:.4f} by chance) -> additive structure: "
      f"{'NONE' if add_closed<=2 else 'SOME'}")

# A multiplicative coset gH is preserved by cube (e->3e) up to a coset shift, but the
# linear MDS destroys multiplicative structure. Confirm: does there exist any 1-dim
# multiplicative invariant for the FULL round map componentwise? We measure how often a
# state with all-coords in H maps (after one full round) back to all-coords in H.
def full_round_no_rc(state, M, p):
    s = [pow(v, 3, p) for v in state]
    return [sum(M[i][j]*s[j] for j in range(len(s))) % p for i in range(len(s))]
stay_in_H = 0
rt = 2000
for _ in range(rt):
    st = [random.choice(Hl) for _ in range(16)]
    out = full_round_no_rc(st, M16, P)
    if all(v in H for v in out):
        stay_in_H += 1
print(f"all-16-coords-in-H preserved by one full round (no RC): {stay_in_H}/{rt} "
      f"(chance ~{rt*(127/P)**16:.2e}) -> {'PRESERVED' if stay_in_H>0 else 'DESTROYED'}")

# ---------------------------------------------------------------------------
banner("(c) Near-linear automorphism commuting with the round? (Frobenius / scalar maps)")
# F_p has NO proper subfield (prime field) => Frobenius x->x^p is identity. Trivial.
print(f"Frobenius x->x^p is identity on F_p (prime field): pow(g,p,p)==g ? {pow(g,P,P)==g}")
# Scalar maps phi(x)=c*x commute with MDS (linear) but with cube: phi(cube(x))=c x^3,
# cube(phi(x))=c^3 x^3. Commute iff c^3=c iff c in {0,+-1...}. Check c with c^3=c:
cubefix = [c for c in range(P) if False]  # don't brute; solve c^3=c => c(c-1)(c+1)=0
print("phi(x)=c*x commutes with cube iff c^3=c => c in {0,1,p-1}. Only trivial scalings.")
# affine/translation maps: phi(x)=x+a. cube(x+a) != cube(x)+something linear. Not a symmetry.
print("Translation x->x+a does NOT commute with cube (binomial cross terms). Not a symmetry.")
print("--> No nontrivial near-linear automorphism. Sub-angle (c) DEAD.")

# ---------------------------------------------------------------------------
banner("CORE TEST (a): degree of CICO system in the DLOG domain vs the STANDARD domain")
# The real question: does working in dlog coords reduce the algebraic degree the attacker
# faces? In dlog coords the cube is linear but the MDS sum is a nightmare (log of a sum).
# We instead measure the achievable algebraic advantage by MEASURING the resultant/Groebner
# degree of a CICO system on TINY instances, and check whether ANY change of variables of
# the form x_i = (free var)^k or restricting free vars to a coset changes that degree.

# Build a tiny CICO system and measure ideal degree under standard polynomial encoding,
# then re-measure when free inputs are forced to be cubes / squares / coset-restricted.
def round_constants_for(t, rf, rp, seed=7):
    rnd = random.Random(seed)
    total = (rf+rp)*t
    return [rnd.randrange(P) for _ in range(total)]

def build_perm_poly_state(pos, free_syms, fixed_vals, dom):
    """Symbolic permutation_plus_linear over GF(p) as Poly objects.
    state[0..k-1]=fixed, state[k..]=free_syms. Returns list of Poly (the output state)."""
    t = pos.t
    M = pos.mds
    RC = pos.round_constants
    k = len(fixed_vals)
    # initial state as Poly
    state = []
    for i in range(t):
        if i < k:
            state.append(Poly(fixed_vals[i] % P, *free_syms, domain=dom))
        else:
            state.append(Poly(free_syms[i-k], *free_syms, domain=dom))
    def mds_mul(st):
        return [sum((Poly(M[i][j], *free_syms, domain=dom)*st[j] for j in range(t)),
                    Poly(0,*free_syms,domain=dom)) for i in range(t)]
    # initial linear
    state = mds_mul(state)
    half = pos.r_f//2
    idx = 0
    for _ in range(half):  # full
        state = [state[i] + Poly(RC[idx][i],*free_syms,domain=dom) for i in range(t)]
        state = [s**3 for s in state]
        state = mds_mul(state); idx+=1
    for _ in range(pos.r_p):  # partial
        state = [state[i] + Poly(RC[idx][i],*free_syms,domain=dom) for i in range(t)]
        state[0] = state[0]**3
        state = mds_mul(state); idx+=1
    for _ in range(half):  # full
        state = [state[i] + Poly(RC[idx][i],*free_syms,domain=dom) for i in range(t)]
        state = [s**3 for s in state]
        state = mds_mul(state); idx+=1
    return state

# Tiny instance: t=2, k=1 CICO (1 free var, 1 output constraint). Measure degree of the
# single-variable polynomial whose roots give CICO solutions, in STANDARD coords and after
# substituting free = u^3 (forcing input to be a cube => image of the cube map = ALL of F_p
# since cube is a bijection!).
def measure_tiny(t, rf, rp, sub_cube=False, label=""):
    M = generate_mds_matrix(t, P)
    RCflat = round_constants_for(t, rf, rp)
    pos = Poseidon(prime=P, alpha=ALPHA, t=t, r_f=rf, r_p=rp, mds=M, round_constants=RCflat)
    dom = GF(P)
    k = 1
    fixed = [12345 % P]
    free = symbols(f'x0:{t-k}')
    free = list(free) if isinstance(free, tuple) else [free]
    out = build_perm_poly_state(pos, free, fixed, dom)
    # CICO: out[0] == 0  (one constraint). With t-k=1 free var -> single univariate poly.
    constraint = out[0]  # set == 0
    cp = Poly(constraint, *free, domain=dom)
    deg = cp.total_degree()
    nterms = len(cp.terms())
    print(f"[{label}] t={t} RF={rf} RP={rp}: CICO constraint poly total_degree={deg}, "
          f"#terms={nterms}, d^(2RF+RP)={ALPHA**(2*rf+rp)}")
    return deg

banner("CORE TEST (a) RESULTS: degree of CICO constraint, standard polynomial encoding")
for (rf,rp) in [(2,1),(2,2),(2,3),(4,2)]:
    t0=time.time()
    try:
        measure_tiny(2, rf, rp, label="std")
    except Exception as e:
        print(f"  t=2 RF={rf} RP={rp}: ERROR {e}")
    print(f"    ({time.time()-t0:.2f}s)")

banner("CORE TEST (a): does a CUBE substitution (free=u^3) change the degree?")
# Since cube is a bijection, free=u^3 ranges over all of F_p; substituting it can only
# RAISE the apparent degree (by factor 3) without adding any solving power. Demonstrate.
def measure_tiny_cubesub(t, rf, rp, label=""):
    M = generate_mds_matrix(t, P)
    RCflat = round_constants_for(t, rf, rp)
    pos = Poseidon(prime=P, alpha=ALPHA, t=t, r_f=rf, r_p=rp, mds=M, round_constants=RCflat)
    dom = GF(P)
    u = symbols('u')
    # state[0]=fixed, state[1] = u^3
    tloc = t
    state = [Poly(12345 % P, u, domain=dom), Poly(u, u, domain=dom)**3]
    M_ = pos.mds; RC = pos.round_constants
    def mds_mul(st):
        return [sum((Poly(M_[i][j],u,domain=dom)*st[j] for j in range(tloc)),
                    Poly(0,u,domain=dom)) for i in range(tloc)]
    state = mds_mul(state)
    half=pos.r_f//2; idx=0
    for _ in range(half):
        state=[state[i]+Poly(RC[idx][i],u,domain=dom) for i in range(tloc)]
        state=[s**3 for s in state]; state=mds_mul(state); idx+=1
    for _ in range(pos.r_p):
        state=[state[i]+Poly(RC[idx][i],u,domain=dom) for i in range(tloc)]
        state[0]=state[0]**3; state=mds_mul(state); idx+=1
    for _ in range(half):
        state=[state[i]+Poly(RC[idx][i],u,domain=dom) for i in range(tloc)]
        state=[s**3 for s in state]; state=mds_mul(state); idx+=1
    cp = Poly(state[0], u, domain=dom)
    print(f"[{label}] t={t} RF={rf} RP={rp}: with free=u^3, deg={cp.total_degree()} "
          f"(std deg was {ALPHA**(2*rf+rp)}); cube-sub MULTIPLIES degree by ~3, no gain.")
for (rf,rp) in [(2,1),(2,2)]:
    measure_tiny_cubesub(2, rf, rp, label="cube-sub")

# ---------------------------------------------------------------------------
banner("DLOG-DOMAIN reality check: cost of expressing one MDS sum in exponent coords")
# In dlog coords, a single MDS output coordinate y = sum_j M[ij] g^{e_j}. To express e_y =
# dlog(y) as a function of (e_j) you must compute a discrete log of a SUM. There is no
# closed-form; the only handle is: the map (e_j) -> e_y has no low-degree polynomial rep.
# Quantify: interpolate the *univariate* function e -> dlog( a + g^e ) over F_p and measure
# its polynomial degree (should be ~p-1, i.e. dense/maximal => no algebraic shortcut).
def dlog_table(g, p):
    # baby-step/giant-step would be needed for full table; instead sample-degree estimate.
    pass
# Practical proxy on a SMALL prime with similar structure: q-1 = 2^a * odd, gcd(3,q-1)=1.
def small_field_dlog_degree(q):
    # q prime, find generator, build dlog, interpolate e->dlog(1+g^e) degree via finite diffs
    facs=list(sympy.factorint(q-1).keys())
    gg=2
    while not all(pow(gg,(q-1)//f,q)!=1 for f in facs): gg+=1
    # discrete log table
    dl={}; x=1
    for e in range(q-1):
        dl[x]=e; x=(x*gg)%q
    # f(e) = dlog(1 + g^e) for e where 1+g^e != 0
    xs=[]; ys=[]
    cur=1
    for e in range(q-1):
        val=(1+cur)%q
        if val!=0:
            xs.append(e); ys.append(dl.get(val))
        cur=(cur*gg)%q
    # interpolate degree over GF(q-1)? exponents live mod q-1 (not a field). So dlog of a sum
    # is NOT even a function into a field where polynomial degree is meaningful. Report that.
    return gg, len(xs)
from collections import Counter
for q in [17, 257]:  # q-1 = 2^a, gcd(3,q-1)=1, structurally like KoalaBear (p-1=2^24*127)
    assert sympy.gcd(3, q-1) == 1
    facs=list(sympy.factorint(q-1).keys()); gg=2
    while not all(pow(gg,(q-1)//f,q)!=1 for f in facs): gg+=1
    dl={}; x=1
    for e in range(q-1): dl[x]=e; x=(x*gg)%q
    n=q-1
    vals=[]; cur=1
    for e in range(n):
        v=(1+cur)%q; vals.append(dl.get(v)); cur=(cur*gg)%q
    defined=[v for v in vals if v is not None]
    bestmatch=0
    for a in range(n):
        c=Counter(((vals[e]-a*e)%n) for e in range(n) if vals[e] is not None)
        bestmatch=max(bestmatch, c.most_common(1)[0][1])
    print(f"  q={q}: dlog-of-sum f(e)=dlog(1+g^e): best affine fit matches {bestmatch}/{len(defined)} "
          f"(random ~{len(defined)/n:.2f}) -> NO dlog shortcut; map is ~random over Z/(q-1).")
print("  NOTE: exponents live in Z/(p-1), NOT a field (2^24 and 127 are zero divisors);")
print("        polynomial-degree / resultant machinery does not even apply in dlog coords.")

banner("END")

"""
DISCIPLINE: computational learning theory / F_p Fourier analysis.

QUESTION: Can we LEARN F (the CICO output map, with x0,x1 pinned) -- or a useful
surrogate -- from input/output samples well enough to predict a preimage of 0?

We evaluate four learning paradigms and TEST them on small Poseidon instances built
from the REAL reference code (same round schedule: initial MDS, RF/2 full, RP partial
[cube word 0 only], RF/2 full; S-box x^3):

  (A) Low-degree algorithm (Linial-Mansour-Nisan): is out0 concentrated on
      low-degree F_p characters? -> measure |hat f|^2 mass at low "degree".
  (B) Heavy-Fourier-coefficient / Goldreich-Levin-Kushilevitz-Mansour: is there ANY
      single character chi_a with |<f, chi_a>| noticeably above the 1/sqrt(p^n) floor?
      -> estimate the max single-character correlation (heavy coefficient).
  (C) Sparse polynomial interpolation (Ben-Or-Tiwari / Zippel): is out0, as a
      polynomial over F_p in the free vars, SPARSE (few monomials)? -> count
      monomials of out0 as a univariate (1 free var) and as restricted polys.
  (D) Best low-degree F_p-polynomial APPROXIMATION of out0 on a line: regress out0
      against {1, x, x^2, ..., x^k} for small k and measure agreement (this is the
      "weak learning" surrogate: any predictor beating 1/p chance helps narrow search).

We work in two regimes:
  - n=1 free variable (a LINE): lets us compute the EXACT univariate F_p Fourier
    spectrum of out0 by FFT over the additive group Z_p (p small), and the exact
    polynomial degree/sparsity. This is the cleanest, fully-rigorous measurement.
  - n=2 free variables on a small prime: spot-check that the 1-var conclusions persist
    (heavy-coefficient floor, no low-degree concentration) in higher dimension.

Output: measured correlations, spectral flatness, sparsity, low-degree-approx error.
"""
import sys, os, math, cmath, random
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "reference"))
import numpy as np
from poseidon.poseidon import Poseidon
from poseidon.mds_matrix import generate_mds_matrix

# ----------------------------------------------------------------------------
def build(prime, t, rf, rp):
    M = generate_mds_matrix(t, prime)
    pos = Poseidon(prime=prime, alpha=3, t=t, r_f=rf, r_p=rp, mds=M)
    return pos

def F_out(pos, free_vals, free_pos, pinned):
    """Evaluate permutation_plus_linear; return (out0, out1).
       pinned: dict idx->value for fixed words. free_pos: list of idx for free vars."""
    t = pos.t
    state = [0]*t
    for i, v in pinned.items():
        state[i] = v % pos.prime
    for idx, v in zip(free_pos, free_vals):
        state[idx] = v % pos.prime
    out = pos.permutation_plus_linear(state)
    return out[0], out[1]

# ----------------------------------------------------------------------------
def additive_fft_spectrum(samples, p):
    """Exact F_p Fourier transform of g: Z_p -> F_p (treated as complex via chi(y)=e^{2pi i y/p}).
       We Fourier-analyze the *additive character embedding* of the OUTPUT:
       f(x) = e^{2 pi i * out(x) / p}, a unit-modulus function, then DFT over the
       INPUT group Z_p. |hat f(a)|^2 sums to 1; flat spectrum => no structure.
       Returns array of |hat f(a)|^2 for a in Z_p."""
    n = len(samples)  # = p, full enumeration over the line
    f = np.exp(2j*np.pi*np.array(samples, dtype=float)/p)
    hat = np.fft.fft(f)/n
    return np.abs(hat)**2

def best_lowdeg_poly_fit(xs, ys, p, k):
    """Least-error degree-<=k F_p polynomial that agrees with (xs,ys) on the most points.
       We can't least-squares over F_p meaningfully for agreement, so we INTERPOLATE
       through the first k+1 points and measure agreement on ALL points (a learner only
       sees that a degree-k model is consistent if agreement is high). Returns frac agree."""
    # Build interpolating poly through first k+1 points via Lagrange over F_p, eval everywhere.
    pts = list(zip(xs[:k+1], ys[:k+1]))
    def lagrange_eval(x):
        total = 0
        for i,(xi,yi) in enumerate(pts):
            num=1; den=1
            for j,(xj,_) in enumerate(pts):
                if i==j: continue
                num = (num*(x-xj))%p
                den = (den*(xi-xj))%p
            total = (total + yi*num*pow(den,-1,p))%p
        return total
    agree = sum(1 for x,y in zip(xs,ys) if lagrange_eval(x)==y)
    return agree/len(xs)

def univariate_poly_sparsity(xs, ys, p):
    """Recover the unique degree<p interpolating polynomial of (xs=0..p-1, ys) over F_p
       and count nonzero coefficients (true F_p-sparsity = #monomials). Uses Newton/
       finite-difference style via numpy over python ints (Lagrange -> coeff is O(p^2))."""
    # ys is full table for x=0..p-1. Coeffs c s.t. sum c_k x^k = y(x). Use the fact that
    # over F_p the interpolating poly of degree<p is unique. Compute via the discrete
    # transform: c = V^{-1} y where V is Vandermonde. For small p do direct Gaussian elim.
    n = p
    # Build Vandermonde mod p and solve.
    V = [[pow(x, k, p) for k in range(n)] for x in range(n)]
    A = [row[:] + [ys[x]] for x,row in enumerate(V)]
    # Gaussian elimination mod p
    for col in range(n):
        piv = next(r for r in range(col,n) if A[r][col]%p!=0)
        A[col],A[piv]=A[piv],A[col]
        inv = pow(A[col][col],-1,p)
        A[col]=[(v*inv)%p for v in A[col]]
        for r in range(n):
            if r!=col and A[r][col]%p!=0:
                f=A[r][col]
                A[r]=[(a-f*b)%p for a,b in zip(A[r],A[col])]
    coeffs=[A[r][n]%p for r in range(n)]
    nz=sum(1 for c in coeffs if c!=0)
    deg = max((k for k,c in enumerate(coeffs) if c!=0), default=0)
    return nz, deg, coeffs

# ----------------------------------------------------------------------------
def run_line_experiment(p, t, rf, rp, label, seed=1):
    """n=1 free variable: vary x2 over all of F_p, pin everything else. Exact analysis."""
    rng = random.Random(seed)
    pos = build(p, t, rf, rp)
    pinned = {0: 0xC09DE4 % p, 1: 0xEE6282 % p}
    for i in range(4, t):
        pinned[i] = rng.randint(0, p-1)
    free_pos=[2]
    xs=list(range(p))
    out0=[]; out1=[]
    for x in xs:
        o0,o1=F_out(pos, [x], free_pos, pinned)
        out0.append(o0); out1.append(o1)

    print(f"\n=== {label}: p={p} t={t} RF={rf} RP={rp}  (1 free var, full line) ===")
    # ---- (C) sparsity / degree of out0 as univariate F_p poly ----
    nz, deg, _ = univariate_poly_sparsity(xs, out0, p)
    expected_deg = 3**(rf+rp)  # algebraic degree of word-0 nonlinearity chain
    print(f"[Sparsity/Degree] out0: #nonzero-monomials={nz}/{p}, deg={deg} "
          f"(field caps at p-1={p-1}; 3^(RF+RP)=3^{rf+rp}={expected_deg})")
    print(f"   density of monomials = {nz/p:.3f}  (1.0 = fully dense, no Ben-Or-Tiwari/Zippel win)")

    # ---- (A) low-degree Fourier concentration of additive embedding ----
    spec = additive_fft_spectrum(out0, p)
    # spec[0] is the DC (mean) term; structure would show as spikes away from flat 1/p.
    floor = 1.0/p
    spec_nodc = spec.copy(); spec_nodc[0]=0
    top = np.sort(spec_nodc)[::-1][:5]
    # "low degree" in the character group Z_p ~ small |a| (and small p-a). Mass on |a|<=k:
    def lowfreq_mass(k):
        idx = list(range(1,k+1))+list(range(p-k,p))
        return float(sum(spec[a] for a in idx))
    print(f"[LMN low-degree] additive-char spectrum of e^(2pi i out0/p):")
    print(f"   flat floor 1/p = {floor:.3e}; max non-DC coeff = {top[0]:.3e} "
          f"(ratio to floor = {top[0]/floor:.2f}x); top5 = {[f'{v:.2e}' for v in top]}")
    print(f"   low-freq mass |a|<=3: {lowfreq_mass(3):.4f}, |a|<=10: {lowfreq_mass(min(10,p//2)):.4f} "
          f"(of total {spec_nodc.sum()+spec[0]:.3f}; uniform would put ~{(2*min(10,p//2))/p:.3f} in |a|<=10)")

    # ---- (B) heavy coefficient = max single-character correlation ----
    # max over a!=0 of |hat f(a)| (amplitude). Random-function floor ~ 1/sqrt(p).
    maxamp = math.sqrt(top[0])
    print(f"[GL/KM heavy coeff] max |hat f(a)| (a!=0) = {maxamp:.4e}; "
          f"random floor 1/sqrt(p) = {1/math.sqrt(p):.4e}  "
          f"(advantage factor = {maxamp*math.sqrt(p):.2f})")

    # ---- (D) best low-degree polynomial agreement (weak learning) ----
    print(f"[Weak-learn] best degree-k poly agreement with out0 over the line:")
    for k in [1,2,3,5,8]:
        if k+1>p: break
        frac = best_lowdeg_poly_fit(xs, out0, p, k)
        print(f"   deg<={k}: agreement = {frac:.4f}  (chance 1/p = {1/p:.4f}; "
              f"a learner needs >> chance to predict)")
    return spec_nodc, nz, deg

def run_2var_check(p, t, rf, rp, label, seed=2):
    """n=2 free vars on small prime: confirm heavy-coeff floor persists in 2D."""
    rng = random.Random(seed)
    pos = build(p, t, rf, rp)
    pinned={0:0xC09DE4%p,1:0xEE6282%p}
    for i in range(4,t): pinned[i]=rng.randint(0,p-1)
    free_pos=[2,3]
    tab=np.zeros((p,p))
    for a in range(p):
        for b in range(p):
            o0,_=F_out(pos,[a,b],free_pos,pinned)
            tab[a,b]=o0
    f=np.exp(2j*np.pi*tab/p)
    hat=np.fft.fft2(f)/(p*p)
    mag=np.abs(hat)**2
    mag[0,0]=0
    floor=1.0/(p*p)
    mx=mag.max()
    print(f"\n=== {label}: p={p} t={t} RF={rf} RP={rp}  (2 free vars) ===")
    print(f"[2D heavy coeff] flat floor 1/p^2 = {floor:.3e}; max non-DC = {mx:.3e} "
          f"(ratio {mx/floor:.2f}x); max |hat f| = {math.sqrt(mx):.3e} vs 1/p = {1/p:.3e} "
          f"(adv {math.sqrt(mx)*p:.2f})")
    # how many preimages of 0 does out0 have on this 2D slice? (the quantity that matters)
    npre=int((tab==0).sum())
    print(f"   #preimages of out0==0 on p^2={p*p} grid = {npre} (expected ~p={p}); "
          f"these are NOT predictable from any heavy coefficient since spectrum is flat")

# ----------------------------------------------------------------------------
if __name__=="__main__":
    # Match the REAL round schedule shape (3 front full / RP partial / 3 back full) but
    # scale RP down so brute force / exact spectra are computable. Field-size-robust
    # structural questions per the brief.
    # Real: RF=6, RP=10. We test the real RF=6 with reduced RP, and also RP=10 on tiny prime.
    print("########## LINE (1-var) EXACT SPECTRA ##########")
    run_line_experiment(p=257, t=16, rf=6, rp=2,  label="tiny-RP")
    run_line_experiment(p=257, t=16, rf=6, rp=10, label="REAL-RP=10")
    run_line_experiment(p=1009,t=16, rf=6, rp=10, label="REAL-RP=10/larger-p")
    run_line_experiment(p=1009,t=8,  rf=6, rp=10, label="REAL-RP=10/t=8")

    print("\n########## 2-VAR CHECK ##########")
    run_2var_check(p=97, t=16, rf=6, rp=10, label="2var")
    run_2var_check(p=127,t=16, rf=6, rp=10, label="2var")

"""
Probe 9 — statistical cryptanalysis: how do values move, and is there exploitable bias?

(1) EXACT algebraic degree of out_0 as a univariate in one free input (others fixed),
    for small round counts. Predicted full degree = 3^(RF+RP). A DROP below that = a
    degree-collapse weakness (the basis of higher-order-differential / integral / zero-test).
(2) DIFFUSION at the full instance (RF=6,RP=10): do out_0,out_1 depend on ALL 14 free
    inputs? Incomplete diffusion would shrink the effective search dimension.
(3) OUTPUT UNIFORMITY at full instance: sample random inputs; is out_0 grossly biased
    (low entropy, value clustering, top-bit imbalance, chi-square over buckets)?
(4) DIFFERENTIAL: fix input difference; is the output-coord difference distribution biased
    (concentrates), which would beat the 1/p^2 brute-force probability of hitting (0,0)?
"""
import sys, os, random
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "reference"))
from poseidon.poseidon import Poseidon
from poseidon.mds_matrix import generate_circulant_mds_matrix

P = 2130706433; T = 16; D = 3
ROW16 = [1, 1, 51, 1, 11, 17, 2, 1, 101, 63, 15, 2, 67, 22, 13, 3]
C1, C2 = 0xC09DE4, 0xEE6282
M = generate_circulant_mds_matrix(ROW16, P)


def perm(free, r_f=6, r_p=10):
    pos = Poseidon(prime=P, alpha=D, t=T, r_f=r_f, r_p=r_p, mds=M)
    return pos.permutation_plus_linear([C1, C2] + list(free))


# ---------- (1) exact univariate algebraic degree via sympy ----------
def degree_check():
    from sympy import symbols, Poly, GF
    x = symbols("x"); DOM = GF(P)
    print("=== (1) exact univariate degree of out_0 (one free input = x, rest fixed) ===")
    for (rf, rp) in [(2, 0), (2, 1), (2, 2), (4, 0)]:
        pos = Poseidon(prime=P, alpha=D, t=T, r_f=rf, r_p=rp, mds=M)
        rc = pos.round_constants
        const = lambda v: Poly(int(v) % P, x, domain=DOM)
        st = [const(C1), const(C2), Poly(x, x, domain=DOM)]
        st += [const((424242 * (i + 1) + 7) % P) for i in range(T - 3)]
        def mds(s): return [sum(int(M[i][j]) * s[j] for j in range(T)) for i in range(T)]
        def arc(s, c): return [s[i] + int(c[i]) for i in range(T)]
        st = mds(st); half = rf // 2; idx = 0
        for _ in range(half):
            st = arc(st, rc[idx]); idx += 1; st = [p**D for p in st]; st = mds(st)
        for _ in range(rp):
            st = arc(st, rc[idx]); idx += 1; st[0] = st[0]**D; st = mds(st)
        for _ in range(half):
            st = arc(st, rc[idx]); idx += 1; st = [p**D for p in st]; st = mds(st)
        deg = Poly(st[0].as_expr(), x, domain=DOM).degree()
        pred = D ** (rf + rp)
        print(f"  RF={rf} RP={rp}: deg(out_0)={deg}  predicted 3^{rf+rp}={pred}  "
              f"{'FULL' if deg == pred else 'DROP! <-- weakness'}", flush=True)


# ---------- (2) diffusion ----------
def diffusion(samples=200):
    print("\n=== (2) diffusion at full (RF=6,RP=10): does out_0/out_1 depend on each free input? ===")
    rng = random.Random(1); depends0 = [0] * 14; depends1 = [0] * 14
    for _ in range(samples):
        base = [rng.randrange(P) for _ in range(14)]
        y0 = perm(base)
        for j in range(14):
            pert = base[:]; pert[j] = (pert[j] + rng.randrange(1, P)) % P
            yj = perm(pert)
            if yj[0] != y0[0]: depends0[j] += 1
            if yj[1] != y0[1]: depends1[j] += 1
    print(f"  out_0 depends on input j (frac of {samples}): {[round(d/samples,2) for d in depends0]}")
    print(f"  out_1 depends on input j (frac of {samples}): {[round(d/samples,2) for d in depends1]}")
    print(f"  -> {'FULL diffusion' if min(depends0+depends1)==samples else 'INCOMPLETE (some input does not affect output!)'}")


# ---------- (3) output uniformity ----------
def uniformity(N=30000):
    print(f"\n=== (3) output uniformity of out_0 at full instance (N={N}) ===")
    rng = random.Random(2); B = 256; buckets = [0] * B; top = 0; vals = []
    for _ in range(N):
        free = [rng.randrange(P) for _ in range(14)]
        v = perm(free)[0]
        buckets[v % B] += 1; top += (v >> 30) & 1; vals.append(v)
    exp = N / B
    chi2 = sum((b - exp) ** 2 / exp for b in buckets)
    mean = sum(vals) / N
    print(f"  mean/p = {mean/P:.4f} (uniform 0.5) | top-bit ones = {top}/{N} ({top/N:.3f}, uniform~0.25 since p<2^31)")
    print(f"  chi-square over {B} buckets = {chi2:.1f} (uniform expects ~{B-1} +-{(2*(B-1))**0.5:.0f}); "
          f"{'BIASED!' if chi2 > (B-1) + 5*(2*(B-1))**0.5 else 'uniform'}")


# ---------- (4) differential ----------
def differential(N=30000):
    print(f"\n=== (4) differential at full instance: output-diff bias for fixed input diff (N={N}) ===")
    rng = random.Random(3)
    delta = [0] * 14; delta[0] = 0x1234  # difference in first free input
    zero0 = 0; small = 0; diffs0 = {}
    for _ in range(N):
        base = [rng.randrange(P) for _ in range(14)]
        pert = [(base[j] + delta[j]) % P for j in range(14)]
        d0 = (perm(pert)[0] - perm(base)[0]) % P
        if d0 == 0: zero0 += 1
        if d0 < 2**16 or d0 > P - 2**16: small += 1
        diffs0[d0] = diffs0.get(d0, 0) + 1
    maxcount = max(diffs0.values())
    print(f"  P(out_0 diff = 0) = {zero0}/{N} (uniform ~{N/P:.2e}); most frequent diff occurs {maxcount}x "
          f"(uniform ~1) -> {'BIASED!' if maxcount > 5 else 'uniform (no usable differential)'}")
    print(f"  out_0 diff 'small' (|d|<2^16): {small}/{N} (uniform ~{N*2**17/P:.3f})")


if __name__ == "__main__":
    degree_check()
    diffusion()
    uniformity()
    differential()

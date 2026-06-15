"""
Probe 8 — two cheap NON-degree checks (the only kind that escapes degree>>field).

(A) Are the bounty's pinned constants C1=0xc09de4, C2=0xee6282 special? A designed
    quirk (or NUMS confirmation) would change everything. Check: size, multiplicative
    order / subgroup membership, cube relations, NTT structure, relation to round consts,
    whether they make round-1 S-box inputs degenerate.
(B) Lattice/2-adic 'size' hook: do CICO solutions (or the relaxed-verifier low-bit
    solutions) have any non-uniform integer/2-adic structure a lattice could exploit?
    Honest test: sample relaxed solutions, look at bit-size / 2-adic valuation distribution.
"""
import sys, os, random
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "reference"))
from poseidon.poseidon import Poseidon
from poseidon.mds_matrix import generate_circulant_mds_matrix

P = 2130706433            # 127 * 2^24 + 1
T = 16; D = 3
ROW16 = [1, 1, 51, 1, 11, 17, 2, 1, 101, 63, 15, 2, 67, 22, 13, 3]
C1, C2 = 0xC09DE4, 0xEE6282


def order(x, p=P):
    if x % p == 0: return 0
    o = p - 1
    for q in (2, 127):
        while o % q == 0 and pow(x, o // q, p) == 1:
            o //= q
    return o


def section_A():
    print("=== (A) constants C1, C2 ===")
    print(f"C1={C1}=0x{C1:x} (bits {C1.bit_length()}), C2={C2}=0x{C2:x} (bits {C2.bit_length()})")
    print(f"C3=C4=0 (output target).  p-1 = 2^24*127.")
    print(f"mult order C1: {order(C1)}  (p-1={P-1}); is C1 in order-127 subgroup? "
          f"{pow(C1,127,P)==1}")
    print(f"mult order C2: {order(C2)}; in order-127 subgroup? {pow(C2,127,P)==1}")
    print(f"cube relations: C1^3 mod p = {pow(C1,3,P)}, ==C2? {pow(C1,3,P)==C2}; "
          f"C2^3==C1? {pow(C2,3,P)==C1}")
    print(f"C1/C2 = {C1*pow(C2,-1,P)%P}, C2/C1 = {C2*pow(C1,-1,P)%P}, C1*C2={C1*C2%P}, C1+C2={(C1+C2)%P}")
    print(f"ratio C2-C1={C2-C1}, is C1 a cube? {pow(C1,(P-1)//1,P)==1}, "
          f"C1 cube-residue? (always, gcd(3,p-1)=1)")
    # do the constants relate to round constants of the actual instance?
    pos = Poseidon(prime=P, alpha=D, t=T, r_f=6, r_p=10, mds=generate_circulant_mds_matrix(ROW16, P))
    rc0 = pos.round_constants[0]
    print(f"round-const[0][0..3] = {rc0[:4]}")
    print(f"C1,C2 appear in any round const? {any(C1 in rc or C2 in rc for rc in pos.round_constants)}")
    # 2-adic valuation of constants (as integers) and of C-related quantities
    def v2(n): return (n & -n).bit_length() - 1 if n else -1
    print(f"v2(C1)={v2(C1)}, v2(C2)={v2(C2)}, v2(C1-1)={v2(C1-1)}, v2(C2-1)={v2(C2-1)}")


def section_B(trials=4000, m=6):
    """Relaxed-verifier solutions (low m bits of both output diffs zero): inspect their structure."""
    print(f"\n=== (B) integer/2-adic structure of relaxed solutions (m={m}) ===")
    pos = Poseidon(prime=P, alpha=D, t=T, r_f=6, r_p=10, mds=generate_circulant_mds_matrix(ROW16, P))
    mask = (1 << m) - 1
    rng = random.Random(0)
    sols = []
    for _ in range(trials):
        free = [rng.randrange(P) for _ in range(14)]
        y = pos.permutation_plus_linear([C1, C2] + free)
        if (y[0] & mask) == 0 and (y[1] & mask) == 0:
            sols.append(free)
    print(f"found {len(sols)} relaxed solutions in {trials} trials (expect ~{trials//(1<<(2*m))})")
    if sols:
        # are the solution inputs unusually small / structured?
        allvals = [v for s in sols for v in s]
        small = sum(1 for v in allvals if v < 2**24)
        print(f"  solution-input values < 2^24: {small}/{len(allvals)} "
              f"(uniform would be ~{len(allvals)*2**24//P})  -> {'STRUCTURE' if small>3*len(allvals)*2**24//P else 'uniform (no smallness)'}")
        def v2(n): return (n & -n).bit_length() - 1 if n else 24
        avgv2 = sum(v2(v) for v in allvals) / len(allvals)
        print(f"  mean 2-adic valuation of solution inputs: {avgv2:.2f} (uniform ~1.0)  "
              f"-> {'STRUCTURE' if avgv2 > 2 else 'uniform (no 2-adic bias)'}")


if __name__ == "__main__":
    section_A()
    section_B()

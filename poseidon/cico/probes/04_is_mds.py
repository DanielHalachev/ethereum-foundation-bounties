"""
Probe 4 — is the Plonky3 KoalaBear circulant matrix actually MDS?

MDS  <=>  every square submatrix is non-singular (all minors != 0).
This decides whether the CheapLunch *non-MDS* round-skip lever (skip >1 round,
reduce D_I by more than d^k) is available for our instance.

We exhaustively check all 1x1 and 2x2 minors, exhaustively 3x3, and heavily
sample 4x4..8x8.  A single singular submatrix proves NON-MDS.
"""
import sys, os, itertools, random
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "reference"))
from poseidon.mds_matrix import generate_circulant_mds_matrix, generate_mds_matrix

P = 2130706433
T = 16
ROW16 = [1, 1, 51, 1, 11, 17, 2, 1, 101, 63, 15, 2, 67, 22, 13, 3]


def det_mod(sub, p):
    """Determinant of a small matrix mod p via fraction-free / Gaussian elimination."""
    n = len(sub)
    M = [[x % p for x in row] for row in sub]
    det = 1
    for col in range(n):
        piv = next((r for r in range(col, n) if M[r][col] % p != 0), None)
        if piv is None:
            return 0
        if piv != col:
            M[col], M[piv] = M[piv], M[col]; det = (-det) % p
        det = (det * M[col][col]) % p
        inv = pow(M[col][col], -1, p)
        for r in range(col + 1, n):
            f = (M[r][col] * inv) % p
            if f:
                M[r] = [(M[r][c] - f * M[col][c]) % p for c in range(n)]
    return det % p


def submatrix(M, rows, cols):
    return [[M[i][j] for j in cols] for i in rows]


def check_mds(M, name, rng):
    print(f"=== {name} ===", flush=True)
    idx = list(range(T))
    worst = None
    # exhaustive 1x1, 2x2, 3x3
    for size in (1, 2, 3):
        sing = 0; total = 0
        for rows in itertools.combinations(idx, size):
            for cols in itertools.combinations(idx, size):
                total += 1
                if det_mod(submatrix(M, rows, cols), P) == 0:
                    sing += 1
                    if worst is None:
                        worst = (size, rows, cols)
        print(f"  {size}x{size}: {sing}/{total} singular", flush=True)
    # sampled 4x4..8x8
    for size in range(4, 9):
        sing = 0; trials = 20000
        for _ in range(trials):
            rows = tuple(rng.sample(idx, size)); cols = tuple(rng.sample(idx, size))
            if det_mod(submatrix(M, rows, cols), P) == 0:
                sing += 1
                if worst is None or size < worst[0]:
                    worst = (size, rows, cols)
        print(f"  {size}x{size}: {sing}/{trials} singular (sampled)", flush=True)
    verdict = "NON-MDS" if worst else "MDS (no singular minor found)"
    print(f"  --> {verdict}", flush=True)
    if worst:
        print(f"      smallest singular minor: size {worst[0]}, rows {worst[1]}, cols {worst[2]}", flush=True)
    return worst is None


def main():
    rng = random.Random(1)
    M = generate_circulant_mds_matrix(ROW16, P)
    check_mds(M, "Plonky3 circulant (operative matrix)", rng)
    print()
    check_mds(generate_mds_matrix(T, P), "Cauchy default (reference, should be MDS)", rng)


if __name__ == "__main__":
    main()

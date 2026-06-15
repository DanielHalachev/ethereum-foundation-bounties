"""
Probe 1 — linear-layer structure of the CICO target matrix.

Question: does the Plonky3 KoalaBear circulant matrix (the website's named linear
layer for the 2026 Poseidon1-31 CICO challenge) have eigen/subspace structure that
the textbook Cauchy MDS lacks, and could that buy extra round-skips?

We compute, over GF(p), p = 2^31 - 2^24 + 1:
  A. Eigenvalues of the circulant (= NTT of its first row), their multiplicative
     orders, distinctness, and any special structure. Cross-checked vs the
     reference char-poly roots.
  B. The Krylov / subspace-trail capacity: rank of {e0, e0·M, ..., e0·M^(l-1)}
     for l=1..t. dim S^(l) = t - rank. l_max = largest l with dim S^(l) >= 1
     bounds how many partial rounds a subspace trail can linearize.
  C. The same for the Cauchy default, for contrast.
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "reference"))

from poseidon.mds_matrix import (
    generate_circulant_mds_matrix, generate_mds_matrix,
    _char_poly, _roots_over_gfp, _row_space_basis, _mat_vec_mul,
)

P = 2130706433                      # KoalaBear, p-1 = 2^24 * 127
T = 16
ROW16 = [1, 1, 51, 1, 11, 17, 2, 1, 101, 63, 15, 2, 67, 22, 13, 3]  # Plonky3 koala-bear/src/mds.rs


def mult_order(x, p):
    """Multiplicative order of x in GF(p)*, p-1 = 2^24 * 127."""
    if x % p == 0:
        return 0
    order = p - 1
    for q in (2, 127):
        while order % q == 0 and pow(x, order // q, p) == 1:
            order //= q
    return order


def find_root_of_unity(n, p):
    """A primitive n-th root of unity in GF(p) (n | p-1)."""
    assert (p - 1) % n == 0
    cof = (p - 1) // n
    a = 2
    while True:
        w = pow(a, cof, p)
        if pow(w, n, p) == 1 and pow(w, n // 2, p) != 1:
            return w
        a += 1


def circulant_eigenvalues(row, p):
    """Eigenvalues of the circulant C[i][j]=row[(j-i)%t]: lambda_j = sum_k row[k] w^{jk}."""
    t = len(row)
    w = find_root_of_unity(t, p)
    eig = []
    for j in range(t):
        s = 0
        for k in range(t):
            s = (s + row[k] * pow(w, (j * k) % t, p)) % p
        eig.append(s)
    return eig, w


def left_krylov_rank_sequence(M, e0, p, t):
    """rank of {e0, e0 M, e0 M^2, ...} (row-vector * matrix) for l=1..t."""
    # left action: r_{j+1} = r_j * M  (row vector times matrix)
    def vec_mat(v):
        return [sum(v[i] * M[i][c] for i in range(t)) % p for c in range(t)]
    rows, seq, r = [], [], e0[:]
    for l in range(1, t + 1):
        rows.append(r)
        rank = len(_row_space_basis([row[:] for row in rows], p))
        seq.append(rank)
        r = vec_mat(r)
    return seq


def report(name, M):
    print(f"\n===== {name} =====")
    e0 = [1 if i == 0 else 0 for i in range(T)]
    seq = left_krylov_rank_sequence(M, e0, P, T)
    print("Krylov rank of {e0, e0 M, ...} for l=1..16:")
    print("  ", seq)
    lmax = max(l for l in range(1, T + 1) if (T - seq[l - 1]) >= 1)
    print(f"  -> dim S^(l) = 16 - rank; subspace nonempty (dim>=1) up to l_max = {lmax}")
    print(f"  -> e0 is {'CYCLIC (full Krylov, dim 16)' if seq[-1] == T else 'NOT cyclic (deficient)'}")


def main():
    M = generate_circulant_mds_matrix(ROW16, P)
    eig, w = circulant_eigenvalues(ROW16, P)
    print("Plonky3 circulant eigenvalues (NTT of first row):")
    print("  ", eig)
    print("  sum(row) = lambda_0 =", sum(ROW16) % P)
    print("  distinct eigenvalues:", len(set(eig)), "of", T)
    cp_roots = set(_roots_over_gfp(_char_poly(M, P), P))
    print("  match reference char-poly GF(p) roots:", set(eig) == cp_roots)
    orders = [mult_order(e, P) for e in eig]
    print("  multiplicative orders:", orders)
    print("  any eigenvalue = 1?", 1 in eig, " | any = 0?", 0 in eig)
    # is the eigenvalue list a geometric progression / related by a common ratio?
    print("  distinct orders:", sorted(set(orders)))

    report("Plonky3 circulant M16", M)
    report("Cauchy default M16", generate_mds_matrix(T, P))


if __name__ == "__main__":
    main()

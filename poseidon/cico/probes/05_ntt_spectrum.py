"""
Probe 5 (Lens 1) — structure of the NTT eigenbasis reformulation.

In the size-16 NTT basis the circulant M is diagonal: s_j <- lambda_j * s_j.
Partial round: s_j <- lambda_j*(s_j + c_j + Delta), Delta = m^3 - m, m = mean.
Moments M_k = sum_j lambda_j^k s_j undergo a SHIFT (M_k -> M_{k+1}) + rank-1 cubic kick.

We hunt for exploitable structure in the spectrum {lambda_j}:
  - distinctness, multiplicative orders (2-adic levels), 2-power-tower structure
  - multiplicative relations: is {lambda_j} closed under x*y or x^2? any lambda_a*lambda_b=lambda_c?
  - additive relations / small subset sums = 0 (would shrink the kick space)
  - the kick direction 1=(1,..,1): its Krylov space dim under D (caps reachable nonlinearity)
  - characteristic polynomial of D (= the moment linear-recurrence); its factorization
A rich, relationless spectrum = negative (no resonance shortcut); any relation = a lead.
"""
import sys, os, itertools
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "reference"))
from poseidon.mds_matrix import generate_circulant_mds_matrix

P = 2130706433             # 2^31 - 2^24 + 1 = 127 * 2^24 + 1
T = 16
ROW16 = [1, 1, 51, 1, 11, 17, 2, 1, 101, 63, 15, 2, 67, 22, 13, 3]


def find_root_of_unity(n, p):
    cof = (p - 1) // n
    a = 2
    while True:
        w = pow(a, cof, p)
        if pow(w, n, p) == 1 and pow(w, n // 2, p) != 1:
            return w
        a += 1


def eigenvalues(row, p):
    w = find_root_of_unity(len(row), p)
    return [sum(row[k] * pow(w, (j * k) % len(row), p) for k in range(len(row))) % p for j in range(len(row))]


def mult_order(x, p):
    if x % p == 0: return 0
    o = p - 1
    for q in (2, 127):
        while o % q == 0 and pow(x, o // q, p) == 1:
            o //= q
    return o


def two_adic_level(x, p):
    """e such that ord(x) = (p-1)/2^e, i.e. x is a 2^e-th power but not 2^{e+1}-th."""
    o = mult_order(x, p)
    e = 0
    while (p - 1) % (2 ** (e + 1)) == 0 and ((p - 1) // o) % (2 ** (e + 1)) == 0:
        e += 1
    return e


def rank_modp(rows, ncols, p):
    rows = [r[:] for r in rows]; rk = 0; pr = 0
    for c in range(ncols):
        piv = next((r for r in range(pr, len(rows)) if rows[r][c] % p), None)
        if piv is None: continue
        rows[pr], rows[piv] = rows[piv], rows[pr]
        inv = pow(rows[pr][c], -1, p)
        rows[pr] = [(x * inv) % p for x in rows[pr]]
        for r in range(len(rows)):
            if r != pr and rows[r][c] % p:
                f = rows[r][c]; rows[r] = [(rows[r][k] - f * rows[pr][k]) % p for k in range(ncols)]
        pr += 1; rk += 1
    return rk


def main():
    lam = eigenvalues(ROW16, P)
    S = set(lam)
    print("lambda_j:", lam)
    print("distinct:", len(S), "of", T, "| contains 0?", 0 in S, "| contains 1?", 1 in S)
    print("mult orders:", [mult_order(x, P) for x in lam])
    print("2-adic levels e (ord=(p-1)/2^e):", [two_adic_level(x, P) for x in lam])

    # multiplicative closure / relations
    prod_in = sum(1 for a in lam for b in lam if (a * b) % P in S)
    sq_in = sum(1 for a in lam if pow(a, 2, P) in S)
    abc = [(i, j, k) for i, a in enumerate(lam) for j, b in enumerate(lam)
           for k, c in enumerate(lam) if (a * b) % P == c and i <= j]
    print(f"\nmult: a*b in set: {prod_in}/{T*T} | a^2 in set: {sq_in}/{T} | a*b=c triples: {len(abc)}")

    # additive: any small subset of lambda summing to 0 (would give a degenerate kick combo)?
    zero_pairs = [(i, j) for i in range(T) for j in range(i + 1, T) if (lam[i] + lam[j]) % P == 0]
    zero_triples = sum(1 for c in itertools.combinations(range(T), 3) if sum(lam[i] for i in c) % P == 0)
    print(f"additive: lambda_i+lambda_j=0 pairs: {len(zero_pairs)} | zero-sum triples: {zero_triples}")

    # Krylov space of the kick direction 1=(1,...,1) under D=diag(lambda):  D^k 1 = (lambda_j^k)_j
    print("\nKrylov dim of all-ones kick direction under D (vs #partial rounds=10):")
    for r in (5, 10, 16):
        rows = [[pow(lam[j], k, P) for j in range(T)] for k in range(r)]
        print(f"  first {r} kicks span dim {rank_modp(rows, T, P)}")

    # characteristic polynomial of D = prod (x - lambda_j): splits into linear (lambda in F_p).
    # The moment sequence M_k = sum lambda_j^k s_j obeys this order-16 linear recurrence.
    print("\nchar poly of D splits completely over F_p (all roots in F_p):", len(S) == T,
          "-> moment recurrence has 16 distinct F_p modes (no order reduction).")

    # any lambda equal to a simple value (2-power, small) that could anchor a shortcut?
    print("small |lambda| (<10000):", sorted(x for x in lam if x < 10000))
    print("lambda that are pure 2-powers mod p:", [x for x in lam if (x & (x - 1)) == 0])


if __name__ == "__main__":
    main()

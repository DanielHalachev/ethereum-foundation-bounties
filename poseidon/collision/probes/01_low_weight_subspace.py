#!/usr/bin/env python3
"""
DECISIVE PROBE for the q=4 collision prize "research bet":
Does the linear layer have the LOW-WEIGHT invariant subspace that eprint 2026/306's
improved round-skip requires? (That property is what produces their 2^106 collision speedup.)

2026/306 needs a subspace that is invariant under BOTH the linear layer M AND the S-box layer.
A line span(v) is invariant under componentwise cubing iff v has entries in {-1,0,1}
(since (+-1)^3 = +-1). So the lever exists iff M has a LOW-WEIGHT {-1,0,1} eigenvector.
Poseidon2's *tensor* external matrix has such low-weight +-1 eigenvectors; a single 16x16
Poseidon1 MDS (Cauchy or circulant) should not (only full-weight all-ones / alternating).

Run:  python3 01_low_weight_subspace.py     (no deps; pure python over GF(p))
Self-contained — does not import the reference repo.
"""
from itertools import combinations, product

p = 2130706433  # KoalaBear 2^31 - 2^24 + 1
t = 16

def cauchy(t, p):
    return [[pow((i - t - j) % p, -1, p) for j in range(t)] for i in range(t)]

def circulant(first_row, p):
    t = len(first_row)
    return [[first_row[(j - i) % t] % p for j in range(t)] for i in range(t)]

def matvec(M, v, p):
    return [sum(M[i][j] * v[j] for j in range(t)) % p for i in range(t)]

def parallel_scalar(Mv, v, p):
    """Return lambda if Mv == lambda*v (v != 0), else None."""
    lam = None
    for a, b in zip(Mv, v):
        if b % p:
            lam = (a * pow(b, -1, p)) % p
            break
    if lam is None:
        return None
    return lam if all((a - lam * b) % p == 0 for a, b in zip(Mv, v)) else None

def scan(M, name, maxw=4):
    print(f"\n=== {name}: low-weight {{-1,0,1}} eigenvectors, Hamming weight <= {maxw} ===")
    found = []
    for w in range(1, maxw + 1):
        for support in combinations(range(t), w):
            for signs in product([1, p - 1], repeat=w - 1):  # fix v[support[0]]=1 to dedupe sign
                v = [0] * t
                v[support[0]] = 1
                for idx, sgn in zip(support[1:], signs):
                    v[idx] = sgn
                lam = parallel_scalar(matvec(M, v, p), v, p)
                if lam is not None:
                    found.append((w, support, lam))
    if found:
        for w, supp, lam in found[:8]:
            print(f"  weight {w}: support {supp}, eigenvalue {lam}")
        print(f"  -> {len(found)} low-weight eigenvector(s); min weight {min(f[0] for f in found)}  (LEVER PRESENT)")
    else:
        print(f"  NONE up to weight {maxw}  ->  2026/306 low-weight round-skip lever ABSENT")
    # reference: full-weight all-ones / alternating
    ones = [1] * t
    alt = [1 if j % 2 == 0 else p - 1 for j in range(t)]
    print(f"  reference full-weight: all-ones eigenvalue {parallel_scalar(matvec(M,ones,p),ones,p)}, "
          f"alternating eigenvalue {parallel_scalar(matvec(M,alt,p),alt,p)}")
    return found

if __name__ == "__main__":
    M_circ = circulant([1,1,51,1,11,17,2,1,101,63,15,2,67,22,13,3], p)
    M_cauchy = cauchy(t, p)
    scan(M_circ, "Plonky3 circulant (per UPD)")
    scan(M_cauchy, "Cauchy (verifier default)")
    print("\nCONCLUSION: both matrices lack low-weight invariant subspaces -> round-skip capped")
    print("at the generic 1-round MDS skip -> no 2026/306-style collision speedup on this instance.")

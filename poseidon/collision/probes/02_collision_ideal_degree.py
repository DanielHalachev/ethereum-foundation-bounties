#!/usr/bin/env python3
"""
Confirms the SECOND load-bearing fact for the q=4 collision bet:
the algebraic-collision ideal degree grows as ~alpha^(2*R_P) -- DOUBLE the preimage rate
(eprint 2026/306, Lemma 5.5: collision scales worse than preimage in the partial rounds).

Measures, on a tiny faithful instance, the resultant (ideal) degree of the 2-sided collision
system vs the single-side (preimage-proxy) degree, for R_P = 0, 1.

Two sides use DISTINCT fixed lanes so the trivial diagonal x=y is not a solution
(otherwise (s-u) divides both equations and the resultant is identically 0).

Run:  ../../.venv/bin/python 02_collision_ideal_degree.py   (needs sympy)

Observed (q=1009, t=4, R_F=2):
  R_P  1-side deg(=3^(RP+2))   collision ideal deg
   0          9                       54
   1         27                      558           ratio 558/54 = 10.3 ~ 3^2.1
-> ~9x (=3^2) per partial round for the collision vs ~3x for the preimage. Doubling confirmed.
Extrapolated to the real full instance (R_F=8, R_P=20, 1-round skip): collision ideal degree
is dominated by alpha^(2*R_P)=3^40 ~ 2^63 from partial rounds alone, so the algebraic solve
(>= ideal degree, typically ^omega) is >> 2^62 generic birthday. Algebra LOSES to brute force.
"""
import sympy as sp, random

q = 1009; t = 4
s, u = sp.symbols('s u')

def cauchy(t, q):
    return [[pow((i - t - j) % q, -1, q) for j in range(t)] for i in range(t)]

M = cauchy(t, q)
random.seed(1)
RC = [[random.randrange(q) for _ in range(t)] for _ in range(80)]
P = lambda c, v: sp.Poly(c, v, modulus=q)

def mds(st, v):
    return [sum((M[i][j] * st[j] for j in range(t)), P(0, v)) for i in range(t)]

def perm(v, RF, RP, consts):
    st = [sp.Poly(v, v, modulus=q)] + [P(c, v) for c in consts]
    rc = 0
    def full(st):
        nonlocal rc
        st = [st[i] + RC[rc][i] for i in range(t)]; rc += 1
        return mds([pp ** 3 for pp in st], v)
    def part(st):
        nonlocal rc
        st = [st[i] + RC[rc][i] for i in range(t)]; rc += 1
        return mds([st[0] ** 3] + st[1:], v)
    for _ in range(RF // 2): st = full(st)
    for _ in range(RP):      st = part(st)
    for _ in range(RF // 2): st = full(st)
    return st

def dgr(poly):
    try: return int(poly.degree())
    except Exception: return -1

if __name__ == "__main__":
    print("RF=2 | 1-side deg=3^(RP+2) ; collision ideal deg ~ alpha^(2*RP)-scaled (square of 1-side)")
    print(f"{'RP':>3} | {'1side':>6} | {'COLL dI':>8}")
    for RP in [0, 1]:  # RP=2 (deg 81) makes sympy's modular resultant very slow; trend is clear at 0,1
        oX = perm(s, 2, RP, (7, 11, 13))
        oY = perm(u, 2, RP, (8, 12, 14))   # distinct fixed lanes -> no x=y diagonal
        one = dgr(oX[0])
        e1 = sp.Poly(oX[0].as_expr() - oY[0].as_expr(), u, s, modulus=q)  # u first -> eliminate u
        e2 = sp.Poly(oX[1].as_expr() - oY[1].as_expr(), u, s, modulus=q)
        dI = dgr(sp.Poly(e1.resultant(e2), s, modulus=q))
        print(f"{RP:>3} | {one:>6} | {dI:>8}")

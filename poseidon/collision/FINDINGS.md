# Poseidon **Collision Prize** (q=4, $64K) — analysis & verdict

> This is a **different bounty** from the CICO challenge in [`../cico/`](../cico/). Do not conflate them.
> | | This (Collision Prize) | CICO (`../cico/`) |
> |---|---|---|
> | Target | hash `H` in **compression mode** (`permutation()` + feed-forward) | the **permutation** `permutation_plus_linear()` (extra initial MDS) |
> | Rounds | **full**: R_F=8, R_P=20 (≈ recommended set) | **round-reduced**: R_F=6, R_P∈{6,8,10} |
> | Goal | two 15-word inputs X≠Y with H(0xc09de4,X)=H(0xc09de4,Y) on first q words | one input hitting pinned output words |
> | q=4 prize | **$64,000** (q=3 claimed 2026-04-06) | $15K (R_P=10 open) |

## Verdict: **NO feasible attack.** Generic birthday 2⁶² is the best known; ~$3.4M ≫ $64K prize.

Three independent lines, each checked here, all land at/above the generic bar:

1. **No algebraic shortcut (round-skip lever absent).** eprint 2026/306's collision speedup needs a linear layer with **low-weight ±1 invariant subspaces** (Poseidon2's *tensor* external matrix has them). Poseidon1 here uses a single 16×16 MDS. [`probes/01_low_weight_subspace.py`](probes/01_low_weight_subspace.py) finds **none** (weight≤4) in either the Cauchy default or the Plonky3 circulant — only full-weight all-ones/alternating. So the round-skip is capped at the generic **1-round** MDS skip. **Matrix-independent.**

2. **Collisions scale *worse* than preimages here.** 2026/306 Lemma 5.5: collision ideal degree ∝ α^(2·R_P) — double the preimage's α^(R_P). [`probes/02_collision_ideal_degree.py`](probes/02_collision_ideal_degree.py) confirms ~3²≈9×/partial round (vs preimage 3×/round). With R_P=20 the algebraic collision is **≥ ~2⁷¹ > 2⁶²** — algebra loses to brute force. (Paper §6: for α≥3, full-round attacks exceed 2¹²⁸.)

3. **Everything else is blocked by design.** Differential/rebound (x³ DP=2/p ⇒ trail floor ≥2⁶⁶⁰), subspace/integral (no-invariant-subspace-trail MDS), Wagner (k=2, inapplicable). Full Gröbner/CICO on this exact permutation ≈ 2¹¹⁰ (eprint 2026/150 Table 3, the team that claimed the EF CICO ladder; CICO-only, never collisions, never full R_P).

## Compute model
Measured **2545 ns/perm** (KAT-verified C, M2 Pro) — [`probes/poseidon_kb_bench2.c`](probes/poseidon_kb_bench2.c). VW parallel collision search at 2⁶²: ~356 GPU-days / ~$3.4M on a 1000-GPU farm. Break-even vs the $64K prize is exponent ≈56; nothing found below 62.

## Known unknowns (what would change this)
- **Exact instance unconfirmed:** spec/UPD says "MDS arbitrary, e.g. Plonky3 circulant"; the reference verifier *defaults* to Cauchy+Grain and the q=3 submission verifies under Cauchy, not Plonky3. **Immaterial to the verdict** — every result above is matrix-independent (probe 01 checks both). Worth confirming with the organizer before any submission.
- **q=3 method unpublished** (report due 2027). Its solution differs only in input words 0–2 ⇒ looks like restricted-domain birthday (~2⁴⁶·⁵), which gives exactly 2⁶² at q=4 — consistent with this verdict, not a hidden shortcut.
- **Verifier gap (disclosure item, not a prize):** `verify_collision_solution` forwards solver-supplied `round_constants` with no validation; if ever allowed as input, collisions are plantable for any q.

Reproduce: `python3 probes/01_low_weight_subspace.py` and `../.venv/bin/python probes/02_collision_ideal_degree.py`.
Papers in [`../docs/`](../docs/): 2026-306 (round-skip+collision), 2026-150 (CICO ladder), 2025-2040 (CheapLunch), 2025-954 (subspace-trail GB), 2024-347 (FreeLunch).

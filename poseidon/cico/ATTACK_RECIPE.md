# CICO Attack Recipe — Poseidon1, R_F=6, R_P=10, $15,000

*Phase-2 synthesis of the modern literature into an implementable plan.*
*Companion to [CICO_THEORY.md](CICO_THEORY.md). Last updated: 2026-06-10.*

Papers read in full (PDFs in [docs/](docs/)):
- **2026/150** — *Ethereum Poseidon bounties* (Bak, Bariant, Jazeron, Bœuf, Hostettler) — **the team that claimed this exact ladder.** ← primary blueprint
- **2024/347** — *The Algebraic FreeLunch* (CRYPTO 2024)
- **2025/2040** — *The Algebraic CheapLunch* (FreeLunch → CICO-2)
- **2025/954** — *Poseidon & Neptune: GB Cryptanalysis Exploiting Subspace Trails*
- **2025/937** — *Attacking Poseidon via Graeffe-Based Root-Finding over NTT-Friendly Fields*

---

## 0. Headline conclusions

1. **The winning method is a resultant-based bivariate solve, not a full Gröbner basis.** Because our problem is **CICO-2** (k=2 output words = 0), the system reduces to two bivariate equations `{P(X,Y)=0, Q(X,Y)=0}`, and eliminating one variable by **resultant** costs `Õ(δ·D_I)` — vastly cheaper than building the full quotient ring (FGLM ≈ `D_I^ω ≈ 2^98`). This is exactly what 2026/150 used for the KoalaBear/Mersenne CICO-2 instances, and what the organizers mean by *"the best attack is a resultant attack."*

2. **Estimated cost for R_P=10: ≈ 2^55 field operations** (with the standard 1-round skip), vs brute force `p² ≈ 2^62`. So the attack *works* — but the margin over brute force is only ~2^7, and **memory is the real wall** (the resultant approach is quadratic in `δ`). That thin margin + memory pressure is precisely why R_P=10 is the lone unclaimed rung.

3. **The bonus/improvement lever is round-skipping.** 2026/150 skipped only **1** round. Theory (BBLP22 / CheapLunch) says an MDS front may allow skipping up to `s = ⌊t/k⌋ − 2 = 6` rounds for our `t=16, k=2`. Each extra skipped round divides `D_I` by `d=3`. Closing the gap from 1→more skipped rounds is the concrete, **bonus-eligible** research target ("new attack ideas might qualify for an additional bonus").

---

## 1. ⚠️ Reconciliation: which instance is the $15k, exactly

There are **two bounty editions**, and they differ. Pin this down before any compute:

| | 2026/150's tables (the *2024/25* EF round) | The repo `bounty2026.tex` + live site (*our* target) |
|---|---|---|
| KoalaBear instance name | **Poseidon2**-31k | **Poseidon1** (+ initial linear layer) |
| t, d, k | 16, 3, 2 (CICO-2) | 16, 3, 2 (CICO-2) |
| Reduced R_P rungs | **{1, 3, 4}** (+ full 8/20) | **{6, 8, 10}** |
| Status | Bak et al. solved up to R_F=6/R_P=4 | R_P=6,8 gone; **R_P=10 open** |

So the EF **escalated** for 2026: bigger `R_P` and a switch of base permutation. Our verifier ([cico_verifier.py:158](reference/bounties/cico_verifier.py#L158)) calls `permutation_plus_linear` = **Poseidon1 rounds + one extra MDS before round 1**. That initial-MDS front is Poseidon2-like, which matters for round-skipping (see §4). The 2026/150 *methodology* transfers directly; its *numbers* do not — we recompute below for our exact instance.

> **Action:** confirm on poseidon-initiative.info at submit time that R_P=10/Poseidon1 is still the open $15k rung, and that the MDS + constants match [reference/poseidon/](reference/poseidon/) (the attack constants are instance-specific).

---

## 2. Polynomial modeling of CICO-2 (2026/150 §4.1)

Write the permutation `𝒫 = F_R ∘ … ∘ F_1` (here R = R_F + R_P = 16 rounds, plus the initial MDS folded into `F_1`). Then:

1. Introduce **k=2 variables** `X, Y`. Set the input to the affine combination
   `input = X·A₁ + Y·A₂ + B`, where `A₁, A₂, B ∈ 𝔽ₚ¹⁶` are fixed so the input already satisfies the CICO **input** constraint (the 2 pinned words = the constants C₁,C₂; equivalently, in the homogenised "= 0" form, the constrained coords are 0). For us the 2 fixed input words are `0xC09DE4, 0xEE6282`; fold them into `B` and let `A₁,A₂` span the free directions.
2. Evaluate `𝒫` **symbolically** on this input — the round maps are affine + `x³`, so they act on polynomials cheaply. After all rounds you get `𝒫(input) = (P₁(X,Y), …, P₁₆(X,Y))`, each of total degree `d^R`.
3. Impose the 2 **output** constraints. The CICO-2 system is the **bivariate pair**:
   ```
   P(X,Y) := 𝒫(input)₀ = 0
   Q(X,Y) := 𝒫(input)₁ = 0
   ```
   (the first two output words forced to 0). 14 input dof − 2 = 12 leftover, so solutions exist generically.

The forward symbolic evaluation never introduces per-S-box variables (the maps are low-degree), so the model stays a clean **2-variable** system — the whole reason the resultant route is viable.

---

## 3. Solving by resultant (2026/150 §3.2) — the main pipeline

Given `{P(X,Y)=0, Q(X,Y)=0}` of degree `δ`:

1. **Eliminate Y:** `R(X) = Res_Y(P, Q)` via **evaluation–interpolation**:
   - Write `P = Σ Pᵢ(X)Yⁱ`, `Q = Σ Qᵢ(X)Yⁱ`.
   - Fast multipoint-evaluate the `Pᵢ, Qᵢ` on `D+1` geometric points (`D` = `deg_X R ≈ D_I`).
   - At each point compute a univariate `Res_Y` via half-gcd, then interpolate `R(X)`.
   - Cost `O(δ·M(D) + D·M(δ)·log δ)`, i.e. `Õ(δ·D_I)`.
2. **Root-find `R(X)` over 𝔽ₚ:** `Q(X)=Xᵖ−X mod R(X)` by double-and-add, then `G=gcd(R,Q)`; roots of `G` are the 𝔽ₚ-roots. (Optionally Graeffe — see §6.)
3. **Recover Y:** for each root `x₀`, `gcd(P(x₀,Y), Q(x₀,Y))` gives `Y`. Back-substitute `(X,Y) → (x₃,…,x₁₆)` and **check against the verifier**.

**Implementation:** NTL for univariate arithmetic; **PML** (Hyun–Neiger–Schost) for the bivariate resultant via Kronecker substitution; parallelised `gcd(Xᵖ−X, R)`. Reference code ships with 2026/150 and with [aurelbof/algebraic-freelunch], [atharva-simulauib/algebraic-cheaplunch].

---

## 4. Round-skipping — the degree-reducer (2026/150 §4.2, from BBLP22; CheapLunch App. D.1)

The ideal/resultant degree is set by the number of *active* rounds. Skipping front rounds is therefore the biggest lever.

- **The trick:** with capacity `c = k = 2` and `t ≥ c(c+2) = 8` (we have `t=16` ✓), choose the input direction vectors so that, after the first MDS, several branches enter the S-box layer as `a·X^{1/d}` and exit as `a^d·X` — **degree 1 instead of degree d**. This makes the first round(s) affine in the unknowns, removing them from the degree growth. Works whenever the relevant submatrix of `M⁻¹` has full rank — **always true for an MDS `M`** (our Cauchy matrix qualifies).
- **Poseidon note:** plain Poseidon1 applies an S-box layer *before* the first matrix, but `0^{1/d}=0`, so zeros in the capacity pre-matrix ≡ zeros in the input — the trick still applies. Our instance additionally has the extra initial MDS, which is the Poseidon2-style front the trick was written for.

**Ideal degree (CICO-2, [BSGL20] Thm 10):** `D_I = d^{k·R_F + R_P}`.
- Skipping **1** round (what 2026/150 actually did): `D_I = d^{2(R_F−1)+R_P}`, and input-poly degree `δ = d^{R_F−1+R_P}`.
- Conjectured max skip for MDS front (CheapLunch): `s = ⌊t/k⌋−2 = 6` rounds → divide `D_I` by `d^s = 3^6 = 729`. **Realising s>1 here is unverified and is the improvement target.**

---

## 5. Complexity table for our instance (Poseidon1, t=16, d=3, R_F=6, k=2)

Using `δ·D_I = d^{3(R_F−1)+2R_P}` (1-round skip) for the resultant solve:

| R_P | D_I = 3^{2(R_F−1)+R_P} | δ = 3^{R_F−1+R_P} | Resultant solve `Õ(δ·D_I)` | Brute `p²` | Status |
|----:|---|---|---|---|---|
| 6 | 3^16 ≈ 2^25.4 | 3^11 ≈ 2^17.4 | 3^27 ≈ **2^42.8** | 2^62 | solved (Apr 2026) |
| 8 | 3^18 ≈ 2^28.5 | 3^13 ≈ 2^20.6 | 3^31 ≈ **2^49.1** | 2^62 | submitted (May 2026) |
| **10** | 3^20 ≈ 2^31.7 | 3^15 ≈ 2^23.8 | 3^35 ≈ **2^55.5** | 2^62 | **OPEN — $15k** |

**Calibration (2026/150 Table 2, real runs):** Poseidon2-31k R_F=6/R_P=4 ≈ 2^46, **2.6 h, 10 TB**; Poseidon-256 R_F=6/R_P=9 ≈ 2^45.6, **8.5 days, 1.3 TB**. Extrapolating, R_P=10 at ≈2^55 is **~weeks of compute and very large memory** with the 1-round skip — feasible but at the frontier, which is consistent with it being unclaimed.

**With an aggressive skip (e.g. s=3):** `δ·D_I = d^{3(R_F−s)+2R_P} = 3^{3·3+20} = 3^29 ≈ 2^46` — back into the comfortably-practical regime. **This is why pushing the skip past 1 round is the whole game for R_P=10.**

---

## 6. Final root-finding: standard vs Graeffe (2025/937)

- The final univariate has degree `D_I ≈ 3^20 ≈ 2^31.7`, which **exceeds KoalaBear's radix-2 NTT ceiling of 2^24** (2-adicity = 24). So pure radix-2 NTTs don't fit; you'd need **mixed-radix / Bluestein** NTTs (exactly the regime 2025/937 used for Goldilocks).
- **Graeffe/tangent-Graeffe** replaces the `O(log p)` modular-exponentiation factor of the standard `gcd(Xᵖ−X, R)` with `O(log r)` iterations + an `O(s log s)` NTT — a log/constant-factor win, not an asymptotic class change.
- **But root-finding is NOT the bottleneck of the resultant attack** — the `Õ(δ·D_I)` resultant computation dominates. On KoalaBear `log p ≈ 31` is already ~8× smaller than the 256-bit fields where Graeffe paid off, so the absolute headroom is small. **Verdict: implement standard CZ/NTL root-finding first; treat Graeffe as a later optimisation only if profiling shows root-finding matters.**

---

## 7. Method selection summary

| Method | When it wins | For us (CICO-2, R_P=10) |
|---|---|---|
| **Resultant** (2026/150 §3.2) | k=2, bivariate, `q ≥ D_I` | ✅ **primary** — `Õ(δ·D_I) ≈ 2^55`, quadratic memory |
| **GCD** (2026/150 §4.4) | resultant memory too high | fallback — linear memory, but `O(q·M(δ)) ≈ 2^31·…` time |
| **FreeLunch** (2024/347) | CICO-**1** only (single x₀) | ✗ inapplicable (needs ≥2 zeros) — but its skip/det ideas inform us |
| **CheapLunch** (2025/2040) | k≥2 via GB, or k≥3 | ✗ for speed (FGLM ≈ 2^98), ✓ for its **round-skip theory** (s=6) |
| **Subspace-trail GB** (2025/954) | margin analysis, large t | reference for why skipping works; not a faster solver here |

---

## 8. Build plan (this is where code starts)

1. **Exact model harness.** Wrap [reference/poseidon/poseidon.py](reference/poseidon/poseidon.py) symbolically (Sage / NTL / `galois`). Reproduce `permutation_plus_linear` over `𝔽ₚ[X,Y]` and confirm a random `(X,Y)` matches the integer reference bit-for-bit.
2. **Toy validation.** Build `{P,Q}` for **R_P=1,2** and solve by resultant; cross-check the root against the relaxed verifier, then the exact verifier. Lock the pipeline at tiny scale.
3. **Round-skip module.** Implement the 1-round skip (§4) — construct `A₁,A₂,B` from `M⁻¹`. Verify `D_I` drops from `3^{2R_F+R_P}` to `3^{2(R_F−1)+R_P}` empirically at small R_P.
4. **Scale R_P → 10** with PML bivariate resultants; profile time **and memory**. Decide resultant vs GCD per memory budget.
5. **Improvement push (bonus):** attempt s>1 round-skip using the CheapLunch/BBLP22 MDS construction for `t=16, k=2`. Each extra skipped round = `÷3` on `D_I`. This is the path from "feasible but expensive" to "practical," and is independently publishable.
6. **Submit** `(x₃,…,x₁₆)`; re-confirm the live rung first.

---

## 9. Open questions / risks
- **How many rounds can actually be skipped** for Poseidon1 + initial-MDS, t=16, k=2? (1 is proven-in-practice; up to 6 is conjectured.) This dominates feasibility.
- **Memory** at R_P=10: resultant storage `~δ²` may be prohibitive; the GCD fallback trades it for `O(q)≈2^31` iterations.
- **Instance drift:** Poseidon1 vs Poseidon2, exact MDS/constants — regenerate `A₁,A₂,B` and `D_I` against the *live* spec, not the 2026/150 tables.
- Confirm `[BSGL20]` degree `D_I = d^{kR_F+R_P}` holds with the extra initial MDS (it should — the MDS is absorbed into `F_1`).

---

## 10. Probe log (empirical)

**Matrix resolved (website > repo authority):** linear layer = **Plonky3 KoalaBear circulant**
`first_row=[1,1,51,1,11,17,2,1,101,63,15,2,67,22,13,3]`, not the repo's Cauchy default.

**Probe 1 — diagonalizability is NOT an attack lever** ([probes/01_matrix_structure.py](probes/01_matrix_structure.py)).
The circulant is fully diagonalizable over GF(p) (16 eigenvalues in 𝔽ₚ; Cauchy has 0), but:
subspace-trail/Krylov skip capacity is identical to Cauchy (`e₀` cyclic, `ℓ_max=15` both);
eigenvalues carry no `x³` structure (not closed under cubing, no `λᵢ=λⱼ³`, all 240 pairwise
ratios distinct). Conclusion: don't chase the matrix; the lever is matrix-independent.

**Probe 2 — `D_I` measured, degree-drop + matrix-independence confirmed**
([probes/02_ideal_degree.py](probes/02_ideal_degree.py), sympy over GF(p), `D_I = deg Res_Y(P,Q)`):

| matrix | R_F | R_P | deg P | measured `D_I` | `d^{2R_F+R_P}` | trivial `d^{2(R_F+R_P)}` |
|---|---|---|---|---|---|---|
| Plonky3 circulant | 2 | 0 | 9 | **81** | 81 | 81 |
| Cauchy default | 2 | 0 | 9 | **81** | 81 | 81 |
| Plonky3 circulant | 2 | 1 | 27 | **243** | 243 | 729 |
| Cauchy default | 2 | 1 | 27 | **243** | 243 | 729 |

→ `D_I = d^{2R_F+R_P}` (the factor-`d^{R_P}` Poseidon drop) is real; `R_P=1` gives 243=3⁵, not
729=3⁶. And **both matrices give identical `D_I`** — empirically sealing Probe 1. So for the
real target `D_I(no skip) = 3^{2·6+10} = 3^{22} ≈ 2^{34.9}`, matrix-independent.

**Probe 3 — BBLP22 first-round skip implemented + validated** ([probes/03_round_skip.py](probes/03_round_skip.py)).
Built the §4.2 construction for the Plonky3 circulant: `a_k` = kernel of the `c×(c+1)` submatrix of
`M⁻¹`, offset `b` from the `c×c` capacity system. **Gated by a round-trip check** — the recovered
input `M⁻¹(s2−C⁽¹⁾)` has its capacity coords exactly 0, so it's a genuine CICO input.

| | roundtrip | no-skip `D_I` | skip-1 `D_I` | predicted `d^{2(R_F−1)+R_P}` |
|---|---|---|---|---|
| RF=2, RP=0 | ✅ | 81=3⁴ | **9=3²** | 9 ✓ |
| RF=2, RP=1 | ✅ | 243=3⁵ | **27=3³** | 27 ✓ |

→ the 1-round skip divides `D_I` by exactly `d²=9`. Confirms for the real target:
**`D_I(skip-1)=3²⁰≈2³¹·⁷`, solve `Õ(δ·D_I)≈2⁵⁵`.**

**`s>1` (multi-round skip) — assessed, NOT a free extension for Poseidon's full rounds.**
The naive push (apply the same construction one round deeper) provably breaks: inverting back
through a 2nd *full* round forces the capacity-zero condition to depend on non-separable `X^{1/d}`
terms, so no linear `a/b` keeps capacity = 0 for all `X,Y`. Literature multi-round bypass (FreeLunch
+ "bypass 3 rounds") is **primitive-specific to Griffin**, and the 2026/150 team used only the
1-round skip for Poseidon in their actual bounty attacks. So `s>1` for Poseidon is the
**bonus-eligible research frontier**, not a quick win. (A 2025 "Improved Resultant Attack" paper may
push further — unread.)

**Feasibility takeaway for R_P=10:** with the validated 1-round skip, the **GCD approach**
(2026/150 §4.4) is the realistic route — `O(q·M(δ)) ≈ 2⁵⁵` time but **linear memory** in `δ=3¹⁵≈10⁷`
(the resultant route's `~δ²` storage is the wall). I.e. R_P=10 is solvable by *engineering the
pipeline at 2⁵⁵ on serious hardware* (cf. 2026/150's 1 TB / 120-thread runs), not blocked on a new skip.

---

## 11. Phase B verdict — R_P=10 on *personal* hardware (the hard truth)

Goal: avoid the 1 TB / 120-core rental via a cheaper attack. After the SOTA (2025/259, CheapLunch
App. D.1) + direct probing, the verdict is **negative with current methods**, and precisely why:

1. **2025/259 (Improved Resultant) does not apply to Poseidon** — it targets *high-degree*-inversion
   ciphers (Griffin/Arion/Rescue/Anemoi, ℓ≥1). Poseidon's `x³` is low-degree (ℓ=0); its machinery is
   moot and it has zero Poseidon results.
2. **Round-skip is hard-capped.** CheapLunch App. D.1: for an **MDS** matrix the skip reduces `D_I`
   by `d^{min(k,⌊t/k⌋−2)} = d²` for us — *no better* than the validated 1-round skip (`D_I=3²⁰`). The
   only way past `d²` is a **non-MDS** matrix (their Poseidon2 weakness, +1 round). But
   [probes/04_is_mds.py](probes/04_is_mds.py) **proves our Plonky3 circulant IS MDS** → lever closed.
3. **Brutal cost** (skip-1, `D_I=3²⁰`, `δ=3¹⁵`): resultant `~2⁵⁵` time / **~800 TB** mem; streamed
   resultant `~14 GB` mem / `~2⁶⁷` time; GCD `~56 MB` mem / `~2⁶²` time. Calibration: 2026/150's R_P=4
   = 2⁴⁶/2.6 h/**10 TB**; +6 R_P ⇒ R_P=10 ≈ `2⁵⁶`/**~7 PB** on their rig. **No route fits a personal box**
   (needs ≤~2⁵⁴ ops *and* ≤~64 GB at once; each route misses an axis by 3–9 orders of magnitude).

**Conclusion:** R_P=10 is unclaimed because it sits past the frontier of cheap methods. On a personal
budget it requires a genuine **new idea** (bonus territory), not a known recipe. Open directions, none
guaranteed: (a) a fast bivariate solve with *sub-quadratic memory*; (b) an MDS skip beating the `d^k`
cap (CheapLunch leaves open); (c) a non-resultant attack exploiting the 10 near-linear partial rounds
or the NTT-friendly field. **Pragmatic alternative:** the **Density / Zero-test** tracks are
"best attack wins" (ranked) — incremental, personal-hardware-sized progress is rewarded, unlike the
all-or-nothing, first-come CICO R_P=10.

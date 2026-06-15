# Poseidon CICO Bounty — Theory & Attack Plan

*Working document. Phase 1: study the theory before writing attack code.*
*Last updated: 2026-06-10.*

> **Phase 2 is done:** all five modern attack papers are read and distilled into an
> implementable plan in **[ATTACK_RECIPE.md](ATTACK_RECIPE.md)** — start there for the
> concrete method (resultant-based bivariate solve), complexity per rung, and build plan.
> This file remains the parameter/permutation/attack-model reference.

---

## 0. TL;DR — what we are attacking and why it's winnable

The **open** CICO bounty is a deliberately **round-reduced, small-field** instance of the
Poseidon1 permutation. It is *not* full production Poseidon. Two harder-than-ours
instances of the same ladder have already fallen in the last two months, so the open
one sits right at the current practical frontier.

| | Value |
|---|---|
| **Problem** | CICO (Constrained-Input-Constrained-Output) on the permutation |
| **Open target** | `R_F = 6`, `R_P = 10` → **$15,000** |
| **Already gone** | `R_P = 6` ($6k, solved 2026-04-10) · `R_P = 8` ($10k, submitted 2026-05-06) |
| **Field** | KoalaBear, `p = 2³¹ − 2²⁴ + 1 = 2 130 706 433` |
| **State width** | `t = 16` |
| **S-box** | `x³` (α = 3) |
| **Expected attack** | *Resultant-based algebraic attack* (stated by the organizers) |

The brute-force success probability is `1/p² ≈ 2⁻⁶²` per random input, so a solution
**must** come from an algebraic/structural attack. That is the whole game.

Sources of truth in this repo:
- Reference implementation: [reference/poseidon/poseidon.py](reference/poseidon/poseidon.py)
- CICO verifier (defines exactly what counts as a win): [reference/bounties/cico_verifier.py](reference/bounties/cico_verifier.py)
- Official bounty spec: [reference/bounties/docs/bounty2026.tex](reference/bounties/docs/bounty2026.tex)
- Poseidon paper: [docs/2019-458.pdf](docs/2019-458.pdf)

> ⚠️ **Version skew to confirm before burning compute.** The live site and
> `bounty2026.tex` list `R_P ∈ {6,8,10}`; the repo `README.md` lists `R_P ∈ {8,10,12}`.
> The live site is authoritative and currently shows **RP=10 as the only open one**.
> Re-check the site the day we submit — the ladder shifts up as instances fall.

---

## 1. Exact problem statement

Let `P⁺ : 𝔽ₚ¹⁶ → 𝔽ₚ¹⁶` be the bounty permutation (defined in §2). Fixed constants:

```
C₁ = 0xC09DE4 = 12 623 332
C₂ = 0xEE6282 = 15 622 786     (both ≈ 2²⁴, i.e. live in the low bits)
C₃ = C₄ = 0
```

**Find** `(x₃, x₄, …, x₁₆) ∈ 𝔽ₚ¹⁴` such that

```
P⁺( C₁, C₂, x₃, x₄, …, x₁₆ )  =  ( 0, 0, *, *, …, * )
```

i.e. the **first two input words are pinned** to `(C₁,C₂)`, the remaining **14 words are
free**, and the **first two output words must be exactly 0** (the other 14 output words
are unconstrained). This is the `k = 2` CICO problem ("CICO-2"): 2 input constraints,
2 output constraints.

Degrees of freedom: 14 free unknowns, 2 polynomial equations ⇒ the solution variety has
dimension ≈ 12. Solutions are abundant; the difficulty is *finding one cheaply* despite
the equations having astronomical degree (`3¹⁶`) when written naïvely in the inputs.

The verifier ([cico_verifier.py:158](reference/bounties/cico_verifier.py#L158)) calls
`permutation_plus_linear` and checks `y[0]==0 and y[1]==0` exactly — no feed-forward, no
hashing wrapper. The CICO challenge lives **purely on the permutation**.

---

## 2. The permutation, precisely

### 2.1 Field
`p = 2³¹ − 2²⁴ + 1`. Then `p − 1 = 2²⁴ · 127`.
- `gcd(3, p−1) = 1` ⇒ `x ↦ x³` is a bijection (valid S-box). Its inverse exponent is
  `3⁻¹ mod (p−1)`, which is huge — so the **forward** direction is the cheap one.
- **2-adicity 24** ⇒ the field is extremely NTT-friendly. This matters a lot: the final
  step of the attack is finding roots of a high-degree univariate over 𝔽ₚ, and NTT-based
  root-finding (Graeffe iteration, eprint 2025/937) is cheap here.

### 2.2 Round schedule (`R_F = 6`, `R_P = 10`, total **16 rounds**)

From [poseidon.py:116-142](reference/poseidon/poseidon.py#L116-L142). With an **extra
initial linear layer** (this is what the `⁺` means — `permutation_plus_linear`):

```
state ← M · state                      # initial linear layer (the "+")
3 ×  full round                        # R_F/2 = 3
10 × partial round                     # R_P
3 ×  full round                        # R_F/2 = 3
```

Each round is `ARC → S-box → M` (constants, then S-box, then MDS):

- **Full round:** add round constants to all 16 words, cube **all 16** words, multiply by `M`.
- **Partial round:** add round constants to all 16 words, cube **only word 0**, multiply by `M`.

So nonlinearity count: `6·16 + 10·1 = 96 + 10 = 106` S-boxes total. The 10 partial-round
S-boxes are the structural soft spot (see §4).

### 2.3 Linear layer `M` — the **Plonky3 circulant**, not Cauchy
**Source authority: website > repo.** The live site's UPD says the matrix is "arbitrary
satisfying no-invariant-subspace-trail conditions" and names the **Plonky3 KoalaBear
circulant** as the example; the repo corroborates with Plonky3 test vectors. So the
operative matrix is the circulant `M[i][j] = first_row[(j−i) mod t]` with
`first_row = [1,1,51,1,11,17,2,1,101,63,15,2,67,22,13,3]` (Plonky3 `koala-bear/src/mds.rs`).
The Cauchy `(i−t−j)⁻¹` in `generate_mds_matrix` is only the repo's **code default** and is
subordinate. Both are MDS and **pass** the repo's anti-subspace-trail audit
(`verify_mds_matrix` = `_algorithm_1/2/3`, [mds_matrix.py:498-573](reference/poseidon/mds_matrix.py#L498-L573)).

**Key structural fact (verified, [probes/01_matrix_structure.py](probes/01_matrix_structure.py)):**
the Plonky3 circulant is **fully diagonalizable over GF(p)** — all 16 eigenvalues lie in
𝔽ₚ (because 16 | p−1), so the 16-point NTT diagonalizes it. The Cauchy matrix has **zero**
GF(p) eigenvalues (its char-poly factors are higher-degree). This is the one genuine
structural difference — but see §4 / the probe log: it does **not** appear to be an attack
lever, because the subspace-trail (Krylov) skip capacity is identical for both matrices
(`e₀` is cyclic either way, `ℓ_max = t−1 = 15`), and the eigenvalues carry no special
structure under the `x³` S-box (no `λᵢ=λⱼ³`, all pairwise ratios distinct).

### 2.4 Round constants
80-bit Grain LFSR ([grain_lfsr.py](reference/poseidon/grain_lfsr.py)), seeded by
`(field type, α, ⌈log₂p⌉, t, R_F, R_P)`. They are deterministic, public, and
nothing-up-my-sleeve. **Important modeling consequence:** the constants depend on `R_P`,
so the `R_P=10` instance has *different* constants than `R_P=8` — a solution does not
transfer between sub-instances; each must be solved fresh.

### 2.5 Folding the "+" away
`M` (initial) and `ARC₁` are known and affine, so the input map
`s ↦ ARC₁(M·s)` is a **known invertible affine map** `A`. Working in the variable
`u = A(s)` (the input to the very first S-box layer) is equivalent and cleaner: the two
pinned input words become **2 linear constraints on `u`**, and `u` ranges over a
**14-dimensional affine subspace** `V ⊂ 𝔽ₚ¹⁶`. This `V` is the freedom we spend in §4.

---

## 3. The algebraic attack model (the standard AO/CICO framing)

Writing the 2 output words directly as polynomials in `(x₃,…,x₁₆)` gives degree `3¹⁶ ≈
4.3·10⁷` — useless. Instead, **introduce one auxiliary variable per S-box output** and
keep every equation **degree 3**. This is the universal arithmetization-oriented (AO)
modeling:

- One variable per active S-box; one cubic equation `v = (lin-form)³` per S-box linking
  it to the (affine) state going in.
- The MDS and constants are absorbed into affine forms — they cost nothing in degree.
- The 2 CICO output equations close the system.

Then solve with **Gröbner basis** (`degrevlex → FGLM → univariate factor`, paper
App. C.2.2, [docs/2019-458.pdf](docs/2019-458.pdf) p.24) or — better here — **resultant
elimination** down to one variable, then root-find over 𝔽ₚ.

The cost is governed by the **degree of regularity** `D_reg` of the system, which scales
with the number of *non-skippable* nonlinear rounds. So every round we can linearize away
is a large multiplicative win. That is what §4 buys.

---

## 4. Where the weakness is — the levers

**Lever 1 — partial rounds are nearly linear.**
In a partial round, 15 of 16 words are an affine function of the previous state; only
word 0 is cubed. So algebraic degree grows *slowly* through the partial section, and the
linear words can be expressed as affine forms for free. All the security in the partial
section comes from those 10 lone S-boxes — that's the entire `RP=8 → RP=10` gap.

**Lever 2 — the "skipping rounds" / subspace-trail trick** (paper Eq. (1), Eq. (8);
[BCD⁺20]; subspace-trail GB eprint 2025/954; hinted explicitly by the organizers).
Define `S⁽ⁱ⁾ =` inputs for which **no S-box is active in the first `i` rounds**. For the
**partial** section, "word-0 S-box input is constant" is *one linear constraint per
round*, so requiring it for the first `i` partial rounds cuts a codim-`i` affine subspace
— and on that subspace those `i` rounds collapse to a **known affine map**. We hold **12
spare linear degrees of freedom** in `V` (14 free − 2 output constraints), so we can
**pin several early partial-round S-box inputs to constants and skip those rounds outright**,
shrinking the effective `R_P` that the algebra has to handle. Combined with the
inside-out / coset start (begin "in the middle" so end rounds are affine on a coset),
this is exactly how the `RP≤8` instances were broken. (We cannot fully skip a *full*
round — that needs all 16 S-box inputs fixed = 16 constraints > 12 spare — so `R_F=6`
stays as the algebraic core. The win is in compressing `R_P`.)

**Lever 3 — small NTT-friendly field.**
After elimination we get a high-degree univariate over a 31-bit field with 2-adicity 24.
Root-finding is cheap (Graeffe/NTT, eprint 2025/937). On a 255-bit field the same step
would dominate; here it nearly free. This is *why the bounty uses KoalaBear* and why the
instances are practically breakable at all.

**Lever 4 — the `⁺` layer and constants add no security.** The initial `M` is invertible
and known (fold it into `A`, §2.5); constants are public affine shifts. Neither raises
`D_reg`. They only fix *which* affine subspace `V` is and *which* concrete univariate we
root-find — relevant for engineering, not for hardness.

---

## 5. State of the art (the attacks we should build on)

| Work | Idea | Relevance to us |
|---|---|---|
| **Poseidon paper App. C.2** ([docs/2019-458.pdf](docs/2019-458.pdf)) | GB attack framing, `D_reg`, the `S⁽ⁱ⁾` subspace skip (Eq. 1, 8) | Baseline model + the skip trick, in the designers' own words |
| **FreeLunch** (Bariant–Boissier–Bouvier–Leurent et al., CRYPTO 2024) | Pick a monomial order so the CICO system is *already* a Gröbner basis ⇒ the costly GB step is **free**; only FGLM + root-find remain | The core efficiency trick; turns GB cost into root-finding cost |
| **CheapLunch** (eprint **2025/2040**) | Extends FreeLunch **beyond CICO-1 to multiple output constraints** | **Directly our case** — we have `k=2` output constraints (CICO-2) |
| **Subspace-trail GB on Poseidon/Neptune** (eprint **2025/954**) | Use subspace trails to linearize partial rounds inside the GB attack | Formalizes Lever 2 for exactly this primitive |
| **"Claiming bounties… resultant-based attacks"** (eprint **2026/150**, Bak–Bariant–Boeuf–Hostettler–Jazeron, Feb 2026) | The team that **actually claimed** the 31k/31m CICO ladder; eliminate to a univariate via **resultants** instead of full GB | The blueprint for the *exact* family we're attacking |
| **Graeffe root-finding over NTT fields** (eprint **2025/937**) | Fast roots of high-degree univariates over NTT-friendly primes | The cheap final step on KoalaBear |

> Update (2026-06-10): all five PDFs (FreeLunch = **2024/347**, CheapLunch 2025/2040,
> subspace-trail 2025/954, resultant 2026/150, Graeffe 2025/937) are now downloaded to
> [docs/](docs/) and read in full. The distilled, implementable method — a **resultant-based
> bivariate solve** with round-skipping — is in **[ATTACK_RECIPE.md](ATTACK_RECIPE.md)**,
> including per-rung complexity and the bonus lever (skip >1 round).

---

## 6. Proposed plan (theory → code bridge)

1. **Read the four modern papers in full** (esp. 2026/150 resultant pipeline + CheapLunch
   CICO-2 setup). Extract their exact variable placement and elimination order.
2. **Stand up an exact reference + tiny oracle.** Wrap [poseidon.py](reference/poseidon/poseidon.py)
   and confirm against the verifier; reproduce the *relaxed* solver (the repo's
   `m`-bit demo) to lock down byte-for-byte agreement on `permutation_plus_linear`.
3. **Build the polynomial model in Sage/`msolve`** (one var per S-box, MDS+constants as
   affine forms). Validate end-to-end on **toy `R_P = 1, 2`** where a full GB is instant.
4. **Implement Lever 2** (skip the first few partial rounds via linear constraints on `V`)
   and measure how far `D_reg` / solving time drops per skipped round.
5. **Scale `R_P` 2 → 10**, switching GB → resultant elimination + NTT root-finding as the
   degree grows. Cross-check intermediate solutions with the relaxed verifier at growing
   `m` before going for exact equality.
6. **Submit** `(x₃,…,x₁₆)` for `R_F=6, R_P=10`; re-confirm the live ladder first.

Even a clean, novel write-up of an improved skipping variant is explicitly eligible for a
**bonus** per the organizers — so partial progress is not wasted.

---

## 7. Open questions / risks
- Exact `R_P` ladder (`{6,8,10}` vs `{8,10,12}`) — confirm on the live site at submit time.
- Does 2026/150 already cover the **Poseidon1** CICO track, or only Poseidon2? The repo
  spec says Poseidon1; the paper abstract says Poseidon2-31k/31m. Resolve which track the
  `$15k` instance belongs to and whether the published attack transfers as-is.
- Memory: full GB at `R_P=10` may blow up; the resultant route is the hedge — verify its
  univariate degree stays within NTT-root-finding reach on KoalaBear.
```

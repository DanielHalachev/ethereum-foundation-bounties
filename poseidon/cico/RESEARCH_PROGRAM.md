# Research Program — a non-generic weakness in Poseidon-31 CICO

*Committing to a genuinely new idea. The goal: a shortcut that beats the `3^r` algebraic-degree
wall by exploiting hidden structure the standard Gröbner/resultant attacks ignore. We deliberately
reach into math not usually applied to symmetric cryptanalysis. Long shots are expected; we rank by
plausibility and always attach a concrete first probe.*

Target recap: `F_p`, `p = 2³¹−2²⁴+1 = 127·2²⁴+1`; state `t=16`; S-box `x³`; linear layer = Plonky3
**circulant MDS** `M` (diagonalizable over `F_p`, eigenvalues = NTT of first row); CICO-2 (2 input
words pinned, 2 output words = 0); `R_F=6`, `R_P=10`. Generic cost ≈ `2⁵⁵`+ and PB memory (see
[ATTACK_RECIPE.md](ATTACK_RECIPE.md) §11). We need structure, not compute.

## The three peculiarities to weaponize
1. **Circulant ⇒ `C₁₆` symmetry.** `M` is multiplication in the group algebra `F_p[C₁₆]=F_p[z]/(z¹⁶−1)`.
   The S-box (cube each coord) is *also* `C₁₆`-equivariant (coordinatewise). So the **constant-free
   permutation commutes with cyclic shift** — the round constants are the *sole* symmetry-breaker.
2. **`x³` (multiplicative) vs MDS (additive).** The round function alternates the two ring structures;
   their incompatibility is exactly what additive-combinatorics / sum-product theory measures.
3. **`p−1 = 2²⁴·127`.** Enormous 2-adicity; pinned constants `C₁,C₂ ≈ 2²⁴`. Smells exploitable by
   lattice (Coppersmith) / 2-adic / valuation methods.

---

## Lens 1 — Harmonic analysis on `F_p[C₁₆]` (NTT eigenbasis) ★ recommended start
**Idea.** Work in the NTT basis where `M` is *diagonal*. Derivation (verified):
- Full round: cube-each-coord becomes a **cyclic triple-convolution** of `ŝ` (a symmetric trilinear form).
- **Partial round collapses beautifully:** `ŝ_j^{(i+1)} = λ_j·(ŝ_j^{(i)} + ĉ_j^{(i)} + Δ_i)`, where
  `Δ_i = m_i³ − m_i` and `m_i = (1/16)Σ_j(ŝ_j^{(i)}+ĉ_j^{(i)})` is the **mean** (DC component). I.e. each
  of the 10 partial rounds is *diagonal scaling + a rank-1 broadcast of one scalar cubic*. The whole
  partial section is affine in `(input, Δ₀..Δ₉)` with a triangular cubic recurrence in the scalars `m_i`.
**Why it might bite.** The instance-specific structure (the `λ_j`, the symmetry) is laid bare here and
invisible to degree-counting. If the `λ_j` spectrum or the `m_i` recurrence has special structure
(e.g. resonances, a closed-form for the `m_i` generating function), rounds could collapse.
**Plausibility:** medium — most concrete, fully computable, exploits the one real peculiarity.
**First probe:** symbolically rewrite the permutation in the NTT basis; print the `m_i` recurrence and
the `λ_j`; quantify how many isotypic components the round constants actually perturb; test whether a
choice of the 14 free inputs can zero the symmetry-breaking in enough components to linearize rounds.

## Lens 2 — Additive combinatorics / sum-product
**Idea.** `x³` structures sets multiplicatively; MDS structures them additively. Sum-product
(Bourgain–Katz–Tao, Elekes) says a set can't be highly structured under both. Conversely, *special*
sets (geometric progressions, multiplicative cosets, subfield-like sets) are where the tension is
weakest — solutions might concentrate there.
**Why it might bite.** A "weak input class" (inputs drawn from a structured set) could shrink the search
from `p²` to something tractable, or expose a low-degree relation.
**Plausibility:** lower for *constructive* solving (sum-product usually yields distinguishers/bounds),
but genuinely under-applied to AO ciphers.
**First probe:** test whether restricting free inputs to a geometric progression / a multiplicative
subgroup coset (note `127 | p−1`, so there's an order-127 subgroup) makes the output-capacity polys
degenerate (lower degree, or factor).

## Lens 3 — Lattice / Coppersmith with the special prime ★ cheap to test
**Idea.** `p = 127·2²⁴+1`; `C₁,C₂ ≈ 2²⁴`. Coppersmith finds small roots of modular polynomials; LLL
finds short integer relations. Formulate CICO as modular polynomial relations and hunt for
small-norm / low-bit-height solutions exploiting the `2²⁴` structure.
**Why it might bite.** If solutions (or intermediate quantities) have non-generic 2-adic size, lattice
reduction finds them far faster than algebra.
**Plausibility:** medium-low (CICO solutions aren't obviously small), but **cheap and independent** —
a good parallel bet.
**First probe:** empirically scan whether any CICO solution exists with small inputs (low bit-length),
or whether the relaxed-verifier solutions cluster in bit-size; build an LLL basis for the linearized
relations and inspect.

## Lens 4 — p-adic / tropical (valuation) geometry
**Idea.** Use the huge 2-adicity: analyze intermediate quantities by 2-adic valuation, or tropicalize
the polynomial system (Newton polytopes / min-plus) to see if its tropical variety has lower effective
dimension/complexity than the naive degree.
**Plausibility:** low-medium, speculative but novel for this setting.
**First probe:** compute 2-adic valuations of the `Δ_i`/`m_i` sequence on random inputs; look for
forced valuation patterns. Tropicalize the 2-variable CICO system; inspect the Newton polytope.

## Lens 5 — Invariant theory / perturbed `C₁₆` symmetry
**Idea.** The naked permutation is `C₁₆`-equivariant; round constants break it. Treat constants as a
perturbation and look for a *relative* invariant or twisted-equivariant structure that survives, or
choose free inputs to restore symmetry on a subspace.
**Plausibility:** low-medium; elegant, but the constants are full-size, not a small perturbation.
**First probe:** decompose the constants into `C₁₆` isotypic components; check if any component is
anomalously small/zero (a design oversight) that a symmetry argument could exploit.

## Lens 6 — Algebraic geometry of the CICO variety (fibration / rational curves) ★ highest upside
**Idea.** The 12-dim solution variety comes from an *iterated* low-degree map — it likely has a tower /
fibration structure. Finding a **rational curve or low-dim rational subvariety** inside it gives a cheap
parametrization of infinitely many solutions.
**Why it might bite.** If even a 1-parameter family of solutions can be written down, we win outright.
**Plausibility:** medium; highest payoff if the structure exists; hardest to find.
**First probe:** at tiny `R_F/R_P`, compute the solution variety's structure (is it rational? unirational?
does the iterated map induce a fibration?) with the symbolic tools we already have.

---

## Working method
- Keep [probes/](probes/) reproducible; log every result (positive or negative) here.
- Bias toward **reformulations** — a new idea usually appears when a change of representation makes a
  hidden structure visible. Lens 1 (NTT) is our computational backbone for that.
- Cross-pollinate: the NTT `m_i` recurrence (L1) is the natural object for the valuation (L4),
  sum-product (L2), and AG (L6) lenses to act on.
- Honesty rule: a lens is "alive" only while it produces falsifiable predictions we can probe. Kill
  dead ones fast (as we killed the diagonalizability lever in Probe 1).

---

# Findings log

## Round 1 (2026-06-11)

### Theory hunt verdict (subagent survey of unconventional areas)
Ruthless prioritization after surveying arithmetic dynamics, Krylov/Sparse-FGLM, tensor
decomposition, sum-product, p-adic/tropical, isogeny:
- **KILLED as constructive paths:** Lens 1 spectral-resonance (see probe 05), arboreal Galois
  (affine conjugation destroys the monomial/PCF structure needed for easy preimages), tensor rank
  (cyclic-cubing tensor is full-rank; rank wouldn't invert anyway), sum-product (obstruction/
  distinguisher only — it actually *proves* the design safe against invariant-coset attacks),
  p-adic (cube root is already free since gcd(3,p-1)=1; 2-adicity merely gifts us the NTT basis).
- **SURVIVING constructive lead = Lead 2 (Krylov / Sparse-FGLM / block-Wiedemann in the eigenbasis).**
  Our moment-shift reformulation is the *dual* of the Sparse-FGLM change-of-order machinery
  (Faugère–Mou, arXiv:1304.1238). Two payoffs:
  (a) **likely breaks the MEMORY wall** — Sparse-FGLM/Wiedemann uses **O(D_I) memory** (matrix-vector
      products + Berlekamp–Massey), not the resultant's δ². For D_I=3²⁰ that's ~tens of GB, not PB.
      2026/150 only weighed resultant (δ² mem) vs GCD (O(q)=2³¹ time) — **they did not list a
      Wiedemann route**; it may give linear memory at ~2⁵⁵ time = *cheap-CPU* feasible (no 1 TB RAM).
  (b) **superfast dream (novel, uncertain):** if the multiplication matrix T_X has low **displacement
      rank** (because the partial-round perturbation is rank-1/round), a structured Hankel/Toeplitz
      solver gives Õ(D_I)≈2³⁵ time. This is the genuine "unexpected theorem" bet (structured-matrix
      theory → AO cryptanalysis); not in the literature for Poseidon.
- **Bonus unexpected angle = Lead 6(b):** view the cubic kicks as *faults* in an otherwise
  linearly-recurrent moment stream → Berlekamp–Massey–Sakata / list-decoding recovery (coding theory).
  Constructive in spirit, shares tooling with Lead 2. (Caveat: 10 kicks aren't obviously "low-weight".)
- **New on-target papers to read:** CheapLunch (2025/2040, *the* CICO-2 paper), subspace-trail GB
  (2025/954), "Revisiting Linear Subspace Trails in Poseidon" (2026/967 — rank-growth tied to MDS
  eigenstructure = our object), "Skipping Class" (2026/306 — weak matrices/modes).

### Probe 05 — NTT spectrum is structureless (negative for spectral resonance)
[probes/05_ntt_spectrum.py]: λ_j = 16 distinct, none 0/1; **no multiplicative relations**
(a·b∈set: 0/256, a²∈set: 0/16, a·b=c: 0), **no additive relations** (0 zero-pairs/triples),
**full Krylov dim** for the all-ones kick (10 kicks→dim 10), char poly splits into 16 distinct F_p
modes (no order reduction). 2-adic levels e∈{0..4}; two small eigenvalues (151, 371) — no use found.
⇒ no cheap spectral collapse. The clean reformulation stands as the substrate for Lead 2.

### Refined focus
The memory wall (not time) is what made R_P=10 unaffordable (10 TB+ RAM). **Lead 2's Wiedemann route
plausibly converts that to O(D_I)≈tens-of-GB memory at ~2⁵⁵ time = cheap many-core CPU**, and the
displacement-structure refinement could cut time to ~2³⁵. Next probe: build the CICO-2 multiplication
matrix at tiny rounds; measure sparsity + displacement rank + Wiedemann linear complexity to test (a)
memory-feasibility and (b) the superfast dream.

### Probe 06 — multiplication-matrix displacement structure (Lead 2 superfast test): NEGATIVE
[probes/06_mult_matrix_structure.py], CICO-2 multiplication matrix T_X in the grevlex basis:
| matrix | D_I | nnz(T_X) | displacement rank (down-shift) |
|---|---|---|---|
| Plonky3 circulant | 81 | 1273 | 18 |
| Plonky3 circulant | 243 | 5798 | 42 |
| Cauchy default | 81 | 1273 | 18 |
| Cauchy default | 243 | 5798 | 42 |

Two conclusions:
1. **MATRIX-INDEPENDENT** — Plonky3 and Cauchy give byte-identical nnz and displacement rank.
   T_X structure does not depend on the linear layer (just like D_I, probe 02). The circulant's
   diagonalizability is irrelevant to solving cost in every metric measured. (Answers the
   "did you test both matrices?" control: yes — no Plonky3 advantage exists.)
2. **dr GROWS with D_I** (18->42 for 81->243, ~D_I^0.77, tracking the staircase border δ_X+δ_Y).
   A displacement-structured solve therefore costs Õ(dr·D_I) ≈ Õ(δ·D_I) = the resultant's ~2^55 —
   it re-derives, not beats, the known complexity. The superfast Õ(D_I)≈2^35 dream is DEAD in this basis.

### Status after Round 1: degree/algebra avenues exhausted
Everything tied to the F_p algebraic-degree structure (D_I, T_X, eigenbasis, skip) is
matrix-independent and bottlenecked at D_I=3^20 / ~2^55. No spectral, displacement, tensor,
arboreal, or sum-product shortcut survived. A genuine break must come from a NON-degree avenue:
- Lens 3 (lattice/Coppersmith — integer/2-adic size, not F_p degree) — UNTESTED, cheap.
- Lens 6 (AG: a rational curve / fibration of the solution variety — bypasses degree) — UNTESTED, highest upside.
- Lens 2 (structured input class / order-127 subgroup) — UNTESTED.
These don't reduce to D_I the way the algebraic attacks do; they are the live frontier.

## Round 2 (2026-06-11) — Algebraic Geometry (user-chosen direction)

### Probe 07 — reducibility / gcd / Ritt decomposition: GENERIC (negative)
[probes/07_ag_reducibility.py], both matrices, (RF,RP) in {(2,0),(2,1)}:
- gcd(P,Q) = 0 (coprime) in ALL cases — no common curve component.
- P, Q irreducible (varying univariate-slice profiles with a large irreducible factor ⇒ irreducible bivariate).
- Resultant R(X) **Ritt-indecomposable** in all cases (no g∘h composition structure).
- R(X) factor profiles generic (one big irreducible + a few small; 1–2 rational roots), as for a random poly.
- Plonky3 vs Cauchy differ only in cosmetic factor degrees — same generic CHARACTER. No matrix weakness.

### AG point-finding literature verdict — fundamentally BLOCKED by degree >> field size
Subagent surveyed Cafure–Matera (constructive F_q-point of a variety), effective Lang–Weil,
unirationality (Kollár), Tsen–Lang/C_i, determinant method (Bombieri–Pila), Weil descent (GHS),
BKK/tropical. **All fail constructively here for one structural reason:** these theorems need the
field large vs the variety degree (Cafure–Matera needs q > 8n²dδ⁴; for δ=3^16..3^20 that's q>2^113..2^139,
violated by 83–108 bits). AO ciphers are deliberately instantiated over fields SMALL relative to the
degrees they generate (3^16 over 2^31) — **exactly the regime where every AG construction theorem dies.**
AG gives existence (~p^12 points) and beautiful structure, but no cheaper CONSTRUCTION than the
existing crypto attacks. No AG tool has ever broken an AO/SPN hash, for this reason.

### Useful correction + pointer from the hunt
- The realistic SOTA is resultant/elimination + fast univariate root-finding. The **tangent-Graeffe
  root-finder (eprint 2025/937) is purpose-built for our EXACT prime p=127·2^24+1** (σ·2^m+1 form) and
  manages memory in the root-find step — relevant if/when we reach the final univariate.
- CAVEAT: the subagent's "~2^40, 14 GB total" applies to CICO-**1** (single output → pure univariate
  root-find). Our problem is CICO-**2** (two outputs, overdetermined) → still needs the bivariate
  elimination at D_I=3^20, ~2^55 time with the memory tradeoff. CICO-2 cost is UNCHANGED.

### Meta-status after Round 2
Every degree/algebra/AG avenue collapses to the same 3^20 / ~2^55 wall, and the COMMON CAUSE is now
explicit: **degree(3^k) >> field(2^31)**. The only escape is an invariant that is NOT the F_p algebraic
degree. Candidates: Lens 3 (lattice/Coppersmith — integer/2-adic SIZE) — UNTESTED, and the principled
next bet since it sidesteps degree>>field; Lens 6b (coding/BMS). And the user's specific (unrevealed) AG theorem.

## Round 3 (2026-06-11) — last non-degree checks (both negative)
[probes/08_constants_and_size.py]:
- (A) Constants C1=0xc09de4, C2=0xee6282 are NUMS: 24-bit hex, C1 full order, neither in order-127
  subgroup, no cube relation, no tie to round constants. (C2 is a 16th power — cosmetic.) No backdoor.
- (B) Lattice/2-adic hook ABSENT: relaxed-verifier solution inputs are uniform integers (0/28 below
  2^24; mean v2 = 0.75 ≈ uniform). No smallness/2-adic bias for LLL/Coppersmith. Lens 3 is vacuous.

## OVERALL VERDICT after 3 rounds / ~9 lenses
No personal-hardware break found, and the reason is now structural and well-understood:
  degree(3^k) >> field(2^31)  +  generic system (irreducible/coprime/indecomposable, matrix-independent)
  +  optimal MDS (no subspace-trail, skip capped at d^k)  +  uniform solutions (no lattice hook).
Every avenue (Gröbner/resultant/FGLM, NTT spectral, displacement/Wiedemann, round-skip, arboreal,
tensor, sum-product, p-adic, AG point-finding, lattice) collapses to the same 3^20 / ~2^55 wall, by
design. This matches reality: world-class teams (Bariant et al.) stopped at R_P=8 for this bounty.
A break needs a genuinely new idea beyond this multi-front search — open research, not a known recipe.
Documented tooling (probes/ + .venv) and the validated facts (skip /d^2, matrix-independence, MDS,
NUMS constants) are reusable. Recommend: (i) write up the structural analysis as a contribution to the
initiative (the bounty explicitly values documented findings), and/or (ii) pivot to the Density /
Zero-test tracks ("best attack wins" — incremental, personal-hardware-sized progress is rewarded).

## Round 4 (2026-06-11) — statistical cryptanalysis (all negative)
[probes/09_statistical.py], full instance RF=6,RP=10 unless noted:
- (1) EXACT algebraic degree of out_0 = 3^(RF+RP) in every reduced case (9,27,81,81) — FULL, no
  degree collapse. Kills integral / higher-order-differential / zero-test degree shortcuts.
- (2) Diffusion COMPLETE — out_0 and out_1 each depend on all 14 free inputs (frac 1.0).
- (3) Output uniformity: mean/p=0.4999, chi^2=298.5 over 256 buckets (expect ~255+-23) — uniform.
- (4) Differential (fixed input diff): P(out_0 diff=0)=0/30000 (~uniform), max diff multiplicity 1 —
  no usable differential.
=> No statistical weakness; matches the design intent of R_F=6 full rounds. Statistical attacks would
   in any case give a DISTINGUISHER, not a constructive CICO solution. Front closed.

# Poseidon-1 CICO ($R_F=6$, $R_P=10$)

Challenge can be found at [poseidon-initiative.info](https://www.poseidon-initiative.info/#h.l22wkxnie0q7).

> **Field**: KoalaBear $2^{31}-2^{24}+1 = 2130706433$. Degree $d$ of power mapping is $3$, the state size $t$ is $16$.
>
> **UPD**: the MDS matrix is arbitrary but satisfying "no-invariant-subspace-trail" conditions.  A example is the [Plonky3 circulant matrix](https://github.com/khovratovich/poseidon-tools/blob/main/tests/test_poseidon.py#L339).
>
> The task is to find a partial preimage of $0$, or, more precisely:
>
> For Poseidon-31 a 62-bit preimage: find $X_1,\dots, X_{14}, Y_1,\dots, Y_{14}$   such that:
> $$
> \text{Perm}(0xc09de4,0xee6282,X_1,...,X_{14})= (0,0,Y_1,...,Y_{14})
> $$
>
> where `Perm` is the inner sponge permutation (bijective mapping) of Poseidon1 with an extra linear layer before the first round.
>
> We encourage cryptanalysts to find an improved attack variant (such as “skipping first rounds” trick)  rather than to find a solution with a brute force. New attack ideas might qualify for a bonus.
>
> Concrete bounties (details here):
>
> - ~~RF=6, RP=6 $6000  Solution verified on 10 April 2026~~
> - ~~RF=6, RP=8 $10000 Solution submitted on 6 May 2026~~
> - RF=6, RP=10 $15000
>
> We expect that the best attack that solves these bounties is a resultant attack. A Groebner basis attack that breaks any of these instances may qualify for an additional bonus.

## 1. Setup and notation

- field $\mathbb{F}_p$
- $p = 2^{31}-2^{24}+1 = 2130706433 \implies p-1 = 2^{24}\cdot 127 \implies \gcd(3,p-1)=1 \implies $x\mapsto x^3$ is a bijection.
- state width $t=16$
- S-box exponent $d=3$
- full rounds $R_F=6$ (split $3$ before and $3$ after the partial layer)
- partial rounds $R_P=10$.
- let $M\in\mathbb{F}_p^{16\times 16}$ be the fixed mixing matrix
- let $c_1,\dots,c_{16}\in\mathbb{F}_p^{16}$ the round constants.
- let $C_1=\texttt{0xc09de4}$, $C_2=\texttt{0xee6282}$

The permutation `permutation_plus_linear` is $P^+ = \Phi_{16}\circ\cdots\circ\Phi_1\circ(s\mapsto Ms)$, where for $\sigma(z)=(z_0^3,\dots,z_{15}^3)$ and $\sigma_0(z)=(z_0^3,z_1,\dots,z_{15})$,
$$
\Phi_i(s)=M\,\sigma(s+c_i)\ \text{(full round)},\qquad \Phi_i(s)=M\,\sigma_0(s+c_i)\ \text{(partial round)}.
$$

**CICO-2 problem.** With , find $(x_2,\dots,x_{15})\in\mathbb{F}_p^{14}$ such that
$$
P^+(C_1,C_2,x_2,\dots,x_{15}) = (0,0,\ast,\dots,\ast).
$$
The solution variety has dimension $14-2=12$, i.e. $\approx p^{12}$ points; the difficulty is exhibiting one. Brute force succeeds with probability $p^{-2}\approx 2^{-62}$.

The verifier is `reference/bounties/cico_verifier.py::verify_cico_solution`. The matrix specified by the website task is the Plonky3 KoalaBear circulant with `first_row = [1,1,51,1,11,17,2,1,101,63,15,2,67,22,13,3]`; the repo verifier additionally *defaults* to a Cauchy matrix $M_{ij}=(i-t-j)^{-1}\bmod p$ as a placeholder and accepts a supplied matrix. Both were tested and give identical degree/structure data (Section 4), so $M$ is immaterial to hardness, but a submission must use the circulant.

## 2. Definitions and theorems

The attack (Section 3) and the obstruction arguments (Sections 4–5) draw on the tools below. The notation fixed here — $D_I$, $\delta$, $d_{\mathrm{reg}}$ — is used throughout.

**MDS matrix.** $M\in\mathbb{F}_p^{t\times t}$ is MDS if every square submatrix is invertible (equivalently, the code $\{(v,Mv)\}$ has minimum distance $t+1$). Consequence: a nonzero input cannot vanish on more than $t-|J|$ output coordinates — so $M$ admits no invariant subspace trail and diffuses densely.

**Resultant / elimination degree.** For $f,g\in K[X][Y]$, the resultant $\operatorname{Res}_Y(f,g)\in K[X]$ eliminates $Y$: it vanishes at $x_0$ iff $f(x_0,Y),g(x_0,Y)$ share a root (or both lead coefficients do). Degree bound $\deg_X\operatorname{Res}_Y(f,g)\le \deg_X f\,\deg_Y g+\deg_Y f\,\deg_X g$. We write $D_I:=\deg_X\operatorname{Res}_Y(P,Q)$ for the two CICO output polynomials $P,Q$. By **Bézout**, two plane curves of degrees $a,b$ with no common component meet in $ab$ points, so generically $D_I\le ab$.

**Gröbner basis / degree of regularity.** For a $0$-dimensional system the solving cost is set by $d_{\mathrm{reg}}$, the top degree reached in an $F_4/F_5$ run: dense linear algebra costs $O\big(\binom{n+d_{\mathrm{reg}}}{d_{\mathrm{reg}}}^{\omega}\big)$, $2\le\omega\le 3$.

**Wiedemann.** For sparse $A\in\mathbb{F}_q^{N\times N}$ with $\omega$ nonzeros, Berlekamp–Massey on the sequence $(u^{\mathsf T}A^i v)_i$ recovers $A$'s minimal polynomial and solves $Ax=b$ in $O(N)$ matrix–vector products — time $O(N\omega+N^2)$, space $O(N+\omega)$. (Block variants cut the iteration count.)

**Displacement rank.** With $Z$ the down-shift matrix, $\nabla A:=A-ZAZ^{\mathsf T}$ and $\operatorname{drank}(A):=\operatorname{rank}\nabla A$. If $\operatorname{drank}(A)=\alpha$ then $Ax=b$ is solvable in $\tilde O(\alpha N)$ (Toeplitz: $\alpha\le 2$) — the only route to a near-linear solve.

**Counting theorems** (used to argue the obstruction, not to solve):
- *Chevalley–Warning:* if $\sum_i\deg f_i<n$ then $\#\{f_i=0\}\equiv 0\pmod p$; if that set is nonempty its size is $\ge q^{\,n-\sum_i\deg f_i}$.
- *Lang–Weil:* for an absolutely irreducible $V\subseteq\mathbb{A}^n$ of dimension $r$ and degree $\delta$, $\big|\#V(\mathbb{F}_q)-q^{r}\big|=O(\delta^2\,q^{\,r-1/2})$ — the error growing as $\delta^2$.
- *Bernstein/BKK:* the number of isolated zeros in $(\mathbb{C}^\ast)^n$ is at most the mixed volume of the Newton polytopes (equality for generic coefficients); homotopy continuation tracks one path per root.

**Belief propagation.** Message passing on a factor graph: on a tree its fixed point gives the exact marginals, but on a graph with cycles it is only a heuristic.

**Hellman time–memory tradeoff.** Inverting $f$ on $N=2^n$ points needs $\Theta(N)$ precomputation (unavoidable), after which online time $T$ and table size $M$ obey $TM^2=N^2$.

## 3. The resultant attack and its complexity

Fix $x_4,\dots,x_{15}$ and set $X=x_2$, $Y=x_3$. Evaluating $P^+$ symbolically turns the two pinned outputs into $P(X,Y),Q(X,Y)\in\mathbb{F}_p[X,Y]$ of total degree $d^{R_F+R_P}=3^{16}$. Eliminating $Y$ gives $R(X)=\operatorname{Res}_Y(P,Q)$; its $\mathbb{F}_p$-roots, lifted back via $\gcd(P(x_0,Y),Q(x_0,Y))$, are the solutions.

**Elimination degree.** The trivial Bézout value is $d^{2(R_F+R_P)}=3^{32}$, but each partial round contributes only degree $d$, giving the Poseidon drop
$$
D_I = d^{\,2R_F+R_P}=3^{22}\approx 2^{34.9}\quad[\text{BSGL20, Thm 10}].
$$
We reproduced this directly with `sympy` over $\mathbb{F}_p$ on both matrices: $(R_F,R_P)=(2,0)\Rightarrow D_I=81=3^4$ and $(2,1)\Rightarrow 243=3^5$ (not the trivial $729$). See `probes/10_ground_truth.py`, `probes/02_ideal_degree.py`.

**Round-skip (one round).** Since $\gcd(d,p-1)=1$, choosing the free inputs as $a\,X^{1/d}$ makes the round-$1$ S-box output $a^d X$ affine in $X$, removing round $1$ from the degree growth:
$$
D_I^{(\mathrm{skip}\text{-}1)}=d^{\,2(R_F-1)+R_P}=3^{20}\approx 2^{31.7},\qquad \delta:=\deg P=\deg Q=d^{\,R_F-1+R_P}=3^{15}.
$$
The skip stops at one round. Entering round $2$ the state is $M\cdot(\text{affine})$, so every coordinate is $\alpha_i+\beta_iX+\gamma_iY$ with $\beta_i,\gamma_i\ne 0$ (MDS makes it dense); its cube has total degree $3$ and is not $a\,W^3$ for any single linear form $W$, so no cube-root substitution linearizes round $2$ — we verified all $16$ round-$2$ S-box inputs are genuinely bivariate (`probes/wf_two_sided_skip.py`). The output side gives nothing either: $\mathrm{out}_0,\mathrm{out}_1$ are linear in $\sigma(s_{15}+c_{16})$, so stopping before the final $M$ trades the degree-$D_I$ polynomial for a degree-$3$ condition on $s_{15}$ — and $\deg_{X,Y}s_{15}=D_I/d$, so the product is again $D_I$, unchanged.

**Cost.** Evaluation–interpolation of $R$ costs $\tilde O(\delta\,D_I)=\tilde O(3^{35})\approx 2^{55.5}$ field operations; streaming keeps memory at $O(D_I)\approx 2^{31.7}$ coefficients ($\sim 15$–$30$ GB), avoiding the $\delta^2$ of a dense Sylvester solve. Root-finding $R$ via $\gcd(X^p-X,R)$ is subdominant. The end-to-end solver `probes/11_resultant_solver.py` returns `verify_cico_solution -> True` for $(R_F,R_P)=(2,0)$ and $(2,1)$ with the real $C_1,C_2$; only the scale-up to $3^{20}$ is missing.

**Guess-and-determine doesn't help.** Fixing $g$ partial-round S-box inputs reduces $R_P\to R_P-g$, but each fixed input is a free $\mathbb{F}_p$ value, so total cost is $p^{g}\cdot 3^{\,3(R_F-1)+2(R_P-g)}=2^{31g}\cdot 3^{35-2g}$. Since $31\ln 2>2\ln 3$ this grows with $g$; the optimum is $g=0$ (`probes/wf_guess_determine.py`).

## 4. Methods attempted and their obstructions

**Structural / algebraic.** Every one collapses to the same $D_I=3^{20}$, and all are matrix-independent — the Plonky3 circulant and the Cauchy default give identical degree and structure data.

- *Subspace trails / eigenstructure.* The circulant is diagonalized by the $16$-point NTT, but its eigenvalues admit no relation $\lambda_i=\lambda_j^3$ and no additive degeneracy; the Krylov reach of $e_0,e_1$ under $M$ is full ($=16$), as for Cauchy. MDS forbids an invariant trail (Section 2). `probes/01_*, 05_*`.
- *Krylov / sparse linear algebra* (Wiedemann, block-Wiedemann, sparse-FGLM). Cuts storage to $O(N+\omega)$ with $N=D_I=3^{20}$ — a genuine RAM win — but time is $O(N\omega)$, and the multiplication matrix's $\omega$ grows $\propto D_I$, giving $\approx 2^{62}$. The memory wall is merely traded for a time wall. `probes/06_*`.
- *Displacement-rank / structured solve* — the only route to $\tilde O(D_I)$. The measured $\operatorname{drank}$ of the Sylvester/multiplication matrix grows with $D_I$ ($18$ at $D_I=81$, $42$ at $243$), exactly as for a random pair; no Toeplitz/Hankel structure survives the partial rounds (the cubed coordinate spreads to all $16$ words at the next $M$ and re-triples every round: $3,9,27,\dots,6561$ measured). `probes/06_*`.
- *Reducibility / decomposition.* $\gcd(P,Q)=1$; $P,Q$ are irreducible; $R(X)$ is Ritt-indecomposable with a generic factor profile; the gcd of its root-exponents is $1$ (no $X\mapsto X^m$ structure for baby-step/giant-step). `probes/07_*`.
- *Constants.* The $256$ Grain-LFSR constants have no zeros, no repeats, none below $2^{16}$, and no alignment with the pinned coordinates $\{0,1\}$ — nothing-up-my-sleeve. `probes/08_*, 10_*`.

**Outside cryptanalysis.**

- *Statistics / learning.* On a $2$-variable slice, $P^+$ is statistically a uniform random function: outputs are uniform ($\chi^2=105.8$ vs $100\pm 20$ at $q=101$), and the largest linear character coefficient is $\max_{a\ne0}|\hat f(a)|=0.0301$, matching the random-function extreme-value prediction $\sqrt{2\ln q}/q=0.0301$. No $\varepsilon$-heavy Fourier coefficient $\Rightarrow$ not weakly learnable (Linial–Mansour–Nisan, Kushilevitz–Mansour). `probes/bm_learning_fourier.py`.
- *Belief / survey propagation.* BP (validated exact on a tree) converges on the Poseidon factor graph to the uniform fixed point — every belief $=1/p$ to machine $\varepsilon$, so decimation is just random guessing. Even the exact single-variable marginal over the solution set is near-uniform (entropy $2.806$ vs max $2.833$): only the joint distribution pins a solution, which BP cannot see. Cause: MDS $\Rightarrow$ a dense, girth-$4$ factor graph with no tree-likeness. `probes/bm_inference_bp.py`.
- *Time–memory / collision / birthday.* A single target makes Hellman precomputation $=\Theta(2^{62})=$ brute force; there is no group structure for rho/BSGS; and the surfaces $\{\mathrm{out}_0=0\},\{\mathrm{out}_1=0\}$ are each $14$-dimensional, so independent samples collide only with probability $\approx p^{-13}$.
- *Analytic number theory.* Lang–Weil's $\asymp\delta^2 q^{r-1/2}$ error (with $\delta\in[3^{16},3^{20}]$) swamps the main term; for a $0$-dimensional slice ($r=0$) it gives only $\#\le\delta$ with no usable lower bound — it counts but cannot locate. Stepanov needs $\deg\ll\sqrt p$, but $\deg\approx 3^{16}=2^{25.4}\gg\sqrt p=2^{15.5}$.
- *Algebraic-geometry point-finding.* Cafure–Matera (constructive $\mathbb{F}_q$-point) requires $q>8\,n^2 d\,\delta^4$, i.e. $q>2^{113}$ here, while $q=2^{31}$. Unirationality and determinant-method bounds fail in the same $\deg\gg q$ regime.
- *Lattices / $p$-adic.* Coppersmith/LLL need low-height roots, but the relaxed-verifier solutions are uniform (mean $2$-adic valuation $0.75$, no bias). Tropical/BKK: the post-elimination Newton polytope is full (support density $\approx 1$), so the mixed volume equals Bézout $=3^{32}$ — worse — and the homotopy path count is $\ge D_I$.
- *Optimization.* Sum-of-squares/Lasserre need an ordered field (absent in $\mathbb{F}_p$); even lifted, the relaxation order $\sim 3^{16}$ blows the moment matrix past $2^{55}$ before exactness.

## 5. Structural obstruction

Three properties hold simultaneously, each disabling a whole family of methods:

1. **Degree $\gg$ field.** $\deg=3^{16}=2^{25.4}$ and $D_I=3^{20}\approx p$. This kills constructive point-counting (Lang–Weil's error $\asymp\delta^2$, Stepanov needs $\deg\ll\sqrt p$, Cafure–Matera needs $q\gg\delta^4$) and collapses any interpolate-vs-scan tradeoff (since $D_I\approx p$).
2. **Statistically random.** The character spectrum is flat to the $\sqrt{\ln q}$ fluctuation, and the single-variable solution marginals are near-uniform. This kills learning, inference (BP/SP), and collision/TMTO methods — none see a per-coordinate signal.
3. **Generic.** $P,Q$ are irreducible and coprime, the displacement rank is that of a random pair, the Newton polytope is full after elimination, and the constants are NUMS. This kills sparsity-, symmetry-, group-, and structure-exploiting solvers.

Chevalley–Warning ($\sum\deg=2\cdot 3^{16}\gg 14$, so inapplicable) and Lang–Weil only confirm that the $\approx p^{12}$ solutions exist; neither constructs one. Every method collapses to the same invariant: $\approx 1$ $\mathbb{F}_p$-solution per slice, locatable only at elimination degree $D_I=3^{20}$, which nothing here reduces.

## 6. Status

No sub-$2^{55}$ method emerged anywhere in the algebraic, structural, and cross-disciplinary search — consistent with published work stopping at $R_P=8$. The only route to a solution is the expected one: the skip-$1$ resultant at $D_I=3^{20}$, about $2^{55}$ operations and $15$–$30$ GB of streamed memory — a rented-compute job of order \$3k–\$15k against a \$15k prize, which we have chosen not to run. The deliverables are this analysis and the verified small-scale solver `probes/11_resultant_solver.py`.

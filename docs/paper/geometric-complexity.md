---
title: "Geometric Complexity in Cyclic Regime Changes"
subtitle: "Full, Fundamental-Subspace, Shared-Orbit and Approximate-Orbit Models under Minimum Description Length"
author: Mac Mayo
---

# Abstract

A categorical regime change can be modelled at several levels of structural
constraint: each segment free in the full parameter space, each confined
independently to an invariant subspace, or both required to arise from one
continuous state differing by a finite group action. At a known boundary these
cost $\tfrac{g-1}{2}\log L$, $\tfrac{d_{\mathrm{fund}}}{2}\log L$ and — because
the shared-orbit alternative introduces no continuous parameter — only $O(1)$.

That $O(1)$ regime is the substance of the paper. Orbit collapse is a *benign
finite-branch singularity*: the relative label is nonidentifiable where the
branches meet, so the union is singular as a statistical model, but every fixed
branch has nonsingular Fisher information and the marginal likelihood is a finite
sum of regular Laplace integrals. The real log canonical threshold is therefore
unchanged and the zero increment survives; what collapse alters is the bounded
term and the limiting law. We derive that law — a maximum over the nonidentity
shifts of differences of Gaussian projection energies — prove it, and validate it
against the detector. Two consequences follow. It supplies asymptotic critical
values, making the statistic operational. And it shows that no
two-part label constant controls the raw rule: the unpenalised gain exceeds zero
with probability exactly $(g-1)/g$, and subtracting $\log(g-1)$ still leaves an
asymptotic firing rate of 0.14 to 0.50.

A 468,000-dataset study over $C_2$ through $C_6$ accompanies this. Calibrated at
a common 5%, the shared-orbit detector needs roughly 31–39% fewer observations
than the unrestricted one on exact-orbit data for $g = 4,5,6$, and loses that
entirely under higher-mode misspecification; since a per-configuration penalty
cancels exactly from a calibrated comparison, these figures measure the targeting
of the likelihood statistics rather than the effect of the codes. Three quantities
measured in nats are kept distinct throughout: the code constant, the testing
threshold, and the calibrated operating point.

**Keywords:** minimum description length; cyclic groups; changepoint models;
representation theory; information geometry; categorical data; model selection;
symmetry.

---

# 1. Introduction

A standard two-segment regime comparison asks whether one distribution explains
an entire sample or whether separate distributions are needed on either side of
a proposed boundary. In an unrestricted categorical model the alternative
assigns one probability vector to the left segment and another to the right. The
additional complexity is governed by the dimension of the entire simplex.

Many regime structures are not arbitrary. Physical phases, reading-frame states,
rotational sectors and cyclic operating modes may be related by symmetry. Once
symmetry is admitted, however, there are at least two fundamentally different
ways to use it.

First, symmetry may identify a low-dimensional invariant subspace in which each
regime is allowed to vary independently. This is an ordinary reduction of the
admissible parameter space. Second, the regimes may be required to be exact
group transforms of one shared state. This is a *relational* constraint: it
removes the independent continuous parameter vector introduced by an ordinary
split and replaces it with a discrete group label.

These possibilities lead to the models

- **A. Full independent.** Each segment receives an unrestricted parameter
  vector.
- **B. Independent fundamental-subspace.** Each segment receives its own
  parameter vector inside a selected invariant subspace.
- **C. Shared exact orbit.** The segments share one continuous parameter vector
  and differ only by a relative group action.
- **D. Approximate orbit.** The segments share one continuous vector up to a
  *shrunk* deviation, interpolating between B and C.

with continuous-dimension increments $g-1$, $d_{\mathrm{fund}}$, $0$ and
(asymptotically) $d_{\mathrm{fund}}$ respectively.

The distinction is not semantic. It determines the leading MDL penalty and the
sample length required for reliable detection. The contributions of this paper
are:

1. the known-boundary MDL penalties for the four model classes, and the limit
   law governing the shared-orbit statistic where that penalty is $O(1)$;
2. a Fisher-orthonormal Fourier parameterisation for direct cyclic categorical
   families, and its block generalisation in which the group acts on phases
   while each phase carries its own alphabet;
3. implementations of all detectors with consistent known-boundary penalties and
   detector-specific population gains;
4. a large Monte Carlo comparison testing both asymptotic score behaviour and
   calibrated practical power;
5. a model-selection procedure that identifies the geometry from data rather
   than assuming it;
6. a negative result: the two-part label constant does not control the raw
   decision rule near orbit collapse. A constant *can* fix a chosen asymptotic
   level — that is what Corollary 2 supplies — but no finite constant makes the
   false-alarm probability vanish there, and the structural constant is a nat or
   more short of even the level-$0.05$ threshold.

The analysis is explicitly an offline **known-boundary** model comparison. It is
not a changepoint-discovery algorithm. Unknown-boundary scanning and online
detection are discussed separately because they add a location or multiplicity
cost to every detector.

## 1.1 Relation to prior work

The ingredients are individually standard. Regular dimension penalties descend
from Schwarz (1978); MDL as a coding principle rather than a fixed penalty
formula is set out by Barron, Rissanen and Yu (1998) and Grünwald (2007).
Information-criterion changepoint selection dates to Yao (1988), with MDL
segmentation developed by Davis, Lee and Rodriguez-Yam (2006). The closest
categorical antecedents are Wang, Zou and Yin (2018) and Truong and Runge
(2024), both of which test unrestricted probability change rather than
invariant-subspace or orbit alternatives. Dimension-reduced changepoints appear
in Yu et al. (2026), where the subspace is learned from sparse structure rather
than fixed by a representation. Group-invariant testing via maximal invariants
and e-statistics is developed by Pérez-Ortiz, Lardy, de Heide and Grünwald
(2024). Singular-model asymptotics, directly relevant at orbit collapse, are due
to Watanabe (2009, 2013) and Drton and Plummer (2017).

The contribution here is their combination: a known-boundary MDL comparison that
separates *independent invariant-subspace change* from *shared-orbit
transition*, and shows the latter introduces no new continuous parameter across
the split on a regular orbit stratum. We are not aware of prior work making that
distinction explicitly, but this is a search-based inference about absence rather
than a proof of it, and the claim should be read accordingly.

The claim is deliberately narrow, and it is worth saying what is *not* claimed.
Shared parameters across a changepoint, transformation alternatives, and
group-invariant testing are all established. What we have not found elsewhere is
the explicit four-way MDL separation of unrestricted change, independent
representation-constrained change, shared finite-group-orbit change, and
shrinkage toward such an orbit, together with the observation that the third of
these carries no continuous-dimension increment at all while the second carries
$d_{\mathrm{fund}}$.

The relation to group-invariant e-testing deserves particular care, because it is
the nearest neighbour. Pérez-Ortiz et al. (2024) show that among all e-statistics
for testing between two group models, the likelihood ratio of the maximal
invariant is growth-rate optimal, and that an anytime-valid test can be built on
it. That is a different question from ours — theirs is a sequential testing
problem with no coding component, ours a fixed-boundary description-length
comparison across a hierarchy of alternatives — but the group being quotiented is
the same one. We do **not** claim Model C's statistic is the maximal-invariant
likelihood ratio: by Wijsman's representation [23] — see also Eaton [12] — that
ratio *averages* the densities over the group with respect to Haar measure,
whereas Model C aligns and profiles, then *maximises* over the nonidentity
elements. Those are different statistics, and establishing a
correspondence would take an argument this paper does not give. Their framework
remains the most promising route to a sequential version of this comparison, one
that would not depend on a known boundary or on two-part boundary coding;
Section 14.7 develops the point.

---

# 2. Statistical setup

## 2.1 Known-boundary two-segment comparison

Let $X_1,\dots,X_L$ be independent categorical observations with a proposed
boundary after $L_1$ observations, $L = L_1 + L_2$. The left and right segments
are denoted $S_1$ and $S_2$. **The boundary is supplied to the detector**; no
search over candidate locations is performed.

For model class $M$, write $H_0^M$ for the one-regime null and $H_1^M$ for the
two-regime alternative. The detector score is

$$ T_M \;=\; \big[\hat\ell_1^M - \hat\ell_0^M\big] \;-\; \mathrm{pen}_M , $$

where $\hat\ell$ are maximised log-likelihoods in nats. The raw MDL rule declares
a change when $T_M > 0$.

**Three objects, kept apart.** Much of what follows depends on not conflating
them, and they are easy to conflate because all three are measured in nats.

1. The **structural code constant** $\mathrm{pen}_M$ — a codelength, fixed by a
   coding convention, and the thing the complexity laws of Section 4 are about.
2. A **testing threshold** — a critical value, fixed by choosing an error rate
   $\alpha$, and derivable for Model C from the limit law of Section 6.4.
3. **Calibrated power** — a comparison at a common nominal false-positive rate,
   which Section 8.3 shows is invariant to (1) and therefore reports on the
   likelihood statistics rather than on the codes.

Section 9.1 shows that (1) does not serve as (2), which is the paper's negative
result; Section 9.2 reports (3) and says what it does and does not license.

## 2.2 Full parameter family

Let $\{P_\theta\}$ be a regular parametric family of continuous dimension
$d_{\mathrm{full}}$. For a direct $g$-category multinomial model,
$d_{\mathrm{full}} = g-1$. The notation is kept general because applications with
several categorical blocks have larger dimensions; Section 12 develops that case.

## 2.3 Notation

One symbol per quantity, throughout:

| Symbol | Meaning |
|---|---|
| $g$ | **group order** — the order of the cyclic group $C_g$, and the number of phase blocks in a block family |
| $a$ | **alphabet size** within one block. In the direct model of Sections 3–11 the group acts on the categories themselves, so $a = g$ |
| $L$, $L_1$, $L_2$ | total and per-segment sample lengths; $\rho = L_1/L$ |
| $d_{\mathrm{full}}$, $d_{\mathrm{fund}}$ | continuous dimensions of the full and fundamental families |
| $\eta$ | fundamental-family coordinate; $\|\eta\|_2$ is a Fisher norm |
| $r$ | relative group element (shift); $R_g$ its action on coordinates |
| $\delta$, $\tau$ | Model D's deviation and its prior scale |

Because the direct model sets $a = g$, statements about "raising $g$" in
Sections 3–11 change the alphabet too — which is exactly the confound Section 12
removes.

Readers moving between the paper and the reference implementation should note
that the code names the direct model's group order `m`, since there it is
simultaneously the alphabet size; the block module uses `g` and `a` as here.

## 2.4 Cyclic group action

Let $C_g = \{e, r, \dots, r^{g-1}\}$ act on the family through transformations
$T_r$. In the direct categorical model the alphabet size equals the group order
and $T_r$ cyclically permutes the category coordinates. In block models the group
permutes phase blocks while leaving the within-block alphabet unchanged.

At a symmetric reference distribution $P_0$ the tangent space carries an induced
real representation of $C_g$. Let $V_{\mathrm{fund}}$ be the selected
fundamental invariant component, of real dimension $d_{\mathrm{fund}}$. For a
direct cyclic categorical model,

$$ d_{\mathrm{fund}} = \begin{cases} 1 & g = 2 \ \text{(the sign representation)} \\ 2 & g \geq 3 .\end{cases} $$

---

# 3. Four model classes

## 3.1 Model A: full independent change

Under the null both segments share one unrestricted parameter $\theta$. Under
the alternative $\theta_1$ and $\theta_2$ are independently fitted, so

$$ \Delta d_A = d_{\mathrm{full}} = g-1 . $$

## 3.2 Model B: independent fundamental-subspace change

Let $\eta \in \mathbb{R}^{d_{\mathrm{fund}}}$ parameterise a smooth
exponential-family chart inside the selected invariant subspace,

$$ p(\eta) \;=\; \mathrm{softmax}\!\big(B\eta\big), $$

where the columns of $B$ span the fundamental component in logit coordinates.
Under the alternative $\eta_1$ and $\eta_2$ are independent, so

$$ \Delta d_B = d_{\mathrm{fund}} . $$

Model B states that both regimes lie in the same invariant subspace. It does
**not** require the right regime to be a group transform of the left.

## 3.3 Model C: shared exact-orbit transition

Under the alternative both segments share one continuous vector $\eta$, and the
right segment is $T_r$ applied to it:

$$ p_1 = p(\eta), \qquad p_2 = T_r\, p(\eta), \qquad r \in \{1,\dots,g-1\}. $$

A simultaneous shift of both labels is observationally redundant, so the left
label is fixed as a gauge.

> **Proposition 1 (zero continuous increment on a regular stratum).** Assume
>
> 1. the group order $g$ is fixed and known, and the segment boundary is known;
> 2. the shared state $\eta$ lies on a **regular orbit** — its stabiliser in
>    $C_g$ is trivial, so the $g$ points $\{R^s\eta\}$ are distinct — and in
>    particular $\eta \neq 0$;
> 3. the chart of Section 5.3 is used, so the model is a regular exponential
>    family in a neighbourhood of $\eta$ with non-singular Fisher information;
> 4. the relative shift $r$ is encoded by the uniform two-part code over the
>    $g-1$ nonidentity elements.
>
> Then Model C introduces no continuous parameter beyond the null, so
> $$ \Delta d_C = 0, $$
> and its total incremental cost over the null is $\log(g-1)$ nats, constant in
> $L$.
>
> Assumption 4 fixes the constant but is not forced by the rest. Model C's
> alternative excludes the identity, so it is not a proper extension of its own
> null, and the cost depends on whether the code is taken conditional on a
> change, includes the identity ($\log g$), carries a model-selection cost, or
> mixes over $r$ rather than maximising. **The robust content of the proposition
> is $\mathrm{pen}_C = O(1)$ with no continuous $\log L$ increment**; the exact
> constant is a convention, and Section 9.1 measures what turns on the choice.

Assumption 2 is the one that bites. At $\eta = 0$ the orbit collapses to a point,
every $R^s$ acts trivially, and the parameter is a fixed point of the whole group;
more generally at any $\eta$ with a nontrivial stabiliser the shift is not
identifiable, so the model is singular there in the sense of Watanabe [7,8] and
Drton–Plummer [9]. Section 6.3 shows the singularity is benign: the alternative is
a finite union of branches each with nonsingular Fisher information, the real log
canonical threshold is unchanged, and the zero increment survives. What collapse
alters is the bounded term, and Section 6.4 derives the nondegenerate $O_p(1)$
law that governs it. Section 9.1 shows what that does to the raw rule in practice. At $g=2$
there is a single nonidentity shift and the label cost is $\log 1 = 0$.

## 3.4 Model D: approximate orbit

Exact symmetry is a strong assumption. Model D relaxes it:

$$ \eta_2 \;=\; R_g^{\,r}\,\eta_1 + \delta, \qquad \delta \sim \mathcal{N}(0,\tau^2 I), $$

with a shrinkage code on $\delta$. Setting $\tau = 0$ pins $\delta$ out and
recovers Model C *exactly*; letting $\tau \to \infty$ leaves $\delta$ free and
recovers Model B's maximised gain.

The cost of $\delta$ follows from the Laplace approximation, but the shared state
$\eta_1$ must be handled with care: it is *estimated jointly with* $\delta$, not
known, and the two are correlated, because shifting $\eta_1$ and compensating
with $\delta$ leaves the right segment's fit unchanged. In Fisher-orthonormal
coordinates the joint information in $(\eta_1, \delta)$ is, per direction,

$$ H \;=\; \begin{pmatrix} L_1 + L_2 & L_2 \\ L_2 & L_2 \end{pmatrix}, $$

so what constrains $\delta$ after $\eta_1$ is profiled away is not $L_2$ but the
Schur complement

$$ J_{\mathrm{eff}} \;=\; L_2 - \frac{L_2^{\,2}}{L_1 + L_2} \;=\; \frac{L_1 L_2}{L_1 + L_2}, $$

giving

$$ \mathrm{pen}_\delta \;=\; \frac{d_{\mathrm{fund}}}{2}\,\log\!\Big(1 + \tau^2\,\frac{L_1 L_2}{L_1 + L_2}\Big) \ \text{nats}. $$

Using $L_2$ in place of $J_{\mathrm{eff}}$ would be correct only if the shared
state were known, or if the left segment were infinitely informative; on a
balanced split it uses $L/2$ where joint estimation gives $L/4$, inflating the
penalty by up to $\tfrac{d_{\mathrm{fund}}}{2}\log 2$. That error is bounded, so
it would leave the leading coefficient alone — but Model D's entire contribution
*is* the bounded term, so it is exactly the term that must be right.

More generally, away from the Fisher reference point the identity-information
approximation is no longer exact, and the correct expression is

$$ \tfrac12 \log\det\!\big(I + \tau^2 J_{\mathrm{eff}}\big), $$

with $J_{\mathrm{eff}}$ the Schur complement of the observed information. Fisher
orthonormality holds *at* $\eta = 0$; it does not make the metric globally the
identity. The isotropic formula above is the reduction of this at the reference
point, and it is what the implementation uses.

The discrepancy is $O(\|\eta\|^2)$, for every group order, and the argument is
short enough to give. Because $TB = BR$, the chart's information is equivariant,
$I(R_\theta\eta) = R_\theta^{\top} I(\eta) R_\theta$. Writing
$I(\eta) = I + A(\eta) + O(\|\eta\|^2)$ with $A$ linear, equivariance forces
$A(R_\theta\eta) = R_\theta^{\top} A(\eta) R_\theta$. Split $A$ into trace and
traceless parts. The trace part is a rotation-invariant *linear* functional on
$\mathbb{R}^2$, hence identically zero for every $g \ge 2$. The traceless part
transforms with weight $2$ while $\eta$ carries weight $1$, so it survives only
when $3\theta \equiv 0$, that is only at $g = 3$.

So the *metric* error is genuinely first order at $g=3$ and second order
elsewhere. The *penalty* error is second order at every $g$, because
$\tfrac12\log\det$ is blind to a traceless first-order perturbation: the Schur
complement's first-order change satisfies
$\operatorname{tr}\delta S = \kappa_1 \operatorname{tr}E_2 + \kappa_2 \operatorname{tr}E_1$
(using $\operatorname{tr}(RER^{\top}) = \operatorname{tr}E$), and
$\operatorname{tr}A \equiv 0$ by the trace argument above. Measured exponents
over a decade in $\|\eta\|$, at $L_1 = L_2 = 1{,}200$ and $\tau = 0.15$: the
metric gap scales as $2.00$ at $g = 2,4,5,6$ and as $1.00$ at $g=3$, while the
penalty gap scales as $1.99$–$2.03$ at every $g$. Section 13 records the
approximation as a limitation regardless, since second order is not zero and
Model D's margins are of order a quarter of a nat.

**A caveat belongs with this.** For any *fixed* $\tau > 0$ the leading
coefficient is $d_{\mathrm{fund}}/2$ — Model B's rate, not an intermediate one.
The interpolation lives entirely in the bounded term, that is, at finite $L$.
That suits the question being asked, since tolerance to imperfect symmetry is a
finite-sample matter, but a genuine interpolation of the *leading* coefficient
would require $\tau$ shrinking with $L$.

**Table 1.** Four levels of geometric constraint at a known boundary.

| Model | Alternative relation between segments | $\Delta d$ | Leading incremental penalty |
|---|---|---|---|
| A. Full independent | arbitrary separate parameters | $g-1$ | $\tfrac{g-1}{2}\log L$ |
| B. Independent fundamental | separate parameters in $V_{\mathrm{fund}}$ | $d_{\mathrm{fund}}$ | $\tfrac{d_{\mathrm{fund}}}{2}\log L$ |
| C. Shared exact orbit | one shared state plus relative $r$ | $0$ | $\log(g-1)$, constant in $L$ |
| D. Approximate orbit | one shared state, shrunk deviation | $d_{\mathrm{fund}}$ (fixed $\tau$) | $\log(g-1) + \tfrac{d_{\mathrm{fund}}}{2}\log\big(1+\tau^2 L_1L_2/(L_1{+}L_2)\big)$ |

---

# 4. MDL complexity laws

## 4.1 Regular codelength expansion

For a regular $d$-dimensional family fitted to $L$ observations, BIC, the regular
Laplace marginal likelihood and standard parametric-complexity expansions share
the leading form $\tfrac{d}{2}\log L + O(1)$, where the $O(1)$ term depends on
the coding convention, prior, Fisher information and parameter-space geometry.

## 4.2 BIC/Laplace split increment

For a known split, the incremental complexity produced by replacing one
$d$-dimensional fit with two is, under the BIC/Laplace expansion of Section 4.1,

$$ \mathrm{pen}_{\mathrm{split}}(d; L_1, L_2) \;=\; \frac{d}{2}\Big[\log L_1 + \log L_2 - \log (L_1+L_2)\Big] . $$

With $L_1 = \rho L$ and $L_2 = (1-\rho)L$ this rearranges exactly to

$$ \mathrm{pen}_{\mathrm{split}} \;=\; \frac{d}{2}\log L \;+\; \frac{d}{2}\log\big[\rho(1-\rho)\big] . $$

Thus the coefficient of $\log L$ is $d/2$, while the split fraction affects only
the bounded term. The expression is *algebraically* exact within the BIC/Laplace
approximation; it is not an exact universal codelength, and the name is chosen to
say so. Appendix B.4 examines the prediction empirically.

## 4.3 Unknown boundaries

If the boundary is unknown, a detector must encode or search over approximately
$L$ candidate locations. A simple two-part location code adds $\log L$ to **every**
model, so the leading coefficients become $\tfrac{d}{2}+1$ in each case. A
location cost cannot be applied to only one detector in a known-boundary
comparison; doing so confounds changepoint multiplicity with model dimension.

## 4.4 Units, and the constant $K^\ast$

In bits, a regular penalty is $\tfrac{\Delta d}{2}\log_2 L = \tfrac{\Delta d}{2\ln 2}\ln L$.
Writing

$$ K^\ast \;=\; \frac{1}{2\ln 2} \;=\; 0.7213475\ldots $$

for the per-dimension rate in bits per e-fold, the leading coefficients are
$(g-1)K^\ast$ for Model A, $d_{\mathrm{fund}}K^\ast$ for Model B and zero for
Model C.

This is a **units conversion and nothing more**. $K^\ast$ is Schwarz's one-half
expressed in base 2, so any quantity that counts half a parameter per e-fold and
reports in bits produces it, and the observation that the coefficients are whole
multiples of it says only that regular BIC counts a whole number of parameters.
The paper draws no conclusion from it. Appendix C records a numerical
coincidence involving the same constant and explains why that coincidence should
not be leaned on either.

---

# 5. Cyclic Fourier geometry

## 5.1 Fisher metric at the uniform distribution

For the direct $g$-category model let $u = (1/g,\dots,1/g)$. The tangent space is
$\{v \in \mathbb{R}^g : \sum_j v_j = 0\}$, and at $u$ the Fisher inner product is
$\langle v, w\rangle_F = g\sum_j v_j w_j$.

## 5.2 Fundamental Fourier basis

For $g \geq 3$ define $\phi_j = 2\pi j/g$ and take logit directions
$\sqrt{2}\cos\phi_j$ and $\sqrt{2}\sin\phi_j$; for $g = 2$ take the sign
direction. Their span is invariant under cyclic permutation, and in coefficient
space a one-step cyclic shift acts as a planar rotation by $2\pi/g$.

The implementation uses **Cartesian** coordinates rather than amplitude–angle
coordinates, which avoids the unidentifiable angular coordinate at zero
amplitude.

## 5.3 Positivity-preserving family

The direct fundamental family is defined through softmax logits,
$p(\eta) = \mathrm{softmax}(B\eta)$, with $B$ scaled so that $\partial p/\partial\eta$
at $\eta = 0$ is the Fisher-orthonormal tangent basis — making $\|\eta\|_2$ the
Fisher norm of the local perturbation. If $R_g$ is the fundamental rotation
matrix, the family satisfies the equivariance identity

$$ p\big(R_g^{\,r}\eta\big) \;=\; T_r\, p(\eta). $$

Automated tests verify Fisher orthonormality and numerical equivariance to
$10^{-13}$ through $g = 8$.

## 5.4 Local Jensen–Shannon geometry

For small perturbations the weighted Jensen–Shannon divergence between $p(\eta)$
and $p(R_g\eta)$ satisfies

$$ \mathrm{JSD} \;\to\; \tfrac{1}{8}\,\big\|R_g\eta - \eta\big\|^2 \;=\; \frac{1-\cos(2\pi/g)}{4}\,\|\eta\|^2 , $$

giving coefficients $1/2$, $3/8$ and $1/4$ for $g = 2, 3, 4$. No extra factor of
$g$ appears when $\eta$ is measured in Fisher norm.

Verified numerically, and the tolerance is not uniform across group orders: at
$\|\eta\| = 10^{-2}$ the relative error is $3.5\times10^{-3}$ at $g = 3$ and
$2.5\times10^{-5}$–$5.0\times10^{-5}$ at every other group order through $g = 8$
— two orders of magnitude smaller. The gap is the same representation-theoretic
fact as in Section 3.4: the chart metric carries a nonzero *linear* term only at
$g = 3$, so the cubic correction to the divergence is an order of magnitude
larger there. Quoting a single tolerance for all $g$ would report the $g = 3$
figure as though it were a general numerical limit.

---

# 6. Population gains and detection boundaries

## 6.1 Detector-specific gain

The relevant signal strength is the expected log-likelihood advantage *within
the model being fitted*. For the full multinomial model with equal segment sizes
the pooled null is the arithmetic mixture and the gain is the weighted
Jensen–Shannon divergence. For a restricted family the pooled null is generally
a KL projection rather than the arithmetic mixture, so **ordinary JSD must not be
reused for Models B, C and D**. The implementation computes each population gain
by optimising the corresponding population log-likelihood.

## 6.2 Local detection laws

If $G_M \approx c_M \epsilon^2$, the regular models have boundaries
$\epsilon^2 \sim \tfrac{d}{2}\,\tfrac{\log L}{L}$. For Model C, *any* fixed code
constant $c$ gives a deterministic crossing at $\epsilon^2 \sim c/L$, so the
scaling is $O(L^{-1})$ rather than $O(L^{-1}\log L)$ — the improvement comes from
*sharing the continuous state across the boundary*, not merely from a smaller
tangent space.

That is a statement about where the mean gain crosses a constant, and it is not a
detection boundary. The raw rule with constant $c$ has an asymptotic false-alarm
rate of $\Pr(W_{g,\rho}(h) > c)$, which by Section 6.4 does not vanish; at $g=2$
the natural constant is $\log 1 = 0$ and the crossing statement degenerates
entirely. Level-$\alpha$ local power for Model C is governed by the quantiles of
$W_{g,\rho}(h)$, and is properly stated there.

## 6.3 Orbit collapse is a benign finite-branch singularity

At $\eta = 0$ all orbit elements coincide and the relative label is
unidentifiable, so the map $(\eta, r) \mapsto P_{\eta,r}$ is not one-to-one
there: the $g-1$ points $(0, r)$ index the same distribution. Under the standard
definition — a one-to-one parameterisation with positive-definite Fisher
information [8] — **Model C is therefore singular at collapse**, and calling it
regular would be wrong.

The singularity is nonetheless benign, in a way that matters for everything that
follows. For each *fixed* nonidentity $r$, the map
$\eta \mapsto \big(p(\eta),\, p(R^r\eta)\big)$ is a smooth
$d_{\mathrm{fund}}$-dimensional family with positive-definite quadratic
Kullback–Leibler expansion at the origin. Nothing about it degenerates. Model C's
parameter space is the *union* of the $g-1$ such branches, all passing through
the same point $(u,u)$, and the null is a further $d_{\mathrm{fund}}$-dimensional
manifold through it. Every branch has nonsingular Fisher information; what fails
is only identifiability *between* branches, and at a single point. The
Kullback–Leibler zero set at collapse is the finite set $\{(0,r)\}$ of isolated,
non-degenerate zeros, so the marginal likelihood is a finite sum of ordinary
Laplace integrals, each of order $L^{-d_{\mathrm{fund}}/2}$, and the branch
multiplicity enters the constant rather than the exponent. The real log canonical
threshold is thus $\lambda = d_{\mathrm{fund}}/2$ on both sides, giving

$$ \Delta\lambda_C \;=\; 0 $$

**at collapse as well as on a regular orbit.** So the singularity does not
disturb the $\log L$ coefficient: Proposition 1's zero increment survives the
point that appeared to threaten it. What collapse changes is the bounded term and
the limiting law of the statistic — and that is not a technicality, because for
Model C the bounded term is the entire penalty. Section 6.4 derives the law, and
it is where the paper's negative result about two-part MDL thresholding comes
from.

The only collapsed state in this family is the origin. A nonidentity element of
$C_g$ acts on $V_{\mathrm{fund}}$ as a rotation by $2\pi r/g \neq 0$, so
$R^r - I$ is invertible and fixes nothing but $0$; numerically
$\min_{r\neq e}|\det(R^r - I)| \geq 0.586$ through $g = 8$. Statements about
"states with nontrivial stabiliser" are therefore vacuous here, though they would
not be for a nonfaithful higher mode or a richer group action.

## 6.4 The collapse law

Write $\rho = L_1/L$ and let the shared state sit at the local value
$\eta_L = h/\sqrt{L}$, so that $h = 0$ is exact collapse and large $\|h\|$ is a
regular orbit state. Let $U_1, U_2 \sim N(0, I_{d_{\mathrm{fund}}})$ be the two
segments' normalised scores, and $w_r = R^r h - h$.

Assume throughout that $L_1, L_2 \to \infty$ with $L_1/L \to \rho \in (0,1)$; that
$g$ and $h$ are fixed; and that the chart of Section 5.3 is used, so the family is
a smooth exponential family with Fisher information continuous at the origin and
$I(0) = I_{d_{\mathrm{fund}}}$.

> **Theorem 1 (local law for the shared-orbit statistic).** Under the null of no
> change with shared state $\eta_L = h/\sqrt{L}$, the shared-orbit raw gain
> converges in distribution to
>
> $$ W_{g,\rho}(h) \;=\; \max_{r \neq e}\Big[ \tfrac12\big\|\sqrt{\rho}\,U_1 + \sqrt{1-\rho}\,R^{-r}U_2 - (1-\rho)R^{-r}w_r\big\|^2 \;-\; \tfrac12\big\|\sqrt{\rho}\,U_1 + \sqrt{1-\rho}\,U_2\big\|^2 \;+\; \sqrt{1-\rho}\,U_2^{\top}w_r \;-\; \tfrac{1-\rho}{2}\|w_r\|^2 \Big]. $$

*Proof.* Write $\rho_L = L_1/L \to \rho$, $\rho_{1,L} = \rho_L$,
$\rho_{2,L} = 1-\rho_L$, and $w_r = R^r h - h$. Let $S_{i,L}$ be segment $i$'s
score at $\eta_L$ and $U_{i,L} = S_{i,L}/\sqrt{L_i}$ its normalisation; the two
segments are independent, $\eta_L \to 0$, and $I$ is continuous with
$I(0) = I_{d_{\mathrm{fund}}}$, so the joint central limit theorem gives
$$ (U_{1,L}, U_{2,L}) \;\Longrightarrow\; (U_1, U_2) \sim N\big(0, I_{2d_{\mathrm{fund}}}\big), \qquad I(\eta_L) \to I_{d_{\mathrm{fund}}}. $$
Reparameterise the shared state as $\eta_1 = (h+t)/\sqrt{L}$, so that under shift
$r$ the two segments sit at local displacements $t/\sqrt{L}$ and
$(w_r + R^r t)/\sqrt{L}$ from the truth.

**Localisation.** Each branch objective $\ell_1(\eta) + \ell_2(R^r\eta)$ is
strictly concave in $\eta$: the chart of Section 5.3 makes each term a regular
exponential-family log-likelihood in the logits, $B$ and $R^r$ are linear, and
the information is positive definite near the origin. Its population maximiser
is the truth and is consistent, so the exact profiled maximiser satisfies
$\sqrt{L}\,\hat\eta_{r,L} = O_P(1)$ for each $r$, and uniformly over the fixed
finite set of shifts. With probability tending to one the profiling may therefore
be restricted to a compact neighbourhood of the origin in $t$, on which the LAN
remainder below is uniform — which is what licenses substituting the quadratic
expansion into the profiling steps rather than merely maximising the expansion.

**Expansion.** By local asymptotic normality, for each segment and uniformly on
compact sets of $\Delta_i$,
$$ \ell_i\big(\eta_L + \Delta_i/\sqrt{L}\big) - \ell_i(\eta_L) \;=\; \sqrt{\rho_{i,L}}\,U_{i,L}^{\top}\Delta_i \;-\; \tfrac{\rho_{i,L}}{2}\,\Delta_i^{\top}I(\eta_L)\,\Delta_i \;+\; o_P(1). $$

**Null branch.** At $r = e$ both displacements equal $t$; the quadratic
coefficients combine as $\rho_L/2 + (1-\rho_L)/2 = 1/2$, and the objective is
$$ \big(\sqrt{\rho_L}\,U_{1,L} + \sqrt{1-\rho_L}\,U_{2,L}\big)^{\top}t \;-\; \tfrac12\|t\|^2 \;+\; o_P(1), $$
whose maximum over the localised region is
$\tfrac12\big\|\sqrt{\rho_L}U_{1,L} + \sqrt{1-\rho_L}U_{2,L}\big\|^2 + o_P(1)$.

**Shifted branches.** For fixed $r \neq e$, substituting $\Delta_2 = w_r + R^r t$
and using $w_r^{\top}R^r t = (R^{-r}w_r)^{\top}t$ gives
$$ \big(\sqrt{\rho_L}\,U_{1,L} + \sqrt{1-\rho_L}\,R^{-r}U_{2,L} - (1-\rho_L)R^{-r}w_r\big)^{\top} t \;-\; \tfrac12\|t\|^2 \;+\; \sqrt{1-\rho_L}\,U_{2,L}^{\top}w_r \;-\; \tfrac{1-\rho_L}{2}\|w_r\|^2 \;+\; o_P(1), $$
with maximum
$\tfrac12\big\|\sqrt{\rho_L}U_{1,L} + \sqrt{1-\rho_L}R^{-r}U_{2,L} - (1-\rho_L)R^{-r}w_r\big\|^2
+ \sqrt{1-\rho_L}\,U_{2,L}^{\top}w_r - \tfrac{1-\rho_L}{2}\|w_r\|^2 + o_P(1)$.

**Conclusion.** The detector's raw gain is the maximum over the $g-1$ nonidentity
shifts of the second display minus the first. That maximum is a continuous
function of $(U_{1,L}, U_{2,L})$ and of $\rho_L$, the index set is finite, and
$\rho_L \to \rho$, so the continuous mapping theorem applied to
$(U_{1,L}, U_{2,L}) \Rightarrow (U_1, U_2)$ delivers the stated limit.
$\square$

> **Corollary 1 (exact collapse).** At $h = 0$ every $w_r$ vanishes and the limit
> is a difference of Gaussian projection energies,
> $$ W_{g,\rho}(0) \;=\; \tfrac12\max_{r\neq e}\Big[\big\|\sqrt{\rho}U_1 + \sqrt{1-\rho}R^{-r}U_2\big\|^2 - \big\|\sqrt{\rho}U_1 + \sqrt{1-\rho}U_2\big\|^2\Big], $$
> the maximum over the nonidentity shifts of the gain from aligning rather than
> pooling. It does not depend on $L$: the statistic is $O_p(1)$ at collapse, and
> the zero-threshold false-positive rate therefore does **not** vanish as data
> accumulates.

**A readable form on a balanced split.** Take $L_1 = L_2 = n$ and the
*per-segment* convention $\eta_n = h_n/\sqrt{n}$. With $Z_1, Z_2$ iid
$N(0, I_{d_{\mathrm{fund}}})$ put $Y_i = h_n + Z_i$ — signal plus noise in each
segment. Then

$$ W_g(h_n) \;=\; \tfrac14 \max_{r \neq e}\Big[\;\big\|Y_1 + R^{-r}Y_2\big\|^2 \;-\; \big\|Y_1 + Y_2\big\|^2\;\Big]: $$

aligned against unaligned, and nothing else. This is Theorem 1 at $\rho = 1/2$
and $h = \sqrt2\,h_n$; the two expressions differ by
$h^{\top}R^{-r}h - h^{\top}R^{r}h$, which vanishes for a rotation, and agree to
machine precision under common random numbers. The convention matters and is
easy to get wrong: the naive extension of the $Y$ form to unbalanced splits,
weighting $Y_i$ by $\sqrt{\rho_i}$, reproduces Corollary 1 exactly and is wrong
everywhere off collapse.

> **Corollary 2 (critical values).** For each $g$, $\rho$ and $\alpha$ the
> $1-\alpha$ quantile of $W_{g,\rho}(0)$ is a pointwise asymptotic level-$\alpha$
> critical value for the raw shared-orbit gain at collapse.

One feature of the collapse law is available in closed form, and it is the one
that does the most work later.

> **Corollary 3 (zero-threshold rate at collapse).** For every $g \geq 2$ and
> every $\rho \in (0,1)$,
> $$ \Pr\big(W_{g,\rho}(0) > 0\big) \;=\; \frac{g-1}{g}. $$

*Proof.* For $r = 0,\dots,g-1$ put
$A_r = \|\sqrt{\rho}U_1 + \sqrt{1-\rho}R^{-r}U_2\|^2$, so that
$W_{g,\rho}(0) = \tfrac12\big(\max_{r\neq e}A_r - A_0\big)$. Expanding,
$$ A_r \;=\; \rho\|U_1\|^2 \;+\; (1-\rho)\|U_2\|^2 \;+\; 2\sqrt{\rho(1-\rho)}\;U_1^{\top}R^{-r}U_2, $$
in which only the cross term depends on $r$. Hence
$\arg\max_r A_r = \arg\max_r U_1^{\top}R^{-r}U_2$ and the event
$\{W_{g,\rho}(0) > 0\} = \{\arg\max_r A_r \neq e\}$ does not depend on $\rho$ at
all. The map $U_2 \mapsto R^{-1}U_2$ preserves the standard Gaussian law and
independence from $U_1$, and carries $A_r$ to $A_{r+1 \bmod g}$, so
$(A_0,\dots,A_{g-1})$ is cyclically exchangeable; ties have probability zero.
Each index is therefore equally likely to be the unique maximiser, giving
$\Pr(\arg\max = e) = 1/g$. $\square$

Equivalently: at collapse the raw gain is positive exactly when some nonidentity
shift aligns the two segments' scores better than the identity does, and by
symmetry that happens $g-1$ times out of $g$. The rate is $1/2$ at $g=2$ and
climbs toward $1$; it is a property of the geometry, not of the sample size.

**Table 2.** The collapse law $W_{g,1/2}(0)$: quantiles from 400,000 draws at the
base seed (`regimeshift/collapse.py`), and the exact rate from Corollary 3.

| $g$ | mean | $q_{0.95}$ | $q_{0.99}$ | $\Pr(W > 0) = \tfrac{g-1}{g}$ |
|---:|---:|---:|---:|---:|
| 2 | $-0.002$ | 1.596 | 2.980 | 0.5000 |
| 3 | 0.434 | 2.463 | 3.889 | 0.6667 |
| 4 | 0.604 | 2.520 | 4.049 | 0.7500 |
| 5 | 0.679 | 2.584 | 4.155 | 0.8000 |
| 6 | 0.716 | 2.610 | 4.189 | 0.8333 |

The simulated rates agree with the exact column to within $0.0011$ at 400,000
draws, which is a check on the implementation rather than on the corollary.

These are **critical values, not codelength constants**. They are defined by a
chosen error rate and move with $\alpha$ — $q_{0.99}$ is more than a nat above
$q_{0.95}$ at every $g$ — so they answer "what threshold makes the raw rule an
asymptotically valid test near collapse", not "what should the relative shift
cost to encode". Section 10 keeps the two apart deliberately.

**Is collapse the least favourable local null?** Corollary 2 is pointwise; using
$q_{1-\alpha}(W_{g,\rho}(0))$ as an operating threshold would need
$\sup_h q_{1-\alpha}(W_{g,\rho}(h)) = q_{1-\alpha}(W_{g,\rho}(0))$. Sweeping
$\|h\| \in [0,6]$ and a full orbit sector at $\alpha = 0.05$, $\rho = 1/2$, with
common random numbers across every $h$:

| $g$ | $q_{0.95}$ at collapse | sweep maximum | excess | $\|h\|$ at maximum | independent re-check |
|---:|---:|---:|---:|---:|---:|
| 2 | 1.603 | 1.605 | 0.002 | 0.05 | $+0.021$ |
| 3 | 2.461 | 2.464 | 0.003 | 0.15 | $+0.000$ |
| 4 | 2.517 | 2.522 | 0.006 | 0.15 | $+0.001$ |
| 5 | 2.585 | 2.586 | 0.001 | 0.10 | $+0.008$ |
| 6 | 2.603 | 2.603 | 0.000 | 0.15 | $+0.009$ |

No point in the sweep exceeds the collapse quantile by more than $0.01$ — the
largest excess anywhere is $0.006$ — the argmax sits at or beside the origin, and
re-evaluating the selected argmax against the origin on *independent* draws gives
differences of at most $0.021$ in either direction, consistent with noise rather
than with a real maximum away from collapse.
Common random numbers are essential to seeing this: estimating a few hundred
quantiles independently and taking the maximum manufactures an excess of about
$0.02$ from Monte Carlo noise alone, the same size as the effect being looked
for, and they still leave a residual selection bias that the independent
re-check is there to expose.

This is a **numerical finding over the evaluated local-null grid**, not uniform
size control. Whether collapse is least favourable for every $h$, and for every
$\rho$ and $\alpha$, is left as a conjecture; nothing in this paper turns on it
beyond the scope of the grid tested.

Theorem 1 is validated against the detector in `tests/test_collapse.py`, at
$h = 0$ for $g = 2,\dots,6$, at $\|h\| = 1.5$ and $3.0$, and at
$\rho = 0.25, 0.1$.

---

# 7. Detector implementation

## 7.1 Full detector

A BIC-scored unrestricted multinomial detector, so that all detectors are
comparable through maximised likelihood plus explicit complexity increments. The
penalty is the BIC/Laplace known-split increment with $d = g-1$. An exact
KT/Dirichlet mixture implementation is *not* used; see Section 14.3.

## 7.2 Fundamental detector

The null fits one coordinate $\eta$ to the combined counts; the alternative fits
$\eta_1$ and $\eta_2$ separately. Optimisation is L-BFGS-B with analytic
gradients in Cartesian Fourier coordinates and a smooth softmax map. The penalty
is the BIC/Laplace known-split increment with $d = d_{\mathrm{fund}}$.

## 7.3 Shared-orbit detector

The null fits one $\eta$ to the combined counts. For each nonidentity shift $r$
the right counts are aligned by $T_{-r}$, pooled with the left counts, and fitted
with one shared $\eta$. The alternative takes the shift with the largest
shared-state likelihood. Its penalty is $\log(g-1)$, with no location cost and no
continuous-dimension increment.

## 7.4 Approximate-orbit detector

For each nonidentity shift, $(\eta,\delta)$ are fitted jointly by maximising the
penalised likelihood $\ell_1(\eta) + \ell_2(R^r\eta + \delta) - \|\delta\|^2/2\tau^2$
with analytic gradients. The penalty is $\log(g-1) + \mathrm{pen}_\delta$, with
$\mathrm{pen}_\delta$ the profiled-information expression of Section 3.4. That
the joint fit estimates $\eta$ rather than conditioning on it is precisely why
the penalty must use $L_1L_2/(L_1+L_2)$ and not $L_2$: the code and the penalty
have to describe the same procedure.

`tests/test_approximate_orbit.py` grounds the closed form against brute-force
numerical marginalisation of the joint likelihood under the Gaussian code. At
$g=2$, $\tau = 0.15$ the two agree to within 0.003 nats at the Fisher reference
point, where the conditional $L_2$ form is off by 0.27.

## 7.5 Numerical notes

Convergence is judged on the first-order condition rather than on the optimiser's
own success flag. Under a deliberately tight `ftol`, L-BFGS-B reports abnormal
termination whenever its line search cannot improve at machine precision, which
happens *at* the optimum: observed gradient norms are around $5\times10^{-8}$,
giving roughly 2.5% false alarms if the flag is trusted. Judging convergence by a
relative gradient criterion gives zero failures across the production grid.

A second numerical property concerns existence of the fundamental MLE, and it is
sharper than a condition on empty cells. The fundamental family is a *linear*
exponential family — the logits are confined to the column span of $B$ — so its
sufficient statistic is $\bar t = \sum_j f_j B_j$ and, by the standard condition,
the MLE exists exactly when $\bar t$ lies in the **interior** of
$\mathrm{conv}\{B_1,\dots,B_g\}$. That hull is the regular $g$-gon on the $B_j$,
and every proper face of a regular polygon is spanned by cyclically *adjacent*
vertices. The criterion is therefore decidable from the support alone, in $O(g)$:

> **Proposition 2 (existence).** For $g = 2$ the MLE exists iff both categories
> have positive count. For $g \ge 3$ it exists iff the empirical support is
> contained in neither a single category nor a pair of cyclically adjacent
> categories.

A zero cell is thus neither necessary nor sufficient once $g \ge 3$. At $g = 4$,
$(0,10,10,10)$ has an ordinary interior optimum at $\|\hat\eta\| = 0.490$, and so
does $(10,0,10,0)$ with *two* empty cells, at $\hat\eta = 0$, because opposite
vertices span no face. What genuinely fails is $(10,10,0,0)$: an adjacent pair
carrying all the mass puts $\bar t$ on an edge, and the likelihood is
asymptotically flat along the escape direction — numerically
$\hat\ell(4\hat\eta) = \hat\ell(\hat\eta)$ to machine precision, which is why
different optimiser starts halt at very different coordinates while agreeing on
the log-likelihood. Only at $g = 2$, where the hull is a segment and both
endpoints are faces, does a zero cell coincide with nonexistence.

Nor are zero counts common on this design. The shortest segment on the
production grid is 100 observations and the smallest cell probability anywhere on
the grid is $0.091$, so the worst configuration's probability of an empty cell is
$7.2\times10^{-5}$; summed over the design with each configuration's own lengths
and draw counts, the expected number of zero-cell segments across all 936,000 the
run scores is $0.2$. And a zero cell would not be the problem even then: at that
same corner a genuine failure of Proposition 2 needs the whole segment inside two
adjacent cells, and no adjacent pair there carries more than $0.43$ of the mass,
bounding it by $10^{-37}$. Detector scores consume only likelihoods and are
start-independent either way, since the supremum is finite even when unattained;
what should not be interpreted is a fitted *coordinate* when Proposition 2 fails.

## 7.6 Structural validation

Automated tests verify Fisher orthonormality of the Fourier basis; exact cyclic
equivariance; the intended continuous-dimension increments; equality of Models A
and B for $g = 2$ and $g = 3$, where the fundamental component spans the full
nontrivial tangent space; recovery of planted relative shifts; and the absence of
sample-length dependence in Model C's penalty. All pass.

---

# 8. Monte Carlo design

## 8.1 Simulation grid

**Table 3.** Production design.

| Component | Values |
|---|---|
| Cyclic groups | $C_2$ … $C_6$ |
| Effect coordinates | 0.08, 0.12, 0.18, 0.25 |
| Segment length per side | 100, 200, 400, 800, 1,600, 3,200 |
| Total length $L$ | 200 … 6,400 |
| Alternative trials | 500 per configuration |
| Null calibration trials | 1,000 per configuration |
| Calibration target | 5% |
| Full detector | multinomial BIC |
| Configurations | 312 |
| Detector-level result rows | 936 |
| Simulated two-segment datasets | 468,000 |
| Base random seed | 20260713 |

Configurations carry deterministic content-derived seeds, so results do not
depend on grid ordering, worker count or completion order. Runs are checkpointed
and resumable.

Making that guarantee hold of every artefact takes one further step. Rows are
written as workers finish, so the results frame would otherwise be ordered by
completion. Point estimates are indifferent to it — every per-configuration seed
is content-derived and the crossover interpolator sorts its own inputs — but the
crossover bootstrap draws one random number per row, so an arbitrary order would
permute which draw lands on which segment length and its intervals would move
between runs of identical data. The runner therefore emits rows in a canonical
order and the bootstrap sorts within each group; a test asserts that shuffling
the input changes no analysis output.

## 8.2 Data-generating scenarios

**Exact orbit.** The left state is generated from a fundamental coordinate with
norm equal to the specified effect; the right state is its one-step cyclic
transform. Matches Model C and is contained in Models A and B.

**Independent fundamental.** Both regimes lie in the fundamental family but are
not related by a cyclic shift. For $g \geq 3$ the right coordinate has radius
0.72 times the left and angle 0.713 radians; for $g = 2$ it is $-0.55$ times the
left. Matches Model B, generally violates Model C. **Section 9.5 shows this
parameterisation is defective and gives a variant with a different confound.**

**Full-space higher-mode change.** A mode-2 Fourier component with amplitude 0.85
times the effect is added with opposite signs on the two sides of the boundary,
*both segments sharing one fundamental coordinate*. For $g = 4$ this is the
one-dimensional sign representation; for $g = 5,6$ a higher two-dimensional
Fourier mode. The change lies entirely outside the fundamental component, which
is what makes it a misspecification test. Confining the change to the mode-2 flip
matters: a scenario that also rotated the fundamental coordinate would carry a
full-strength exact-orbit component alongside the mode, leaving Model C half the
population gain and competitive under what is supposed to be misspecification.

**Approximate orbit.** A one-step orbit displaced perpendicular to the rotated
state by a controllable multiple of the effect. Zero deviation is the exact
orbit; growing it sweeps continuously toward an independent subspace change.

## 8.3 Two distinct analyses

**Raw MDL score analysis.** For each matched detector and scenario the mean
uncalibrated score is regressed on $L\,G_M$, $\log L$ and effect indicators. The
empirical penalty slope is minus the $\log L$ coefficient; its expected value is
$\Delta d/2$. Each group-level regression uses 24 design points (four effects by
six lengths), and because those aggregate means have unequal Monte Carlo
precision, both ordinary and variance-weighted fits are reported.

> **What this analysis can and cannot establish.** The detector score is
> constructed as $T_M = \widehat{G}_M - \mathrm{pen}_M$ with the proposed penalty
> *inserted*. Regressing the mean score on $L\,G_M$ and $\log L$ will therefore
> tend to return the inserted coefficient unless the maximised-likelihood bias
> itself carries substantial $\log L$ dependence. This is a finite-sample
> compatibility check on the implementation, **not an independent empirical
> discovery of the complexity law.** The decomposition in Appendix B makes the
> point exactly: the fitted slope is $\Delta d/2 - s$, where $s$ is the
> gain-residual slope, so the regression measures $s$ and nothing else. What the
> results below establish is that $s$ is small — that finite-sample corrections
> do not materially obscure the prescribed regular coefficient over the tested
> range. An independent validation would compare a fully specified universal code
> — normalised maximum likelihood, a prequential code, or a Bayesian marginal
> likelihood with a stated prior — against the BIC approximation. That is not
> done here.

**Calibrated power analysis.** At every configuration an additive critical value
is set to the empirical 95th percentile of 1,000 null scores, and power is
estimated from 500 alternative samples.

> **The penalty cancels from this analysis exactly.** $\mathrm{pen}_M$ is a
> deterministic function of $(d, L_1, L_2, g)$ — constant across draws within a
> configuration, and identical under the null and the alternative. So
> $$ q_{1-\alpha}\big(T_M \mid H_0\big) \;=\; q_{1-\alpha}\big(\widehat{G}_M \mid H_0\big) - \mathrm{pen}_M, $$
> and therefore
> $$ T_M > q_{1-\alpha}(T_M \mid H_0) \iff \widehat{G}_M > q_{1-\alpha}(\widehat{G}_M \mid H_0). $$
> Calibrated power is a functional of the maximised-likelihood gain alone. Every
> figure in Section 9.2 is unchanged if the penalties are set to zero, and the
> implementation confirms it: recomputing the grid's calibrated power from raw
> gains reproduces it bit for bit.
>
> The consequence is a scope limit, not a defect. **The calibrated comparison
> shows which likelihood statistic is better targeted at orbit-structured change;
> it is not evidence that the MDL penalties cause the sample-size reduction, and
> it cannot be.** A design that did test the penalties would have to use the raw
> rule — which for Model C fails for the reason Section 9.1 gives.

These curves support fair practical comparisons at a common nominal
false-positive target, read in those terms.

---

# 9. Results

All numbers below come from the committed production run, whose manifest records
the commit, environment, package versions and a SHA-256 per file.

The section is organised around the distinction Section 2.1 draws. Section 9.1
concerns the **structural code constant** and finds it wanting. Section 9.2
concerns **calibrated power**, which — as Section 8.3 shows — is a property of
the likelihood statistics and carries no information about the penalties at all.
The regressions that check the implemented penalty coefficients against their
prescribed values are an implementation matter and sit in Appendix B.

## 9.1 The two-part label constant does not control the raw rule

This is the paper's negative result, and it is the one place where the
complexity law meets data with something at stake.

Model C's penalty does not grow with $L$, and at $g = 2$ it is exactly zero, so
its raw zero-threshold rule gains no protection as the sample grows. Theorem 1
says why: at collapse the statistic converges to $W_{g,\rho}(0)$, an $O_p(1)$
law with $\Pr(W > 0)$ between $0.50$ and $0.83$. No constant subtracted from an
$O_p(1)$ statistic makes it a test at a nominal level unless the constant is
chosen from that statistic's own distribution — and $\log(g-1)$ is not.

Against the limit law, the two candidate structural constants fare like this:

**Table 4.** Asymptotic firing rate at collapse against each candidate label
cost, from 400,000 draws of $W_{g,1/2}(0)$, with the level-$0.05$ critical value
of Table 2 for comparison.

| $g$ | $\log(g-1)$ | $\Pr(W > \log(g-1))$ | $\log g$ | $\Pr(W > \log g)$ | $q_{0.95}$ |
|---:|---:|---:|---:|---:|---:|
| 2 | 0.000 | 0.500 | 0.693 | 0.156 | 1.596 |
| 3 | 0.693 | 0.337 | 1.099 | 0.220 | 2.463 |
| 4 | 1.099 | 0.234 | 1.386 | 0.171 | 2.520 |
| 5 | 1.386 | 0.175 | 1.609 | 0.138 | 2.584 |
| 6 | 1.609 | 0.141 | 1.792 | 0.117 | 2.610 |

Neither convention is close. Including the identity in the label code helps and
does not rescue: it moves $g=2$ from $0.50$ to $0.16$ and $g=6$ from $0.14$ to
$0.12$, still two to three times the chosen 5% level, and the level-$0.05$
threshold sits between $1.00$ and $1.77$ nats above $\log(g-1)$. The gap is
asymptotic, not a small-sample artefact.

The point is not that some other constant would do better. Corollary 3 fixes the
*unpenalised* gain's zero-threshold rate at exactly $(g-1)/g$, and subtracting any
constant $c$ from an $O_p(1)$ statistic leaves a false-alarm probability
$\Pr(W > c)$ that is positive and does not shrink with $L$ — the column above is
that probability at $c = \log(g-1)$ and at $c = \log g$. A constant *can* pin a chosen asymptotic level —
that is Corollary 2 — but it has to be read off the statistic's own distribution,
and $\log(g-1)$ is a codelength read off the model instead.

**A remark on mixing rather than maximising.** Writing $K = g-1$ and $a_r$ for
the profiled gain at shift $r$, the implemented two-part statistic is
$S_{\max} = \max_r a_r - \log K$ and the normalised mixture is
$S_{\mathrm{mix}} = \log\big(K^{-1}\sum_r e^{a_r}\big)$. These satisfy

$$ \max_r a_r - \log K \;\leq\; S_{\mathrm{mix}} \;\leq\; \max_r a_r, $$

so the mixture is never above the *unpenalised* maximum but is always at least
the penalised one — which is why replacing the two-part code by a mixture fires
*more* often at a zero threshold, not less, and does not help here. It is also
not a marginal code: $\eta$ is still profiled rather than integrated, so
averaging over $r$ alone produces neither a Bayesian marginal likelihood nor the
group average a maximal-invariant statistic would require.

The production grid shows the same thing at finite $L$, since its null states
sit at local values $\sqrt{L}\,\eta$ of order 1 to 20 rather than at collapse
exactly:

| detector | mean | worst |
|---|---:|---:|
| full | 0.0053 | 0.063 |
| fundamental | 0.0090 | 0.063 |
| shared orbit | **0.0531** | **0.281** |

The worst cases are at short segments and weak effects — the corner nearest
collapse — and are *not* concentrated at $g = 2$. The worst single row is
$g = 3$, `independent_fundamental`, effect $0.08$, $L = 200$, and $g = 2$ has the
lowest mean of any group order:

**Table 5.** Zero-threshold null rate for the shared-orbit detector, by group
order.

| $g$ | mean | worst |
|---:|---:|---:|
| 2 | 0.0384 | 0.275 |
| 3 | **0.0590** | **0.281** |
| 4 | 0.0570 | 0.208 |
| 5 | 0.0550 | 0.177 |
| 6 | 0.0532 | 0.139 |

Reading this off the limit law rather than off $g$ alone: $\Pr(W > \log(g-1))$
*decreases* in $g$ because the label cost grows faster than the statistic does,
which is exactly the ordering the grid shows once the local states are accounted
for.

Three conclusions follow, and they are different from each other.

*The raw MDL rule is not a usable test for Model C.* Every power comparison in
this paper is therefore made at a common calibrated 5%, and any performance claim
stated on raw scores would be reading threshold generosity as detection ability.

*The structural constant is a coding convention, not a testing threshold.*
$\log(g-1)$ remains a defensible two-part codelength; what Table 4 shows is that
a codelength is not a critical value and was never going to serve as one.

*A valid threshold exists and is derivable.* Corollary 2 supplies it, at the cost
of fixing an error rate — which is a different kind of object, obtained in a
different way, and Section 10 keeps it separate from the code lengths.

## 9.2 Calibrated comparison of structurally targeted statistics

Section 8.3 establishes that this comparison is invariant to the penalties: the
calibrated decision is $\widehat{G}_M > q_{0.95}(\widehat{G}_M \mid H_0)$
whatever $\mathrm{pen}_M$ is. What follows is therefore a comparison of how well
each *likelihood statistic* is targeted at orbit-structured change, and should be
read that way. It is a real and favourable result for the constrained
detectors — it is not evidence about the complexity law.


**Table 6.** Median calibrated 50%-power crossover-length ratios, exact-orbit
data. A ratio below one favours the numerator. $n$ is the number of the four
effect levels at which *both* detectors cross inside the grid, and so the number
of points the median is taken over.

Intervals are 95% bootstrap percentile intervals over 500 replications of the
whole pipeline — power draw, monotone stabilisation, interpolation and median
across effects — with the three detectors resampled **jointly**, as they are
scored.

| $g$ | shared / full | $n$ | shared / fundamental | $n$ | fundamental / full | $n$ |
|---:|---:|---:|---:|---:|---:|---:|
| 2 | 0.758 [0.707, 0.796] | 2 | 0.758 [0.707, 0.796] | 2 | 1.000 (exact) | 2 |
| 3 | 0.788 [0.752, 0.824] | 2 | 0.788 [0.752, 0.824] | 2 | 1.000 (exact) | 3 |
| 4 | 0.687 [0.657, 0.742] | 3 | 0.842 [0.805, 0.880] | 3 | 0.819 [0.786, 0.868] | 3 |
| 5 | 0.630 [0.567, 0.698] | 3 | 0.838 [0.790, 0.879] | 3 | 0.745 [0.715, 0.798] | 4 |
| 6 | 0.608 [0.563, 0.644] | 4 | 0.852 [0.808, 0.876] | 4 | 0.708 [0.666, 0.759] | 4 |

For exact-orbit data the shared detector required approximately 31%, 37% and 39%
fewer observations than the full detector for $g = 4,5,6$. Stated with the
uncertainty rather than as point estimates, the reduction against the full
detector is 26–34% at $g=4$, 30–43% at $g=5$ and 36–44% at $g=6$; against the
fundamental detector it is 12–19%, 12–21% and 12–19%. The shared-vs-fundamental
intervals exclude one at every $g \geq 4$, so that advantage is not a resampling
artefact; the fundamental-vs-full ratio is *exactly* one at $g = 2, 3$, as it
must be when the two models coincide.

**Four properties of these intervals**, stated rather than hidden.

*Detectors are resampled jointly, and must be.* All three score the same
alternative datasets, so their calibrated outcomes are strongly positively
correlated, and drawing each detector's power from its own marginal binomial
would discard that correlation and inflate every interval on a *ratio*. The grid
supplies its own calibration of how much. At $g = 2$ and $g = 3$ the fundamental
component is the whole nontrivial tangent space, so `full` and `fundamental` are
one detector scoring one set of datasets: their calibrated power is *equal* at
every design point and their mean scores agree to $10^{-11}$, making the
`fundamental/full` ratio identically $1$ with no variance. Independent
resampling returns roughly $\pm 9\%$ on that constant, and returns two different
intervals for `shared/full` and `shared/fundamental`, which are there the same
number.

The run therefore retains, per configuration, the joint calibrated-detection
pattern: the eight counts of outcome triples over the three detectors, which is
the sufficient statistic for resampling them together. Each replicate draws one
multinomial over those eight outcomes and reads all three powers off it, so a
dataset caught by every detector is redrawn as a dataset caught by every
detector. Degenerate rows then come back exact, `shared/full` and
`shared/fundamental` agree exactly where they are the same number, and every
interval that carries content is 20–45% narrower than the independent draw would
give — at $g=6$, widths of $0.081$, $0.068$ and $0.094$ against $0.087$, $0.108$
and $0.116$.

*Every replicate estimates the same quantity.* A ratio can only be formed at an
effect where both detectors cross inside the grid. The subset is frozen to the
one the point estimate uses, and replicates that cannot fill it are discarded
rather than re-medianed over whatever survived; the surviving count is reported.
At $g=3$ only 259 of 500 replicates can fill the three effects `fundamental/full`
uses — itself a signal of how thin that row is.

*Critical values are held fixed.* The resampling covers binomial power noise but
not the variability of the empirical 95th-percentile critical value, which is
re-estimated from 1,000 null draws at every configuration and carries its own
error. The intervals are therefore too narrow in that respect.

*A substantial minority of crossovers fall outside the grid.* Of 156 (group,
effect, detector) crossover estimates, 112 are interior, 25 lie above the longest
simulated length and 19 below the shortest. Only interior estimates enter the
medians. The excluded cases are not missing at random — they concentrate in the
weakest effects (below grid) and in the constrained detectors at large $g$ under
strong effects (above grid) — so the medians describe the interior of the design,
not the whole of it.

The $n$ column of Table 6 shows what that costs: at $g = 2$ and $g = 3$ a "median
across effects" is the midpoint of two numbers. And because each column keeps its
own surviving subset, the columns do not compose — at $g=5$ the reported
shared/full is $0.630$ while (shared/fundamental) $\times$ (fundamental/full) is
$0.624$, and similarly $0.608$ against $0.603$ at $g=6$. The three columns are
each defensible; they are not three views of one consistent set of crossovers.

The cumulative-maximum stabilisation of the power curves is defensible, since
true power is monotone in length, but it is one choice among several; isotonic
regression or a fitted parametric power curve would be a useful sensitivity
check, and is not performed here.

## 9.3 Misspecification


**Table 7.** Mean calibrated power under higher-mode full-space changes at
$L = 6{,}400$.

| $g$ | full | fundamental | shared orbit |
|---:|---:|---:|---:|
| 4 | 1.000 | 0.398 | 0.451 |
| 5 | 0.978 | 0.196 | 0.225 |
| 6 | 0.956 | 0.190 | 0.163 |

The full detector retains high power because it can represent the added mode.
The constrained detectors capture only what lies in their selected subspace. At
the population level the fundamental family retains 5.9% of the full gain at
$g=4$ and about 3% at $g = 5, 6$, and **the shared-orbit gain is negative** at
every group order — aligned pooling is worse than not aligning at all — which is
what drives the constrained detectors toward and below the nominal rate.

This result is important: the sample-efficiency advantage of the constrained
detectors is not a consequence of generally lower thresholds. **It is conditional
on structural correctness.**

## 9.4 Relative-shift recovery


**Table 8.** Relative-shift recovery on exact-orbit data at $L = 6{,}400$.

| $g$ | mean accuracy | minimum across effects |
|---:|---:|---:|
| 2 | 1.0000 | 1.000 |
| 3 | 0.9995 | 0.998 |
| 4 | 0.9990 | 0.996 |
| 5 | 0.9920 | 0.974 |
| 6 | 0.9875 | 0.958 |

The slight decline with $g$ is expected: the detector maximises over more
candidate shifts while adjacent orbit states become geometrically closer.

---

## 9.5 A defect in the independent-fundamental scenario


Section 8.2 fixes the angular offset at 0.713 radians while the one-step rotation
$2\pi/g$ *shrinks* as $g$ grows. The scenario therefore slides toward being an
orbit:

| $g$ | one-step rotation | distance from nearest orbit | best description |
|---:|---:|---:|---|
| 3 | 2.094 rad | 1.12 effects | fundamental |
| 4 | 1.571 rad | 0.76 effects | fundamental |
| 5 | 1.257 rad | 0.53 effects | approximate orbit |
| 6 | 1.047 rad | 0.40 effects | approximate orbit |

Since that distance is exactly the signal a shared-orbit fit cannot capture, the
scenario becomes progressively easier for Model C as $g$ grows, and **any
$g$-dependence in independent-fundamental results is partly an artifact of the
drift rather than a property of the detectors**. It also explains the otherwise
puzzling competitiveness of Model C on this scenario at $g = 6$.

A variant holds the orbit distance fixed at 1.5 effects for every $g$, placing
the right coordinate at the angular midpoint between adjacent orbit points and
solving the radius. On it the selector of Section 10 recovers the fundamental
model 100% of the time at every group order.

**That variant is not confound-free either, and the replacement confound is
larger.** With one free radius there is only one constraint to spend, so fixing
the distance to the orbit lets the size of the change run: at effect $0.25$ the
full population gain ranges from $0.0019$ at $g=2$ to $0.0573$ at $g=4$, a factor
of 30, against roughly 5 for the drift the variant was introduced to remove. At
$g \ge 3$ the change is about fifteen times stronger than in the original
scenario, which is most of why the selector recovers the fundamental model
without error; at $g = 2$ it runs the other way, the right state landing $0.50$
effects from the left one rather than $1.55$. The original scenario has the
mirror-image problem, its change size being near-constant at $0.655$ effects for
$g \ge 3$ but $1.55$ at $g = 2$. Holding both quantities fixed across $g$ needs a
second free parameter, which neither scenario has, so **neither supports a
cross-$g$ comparison** and the fixed-distance variant should be read as a check
that the selector is not fooled at a fixed distance, not as a $g$-comparison.

Distances here and in the table above are to the nearest **nonidentity** orbit
point, which is the signal a shared-orbit fit cannot capture and so the quantity
that matters. At $g \ge 3$ the identity is never the nearest point; at $g = 2$ it
is.

Holding the distance constant necessarily takes the right coordinate off the left
one's radius, and increasingly so with $g$. That is forced rather than chosen,
and two distinct bounds govern it.

On the *same circle* as the left coordinate, the furthest any point can sit from
every orbit point is the chord to the angular midpoint, $2\sin(\pi/2g)$ — not half
the adjacent-vertex chord $2\sin(\pi/g)$, since a chord's midpoint is not itself
on the circle. That is $1.000$ effects at $g=3$ and $0.518$ at $g=6$; the two
readings differ by 13% at $g=3$ and only 3% at $g=6$. Either way a unit radius
cannot reach $1.5$, so the right coordinate must leave that circle. And **no
fixed radius ratio could hold the distance constant**, because the radius solving
$r^2 - 2r\cos(\pi/g) + 1 = d^2$ depends on $g$.

Separately, $\sin(\pi/g)$ *is* a real quantity in this construction: minimising
over radius along the midpoint ray rather than maximising over the circle, it is
the smallest distance that ray can achieve, attained at $r = \cos(\pi/g)$. It is
the feasibility floor the implementation enforces — a target below it makes the
discriminant negative and the scenario undefined.

# 10. Choosing the geometry without an oracle

Every efficiency figure above is an **oracle** figure — the detector matching the
generating geometry is chosen in advance. An analyst does not know whether a
change is unrestricted, confined to the invariant subspace, or an exact orbit.
Deciding that is part of the problem.

Selection cannot use the detector scores. A score is measured against *its own*
null, and those nulls differ: Model A pools an unrestricted multinomial while
Models B, C and D pool a fundamental coordinate. What is comparable is the total
description length of the same data under each hypothesis. Writing $L(\cdot)$ for
that quantity, the detector scores are recovered exactly as differences,

$$ T_A = L(\text{null}_{\mathrm{full}}) - L(\text{full}), \qquad T_B = L(\text{null}_{\mathrm{fund}}) - L(\text{fund}), \qquad T_C = L(\text{null}_{\mathrm{fund}}) - L(\text{orbit}), $$

verified to $10^{-9}$. Selecting the shortest code over six candidates — the two
nulls and the four alternatives — answers whether a change occurred and what kind
in a single step.

**The convention, stated.** "Description length" here is not a complete universal
code, and calling these lengths *absolute* would overstate what is implemented.
Each $L(M)$ is a BIC/Laplace expression, and because Sections 10 and 11 turn
entirely on constants the six are worth writing out rather than compressing into
a single schema. A parameter block fitted from $k$ observations costs
$\tfrac{d}{2}\log k$, so a model that fits a block *per segment* pays twice, at
$L_1$ and $L_2$ separately:

$$
\begin{aligned}
L(\text{null}_{\mathrm{full}}) &= -\hat\ell_{\mathrm{full}}(\text{pooled}) + \tfrac{g-1}{2}\log L \\
L(\text{null}_{\mathrm{fund}}) &= -\hat\ell_{\mathrm{fund}}(\text{pooled}) + \tfrac{d_{\mathrm{fund}}}{2}\log L \\
L(\text{full}) &= -\hat\ell_{\mathrm{full}}(S_1) - \hat\ell_{\mathrm{full}}(S_2) + \tfrac{g-1}{2}\big[\log L_1 + \log L_2\big] \\
L(\text{fund}) &= -\hat\ell_{\mathrm{fund}}(S_1) - \hat\ell_{\mathrm{fund}}(S_2) + \tfrac{d_{\mathrm{fund}}}{2}\big[\log L_1 + \log L_2\big] \\
L(\text{orbit}) &= -\max_{r\neq e}\hat\ell_{\mathrm{fund}}(S_1 \oplus T_{-r}S_2) + \tfrac{d_{\mathrm{fund}}}{2}\log L + \log(g-1) \\
L(\text{approx}) &= -\max_{r\neq e}\hat\ell^{\,\mathrm{pen}}_{r} + \tfrac{d_{\mathrm{fund}}}{2}\log L + \log(g-1) + \mathrm{pen}_\delta
\end{aligned}
$$

The $\rho$-dependence of Section 4.2 is now visible rather than hidden: the
split models carry $\log L_1 + \log L_2 = \log L + \log[\rho(1-\rho)] + \log L$,
so differencing against the corresponding null returns exactly
$\mathrm{pen}_{\mathrm{split}}$, and the detector identities above hold to
$10^{-9}$. Models C and D share the null's single $\tfrac{d_{\mathrm{fund}}}{2}\log L$
because they fit one shared state, which is Proposition 1 in codelength form.

What all six omit is the $O(1)$ terms common to a fully specified code — the
Fisher-volume (Jeffreys) term $\log\int\sqrt{\det I(\theta)}\,d\theta$, the
parameter-space truncation, and the coding convention for $L$ itself.

**These are code lengths, and nothing here is a testing threshold.** The
critical values of Section 6.4 answer a different question and are not
interchangeable with the constants above; Section 9.1 is precisely the
observation that substituting one for the other fails.

Those omissions are not innocuous here, and we flag rather than hide the
consequence. Model C's incremental cost is *already* $O(1)$, and Model D is
separated from B and C purely in bounded terms, so an omitted constant of a nat
or two is the same size as the effects Section 11 measures. Comparisons between
A and B, whose separation grows like $\tfrac{g-1-d_{\mathrm{fund}}}{2}\log L$,
are safe from this; comparisons among C and D at fixed $L$ are not, and the
margins reported below should be read as provisional on the convention. A version
of this procedure with genuinely absolute lengths would need normalised maximum
likelihood, an explicit prequential code, or Bayesian marginal likelihoods with
stated priors. That is the single most valuable next step for Section 10 and is
listed in Section 14.

For this reason Section 11 reports differences $L(\text{null}_{\mathrm{full}}) -
L(M)$ against one fixed reference rather than raw lengths: the reference cancels,
and what remains is comparable under the stated convention even though the
individual lengths are not absolute.

At $g = 6$, effect 0.25, 200 trials per cell, seed 20260713
(`scripts/regenerate_selection_tables.py`):

| generated from | 200/side | 800/side | 3,200/side |
|---|---:|---:|---:|
| exact orbit → recovers shared orbit | 58% | 77% | 90% |
| higher mode → recovers full | 1% | 21% | 100% |
| no change → false-change rate | 4% | 2% | 0% |

The procedure also reports the **margin** over the runner-up. A small margin means
the data does not distinguish the geometries — information the oracle comparison
discards entirely.

**A degeneracy worth naming.** At $g = 2$ and $g = 3$ the fundamental component is
the whole nontrivial tangent space, so "full" and "fundamental" are the same
hypothesis and their code lengths agree to $\sim10^{-12}$. Selecting between them
reads floating-point noise. The procedure reports such ties explicitly and breaks
them toward the *less* structured candidate, so a tie never becomes a claim of
structure the data cannot support.

---

# 11. How much asymmetry can be tolerated

Model D answers the question Model C's all-or-nothing hypothesis cannot. Sweeping
the true deviation from an exact orbit at $g = 6$, effect 0.25, 1,200 observations
per side, with $\tau = 0.05$ (250 trials).

Following Section 10, the table reports **description-length savings against one
common reference**, $L(\text{null}_{\mathrm{full}}) - L(M)$ in nats, not detector
scores. Every column is measured against the same no-change code, so the columns
may legitimately be compared and the largest entry is the shortest total code.

| deviation | A full | B fundamental | C exact orbit | D approximate | best |
|---:|---:|---:|---:|---:|---|
| 0.00 | 5.54 | 23.78 | **27.55** | 27.24 | C |
| 0.25 | 6.12 | 24.54 | 26.15 | **27.12** | D |
| 0.50 | 15.88 | 34.26 | 32.29 | **35.34** | D |
| 0.75 | 24.61 | **42.78** | 34.27 | 41.06 | B |
| 1.00 | 40.14 | **58.38** | 38.64 | 51.64 | B |
| 1.50 | 73.24 | **91.58** | 55.00 | 77.02 | B |
| 3.00 | 190.98 | **209.11** | 80.72 | 147.29 | B |

Three regimes are visible. The rigid orbit code wins at an exact orbit; some
relational code wins out to a deviation of about half the effect size; beyond
that the relation is not worth encoding and the unconstrained subspace model
takes over. The relational advantage degrades *gradually* rather than at a
cliff, which is the qualitative point.

**How large is Model D's band, really?** Comparing means across a table row
understates the uncertainty, because the four models are evaluated on the same
datasets. The paired difference $\min_{A,B,C} L - L_D$ per dataset, in nats:

| deviation | mean paired advantage | s.e. | datasets where D is shortest |
|---:|---:|---:|---:|
| 0.00 | $-0.31$ | 0.04 | 23% |
| 0.25 | $+0.28$ | 0.06 | 58% |
| 0.50 | $-0.18$ | 0.12 | 51% |
| 0.75 | $-2.55$ | 0.37 | 43% |
| 1.00 | $-7.38$ | 0.79 | 34% |

Read honestly, this is a much weaker claim than a "middle band" suggests. Model
D's advantage is statistically distinguishable from zero at only one sweep point,
$0.25$, and there it is worth about a quarter of a nat — a rounding error next to
the tens of nats separating the model classes. At $0.50$ the mean advantage is
not distinguishable from zero, and D is shortest on barely half the datasets.
The paired statistic is also conservative by construction: taking a per-dataset
minimum over three competing models exceeds the minimum of their means, so the
comparison is against the best of three noisy quantities rather than against any
single rival. It is nonetheless the operationally relevant number, since a
selection procedure does take that minimum.

What the sweep supports is therefore narrow, and worth stating plainly:

> A shrinkage code on the deviation is never much worse than the better of the
> two endpoints, and is slightly better in a narrow neighbourhood of a nearly
> exact orbit. Its value is as a robust default when the analyst does not know
> whether the orbit is exact, not as a model that wins outright over a wide
> range.

Whether that band widens under a $\tau$ chosen by the data, rather than fixed at
$0.05$, is open; see Section 14.

---

# 12. Block families: separating group order from alphabet size

The direct model identifies the group order with the alphabet size. Raising $g$
therefore changes the number of candidate shifts, the simplex dimension, the
category sparsity and the geometric separation of neighbouring orbit states
simultaneously. A block family separates them.

## 12.1 The codon-phase model

A codon-phase model has $g = 3$ reading-frame phases and $a = 4$ nucleotide
symbols. **These quantities are not interchangeable.** Each phase carries its own
distribution over the alphabet, so a regime is a tuple of $g$ distributions and
one segment's data is a $g \times a$ count array. The group permutes phase blocks
and leaves the alphabet untouched.

$$ d_{\mathrm{full}} = g(a-1) = 3 \times 3 = 9 . $$

## 12.2 Fundamental isotypic dimension

The phase representation decomposes as
$\mathbb{R}^g_{\mathrm{phase}} = V_{\mathrm{triv}} \oplus V_{\mathrm{fund}} \oplus \cdots$
with $\dim V_{\mathrm{fund}} = 2$ for $g \geq 3$. Each phase coordinate carries an
$(a-1)$-dimensional alphabet-contrast tangent space, so

$$ d_{\text{phase-fund}} \;=\; \dim V_{\mathrm{fund}} \times (a-1) \;=\; 2 \times 3 \;=\; 6 . $$

The geometry is the direct model's **tensored** with the contrast space: the group
acts on the phase-mode index only, identically in every contrast direction, so
the rotation is $R \otimes I$. That is why the dimensions multiply, and why a
phase shift remains exactly a coordinate rotation (verified to $10^{-13}$).

## 12.3 The separation, measured

Under the null a regular $d$-dimensional split has $2 \times \text{gain} \sim \chi^2_d$,
so the mean raw gain is $d/2$. This reads the dimension off simulated data with
no scenario machinery:

| $g$ | $a$ | full: observed / $d/2$ | fundamental: observed / $d/2$ |
|---:|---:|---|---|
| 3 | 4 | 4.535 / 4.5 | 3.050 / 3.0 |
| 6 | 4 | 9.182 / 9.0 | **2.996 / 3.0** |
| 3 | 8 | 10.569 / 10.5 | 7.033 / 7.0 |

Holding the alphabet at four and doubling the group order from three to six
doubles the full model's effective dimension while leaving the fundamental
model's exactly where it was. **In the direct model that comparison cannot be
made, because the two quantities are the same number.**

The shared-orbit detector recovers a planted phase shift — a frameshift, in the
codon reading — for every nonidentity shift at $g = 3, 4, 6$.

## 12.4 Proposed experiment

A frameshift study would estimate phase-specific nucleotide parameters from
verified coding sequences; generate or curate sequences containing known $+1$ and
$+2$ frameshifts; compare the detectors at known boundaries; repeat with hidden
boundaries using the same location treatment for all models; and measure
false-positive rate, detection probability, localisation error and sample length.
Frameshift methods based on phase-specific emissions are established; the
contribution proposed here is not frameshift detection but a
representation-theoretic MDL decomposition of phase constraint. **This experiment
has not been performed.**

---

# 13. Limitations

1. **Independent categorical observations.** Markov, hidden-state and
   continuous-observation extensions may change both the population gain and the
   effective sample size.
2. **Known boundary.** A scanning or online procedure requires multiplicity
   correction, a stopping rule, false-alarm control over time and detection-delay
   analysis. This is *not* a changepoint-discovery method.
3. **Regular approximations, not exact codes.** All penalties are BIC/Laplace
   known-split increments. KT/Dirichlet mixtures and normalised maximum
   likelihood are not implemented, so no claim of exact codelength optimality is
   made. In particular the cross-model lengths of Section 10 omit the
   model-dependent $O(1)$ terms of a complete code, which is the same order as
   the effects Section 11 reports.
4. **Identity-information approximation in Model D.** The deviation penalty uses
   the isotropic reduction $\tfrac{d_{\mathrm{fund}}}{2}\log(1 + \tau^2 J_{\mathrm{eff}})$
   rather than $\tfrac12\log\det(I + \tau^2 J_{\mathrm{eff}})$ with the observed
   information. Fisher orthonormality holds at $\eta = 0$ only, so the two differ
   away from it; Section 3.4 derives the order as $O(\|\eta\|^2)$ at every $g$,
   which needs the observation that the chart metric's linear term is traceless
   (and, at $g \neq 3$, zero). Measured against brute-force marginalisation at
   $g=2$, $\tau = 0.15$, the isotropic form is within 0.003 nats at $\eta = 0$ and
   about 0.03 nats out at $\|\eta\| = 0.3$; the observed-information form is
   within 0.002 nats at both. Second order is not zero, and since Model D's
   conclusions live in bounded terms of order a quarter of a nat, this is not
   negligible.
5. **Size control at collapse is pointwise, not uniform.** Corollary 2 gives a
   valid asymptotic level at $h = 0$. Using it as an operating threshold needs
   collapse to be least favourable over all local nulls; Section 6.4 provides
   numerical evidence over the grid it evaluates — $\|h\| \leq 6$, one orbit
   sector, $\alpha = 0.05$, $\rho = 1/2$ — and states the general case as a
   conjecture. Nothing here establishes uniform size control, and the wording
   throughout is "over the evaluated local-null grid" for that reason.
6. **A label code, not a derived constant.** Model C's alternative excludes the
   identity, so it is not a proper extension of its own null, and $\log(g-1)$ is
   one defensible convention among several. Section 9.1 measures what the
   alternatives do; only a complete cross-model code would settle the constant
   rather than choosing it. Separately, a structural codelength constant is not
   thereby a testing threshold: a level-$\alpha$ threshold *can* be constant for
   fixed $g$ and $\rho$ — Corollary 2 supplies one — but it has to be derived
   from the statistic's null distribution rather than from the label code.
7. **Calibrated power is invariant to the penalties.** The comparison of Section
   9.2 is a comparison of likelihood statistics. It supports no claim that the
   complexity increments cause the sample-size reduction, and the design cannot
   be repaired to support one without abandoning calibration.
8. **Configuration-specific calibration.** Critical values are estimated
   separately for every configuration, which is appropriate for practical power
   comparison but means calibrated crossover curves must not be used to infer the
   raw MDL coefficient.
9. **A simulation study, and this is the largest gap.** Every result comes from
   synthetic data generated by the models under test. No real dataset is
   analysed, and the cyclic-orbit assumption has not been validated empirically
   on one. The selection procedure of Section 10 has therefore never been shown
   to produce interpretable margins outside its own simulator, which is the
   condition under which its output would mean anything to a practitioner. A
   controlled categorical phase dataset would suffice — it need not be the
   frameshift application of Section 12.4 — and would do more for the argument
   than any further simulation. We regard this as the necessary next step rather
   than as future work.
10. **Narrow exploration.** One fundamental Fourier mode and one higher-mode
   misspecification. Non-Abelian groups, representation multiplicities,
   stabilizers and approximate orbit relations beyond the single interpolation of
   Section 11 remain open. Neither independent-fundamental scenario supports a
   cross-$g$ comparison, for the reasons of Section 9.5.
11. **Narrow novelty claim.** Symmetry reduction, MDL dimension penalties and
    constrained changepoint models are established individually. The contribution
    is their explicit organisation into full independent, independent
    invariant-subspace, shared exact-orbit and approximate-orbit hypotheses;
    the collapse law of Section 6.4; and the negative result that follows from
    it.

---

# 14. Extensions

**14.1 Stabilizer-adaptive labels, outside this family.** If $\eta$ is invariant
under a nontrivial subgroup, the orbit contains fewer than $g$ distinct states and
the effective label cost should depend on the stabilizer order. Section 6.3 shows
this cannot arise in the direct model — the fundamental representation is
faithful, so the origin is the only fixed point — so the extension is motivated
only for nonfaithful higher modes, representation multiplicities, or non-Abelian
groups, where a stabilizer-adaptive code would have something to adapt to.

**14.2 Shrinking-effect asymptotics.** The present grid uses fixed effects. A
triangular array with $\epsilon_L = c/\sqrt{L}$ or $c\sqrt{\log L/L}$ answers a
different question, and the shared-orbit residual is most naturally studied there.

**14.3 Exact universal codes.** The unrestricted multinomial family admits the KT
mixture. Comparable exact or numerically integrated universal codes for the
fundamental and shared-orbit families would eliminate reliance on regular BIC
constants and permit direct finite-sample regret comparisons.

**14.4 Sequential extension.** An online extension would replace the fixed split
with a stopping time, raising average run length under the null, expected
detection delay under each geometric alternative, and the cost of repeatedly
searching over both boundary locations and relative group elements.

**14.5 Uniform size control at collapse.** *The highest-value extension.*
Section 6.4 gives the limit law and shows numerically that exact collapse is the
least favourable local null over the grid evaluated there. Proving
$$ \sup_h\, q_{1-\alpha}\big(W_{g,\rho}(h)\big) \;=\; q_{1-\alpha}\big(W_{g,\rho}(0)\big) $$
— for all $h$, and ideally for all $\rho$ and $\alpha$ — would upgrade Corollary 2
from a pointwise statement to genuine size control, and turn the shared-orbit
detector into a test with a guarantee rather than a threshold with evidence.

Two further pieces would follow naturally: the power of that test against
orbit-structured local alternatives, which is the quantity Section 9.2 currently
approaches by simulation; and the same analysis for Model D, whose deviation
prior interacts with the collapse geometry in a way the present treatment does
not describe.

**14.6 A complete cross-model code.** Section 10's selection procedure is only as
sound as its $O(1)$ convention. Replacing the BIC/Laplace lengths with normalised
maximum likelihood where it is well defined, an explicit prequential code, or
Bayesian marginal likelihoods under stated priors would make the cross-model
comparison genuinely absolute and would put the Model C and Model D margins — the
ones that live entirely in bounded terms — on a firm footing. This is the
highest-value extension *for the selection results specifically*.

Model C's label code is the concrete case. Its alternative excludes the identity,
so it is not a proper extension of its own null and $\log(g-1)$ is one convention
among several; Section 9.1 measures what switching to $\log g$ does to the raw
null rate, and the answer is a large fraction of the excess at small $g$. A
complete code would settle the constant rather than choosing it.

**14.7 Relation to group-invariant e-testing.** The e-statistics of
Pérez-Ortiz et al. [18] ask a different question — anytime-valid testing between
two group models, with the likelihood ratio of the maximal invariant as the
growth-rate-optimal statistic. That machinery is a plausible route to a sequential
version of the present comparison (Extension 14.4) that would not rely on a known
boundary or on two-part boundary coding at all.

Relating the two would first require replacing Model C's maximum over shifts with
the group average of Wijsman's representation [23], which is a different statistic and not obviously a
better one: Section 9.1's mixture remark shows the normalised average fires *more*
often at collapse than the penalised maximum, and profiling $\eta$ out before
averaging does not produce a marginal code in any case. The correspondence is
worth establishing; it is not available for free.

**14.8 Critical-value uncertainty in the crossover intervals.** The joint
resampling of Section 9.4 accounts for the correlation between detectors, and the
frozen effect subset keeps every replicate on one estimand, but the intervals
still hold each configuration's empirical 95th-percentile threshold fixed. That
threshold is itself estimated, from 1,000 null draws with about $\alpha n$
observations in the relevant tail, so the reported intervals remain too narrow.
Retaining the null scores, or the joint null pattern alongside the alternative
one, would let the calibration be resampled with everything else.

**14.9 A data-chosen $\tau$.** Model D's prior width is fixed at $0.05$
throughout. Estimating $\tau$ from the data, at the cost of encoding it, is the
natural way to ask whether Model D's narrow advantage (Section 11) widens into
something worth the extra model class.


---

# 15. Conclusion

Cyclic structure constrains a regime-change problem in two different ways. It can
restrict each regime independently to a low-dimensional invariant subspace, or it
can require the regimes to be exact transforms of one shared state. At a known
boundary these produce different leading MDL structures — $\tfrac{g-1}{2}\log L$,
$\tfrac{d_{\mathrm{fund}}}{2}\log L$, and a constant.

The constant is the interesting case, and it is not the degenerate one it first
appears. Orbit collapse is a benign finite-branch singularity: the label is
nonidentifiable where the $g-1$ branches meet, so the union is singular, but each
branch is regular and the real log canonical threshold is unchanged, so the zero
increment holds there too. What changes is that the profiled gain acquires a
nondegenerate $O_p(1)$ law, a maximum over the nonidentity shifts of differences
of Gaussian projection energies. That law is derived, proved, and validated
against the detector, and it does two things at once. It yields asymptotic
critical values, making the shared-orbit statistic operational. And it shows that
the two-part constant $\log(g-1)$ cannot control the raw rule. The unpenalised
gain's zero-threshold rate is exactly $(g-1)/g$; after subtracting the proposed
label cost the asymptotic firing rate is still $0.14$ to $0.50$ over
$g = 2,\dots,6$, and the level-$0.05$ threshold sits $1.0$ to $1.8$ nats above
the codelength. A codelength is not a critical value.

The simulation supports this rather than carrying it. Under matched geometry the
constrained detectors need materially less data, and under higher-mode
misspecification that advantage vanishes; but a per-configuration penalty cancels
from a calibrated comparison exactly, so those figures speak to the targeting of
the likelihood statistics and not to the codes.

What the framework establishes, then, is three separate things: that parameter
dimension fixes the regular $\log L$ rate; that a shared finite-group orbit
removes the incremental rate entirely; and that the bounded term left behind has
a describable distribution which the natural two-part constant fails to
accommodate. **Different levels of geometric knowledge define different
statistical questions**, and the answers differ not only in what they cost to
state but in whether the stated cost is enough to act on.

---

# Code and data availability

Implementation, tests and one complete production run are available at
<https://github.com/MacMayo1993/RegimeShift> under the MIT licence.

| Component | Location |
|---|---|
| Fourier geometry | `regimeshift/fourier.py` |
| Detectors A–D | `regimeshift/detectors.py` |
| Population gains | `regimeshift/gains.py` |
| Scenarios | `regimeshift/scenarios.py` |
| Monte Carlo engine | `regimeshift/simulation.py`, `runner.py` |
| Analysis | `regimeshift/analysis.py` |
| Model selection | `regimeshift/selection.py` |
| Sections 10 and 11 tables | `scripts/regenerate_selection_tables.py` |
| Block families | `regimeshift/blocks.py` |
| Run provenance | `regimeshift/manifest.py` |
| Collapse limit law | `regimeshift/collapse.py` |
| Committed production run | `results/v3-production/` |

The committed run ships a `run_manifest.json` recording the commit, environment,
package versions, timing and a SHA-256 for every file; tests verify each file
against its checksum. The base seed is 20260713. Over 700 automated tests cover
the structural properties and the statistical claims. A Lean development
(`lean/`) machine-checks the equivariance identity, the freeness of regular
orbits, the dimension of the fundamental component, and the $\rho(1-\rho)$
rearrangement of Section 4.2; the asymptotic coding law itself is assumed there
exactly where this paper assumes it.

---

# References

1. Schwarz, G. (1978). Estimating the dimension of a model. *The Annals of Statistics*, 6(2), 461–464.
2. Krichevsky, R. E., and Trofimov, V. K. (1981). The performance of universal encoding. *IEEE Transactions on Information Theory*, 27(2), 199–207.
3. Barron, A. R., Rissanen, J., and Yu, B. (1998). The minimum description length principle in coding and modeling. *IEEE Transactions on Information Theory*, 44(6), 2743–2760.
4. Grünwald, P. D. (2007). *The Minimum Description Length Principle*. MIT Press.
5. Rissanen, J. (2007). *Information and Complexity in Statistical Modeling*. Springer.
6. Grünwald, P. D., and Roos, T. (2019). Minimum description length revisited. *International Journal of Mathematics for Industry*, 11(1), 1930001. arXiv:1908.08484.
7. Watanabe, S. (2009). *Algebraic Geometry and Statistical Learning Theory*. Cambridge University Press.
8. Watanabe, S. (2013). A widely applicable Bayesian information criterion. *JMLR*, 14, 867–897.
9. Drton, M., and Plummer, M. (2017). A Bayesian information criterion for singular models. *JRSS B*, 79(2), 323–380.
10. Amari, S., and Nagaoka, H. (2000). *Methods of Information Geometry*. AMS and Oxford University Press.
11. Serre, J.-P. (1977). *Linear Representations of Finite Groups*. Springer.
12. Eaton, M. L. (1989). *Group Invariance Applications in Statistics*. IMS.
13. Yao, Y.-C. (1988). Estimating the number of change-points via Schwarz' criterion. *Statistics & Probability Letters*, 6(3), 181–189.
14. Davis, R. A., Lee, T. C. M., and Rodriguez-Yam, G. A. (2006). Structural break estimation for nonstationary time series models. *JASA*, 101(473), 223–239.
15. Niu, Y. S., Hao, N., and Zhang, H. (2016). Multiple change-point detection: a selective overview. *Statistical Science*, 31(4), 611–623.
16. Wang, G., Zou, C., and Yin, G. (2018). Change-point detection in multinomial data with a large number of categories. *The Annals of Statistics*, 46(5), 2020–2044.
17. Truong, C., and Runge, V. (2024). An efficient algorithm for exact segmentation of large compositional and categorical time series. *Stat*, 13(4), e70012.
18. Pérez-Ortiz, M. F., Lardy, T., de Heide, R., and Grünwald, P. D. (2024). E-statistics, group invariance and anytime-valid testing. *The Annals of Statistics*, 52(4), 1410–1432.
19. Yu, L., Zhao, R., Huang, J., Zhu, L., and Zhu, X. (2026). A sparse dimension-reduced subspace-based approach for detecting multiple change points in high-dimensional data. *Journal of Multivariate Analysis*, 213, 105594.
20. Csiszár, I., and Shields, P. C. (2004). Information theory and statistics: a tutorial. *Foundations and Trends in Communications and Information Theory*, 1(4), 417–528.
21. Rho, M., Tang, H., and Ye, Y. (2010). FragGeneScan: predicting genes in short and error-prone reads. *Nucleic Acids Research*, 38(20), e191.
22. Cancrini, N., Martinelli, F., Roberto, C., and Toninelli, C. (2008). Kinetically constrained spin models. *Probability Theory and Related Fields*, 140, 459–504.
23. Wijsman, R. A. (1967). Cross-sections of orbits and their application to densities of maximal invariants. *Proceedings of the Fifth Berkeley Symposium on Mathematical Statistics and Probability*, 1, 389–400.

> **Note on references.** Authors, titles and venues have been checked against
> the published record.

---

# Appendix A. Reproducibility

For each configuration, null samples are drawn from a scenario-specific no-change
distribution. The empirical 95th percentile of each detector's null scores is
computed with the higher-quantile convention. Alternative samples are then scored
against those fixed thresholds. The raw score, explicit penalty, population gain,
zero-threshold detection indicator, calibrated detection indicator, joint
calibrated-detection pattern, optimiser failure count and selected orbit shift
are recorded for every detector.

The joint pattern is the eight counts of calibrated-detection outcome triples
over the three detectors. It is the sufficient statistic that lets Section 9.4's
bootstrap resample the detectors together; summing the four patterns in which a
given detector fires returns that detector's calibrated power exactly, and the
analysis checks that identity before using them.

Population gains are computed at the distribution level. The full gain is
weighted JSD. The fundamental gain is the difference between the best two-segment
and pooled expected log-likelihoods inside the fundamental family. The
shared-orbit gain maximises the aligned shared-state expected likelihood over
nonidentity shifts.

The crossover estimator first replaces the empirical power sequence by its
cumulative maximum to suppress Monte Carlo reversals. When power crosses 0.5
within the tested grid, interpolation is performed linearly in log total length.
Crossovers below or above the grid are flagged and excluded from slope regressions
requiring internal estimates.

# Appendix B. Implementation validation: the penalty slopes

These regressions check that the *implemented* penalties are the prescribed ones
and that finite-sample corrections do not obscure them over the tested range.
They are an implementation matter rather than a result, for the reason Section
8.3 gives: the penalty is inserted into the score, so the fitted slope is
algebraically $\Delta d/2 - s$ and the regression measures the gain residual $s$.
They belong here rather than in Section 9, and no claim in the paper rests on
them.

## B.1 Full-model penalty


**Table 9.** Raw-score regression estimates for the full model.

| $g$ | $\hat\beta_{\mathrm{gain}}$ | slope (OLS) | slope (WLS) | predicted | $R^2$ |
|---:|---:|---:|---:|---:|---:|
| 4 | 0.995 | 1.488 | 1.504 | 1.500 | 0.99996 |
| 5 | 0.983 | 1.967 | 1.994 | 2.000 | 0.99989 |
| 6 | 0.996 | 2.488 | 2.518 | 2.500 | 0.99990 |

Gain coefficients are within 1.7% of one and the observed penalty slopes follow
the predicted increase with group order, agreeing to within 0.033 at every $g$.
Under variance weighting the agreement is within 0.018.

## B.2 Fundamental-subspace penalty


**Table 10.** Raw-score regression estimates for the fundamental-subspace model.

| $g$ | $\hat\beta_{\mathrm{gain}}$ | slope (OLS) | slope (WLS) | predicted | $R^2$ |
|---:|---:|---:|---:|---:|---:|
| 2 | 1.000 | 0.522 | 0.474 | 0.500 | 0.99994 |
| 3 | 0.994 | 0.996 | 0.979 | 1.000 | 0.99959 |
| 4 | 0.999 | 0.971 | 0.961 | 1.000 | 0.99960 |
| 5 | 0.981 | 0.964 | 0.988 | 1.000 | 0.99931 |
| 6 | 1.010 | 1.036 | 1.017 | 1.000 | 0.99973 |

The coefficients remain approximately constant from $g=3$ through $g=6$ even
though the full-simplex dimension grows from two to five. Read this in the terms
Section 8.3 sets: the penalty is *inserted* into the score, so the fitted slope
is $\Delta d/2 - s$ and the regression measures the gain residual $s$. **What the
constancy establishes is that $s$ stays small for the fundamental family as the
full-simplex dimension grows** — the two families are prescribed different
complexity laws, and the finite-sample corrections do not obscure either one over
the tested range. It is a compatibility check on the implementation, not
independent evidence that the laws differ; that would need a fully specified
universal code, which Section 14.3 lists and this paper does not implement. For
$g=2$ and $g=3$ the model spaces coincide, since $d_{\mathrm{fund}} = g-1$ there;
the informative separation begins at $g=4$.

## B.3 Shared exact-orbit penalty


**Table 11.** Residual raw-score log-length slopes for the shared-orbit model.

| $g$ | $\hat\beta_{\mathrm{gain}}$ | residual slope (OLS) | (WLS) | structural prediction |
|---:|---:|---:|---:|---:|
| 2 | 0.991 | $-0.228$ | $-0.074$ | 0 |
| 3 | 1.001 | 0.041 | 0.090 | 0 |
| 4 | 1.002 | 0.083 | 0.165 | 0 |
| 5 | 1.006 | 0.139 | 0.142 | 0 |
| 6 | 1.004 | 0.122 | 0.179 | 0 |

These residuals are far below the fundamental and full coefficients but are not
uniformly zero, and by the run's own standard errors several are not close to
zero either. Under variance weighting, which Section 8.3 prefers on this design,
the distances from the structural prediction are $1.6$, $1.4$, $3.7$, $4.4$ and
$8.0$ standard errors for $g = 2$ through $6$ — $0.179 \pm 0.022$ at $g=6$. The
ordinary fit ranks them differently ($2.7$ at $g=2$, $3.2$ at $g=6$) but agrees
that the largest are real. To call this small drift would understate it. The
decomposition below is what makes the significance informative rather than
alarming — these residuals are the maximised-likelihood gain's finite-sample
behaviour and nothing else — but a reproducible eight-sigma departure is a
quantified finite-sample deviation from the leading-order law, not noise to be
set aside. The evidence supports the qualified statement:

> The shared exact-orbit detector has a near-zero leading logarithmic coefficient
> relative to the regular split models, while finite-sample score behaviour
> retains small group-dependent drift.

**Where the residual comes from.** Because this implementation subtracts a
penalty it computes *exactly*, the gain and penalty coefficients need not be
estimated jointly. Writing the mean raw gain as $L G + a + s\log L$ and the exact
penalty as $\tfrac{\Delta d}{2}\log L + c$, the fitted slope satisfies the
algebraic identity

$$ \text{penalty slope} \;=\; \frac{\Delta d}{2} \;-\; s , $$

verified to $10^{-8}$. For Model C, $\Delta d = 0$, so **the residual slopes of
Table 11 are $-s$ and nothing else** — a property of the maximised likelihood gain
(shift maximisation and finite-sample MLE bias), not evidence of a hidden
continuous-dimension penalty.

## B.4 The split fraction affects only the bounded term


Section 4.2 predicts that $\rho$ moves the bounded term
$\tfrac{d}{2}\log[\rho(1-\rho)]$ while leaving the $\log L$ coefficient at $d/2$.
Running the $g = 4$ full-detector regression at a balanced split and at
$\rho = 0.25$ gives slopes of 1.404 and 1.519 against a prediction of 1.5.

Two qualifications, because this is the one section whose numbers are not from
the committed run. The production grid holds $\rho = 0.5$ throughout; these
figures come from
`tests/test_statistical_validation.py::test_penalty_slope_is_invariant_to_the_split_fraction`,
on exact-orbit rather than higher-mode data, over five segment lengths rather
than six, at 300 trials per cell rather than 500 and 1,000, and with an
acceptance band of $\pm 0.3$. That design difference, not a drift in $\rho$, is
why 1.404 sits further from 1.5 than Table 9's 1.488 for nominally the same
detector and group order. Adding $\rho$ to the production grid would remove the
discrepancy.

And note what is being checked. The $\rho(1-\rho)$ decomposition of Section 4.2 is
algebraically exact — it is a rearrangement, machine-checked in `Penalty.lean` —
so no simulation can falsify it. What the two regressions test is that the gain
residual $s$ of Section B.3 picks up no $\rho$-dependent $\log L$ structure of its
own. That is a real check, and a narrower one than confirming the prediction.

## B.5 Reading the slopes

Regressing the mean score on $L G_M$, $\log L$ and effect indicators should
recover a coefficient of one on the first and $-\Delta d/2$ on the second. The
very high $R^2$ values indicate that this affine approximation describes the
simulated mean scores over the tested grid.

Each group-level regression uses only 24 aggregate design points of unequal Monte
Carlo precision, so both ordinary and inverse-variance-weighted fits are
reported. On this grid the two agree to within 0.05 for the constrained
detectors; the largest disagreement is 0.03 for the full detector, and weighting
moves the estimate *toward* prediction in every case.

The calibrated-power crossover slopes of Section 9.2 answer a different question
again. Configuration-specific null calibration adds an empirical threshold that
may itself vary with length, so those crossovers are appropriate for comparing
practical sample requirements at a common false-positive target but their slopes
are not expected to equal the raw MDL coefficient.

# Appendix C. Relation to the broader geometric framework

The broader programme distinguishes three conceptual levels. *Topology*: a flat
bundle can carry nontrivial global holonomy despite vanishing local curvature.
*Dynamics*: a twisted or constrained generator can exhibit a spectral response to
the global sector. *Information*: a statistical observer faces model-dependent
description lengths determined by parameter dimension, invariant subspaces and
shared group relations.

**The present paper establishes a result only at the third level.** It does not
prove that topological holonomy, dynamical spectral gaps and MDL coefficients are
numerically identical.

A motivating numerical coincidence remains: the East model contains an
inverse-gap asymptotic involving $1/(2\ln 2)$, while a one-dimensional regular BIC
increment has coefficient $K^\ast = 1/(2\ln 2)$ when expressed in bits per $\ln L$.

This resemblance should be treated with caution. **$K^\ast$ arising as the
per-dimension BIC rate in bits is definitional, not a discovery** — it is
Schwarz's one-half expressed in base 2, and any quantity that counts half a
parameter per e-fold and reports in bits produces it. The resemblance therefore
carries weight only if the dynamical occurrence is *not* likewise a units
artefact — that is, only if $1/(2\ln 2)$ enters the East-model asymptotic from
the structure of the constrained generator rather than from a choice of base.
Settling that is the first question any bridge between the levels must answer,
and it is sharper than asking whether the numbers agree.

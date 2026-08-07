---
title: "Geometric Complexity in Cyclic Regime Changes"
subtitle: "An Introduction: Symmetry, Description Length, and the Cost of Crossing a Boundary"
author: Mac Mayo
date: Version 1.0
---

# Abstract

When a categorical process changes regime, the two regimes may be related by a
symmetry. This document develops what that relation costs, in the precise sense
of minimum description length, and shows that the answer separates into a
hierarchy of four hypotheses with different complexity increments.

For a known boundary and a regular full model of continuous dimension $d$,
fitting both segments independently costs a leading $\tfrac{d}{2}\log L$. If
both segments are confined to a $d_{\mathrm{fund}}$-dimensional invariant
subspace but vary freely within it, the leading cost is
$\tfrac{d_{\mathrm{fund}}}{2}\log L$. If the two segments are exact transforms of
one shared state under a finite cyclic group $C_g$, the continuous-dimension
increment is **zero**: a two-part code pays only the discrete label $\log(g-1)$,
constant in $L$. Between the last two sits a shrinkage model that encodes a small
deviation from exactness.

We derive these laws, construct a Fisher-orthonormal Fourier parameterisation in
which the group acts as a planar rotation, prove the zero-increment result on
regular orbit strata, implement all four detectors, and evaluate them in a
468,000-dataset Monte Carlo study over $C_2$ through $C_6$. Empirical penalty
slopes reproduce the predicted $\Delta d/2$ to within 0.033 across the full
model, hold constant at $d_{\mathrm{fund}}$ for the subspace model while the
simplex dimension grows, and remain near zero for the shared-orbit model. Under a
common 5% calibration the shared-orbit detector needs roughly a third fewer
observations than the unrestricted detector on exact-orbit data. Core geometric
and dimension-counting results are additionally machine-checked in Lean 4.

**Keywords:** minimum description length; cyclic groups; changepoint models;
representation theory; information geometry; categorical data; symmetry.

---

# 1. The question

Consider a sequence of categorical observations with a proposed boundary. The
standard question is whether one distribution explains the whole sequence or
whether the two sides need separate distributions. In an unrestricted model the
alternative assigns one probability vector to the left and another to the right,
and the extra complexity is governed by the dimension of the simplex.

Many regime structures are not arbitrary. Physical phases, reading-frame states,
rotational sectors and cyclic operating modes are related by symmetry. Once
symmetry is admitted there are two genuinely different ways to use it, and the
distinction is the subject of this document.

**First**, symmetry identifies a low-dimensional invariant subspace in which each
regime varies independently. This is an ordinary reduction of the admissible
parameter space: fewer parameters, same kind of hypothesis.

**Second**, the regimes may be required to be exact group transforms of *one
shared state*. This is a **relational** constraint, and it is different in kind.
It does not shrink the parameter space of a segment; it removes the independent
continuous parameter vector that a split ordinarily introduces, and replaces it
with a discrete group label.

The second possibility is what makes the problem interesting. A model that says
"the right regime is the left regime, rotated" introduces no new continuous
parameter at all. Whatever the left segment's state is, the right segment's state
is determined by it, up to a choice from a finite set. Under minimum description
length that costs $\log(g-1)$ nats — a constant — rather than a term growing with
$\log L$.

This yields four hypotheses:

- **A. Full independent.** Each segment receives an unrestricted parameter.
- **B. Independent fundamental-subspace.** Each segment receives its own
  parameter inside a selected invariant subspace.
- **C. Shared exact orbit.** The segments share one continuous parameter and
  differ by a relative group action.
- **D. Approximate orbit.** The segments share one continuous parameter up to a
  shrunk deviation, interpolating between B and C.

with continuous-dimension increments $g-1$, $d_{\mathrm{fund}}$, $0$ and
(asymptotically) $d_{\mathrm{fund}}$.

The analysis throughout is an offline **known-boundary** model comparison. It is
not a changepoint-discovery algorithm; unknown-boundary scanning adds a location
or multiplicity cost to every detector and is discussed separately in §12.

---

# 2. Setup

## 2.1 The comparison

Let $X_1,\dots,X_L$ be independent categorical observations with a proposed
boundary after $L_1$ of them, $L = L_1 + L_2$, writing $S_1$ and $S_2$ for the
two segments and $\rho = L_1/L$. **The boundary is supplied to the detector**; no
search over candidate locations is performed.

For a model class $M$, write $H_0^M$ for the one-regime null and $H_1^M$ for the
two-regime alternative. The detector score is

$$ T_M \;=\; \big[\hat\ell_1^M - \hat\ell_0^M\big] \;-\; \mathrm{pen}_M , $$

with $\hat\ell$ maximised log-likelihoods in nats. The raw rule declares a change
when $T_M > 0$; §9 explains why that rule is not comparable across models and why
reported power is calibrated instead.

## 2.2 Notation

| Symbol | Meaning |
|---|---|
| $g$ | group order of $C_g$; also the number of phase blocks in a block family |
| $a$ | alphabet size within one block; in the direct model $a = g$ |
| $L$, $L_1$, $L_2$ | total and per-segment lengths; $\rho = L_1/L$ |
| $d_{\mathrm{full}}$, $d_{\mathrm{fund}}$ | continuous dimensions of the full and fundamental families |
| $\eta$ | fundamental-family coordinate; $\|\eta\|_2$ is a Fisher norm |
| $r$ | relative group element; $R_g$ its action on coordinates, $T_r$ on distributions |

## 2.3 The group action

Let $C_g = \{e, r, \dots, r^{g-1}\}$ act through transformations $T_r$. In the
**direct** model the alphabet size equals the group order and $T_r$ cyclically
permutes category coordinates. In **block** models (§11) the group permutes phase
blocks while leaving the within-block alphabet untouched.

At a symmetric reference distribution $P_0$ the tangent space carries an induced
real representation of $C_g$. Let $V_{\mathrm{fund}}$ be the selected fundamental
invariant component, of real dimension $d_{\mathrm{fund}}$. For the direct model,

$$ d_{\mathrm{fund}} = \begin{cases} 1 & g = 2 \ \text{(the sign representation)} \\ 2 & g \geq 3 .\end{cases} $$

§10.2 records a machine-checked proof of this dimension count.

---

# 3. The four hypotheses

## 3.1 Model A: full independent change

Under the null both segments share one unrestricted $\theta$; under the
alternative $\theta_1,\theta_2$ are fitted independently, so for a direct
$g$-category multinomial

$$ \Delta d_A = d_{\mathrm{full}} = g-1 . $$

## 3.2 Model B: independent fundamental-subspace change

Let $\eta \in \mathbb{R}^{d_{\mathrm{fund}}}$ parameterise a smooth
exponential-family chart inside the invariant subspace,

$$ p(\eta) \;=\; \mathrm{softmax}\!\big(B\eta\big), $$

with the columns of $B$ spanning the fundamental component in logit coordinates.
Under the alternative $\eta_1,\eta_2$ are independent, so $\Delta d_B =
d_{\mathrm{fund}}$.

Model B asserts that both regimes lie in the same invariant subspace. It does
**not** require the right regime to be a group transform of the left.

## 3.3 Model C: shared exact-orbit transition

Under the alternative both segments share one $\eta$, and the right segment is
$T_r$ applied to it:

$$ p_1 = p(\eta), \qquad p_2 = T_r\, p(\eta), \qquad r \in \{1,\dots,g-1\}. $$

A simultaneous shift of both labels is observationally redundant, so the left
label is fixed as a gauge.

> **Proposition 1 (zero continuous increment on a regular stratum).** Assume
>
> 1. the group order $g$ is fixed and known, and the boundary is known;
> 2. the shared state $\eta$ lies on a **regular orbit** — its stabiliser in
>    $C_g$ is trivial, so the $g$ points $\{R^s\eta\}$ are distinct, and in
>    particular $\eta \neq 0$;
> 3. the chart of §4.3 is used, so the model is a regular exponential family in a
>    neighbourhood of $\eta$ with non-singular Fisher information;
> 4. the relative shift $r$ is encoded by the uniform two-part code over the
>    $g-1$ nonidentity elements.
>
> Then Model C introduces no continuous parameter beyond the null,
> $$ \Delta d_C = 0, $$
> and its total incremental cost over the null is $\log(g-1)$ nats, constant in
> $L$.

*Proof sketch.* The null's parameter space is the chart domain
$U \subseteq \mathbb{R}^{d_{\mathrm{fund}}}$, of dimension $d_{\mathrm{fund}}$.
The alternative's parameter space is $U \times \{1,\dots,g-1\}$: for each choice
of the discrete label $r$, the pair $(p_1,p_2)$ is the image of the single
continuous coordinate $\eta$ under $\eta \mapsto (p(\eta), T_r p(\eta))$. This is
a disjoint union of $g-1$ copies of a manifold of the *same* dimension as the
null's, so the continuous increment is zero and the only additional description
burden is the label.

Identifiability of that label is exactly assumption 2. If $\eta$ has trivial
stabiliser then $R^s\eta = R^{s'}\eta$ forces $s = s'$, so distinct labels give
distinct alternatives and the two-part code is decodable. $\square$

Assumption 2 is the one that bites. At $\eta = 0$ the orbit collapses to a point,
every $R^s$ acts trivially, and the parameter is a fixed point of the whole
group; more generally at any $\eta$ with nontrivial stabiliser the shift is not
identifiable. At such points ordinary dimension counting is not the correct
marginal-likelihood theory, and the singular framework of Watanabe and of
Drton–Plummer applies instead (§5.3). At $g=2$ there is a single nonidentity
shift and the label cost is $\log 1 = 0$.

## 3.4 Model D: approximate orbit

Exact symmetry is a strong assumption. Model D relaxes it:

$$ \eta_2 \;=\; R_g^{\,r}\,\eta_1 + \delta, \qquad \delta \sim \mathcal{N}(0,\tau^2 I), $$

with a shrinkage code on $\delta$. Setting $\tau = 0$ pins $\delta$ out and
recovers Model C exactly; $\tau \to \infty$ leaves $\delta$ free and recovers
Model B's maximised gain.

The cost of $\delta$ follows from Laplace's approximation, but the shared state
must be handled with care: $\eta_1$ is *estimated jointly with* $\delta$, not
conditioned on. The two are correlated, because moving $\eta_1$ and compensating
with $\delta$ leaves the right segment's fit unchanged. In Fisher-orthonormal
coordinates the joint information in $(\eta_1,\delta)$ is, per direction,

$$ \begin{pmatrix} L_1 + L_2 & L_2 \\ L_2 & L_2 \end{pmatrix}, $$

so the information that actually constrains $\delta$ is the Schur complement —
the profile information remaining after $\eta_1$ is optimised away:

$$ J_{\mathrm{eff}} \;=\; L_2 - \frac{L_2^2}{L_1+L_2} \;=\; \frac{L_1 L_2}{L_1+L_2} . $$

The deviation penalty is therefore

$$ \mathrm{pen}_\delta \;=\; \frac{d_{\mathrm{fund}}}{2}\,\log\!\Big(1 + \tau^2\,\frac{L_1L_2}{L_1+L_2}\Big), $$

the isotropic reduction of $\tfrac12\log\det(I + \tau^2 J_{\mathrm{eff}})$. Using
$L_2$ alone in place of $J_{\mathrm{eff}}$ would be correct only if the shared
state were known; on a balanced split it overstates the effective information by
a factor of two and inflates the penalty by $\tfrac{d}{2}\log 2$. That is a
bounded error, which leaves the leading coefficient alone — but Model D's entire
contribution lives in the bounded term, so it is exactly the term that must be
right.

Note what this does *not* do. For any **fixed** $\tau > 0$ the leading
coefficient is $d_{\mathrm{fund}}/2$, i.e. Model B's, not something in between.
The interpolation lives in the bounded term — the finite-sample regime, which is
where "how much deviation can be tolerated" is a meaningful question. A genuine
interpolation of the *leading* coefficient requires $\tau$ shrinking with $L$.

---

# 4. Complexity laws

## 4.1 The regular expansion

For a regular $d$-dimensional family fitted to $L$ observations, BIC, the regular
Laplace marginal likelihood and standard parametric-complexity expansions share
the leading form

$$ \tfrac{d}{2}\log L + O(1), $$

where the $O(1)$ term depends on the coding convention, prior, Fisher information
and parameter-space geometry.

## 4.2 The split increment

For a known split, replacing one $d$-dimensional fit by two costs, under that
expansion,

$$ \mathrm{pen}_{\mathrm{split}}(d; L_1, L_2) \;=\; \frac{d}{2}\Big[\log L_1 + \log L_2 - \log (L_1+L_2)\Big] . $$

*Derivation.* The alternative pays $\tfrac{d}{2}\log L_1 + \tfrac{d}{2}\log L_2$
for two independent $d$-dimensional fits at their respective sample sizes; the
null pays $\tfrac{d}{2}\log(L_1+L_2)$ for one. The increment is the difference.

Substituting $L_1 = \rho L$, $L_2 = (1-\rho)L$,

$$ \mathrm{pen}_{\mathrm{split}} \;=\; \frac{d}{2}\log L \;+\; \frac{d}{2}\log\big[\rho(1-\rho)\big], $$

so the coefficient of $\log L$ is $d/2$ while the split fraction affects only the
bounded term. The identity is algebraically exact within the approximation; it is
not an exact universal codelength.

This expression should not be read literally at very small sample sizes: nothing
constrains its sign, and when $L_1 L_2 < L$ it is negative. At
$L_1 = L_2 = 1$ it equals $-\tfrac{d}{2}\log 2$. Such sizes are far outside the
regime where the expansion approximates a codelength.

## 4.3 Unknown boundaries

If the boundary is unknown, a detector must encode or search over approximately
$L$ candidate locations. A two-part location code adds $\log L$ to **every**
model, so the leading coefficients become $\tfrac{d}{2}+1$ in each case. A
location cost cannot be applied to one detector only in a known-boundary
comparison; doing so confounds changepoint multiplicity with model dimension.

## 4.4 Units, and the constant $K^\ast$

In bits a regular penalty is
$\tfrac{\Delta d}{2}\log_2 L = \tfrac{\Delta d}{2\ln 2}\ln L$. Writing

$$ K^\ast \;=\; \frac{1}{2\ln 2} \;=\; 0.7213475\ldots $$

for the per-dimension penalty rate in bits per e-fold: **on regular strata and
under the coding convention of §8, every leading coefficient in this framework is
an integer multiple of $K^\ast$.** Model A pays $(g-1)K^\ast$, Model B pays
$d_{\mathrm{fund}}K^\ast$, and Model C pays zero. The hierarchy is a counting
statement — how many $K^\ast$ a model spends to cross the boundary.

$K^\ast$ is *definitional*, not empirical: it is Schwarz's one-half expressed in
base 2, and any quantity counting half a parameter per e-fold in bits produces
it.

The integer-multiple statement is conditional, and the qualification is not
cosmetic. It holds because regular counting charges an integer number of
parameters at half a $\log L$ each. Under singular learning theory the leading
coefficient is a real log canonical threshold, which need not be a half-integer
and so need not be an integer multiple of $K^\ast$ at all — and orbit collapse is
exactly such a singularity. The counting picture describes the regular part of
this problem, not all of it.

---

# 5. Cyclic Fourier geometry

## 5.1 The Fisher metric at the uniform distribution

For the direct $g$-category model let $u = (1/g,\dots,1/g)$. The tangent space is
$\{v \in \mathbb{R}^g : \sum_j v_j = 0\}$, and at $u$ the Fisher inner product is

$$ \langle v, w\rangle_F \;=\; \sum_j \frac{v_j w_j}{u_j} \;=\; g\sum_j v_j w_j . $$

## 5.2 The fundamental basis

For $g \geq 3$ set $\phi_j = 2\pi j/g$ and take the logit directions
$\sqrt{2}\cos\phi_j$ and $\sqrt{2}\sin\phi_j$; for $g = 2$ take the sign
direction. Their span is invariant under cyclic permutation, and in coefficient
space a one-step cyclic shift acts as a planar rotation by $2\pi/g$: writing
$\eta = (a,b)$ and using the angle-addition identities,

$$ a\cos(\phi_j - \theta) + b\sin(\phi_j-\theta) = (a\cos\theta - b\sin\theta)\cos\phi_j + (a\sin\theta + b\cos\theta)\sin\phi_j , $$

which is precisely rotation of $\eta$ by $\theta = 2\pi r/g$ composed with the
index shift $j \mapsto j - r$.

The implementation uses **Cartesian** coordinates rather than amplitude–angle
coordinates, which avoids the unidentifiable angular coordinate at zero
amplitude.

## 5.3 The positivity-preserving family

The direct fundamental family is defined through softmax logits,
$p(\eta) = \mathrm{softmax}(B\eta)$, with $B$ scaled so that
$\partial p/\partial\eta$ at $\eta = 0$ is the Fisher-orthonormal tangent basis —
making $\|\eta\|_2$ the Fisher norm of the local perturbation. If $R_g$ is the
fundamental rotation matrix, the family satisfies the **equivariance identity**

$$ p\big(R_g^{\,r}\eta\big) \;=\; T_r\, p(\eta). $$

*Why it holds.* Softmax commutes with any relabelling of categories, because the
normalising constant is a sum over all categories and is therefore invariant
under a bijection of the index set. So it suffices that the *logits* are
equivariant, which is §5.2. This identity is verified numerically to $10^{-13}$
through $g=8$ and proved for all $g$ in Lean (§10.2).

## 5.4 Local Jensen–Shannon geometry

For small perturbations the weighted Jensen–Shannon divergence between $p(\eta)$
and $p(R_g\eta)$ satisfies

$$ \mathrm{JSD} \;\to\; \tfrac{1}{8}\,\big\|R_g\eta - \eta\big\|^2 \;=\; \frac{1-\cos(2\pi/g)}{4}\,\|\eta\|^2 , $$

using $\|R\eta - \eta\|^2 = 2(1-\cos(2\pi/g))\|\eta\|^2$ for a planar rotation.
This gives coefficients $1/2$, $3/8$ and $1/4$ for $g = 2,3,4$. No extra factor
of $g$ appears when $\eta$ is measured in Fisher norm. Verified numerically to a
relative $3\times10^{-3}$.

---

# 6. Population gains and detection boundaries

## 6.1 Detector-specific gain

The relevant signal strength is the expected log-likelihood advantage *within the
model being fitted*. For the full multinomial model with equal segment sizes the
pooled null is the arithmetic mixture and the gain is the weighted
Jensen–Shannon divergence. For a restricted family the pooled null is generally
a KL projection rather than the arithmetic mixture, so **ordinary JSD must not be
reused for Models B, C and D**. Each population gain is computed by optimising
the corresponding population log-likelihood.

## 6.2 Local detection laws

If $G_M \approx c_M \epsilon^2$, the regular models have detection boundaries

$$ \epsilon^2 \sim \frac{d}{2}\,\frac{\log L}{L}, $$

while for a regular shared-orbit stratum with fixed $g$ and an explicit label
code,

$$ \epsilon^2 \sim \frac{\log(g-1)}{L} . $$

The stronger $L^{-1}$ scaling of Model C arises not merely from a smaller tangent
space but from *sharing the continuous state across the boundary*.

## 6.3 The singular qualification

At $\eta = 0$ all orbit elements coincide, the relative label is unidentifiable,
and the shared-orbit model is singular. The two-part code still adds no
independent continuous parameter vector, but exact Bayesian asymptotics near
orbit collapse may contain nonregular corrections. The empirical study therefore
treats a zero leading coefficient as a **structural prediction** while allowing
finite-sample residual length dependence.

---

# 7. The detectors

**Full (A).** A BIC-scored unrestricted multinomial detector, so all detectors
are comparable through maximised likelihood plus explicit complexity increments.
The penalty is the exact known-split increment with $d = g-1$.

**Fundamental (B).** The null fits one $\eta$ to combined counts; the alternative
fits $\eta_1,\eta_2$ separately. Optimisation is L-BFGS-B with analytic gradients
in Cartesian Fourier coordinates through a smooth softmax map. The penalty is the
known-split increment with $d = d_{\mathrm{fund}}$.

**Shared orbit (C).** The null fits one $\eta$ to combined counts. For each
nonidentity shift $r$ the right counts are aligned by $T_{-r}$, pooled with the
left counts, and fitted with one shared $\eta$; the alternative takes the shift
with the largest shared-state likelihood. The penalty is $\log(g-1)$ — no
location cost, no continuous-dimension increment.

**Approximate orbit (D).** For each nonidentity shift, $(\eta,\delta)$ are fitted
jointly by maximising
$\ell_1(\eta) + \ell_2(R^r\eta + \delta) - \|\delta\|^2/2\tau^2$
with analytic gradients. The penalty is $\log(g-1) + \mathrm{pen}_\delta$. That
the joint fit *estimates* $\eta$ rather than conditioning on it is precisely why
the penalty must use $L_1L_2/(L_1+L_2)$: the code and the penalty must describe
the same procedure. The closed form is grounded against brute-force numerical
marginalisation of the joint likelihood, agreeing to within 0.003 nats at $g=2$,
$\tau=0.15$ at the Fisher reference point.

**Numerical notes.** Convergence is judged on the first-order condition rather
than the optimiser's success flag: under a deliberately tight `ftol`, L-BFGS-B
reports abnormal termination whenever its line search cannot improve at machine
precision, which happens *at* the optimum. Observed gradient norms are around
$5\times10^{-8}$; trusting the flag would give roughly 2.5% false alarms, while a
relative gradient criterion gives zero failures across the production grid.

When a category has zero count — routine on short segments — the fundamental MLE
does not exist: the likelihood rises toward the simplex boundary and is
asymptotically flat along it. Different starts then halt at very different
coordinates (observed $\|\eta\|$ of 8 versus 15) while agreeing on the
log-likelihood to six decimals. Detector scores consume only likelihoods and are
therefore start-independent, but a fitted *coordinate* from a short segment
should not be interpreted.

---

# 8. Choosing the geometry without an oracle

Every efficiency figure in §9 is an **oracle** figure: the detector matching the
generating geometry is chosen in advance. An analyst does not know whether a
change is unrestricted, confined to the invariant subspace, or an exact orbit.
Deciding that is part of the problem.

Selection cannot use the detector scores, because a score is measured against
*its own* null and those nulls differ: Model A pools an unrestricted multinomial
while B, C and D pool a fundamental coordinate. What is comparable is the total
description length of the same data under each hypothesis. Writing $L(\cdot)$ for
that quantity, the scores are recovered exactly as differences,

$$ T_A = L(\text{null}_{\mathrm{full}}) - L(\text{full}), \quad T_B = L(\text{null}_{\mathrm{fund}}) - L(\text{fund}), \quad T_C = L(\text{null}_{\mathrm{fund}}) - L(\text{orbit}), $$

verified to $10^{-9}$. Selecting the shortest code over six candidates — two
nulls and four alternatives — answers whether a change occurred and what kind in
one step.

**The convention, stated.** "Description length" here is not a complete universal
code. Each $L(M)$ is

$$ L(M) \;=\; -\hat\ell_M \;+\; \tfrac{d_M}{2}\log L \;+\; c_M, $$

with structural constants $c_M$ the models do not share: $\log(g-1)$ for the
relative shift in C and D, and $\mathrm{pen}_\delta$ for D. What it omits is the
$O(1)$ terms of a fully specified code — the Fisher-volume (Jeffreys) term
$\log\int\sqrt{\det I(\theta)}\,d\theta$, parameter-space truncation, and the
coding convention for $L$ itself.

Those omissions are consequential and we flag rather than hide it. Model C's
incremental cost is *already* $O(1)$, and Model D is separated from B and C purely
in bounded terms, so an omitted constant of a nat or two is the same size as the
effects §9.5 measures. Comparisons between A and B, whose separation grows like
$\tfrac{g-1-d_{\mathrm{fund}}}{2}\log L$, are safe from this; comparisons among C
and D at fixed $L$ are not, and those margins are provisional on the convention.
Genuinely absolute lengths would require normalised maximum likelihood, an
explicit prequential code, or Bayesian marginal likelihoods with stated priors.

Accordingly §9.5 reports differences against one fixed reference rather than raw
lengths: the reference cancels, and what remains is comparable under the stated
convention even though the individual lengths are not absolute.

At $g = 6$, effect 0.25, 200 trials per cell:

| generated from | 200/side | 800/side | 3,200/side |
|---|---:|---:|---:|
| exact orbit → recovers shared orbit | 58% | 77% | 90% |
| higher mode → recovers full | 1% | 21% | 100% |
| no change → false-change rate | 4% | 2% | 0% |

The procedure also reports the **margin** over the runner-up: a small margin means
the data does not distinguish the geometries, information the oracle comparison
discards entirely.

**A degeneracy worth naming.** At $g = 2$ and $g = 3$ the fundamental component is
the whole nontrivial tangent space, so "full" and "fundamental" are the same
hypothesis and their code lengths agree to $\sim10^{-12}$. Selecting between them
reads floating-point noise. Such ties are reported explicitly and broken toward
the *less* structured candidate, so a tie never becomes a claim of structure the
data cannot support.

---

# 9. Empirical study

## 9.1 Design

| Component | Values |
|---|---|
| Cyclic groups | $C_2$ … $C_6$ |
| Effect coordinates | 0.08, 0.12, 0.18, 0.25 |
| Segment length per side | 100, 200, 400, 800, 1,600, 3,200 |
| Total length $L$ | 200 … 6,400 |
| Alternative trials | 500 per configuration |
| Null calibration trials | 1,000 per configuration |
| Calibration target | 5% |
| Configurations | 312 |
| Detector-level result rows | 936 |
| Simulated two-segment datasets | 468,000 |

Configurations carry deterministic content-derived seeds, so results do not
depend on grid ordering, worker count or completion order; runs are checkpointed
and resumable.

**Generating scenarios.** *Exact orbit*: the left state has fundamental
coordinate of norm equal to the effect, the right is its one-step cyclic
transform. *Independent fundamental*: both regimes lie in the fundamental family,
with the right coordinate placed at the angular midpoint between adjacent orbit
points and its radius solved so the distance from the nearest orbit point is held
at 1.5 effects for every $g$ — the quantity a shared-orbit fit cannot capture is
then constant across group orders, so comparisons across $g$ are not confounded
by it. *Full-space higher mode*: a mode-2 Fourier component at 0.85 times the
effect added with opposite signs on the two sides, both segments sharing one
fundamental coordinate, so the change lies entirely outside the fundamental
component. *Approximate orbit*: a one-step orbit displaced perpendicular to the
rotated state by a controllable multiple of the effect.

**What the score regression can establish.** The score is constructed as
$T_M = \widehat{G}_M - \mathrm{pen}_M$ with the proposed penalty *inserted*, so
regressing mean score on $L\,G_M$ and $\log L$ tends to return the inserted
coefficient unless the maximised-likelihood bias itself carries substantial
$\log L$ dependence. Writing the mean raw gain as $LG + a + s\log L$ and the exact
penalty as $\tfrac{\Delta d}{2}\log L + c$, the fitted slope satisfies

$$ \text{penalty slope} \;=\; \frac{\Delta d}{2} \;-\; s , $$

verified to $10^{-8}$. This is a finite-sample **compatibility check** on a
prescribed complexity law, not an independent empirical determination of it: the
regression measures $s$, the gain-residual slope, and nothing else. What follows
establishes that $s$ is small.

## 9.2 Full model

| $g$ | $\hat\beta_{\mathrm{gain}}$ | slope (OLS) | slope (WLS) | predicted | $R^2$ |
|---:|---:|---:|---:|---:|---:|
| 4 | 0.995 | 1.488 | 1.504 | 1.500 | 0.99996 |
| 5 | 0.983 | 1.967 | 1.994 | 2.000 | 0.99989 |
| 6 | 0.996 | 2.488 | 2.518 | 2.500 | 0.99990 |

Gain coefficients are within 1.7% of one; observed penalty slopes follow the
predicted increase with group order, agreeing to within 0.033 at every $g$, and
within 0.018 under variance weighting.

## 9.3 Fundamental subspace

| $g$ | $\hat\beta_{\mathrm{gain}}$ | slope (OLS) | slope (WLS) | predicted | $R^2$ |
|---:|---:|---:|---:|---:|---:|
| 2 | 1.000 | 0.522 | 0.474 | 0.500 | 0.99994 |
| 3 | 0.994 | 0.996 | 0.979 | 1.000 | 0.99959 |
| 4 | 0.999 | 0.971 | 0.961 | 1.000 | 0.99960 |
| 5 | 0.981 | 0.964 | 0.988 | 1.000 | 0.99931 |
| 6 | 1.010 | 1.036 | 1.017 | 1.000 | 0.99973 |

The coefficients stay approximately constant from $g=3$ to $g=6$ even though the
full-simplex dimension grows from two to five. **This is the clearest direct
evidence that an independently fitted fundamental family obeys a different
complexity law from the unrestricted multinomial family.** At $g=2,3$ the model
spaces coincide, since $d_{\mathrm{fund}} = g-1$ there; the informative
separation begins at $g=4$.

## 9.4 Shared exact orbit

| $g$ | $\hat\beta_{\mathrm{gain}}$ | residual slope (OLS) | (WLS) | structural prediction |
|---:|---:|---:|---:|---:|
| 2 | 0.991 | $-0.228$ | $-0.074$ | 0 |
| 3 | 1.001 | 0.041 | 0.090 | 0 |
| 4 | 1.002 | 0.083 | 0.165 | 0 |
| 5 | 1.006 | 0.139 | 0.142 | 0 |
| 6 | 1.004 | 0.122 | 0.179 | 0 |

These residuals are far below the fundamental and full coefficients but are not
uniformly zero, supporting the qualified statement:

> The shared exact-orbit detector has a near-zero leading logarithmic coefficient
> relative to the regular split models, while finite-sample score behaviour
> retains small group-dependent drift.

Because the implementation subtracts a penalty it computes *exactly*, the
identity above gives $\Delta d_C = 0$, so **these residual slopes are $-s$ and
nothing else** — a property of the maximised likelihood gain (shift maximisation
and finite-sample MLE bias), not evidence of a hidden continuous-dimension
penalty.

**Relative-shift recovery** on exact-orbit data at $L = 6{,}400$: mean accuracy
1.0000, 0.9995, 0.9990, 0.9920, 0.9875 for $g = 2,\dots,6$. The slight decline is
expected — the detector maximises over more candidate shifts while adjacent orbit
states become geometrically closer.

## 9.5 Calibrated advantage, and its limits

Median calibrated 50%-power crossover-length ratios on exact-orbit data, with 95%
bootstrap percentile intervals over 500 replications of the whole pipeline. A
ratio below one favours the numerator.

| $g$ | shared / full | shared / fundamental | fundamental / full |
|---:|---:|---:|---:|
| 2 | 0.758 [0.695, 0.814] | 0.758 [0.701, 0.821] | 1.000 [0.916, 1.098] |
| 3 | 0.788 [0.725, 0.856] | 0.788 [0.725, 0.853] | 1.000 [0.926, 1.094] |
| 4 | 0.687 [0.646, 0.764] | 0.842 [0.782, 0.911] | 0.819 [0.758, 0.909] |
| 5 | 0.630 [0.560, 0.726] | 0.838 [0.772, 0.903] | 0.745 [0.687, 0.826] |
| 6 | 0.608 [0.560, 0.647] | 0.852 [0.787, 0.896] | 0.708 [0.656, 0.772] |

The shared detector required approximately 31%, 37% and 39% fewer observations
than the full detector for $g = 4,5,6$. With uncertainty rather than point
estimates: 24–35% at $g=4$, 27–44% at $g=5$, 35–44% at $g=6$; against the
fundamental detector, 9–22%, 10–23% and 10–21%. The shared-vs-fundamental
intervals exclude one at every $g \geq 4$; the fundamental-vs-full intervals
include one at $g = 2,3$, as they must when the two models coincide.

**Three caveats, stated rather than hidden.** *Critical values are held fixed*:
the resampling covers binomial power noise but not the variability of the
empirical 95th-percentile critical value, so the intervals are too narrow in that
respect. *Detectors are resampled independently* although they score the same
datasets and are positively correlated; for a ratio of positively correlated
quantities this makes the intervals conservative. *A substantial minority of
crossovers fall outside the grid*: of 156 estimates, 112 are interior, 25 lie
above the longest simulated length and 19 below the shortest, and only interior
estimates enter the medians. The excluded cases are not missing at random — they
concentrate in the weakest effects and in the constrained detectors at large $g$
under strong effects — so the medians describe the interior of the design, not
the whole of it.

## 9.6 Misspecification

Mean calibrated power under higher-mode full-space changes at $L = 6{,}400$:

| $g$ | full | fundamental | shared orbit |
|---:|---:|---:|---:|
| 4 | 1.000 | 0.398 | 0.451 |
| 5 | 0.978 | 0.196 | 0.225 |
| 6 | 0.956 | 0.190 | 0.163 |

The constrained detectors lose most of their power when the change lies outside
the fundamental component, while the unrestricted detector retains it. This is
the expected cost of the constraint and the reason §8's selection step matters:
the efficiency gains of §9.5 are contingent on the geometry actually holding.

## 9.7 The raw rule is not comparable across models

The raw MDL rule $T_M > 0$ uses a different implicit false-positive rate for each
model, because the null distribution of $T_M$ differs across models. All power
figures above are therefore calibrated to a common 5% target by an additive
critical value estimated from null draws at each configuration. This is
appropriate for practical comparison, but calibrated curves must not be used to
infer the raw MDL coefficient.

## 9.8 How much asymmetry can be tolerated

Sweeping the true deviation from an exact orbit at $g = 6$, effect 0.25, 1,200
observations per side, $\tau = 0.05$, 250 trials. The table reports
description-length savings against one common reference,
$L(\text{null}_{\mathrm{full}}) - L(M)$ in nats, so columns are comparable and the
largest entry is the shortest total code.

| deviation | A full | B fundamental | C exact orbit | D approximate | best |
|---:|---:|---:|---:|---:|---|
| 0.00 | 5.54 | 23.78 | **27.55** | 27.24 | C |
| 0.25 | 6.12 | 24.54 | 26.15 | **27.12** | D |
| 0.50 | 15.88 | 34.26 | 32.29 | **35.34** | D |
| 0.75 | 24.61 | **42.78** | 34.27 | 41.06 | B |
| 1.00 | 40.14 | **58.38** | 38.64 | 51.64 | B |
| 1.50 | 73.24 | **91.58** | 55.00 | 77.02 | B |
| 3.00 | 190.98 | **209.11** | 80.72 | 147.29 | B |

Three regimes are visible: the rigid orbit code wins at an exact orbit; some
relational code wins out to a deviation of about half the effect size; beyond
that the relation is not worth encoding and the unconstrained subspace model
takes over. The relational advantage degrades *gradually* rather than at a cliff.

Comparing means across a row understates uncertainty, because the four models are
evaluated on the same datasets. The paired difference
$\min_{A,B,C} L - L_D$ per dataset:

| deviation | mean paired advantage | s.e. | datasets where D is shortest |
|---:|---:|---:|---:|
| 0.00 | $-0.31$ | 0.04 | 23% |
| 0.25 | $+0.28$ | 0.06 | 58% |
| 0.50 | $-0.18$ | 0.12 | 51% |
| 0.75 | $-2.55$ | 0.37 | 43% |
| 1.00 | $-7.38$ | 0.79 | 34% |

Read carefully, this supports a narrow claim. Model D's advantage is
statistically distinguishable from zero at only one sweep point, $0.25$, and
there it is worth about a quarter of a nat — a rounding error next to the tens of
nats separating the model classes. At $0.50$ the mean advantage is not
distinguishable from zero and D is shortest on barely half the datasets. What the
evidence supports is:

> A shrinkage code on the deviation is never much worse than the better of the
> two endpoints, and is slightly better in a narrow neighbourhood of a nearly
> exact orbit. Its value is as a robust default when the analyst does not know
> whether the orbit is exact, not as a model that wins outright over a wide
> range.

---

# 10. Machine-checked results

The geometric and dimension-counting content is formalised in Lean 4 against
mathlib, under `lean/`, and built in continuous integration with a check that no
proof is left incomplete.

## 10.1 What is proved

- **Equivariance.** `p(R_g^r η) = T_r p(η)` for all $g$, for the concrete basis
  $e_j = \exp(2\pi i j/g)$ — the identity of §5.3, proved rather than sampled.
- **Chart properties.** The softmax chart is positive and sums to one; the shift
  is an action of $C_g$ on categories.
- **Regular orbits.** For a faithful character and $\eta \neq 0$, the action is
  free: the $g$ orbit points are distinct, the stabiliser is trivial, and the
  relative shift is identified. This is the identifiability content of
  Proposition 1's assumption 2. The collapse at $\eta = 0$ is proved as well.
- **Complexity algebra.** The split increment and its $\rho(1-\rho)$
  reformulation; the statement that leading coefficients are integer multiples of
  $K^\ast$.

## 10.2 The dimension of the fundamental component

The two logit directions of §5.2 are linearly independent exactly when
$g \geq 3$ — the proof turns on $\sin(2\pi/g) \neq 0$, which holds precisely
because $0 < 2\pi/g < \pi$ there — while at $g = 2$ the sine direction vanishes
identically, leaving the one-dimensional sign direction. Hence the span has
dimension $d_{\mathrm{fund}}$ in both regimes, so the value used throughout is a
computed quantity rather than a stipulated numeral.

## 10.3 What is assumed

The step from "the continuous increment has dimension $d$" to "the penalty is
$\tfrac{d}{2}\log L$" requires the BIC/Laplace expansion, which mathlib does not
provide. The split penalty is therefore a *definition* in the formalisation,
matching §4.2, with only the algebra downstream of it proved. Identifying the
span of §10.2 with the fundamental isotypic component of the permutation
representation is likewise not formalised, nor are Fisher orthonormality, the
local JSD coefficient, or any singular asymptotics. The formalisation makes the
boundary between proved and assumed explicit; it does not move it.

---

# 11. Block families: separating group order from alphabet size

The direct model identifies group order with alphabet size, so raising $g$
changes the number of candidate shifts, the simplex dimension, the category
sparsity and the geometric separation of neighbouring orbit states all at once. A
block family separates them.

## 11.1 The codon-phase model

A codon-phase model has $g = 3$ reading-frame phases and $a = 4$ nucleotide
symbols. **These are not interchangeable.** Each phase carries its own
distribution over the alphabet, so a regime is a tuple of $g$ distributions and
one segment's data is a $g \times a$ count array. The group permutes phase blocks
and leaves the alphabet untouched, giving

$$ d_{\mathrm{full}} = g(a-1) = 3 \times 3 = 9 . $$

## 11.2 The fundamental isotypic dimension

The phase representation decomposes as
$\mathbb{R}^g_{\mathrm{phase}} = V_{\mathrm{triv}} \oplus V_{\mathrm{fund}} \oplus \cdots$
with $\dim V_{\mathrm{fund}} = 2$ for $g \geq 3$. Each phase coordinate carries an
$(a-1)$-dimensional alphabet-contrast tangent space, so

$$ d_{\text{phase-fund}} \;=\; \dim V_{\mathrm{fund}} \times (a-1) \;=\; 2 \times 3 \;=\; 6 . $$

The geometry is the direct model's **tensored** with the contrast space: the
group acts on the phase-mode index only, identically in every contrast direction,
so the rotation is $R \otimes I$. That is why the dimensions multiply, and why a
phase shift remains exactly a coordinate rotation (verified to $10^{-13}$).

## 11.3 The separation, measured

Under the null a regular $d$-dimensional split has $2 \times \text{gain} \sim
\chi^2_d$, so the mean raw gain is $d/2$. This reads the dimension off simulated
data with no scenario machinery:

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
codon reading — for every nonidentity shift at $g = 3,4,6$.

---

# 12. Scope

1. **Independent categorical observations.** Markov, hidden-state and
   continuous-observation extensions may change both the population gain and the
   effective sample size.
2. **Known boundary.** A scanning or online procedure requires multiplicity
   correction, a stopping rule, false-alarm control over time and detection-delay
   analysis. This is *not* a changepoint-discovery method.
3. **Regular approximations, not exact codes.** All penalties are BIC/Laplace
   known-split increments. KT/Dirichlet mixtures and normalised maximum
   likelihood are not implemented, so no claim of exact codelength optimality is
   made; the cross-model lengths of §8 omit the model-dependent $O(1)$ terms of a
   complete code, which is the same order as the effects §9.8 reports.
4. **Identity-information approximation in Model D.** The deviation penalty uses
   the isotropic reduction rather than
   $\tfrac12\log\det(I + \tau^2 J_{\mathrm{eff}})$ with observed information.
   Fisher orthonormality holds at $\eta = 0$ only, so the two differ at
   $O(\|\eta\|^2)$: measured against brute-force marginalisation at $g=2$,
   $\tau=0.15$, the isotropic form is within 0.003 nats at $\eta=0$ and about
   0.03 nats out at $\|\eta\| = 0.3$. Since Model D's conclusions live in bounded
   terms of order a quarter of a nat, this is not negligible.
5. **Singularity at orbit collapse.** The two-part code establishes the absence
   of an added continuous vector but does not derive the exact singular
   marginal-likelihood expansion. Characterising the real log canonical threshold
   there would replace a caveat with a result.
6. **Configuration-specific calibration.** Critical values are estimated
   separately for every configuration, appropriate for practical power comparison
   but meaning calibrated crossover curves must not be used to infer the raw MDL
   coefficient.
7. **A simulation study, and this is the largest gap.** Every result comes from
   synthetic data generated by the models under test. No real dataset is analysed
   anywhere, and the cyclic-orbit assumption has not been validated empirically on
   one. The selection procedure of §8 has therefore never been shown to produce
   interpretable margins outside its own simulator, which is the condition under
   which its output would mean anything to a practitioner. A controlled
   categorical phase dataset would suffice, and would do more for the argument
   than any further simulation.
8. **Narrow exploration.** One fundamental Fourier mode and one higher-mode
   misspecification. Non-Abelian groups, representation multiplicities,
   stabilisers, and approximate orbit relations beyond the single interpolation
   of §9.8 remain open.
9. **Narrow novelty claim.** Symmetry reduction, MDL dimension penalties and
   constrained changepoint models are individually established. The contribution
   is their explicit organisation into full independent, independent
   invariant-subspace, shared exact-orbit and approximate-orbit hypotheses,
   together with a matched empirical comparison and the observation that the
   third carries no continuous-dimension increment while the second carries
   $d_{\mathrm{fund}}$.

---

# 13. Relation to prior work

The ingredients are individually standard. Regular dimension penalties descend
from Schwarz (1978); MDL as a coding principle rather than a fixed penalty
formula is set out by Barron, Rissanen and Yu (1998) and Grünwald (2007).
Information-criterion changepoint selection dates to Yao (1988), with MDL
segmentation developed by Davis, Lee and Rodriguez-Yam (2006). The closest
categorical antecedents are Wang, Zou and Yin (2018) and Truong and Runge (2024),
both testing unrestricted probability change rather than invariant-subspace or
orbit alternatives. Dimension-reduced changepoints appear in Yu et al. (2026),
where the subspace is learned from sparse structure rather than fixed by a
representation. Group-invariant testing via maximal invariants and e-statistics
is developed by Pérez-Ortiz, Lardy, de Heide and Grünwald (2024). Singular-model
asymptotics, directly relevant at orbit collapse, are due to Watanabe (2009,
2013) and Drton and Plummer (2017).

The relation to group-invariant e-testing deserves particular care, being the
nearest neighbour. Pérez-Ortiz et al. show that among all e-statistics for
testing between two group models, the likelihood ratio of the maximal invariant
is growth-rate optimal, and that an anytime-valid test can be built on it. That
is a different question from this one — theirs a sequential testing problem with
no coding component, this a fixed-boundary description-length comparison across a
hierarchy — but the objects overlap: the maximal invariant of $C_g$ acting on the
fundamental subspace is exactly what Model C conditions on. Their framework is
the most promising route to a sequential version of this comparison, one
depending on neither a known boundary nor two-part boundary coding.

What we have not found elsewhere is the explicit four-way MDL separation of
unrestricted change, independent representation-constrained change, shared
finite-group-orbit change, and shrinkage toward such an orbit. This is a
search-based inference about absence rather than a proof of it, and should be
read accordingly.

---

# 14. Reproducing

```bash
pip install -e .

# the full design: 312 configurations, 936 detector rows, 468,000 datasets
python -m regimeshift run --grid production --out results/production --workers 4

# the machine-checked development
cd lean && lake exe cache get && lake build
```

Every run writes a manifest recording the commit, environment, package versions,
timing and a SHA-256 for each output file, so a result set carries its own
provenance. Configuration seeds are content-derived, so the numbers do not depend
on worker count or completion order.

---

# 15. Summary

A categorical regime change can be modelled at several levels of structural
constraint, and the levels are not stylistic variants of one another: they are
different statistical hypotheses with different description-length costs.

An unrestricted split pays $\tfrac{g-1}{2}\log L$. Confining both segments to the
fundamental invariant subspace pays $\tfrac{d_{\mathrm{fund}}}{2}\log L$ —
constant in $g$ once $g \geq 3$, even as the simplex grows. Requiring the
segments to be exact group transforms of one shared state pays no
continuous-dimension increment at all, only a $\log(g-1)$ label, because the
relational constraint removes the parameter vector rather than shrinking it. A
shrinkage code on the deviation sits between the last two, and its value is
robustness rather than dominance.

The empirical study reproduces these coefficients to within finite-sample noise
across $C_2$ through $C_6$, and shows the resulting practical advantage: roughly
a third fewer observations for the shared-orbit detector on exact-orbit data at a
common calibration. That advantage is contingent on the geometry actually
holding, which is why the selection procedure — choosing the geometry from data
rather than assuming it — is part of the method rather than an addendum to it.

The zero-increment result is stated for regular orbit strata. At orbit collapse
the model is singular, ordinary dimension counting does not apply, and what
replaces it is the most interesting open question here.

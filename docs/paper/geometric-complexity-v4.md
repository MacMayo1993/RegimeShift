---
title: "Geometric Complexity in Cyclic Regime Changes"
subtitle: "Full, Fundamental-Subspace, Shared-Orbit and Approximate-Orbit Models under Minimum Description Length"
author: Mac Mayo
date: Version 4.0
---

# Abstract

A categorical regime change can be modelled at several levels of structural
constraint. An unrestricted detector allows each segment to occupy the full
parameter space. A representation-constrained detector restricts each segment
independently to a selected invariant subspace. A shared-orbit detector imposes
the stronger requirement that both segment distributions arise from one
continuous state and differ only by a finite group action. These are different
statistical hypotheses with different minimum-description-length (MDL)
complexity increments.

For a known boundary and a regular full model of dimension $d$, independently
fitting both segments introduces the leading penalty $\tfrac{d}{2}\log L$.
Restricting both segments independently to a $d_{\mathrm{fund}}$-dimensional
fundamental representation gives $\tfrac{d_{\mathrm{fund}}}{2}\log L$. If the
segments share the same continuous orbit parameter and differ only by a relative
element of a fixed cyclic group $C_g$, the continuous-dimension increment is
zero; a two-part code pays only the discrete relative-label cost $\log(g-1)$,
subject to finite-sample and singular corrections near orbit-collapse points.

We derive these laws, construct a Fisher-orthonormal Fourier family for direct
cyclic categorical models, implement the detectors, and evaluate them in a
468,000-dataset Monte Carlo study over $C_2$ through $C_6$. Empirical full-model
penalty slopes were 1.488, 1.967 and 2.488 for $g=4,5,6$ against predictions
1.5, 2.0 and 2.5. Fundamental-subspace slopes were 0.522 for $g=2$ and
0.964–1.036 for $g=3$ through $g=6$, against predictions 0.5 and 1.0.
Shared-orbit residual slopes ranged from $-0.228$ to $0.139$ against a
structural prediction of zero. Under a common 5% null calibration the
shared-orbit detector required approximately 31%, 37% and 39% fewer observations
than the full detector on exact-orbit data for $g=4,5,6$. Under higher-mode
misspecification the full detector retained high power while the constrained
detectors did not.

Three extensions are new in this version. An **approximate-orbit** model
interpolates between the subspace and shared-orbit hypotheses with a shrinkage
code on the deviation, quantifying how far from exact symmetry a change can
drift before the relational code stops paying. A **model-selection** procedure
chooses the geometry from data rather than assuming it, on absolute code
lengths. A **block family** separates group order from alphabet size, removing a
confound present in the direct model.

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

1. the known-boundary MDL penalties for the four model classes;
2. a Fisher-orthonormal Fourier parameterisation for direct cyclic categorical
   families, and its block generalisation in which the group acts on phases
   while each phase carries its own alphabet;
3. implementations of all detectors with consistent known-boundary penalties and
   detector-specific population gains;
4. a large Monte Carlo comparison testing both asymptotic score behaviour and
   calibrated practical power;
5. a model-selection procedure that identifies the geometry from data rather
   than assuming it.

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
the split. We are not aware of prior work making that distinction explicitly,
but this is a search-based inference about absence rather than a proof of it,
and the claim should be read accordingly.

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
a change when $T_M > 0$. Section 9.6 shows why that rule is not comparable
across models and why all reported power is calibrated instead.

## 2.2 Full parameter family

Let $\{P_\theta\}$ be a regular parametric family of continuous dimension
$d_{\mathrm{full}}$. For a direct $g$-category multinomial model,
$d_{\mathrm{full}} = g-1$. The notation is kept general because applications with
several categorical blocks have larger dimensions; Section 12 develops that case.

## 2.3 Cyclic group action

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
label is fixed as a gauge. Therefore

$$ \Delta d_C = 0 . $$

The alternative introduces only the discrete relative shift. Under a uniform
two-part label code its cost is $\log(g-1)$ nats, which for fixed $g$ is constant
in $L$. At $g=2$ there is a single nonidentity shift and the cost is $\log 1 = 0$
— a fact with consequences developed in Section 9.6.

## 3.4 Model D: approximate orbit

Exact symmetry is a strong assumption. Model D relaxes it:

$$ \eta_2 \;=\; R_g^{\,r}\,\eta_1 + \delta, \qquad \delta \sim \mathcal{N}(0,\tau^2 I), $$

with a shrinkage code on $\delta$. Setting $\tau = 0$ pins $\delta$ out and
recovers Model C *exactly*; letting $\tau \to \infty$ leaves $\delta$ free and
recovers Model B's maximised gain. Under the Laplace approximation, and because
the coordinates are Fisher-orthonormal so the right segment contributes
information $L_2$ per unit in each direction, the deviation costs

$$ \mathrm{pen}_\delta \;=\; \frac{d_{\mathrm{fund}}}{2}\,\log\!\big(1 + L_2\tau^2\big) \ \text{nats}. $$

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
| D. Approximate orbit | one shared state, shrunk deviation | $d_{\mathrm{fund}}$ (fixed $\tau$) | $\log(g-1) + \tfrac{d_{\mathrm{fund}}}{2}\log(1+L_2\tau^2)$ |

---

# 4. MDL complexity laws

## 4.1 Regular codelength expansion

For a regular $d$-dimensional family fitted to $L$ observations, BIC, the regular
Laplace marginal likelihood and standard parametric-complexity expansions share
the leading form $\tfrac{d}{2}\log L + O(1)$, where the $O(1)$ term depends on
the coding convention, prior, Fisher information and parameter-space geometry.

## 4.2 Exact split increment

For a known split, the incremental regular complexity produced by replacing one
$d$-dimensional fit with two is

$$ \mathrm{pen}_{\mathrm{split}}(d; L_1, L_2) \;=\; \frac{d}{2}\Big[\log L_1 + \log L_2 - \log (L_1+L_2)\Big] . $$

When $L_1/L \to \rho$,

$$ \mathrm{pen}_{\mathrm{split}} \;=\; \frac{d}{2}\log L \;+\; \frac{d}{2}\log\big[\rho(1-\rho)\big] \;+\; O(1). $$

Thus the coefficient of $\log L$ is $d/2$, while the split fraction affects only
the bounded term. Section 9.7 verifies this prediction empirically — the first
time it has been checked in this programme.

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

for the per-dimension penalty rate in bits per e-fold, **every leading
coefficient in this framework is an integer multiple of $K^\ast$**: Model A pays
$(g-1)K^\ast$, Model B pays $d_{\mathrm{fund}}K^\ast$, and Model C pays zero. The
hierarchy is a counting statement — how many $K^\ast$ a model spends to cross the
boundary.

$K^\ast$ is *definitional*, not empirical: it is Schwarz's one-half expressed in
base 2, and any quantity counting half a parameter per e-fold in bits produces
it. Section 13 returns to this.

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
$g$ appears when $\eta$ is measured in Fisher norm. This is verified numerically
to a relative $3\times10^{-3}$.

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
$\epsilon^2 \sim \tfrac{d}{2}\,\tfrac{\log L}{L}$, while for a regular
shared-orbit stratum with fixed $g$ and an explicit label code
$\epsilon^2 \sim \tfrac{\log(g-1)}{L}$. The stronger $L^{-1}$ scaling of Model C
arises not merely from a smaller tangent space, but from *sharing the continuous
state across the boundary*.

## 6.3 Singular qualification

At $\eta = 0$ all orbit elements coincide, the relative label is unidentifiable,
and the shared-orbit model is singular. The two-part code still adds no
independent continuous parameter vector, but exact Bayesian asymptotics near
orbit collapse may contain nonregular corrections. The empirical study therefore
treats a zero leading coefficient as a **structural prediction**, while allowing
finite-sample residual length dependence.

---

# 7. Detector implementation

## 7.1 Full detector

A BIC-scored unrestricted multinomial detector, so that all detectors are
comparable through maximised likelihood plus explicit complexity increments. The
penalty is the exact known-split increment with $d = g-1$. An exact
KT/Dirichlet mixture implementation is *not* used; see Section 15.3.

## 7.2 Fundamental detector

The null fits one coordinate $\eta$ to the combined counts; the alternative fits
$\eta_1$ and $\eta_2$ separately. Optimisation is L-BFGS-B with analytic
gradients in Cartesian Fourier coordinates and a smooth softmax map. The penalty
is the exact known-split increment with $d = d_{\mathrm{fund}}$.

## 7.3 Shared-orbit detector

The null fits one $\eta$ to the combined counts. For each nonidentity shift $r$
the right counts are aligned by $T_{-r}$, pooled with the left counts, and fitted
with one shared $\eta$. The alternative takes the shift with the largest
shared-state likelihood. Its penalty is $\log(g-1)$, with no location cost and no
continuous-dimension increment.

## 7.4 Approximate-orbit detector

For each nonidentity shift, $(\eta,\delta)$ are fitted jointly by maximising the
penalised likelihood $\ell_1(\eta) + \ell_2(R^r\eta + \delta) - \|\delta\|^2/2\tau^2$
with analytic gradients. The penalty is $\log(g-1) + \mathrm{pen}_\delta$.

## 7.5 Numerical notes

Convergence is judged on the first-order condition rather than on the optimiser's
own success flag. Under a deliberately tight `ftol`, L-BFGS-B reports abnormal
termination whenever its line search cannot improve at machine precision, which
happens *at* the optimum: observed gradient norms are around $5\times10^{-8}$,
giving roughly 2.5% false alarms if the flag is trusted. Judging convergence by a
relative gradient criterion gives zero failures across the production grid.

A second numerical property is worth recording. When a category has zero count —
routine on short segments — the fundamental MLE does not exist: the likelihood
rises toward the simplex boundary and is asymptotically flat along that
direction. Different optimiser starts then halt at very different coordinates
(observed $\|\eta\|$ of 8 versus 15) while agreeing on the log-likelihood to six
decimals. Detector scores consume only likelihoods and are therefore
start-independent, but a fitted *coordinate* from a short segment should not be
interpreted.

## 7.6 Structural validation

Automated tests verify Fisher orthonormality of the Fourier basis; exact cyclic
equivariance; the intended continuous-dimension increments; equality of Models A
and B for $g = 2$ and $g = 3$, where the fundamental component spans the full
nontrivial tangent space; recovery of planted relative shifts; and the absence of
sample-length dependence in Model C's penalty. All pass.

---

# 8. Monte Carlo design

## 8.1 Simulation grid

**Table 2.** Production design.

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
and resumable. The committed run took 17.4 minutes on four workers.

## 8.2 Data-generating scenarios

**Exact orbit.** The left state is generated from a fundamental coordinate with
norm equal to the specified effect; the right state is its one-step cyclic
transform. Matches Model C and is contained in Models A and B.

**Independent fundamental.** Both regimes lie in the fundamental family but are
not related by a cyclic shift. For $g \geq 3$ the right coordinate has radius
0.72 times the left and angle 0.713 radians; for $g = 2$ it is $-0.55$ times the
left. Matches Model B, generally violates Model C. **Section 9.8 shows this
parameterisation is defective and gives a corrected variant.**

**Full-space higher-mode change.** A mode-2 Fourier component with amplitude 0.85
times the effect is added with opposite signs on the two sides of the boundary,
*both segments sharing one fundamental coordinate*. For $g = 4$ this is the
one-dimensional sign representation; for $g = 5,6$ a higher two-dimensional
Fourier mode. The change lies entirely outside the fundamental component, which
is what makes it a misspecification test.

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

**Calibrated power analysis.** At every configuration an additive critical value
is set to the empirical 95th percentile of 1,000 null scores, and power is
estimated from 500 alternative samples. These curves support fair practical
comparisons at a common nominal false-positive target but do not expose the raw
MDL penalty coefficient.

---

# 9. Results

All numbers below come from the committed production run, whose manifest records
the commit, environment, package versions and a SHA-256 per file.

## 9.1 Full-model penalty

**Table 3.** Raw-score regression estimates for the full model.

| $g$ | $\hat\beta_{\mathrm{gain}}$ | slope (OLS) | slope (WLS) | predicted | $R^2$ |
|---:|---:|---:|---:|---:|---:|
| 4 | 0.995 | 1.488 | 1.504 | 1.500 | 0.99996 |
| 5 | 0.983 | 1.967 | 1.994 | 2.000 | 0.99989 |
| 6 | 0.996 | 2.488 | 2.518 | 2.500 | 0.99990 |

Gain coefficients are within 1.7% of one and the observed penalty slopes follow
the predicted increase with group order, agreeing to within 0.033 at every $g$.
Under variance weighting the agreement is within 0.018.

> **Correction to v3.1.** The previous version reported 1.515, 2.119 and 2.468,
> and Section 9.1 discussed the $g=5$ value as exceeding 2.0 "by about 0.119, a
> finite-sample deviation of roughly two standard errors". At the same design and
> seed convention this run obtains 1.967, and the gain-residual slope there is
> 0.033, leaving no drift to explain. **The anomaly does not reproduce and the
> discussion of it has been removed.** It appears to have been Monte Carlo noise
> in the earlier run.

## 9.2 Fundamental-subspace penalty

**Table 4.** Raw-score regression estimates for the fundamental-subspace model.

| $g$ | $\hat\beta_{\mathrm{gain}}$ | slope (OLS) | slope (WLS) | predicted | $R^2$ |
|---:|---:|---:|---:|---:|---:|
| 2 | 1.000 | 0.522 | 0.474 | 0.500 | 0.99994 |
| 3 | 0.994 | 0.996 | 0.979 | 1.000 | 0.99959 |
| 4 | 0.999 | 0.971 | 0.961 | 1.000 | 0.99960 |
| 5 | 0.981 | 0.964 | 0.988 | 1.000 | 0.99931 |
| 6 | 1.010 | 1.036 | 1.017 | 1.000 | 0.99973 |

The coefficients remain approximately constant from $g=3$ through $g=6$ even
though the full-simplex dimension grows from two to five. **This is the clearest
direct evidence that an independently fitted fundamental family obeys a
different complexity law from the unrestricted multinomial family.** For $g=2$
and $g=3$ the model spaces coincide, since $d_{\mathrm{fund}} = g-1$ there; the
informative separation begins at $g=4$.

## 9.3 Shared exact-orbit penalty

**Table 5.** Residual raw-score log-length slopes for the shared-orbit model.

| $g$ | $\hat\beta_{\mathrm{gain}}$ | residual slope (OLS) | (WLS) | structural prediction |
|---:|---:|---:|---:|---:|
| 2 | 0.991 | $-0.228$ | $-0.074$ | 0 |
| 3 | 1.001 | 0.041 | 0.090 | 0 |
| 4 | 1.002 | 0.083 | 0.165 | 0 |
| 5 | 1.006 | 0.139 | 0.142 | 0 |
| 6 | 1.004 | 0.122 | 0.179 | 0 |

These residuals are far below the fundamental and full coefficients but are not
uniformly zero. The evidence supports the qualified statement:

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
Table 5 are $-s$ and nothing else** — a property of the maximised likelihood gain
(shift maximisation and finite-sample MLE bias), not evidence of a hidden
continuous-dimension penalty. This settles a question left open in v3.1.

## 9.4 Calibrated crossover advantage

**Table 6.** Median calibrated 50%-power crossover-length ratios, exact-orbit
data. A ratio below one favours the numerator.

| $g$ | shared / full | shared / fundamental | fundamental / full |
|---:|---:|---:|---:|
| 2 | 0.758 | 0.758 | 1.000 |
| 3 | 0.788 | 0.788 | 1.000 |
| 4 | 0.687 | 0.842 | 0.819 |
| 5 | 0.630 | 0.838 | 0.745 |
| 6 | 0.608 | 0.852 | 0.708 |

For exact-orbit data the shared detector required approximately 31%, 37% and 39%
fewer observations than the full detector for $g = 4,5,6$, and 15–16% fewer than
the fundamental detector. At $g = 2$ and $g = 3$ the fundamental/full ratio is
exactly 1.000, as it must be when the two models coincide.

Bootstrap intervals for these ratios are computed by resampling the whole
pipeline — power draw, monotone stabilisation, interpolation and median across
effects. Two limits are stated rather than hidden: the resampling covers binomial
power noise but not the variability of the empirical 95th-percentile critical
value, and detectors are resampled independently although they score the same
datasets and are positively correlated, which makes the ratio intervals
conservative.

## 9.5 Misspecification

**Table 7.** Mean calibrated power under higher-mode full-space changes at
$L = 6{,}400$.

| $g$ | full | fundamental | shared orbit |
|---:|---:|---:|---:|
| 4 | 1.000 | 0.398 | 0.451 |
| 5 | 0.978 | 0.196 | 0.225 |
| 6 | 0.956 | 0.190 | 0.163 |

The full detector retains high power because it can represent the added mode.
The constrained detectors capture only what lies in their selected subspace. At
the population level the fundamental family retains about 3% of the full gain and
**the shared-orbit gain is negative** — aligned pooling is worse than not
aligning at all — which is what drives the constrained detectors toward and below
the nominal rate.

This result is important: the sample-efficiency advantage of the constrained
detectors is not a consequence of generally lower thresholds. **It is conditional
on structural correctness.**

> **Correction to v3.1.** The previous version's higher-mode scenario rotated the
> fundamental coordinate *and* added the mode, so the change carried a
> full-strength exact-orbit component alongside it. Model C then retained 51–56%
> of the population gain and remained competitive under what was supposed to be
> misspecification. Confining the change to the mode-2 flip — both segments
> sharing one fundamental coordinate — is what reproduces the reported pattern.

## 9.6 The raw rule is not comparable across models

Model C's penalty does not grow with $L$, and at $g = 2$ it is exactly zero.
Consequently its raw zero-threshold rule provides no increasing protection under
the null. Across the production grid the mean and worst zero-threshold null rates
were:

| detector | mean | worst |
|---|---:|---:|
| full | 0.0053 | 0.063 |
| fundamental | 0.0090 | 0.063 |
| shared orbit | **0.0531** | **0.281** |

with the worst cases at $g = 2$, short segments and weak effects. This is the
same singular qualification as Section 6.3 — the relative label is unidentifiable
at orbit collapse — and it is precisely why every power comparison in this paper
is made at a **common calibrated 5%** rather than by raw rule. Any performance
claim stated on raw scores would be reading threshold generosity as detection
ability.

## 9.7 The split fraction affects only the bounded term

Section 4.2 predicts that $\rho$ moves the bounded term $\tfrac{d}{2}\log[\rho(1-\rho)]$
while leaving the $\log L$ coefficient at $d/2$. Repeating the $g = 4$ full-detector
regression at a balanced split and at $\rho = 0.25$ gives slopes of 1.404 and
1.519 against a prediction of 1.5. The prediction holds; this is the first time
it has been checked.

## 9.8 A defect in the independent-fundamental scenario

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

A corrected variant holds the orbit distance fixed at 1.5 effects for every $g$,
placing the right coordinate at the angular midpoint between adjacent orbit
points and solving the radius. On it the selector of Section 10 recovers the
fundamental model 100% of the time at every group order.

Holding the distance constant necessarily takes the right coordinate off the left
one's radius, and increasingly so with $g$. That is forced rather than chosen:
adjacent orbit points sit $2\sin(\pi/g)$ apart, so on a circle of the same radius
no point can be more than $\sin(\pi/g)$ from all of them — only 0.5 effects at
$g = 6$. **No fixed radius ratio in the original parameterisation could have held
the distance constant.**

## 9.9 Relative-shift recovery

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

At $g = 6$, effect 0.25:

| generated from | 200/side | 800/side | 3,200/side |
|---|---:|---:|---:|
| exact orbit → recovers shared orbit | 66% | 90% | 94% |
| higher mode → recovers full | 1% | 25% | 100% |
| no change → false-change rate | 4.5% | 0% | 0% |

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
per side, with $\tau = 0.05$ (mean MDL score in nats, 250 trials):

| deviation | A full | B fundamental | C exact orbit | D approximate | best |
|---:|---:|---:|---:|---:|---|
| 0.00 | 5.66 | 13.80 | **17.53** | 16.77 | C |
| 0.25 | 13.49 | 21.72 | **24.47** | 24.29 | C |
| 0.50 | 24.93 | 33.10 | 33.77 | **34.83** | D |
| 0.75 | 38.08 | 46.39 | 46.13 | **47.73** | D |
| 1.00 | 53.52 | 62.13 | 61.31 | **63.22** | D |
| 1.50 | 89.13 | **98.55** | 90.98 | 96.74 | B |
| 3.00 | 226.66 | **241.78** | 171.69 | 211.68 | B |

Three regimes. The rigid orbit code holds its edge out to a deviation of roughly
a quarter of the effect size. Some relational code keeps winning out to about
1.5. Beyond that the relation is not worth encoding. **Model D occupies a genuine
middle band from roughly 0.5 to 1.0 where it beats both endpoints** — Model C is
too rigid to fit the drift and Model B pays full dimension for it.

The relational advantage therefore degrades *gradually* rather than at a cliff,
which is what makes an approximate-orbit code worth having.

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

# 13. Relation to the broader geometric framework

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

This resemblance should be treated with more caution than it has been.
**$K^\ast$ arising as the per-dimension BIC rate in bits is definitional, not a
discovery** — it is Schwarz's one-half expressed in base 2, and any quantity that
counts half a parameter per e-fold and reports in bits produces it. The
resemblance therefore carries weight only if the dynamical occurrence is *not*
likewise a units artefact — that is, only if $1/(2\ln 2)$ enters the East-model
asymptotic from the structure of the constrained generator rather than from a
choice of base. Settling that is the first question any bridge between the levels
must answer, and it is sharper than asking whether the numbers agree.

---

# 14. Limitations

1. **Independent categorical observations.** Markov, hidden-state and
   continuous-observation extensions may change both the population gain and the
   effective sample size.
2. **Known boundary.** A scanning or online procedure requires multiplicity
   correction, a stopping rule, false-alarm control over time and detection-delay
   analysis. This is *not* a changepoint-discovery method.
3. **Regular approximations, not exact codes.** All penalties are BIC-style
   known-split increments. KT/Dirichlet mixtures and normalised maximum
   likelihood are not implemented, so no claim of exact codelength optimality is
   made.
4. **Singularity at orbit collapse.** The two-part code establishes the absence
   of an added continuous vector but does not derive the exact singular
   marginal-likelihood expansion.
5. **Configuration-specific calibration.** Critical values are estimated
   separately for every configuration, which is appropriate for practical power
   comparison but means calibrated crossover curves must not be used to infer the
   raw MDL coefficient.
6. **A simulation study.** Every result comes from synthetic data generated by
   the models under test. No real dataset is analysed, and the cyclic-orbit
   assumption has not been validated empirically on one.
7. **Narrow exploration.** One fundamental Fourier mode and one higher-mode
   misspecification. Non-Abelian groups, representation multiplicities,
   stabilizers and approximate orbit relations beyond the single interpolation of
   Section 11 remain open.
8. **Narrow novelty claim.** Symmetry reduction, MDL dimension penalties and
   constrained changepoint models are established individually. The contribution
   is their explicit organisation into full independent, independent
   invariant-subspace, shared exact-orbit and approximate-orbit hypotheses,
   together with a matched empirical comparison.

---

# 15. Extensions

**15.1 Stabilizer-adaptive labels.** If $\eta$ is invariant under a nontrivial
subgroup, the orbit contains fewer than $g$ distinct states and the effective
label cost should depend on the stabilizer order. A stabilizer-adaptive code could
improve finite-sample behaviour near symmetric strata.

**15.2 Shrinking-effect asymptotics.** The present grid uses fixed effects. A
triangular array with $\epsilon_L = c/\sqrt{L}$ or $c\sqrt{\log L/L}$ answers a
different question, and the shared-orbit residual is most naturally studied there.

**15.3 Exact universal codes.** The unrestricted multinomial family admits the KT
mixture. Comparable exact or numerically integrated universal codes for the
fundamental and shared-orbit families would eliminate reliance on regular BIC
constants and permit direct finite-sample regret comparisons.

**15.4 Sequential extension.** An online extension would replace the fixed split
with a stopping time, raising average run length under the null, expected
detection delay under each geometric alternative, and the cost of repeatedly
searching over both boundary locations and relative group elements.

**15.5 Singular analysis at orbit collapse.** Computing or bounding the relevant
real log canonical threshold would replace the qualified statement of Section 9.3
with a theorem.

---

# 16. Conclusion

Cyclic structure constrains a regime-change problem in two different ways. It can
restrict each regime independently to a low-dimensional invariant subspace, or it
can require the regimes to be exact transforms of one shared state. For a known
boundary these produce different leading MDL structures — $\tfrac{g-1}{2}\log L$,
$\tfrac{d_{\mathrm{fund}}}{2}\log L$, and a constant — with the last understood as
a two-part structural code on a regular orbit stratum.

The Monte Carlo results validate the full and independent-fundamental
coefficients closely, and support a near-zero leading coefficient for the shared
exact-orbit detector while revealing finite-sample residuals that an exact or
singular coding analysis would be needed to explain. The decomposition
$\text{slope} = \Delta d/2 - s$ localises those residuals in the likelihood gain
rather than in any hidden penalty. Under matched geometry the constrained
detectors reduce required sample length; under higher-mode misspecification the
advantage disappears and the full detector dominates.

The essential conclusion is not that cyclic groups possess one universal detection
threshold. It is that **different levels of geometric knowledge define different
statistical questions**, whose answers carry different complexity increments — and,
in bits, differ by whole multiples of $K^\ast$. That distinction is the correct
foundation for future work on cyclic changepoints, phase transitions, biological
reading frames and symmetry-aware sequential detection.

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
| Block families | `regimeshift/blocks.py` |
| Run provenance | `regimeshift/manifest.py` |
| Committed production run | `results/v3-production/` |

The committed run ships a `run_manifest.json` recording the commit, environment,
package versions, timing and a SHA-256 for every file; tests verify each file
against its checksum. The base seed is 20260713. Over 650 automated tests cover
the structural properties and the statistical claims.

---

# References

1. Schwarz, G. (1978). Estimating the dimension of a model. *The Annals of Statistics*, 6(2), 461–464.
2. Krichevsky, R. E., and Trofimov, V. K. (1981). The performance of universal encoding. *IEEE Transactions on Information Theory*, 27(2), 199–207.
3. Barron, A. R., Rissanen, J., and Yu, B. (1998). The minimum description length principle in coding and modeling. *IEEE Transactions on Information Theory*, 44(6), 2743–2760.
4. Grünwald, P. D. (2007). *The Minimum Description Length Principle*. MIT Press.
5. Rissanen, J. (2007). *Information and Complexity in Statistical Modeling*. Springer.
6. Grünwald, P. D., and Roos, T. (2019). Minimum description length revisited. arXiv:1908.08484.
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
17. Truong, C., and Runge, V. (2024). An efficient algorithm for exact segmentation of large compositional and categorical time series. *Stat*, 13(1).
18. Pérez-Ortiz, M. F., Lardy, T., de Heide, R., and Grünwald, P. D. (2024). E-statistics, group invariance and anytime-valid testing. *The Annals of Statistics*, 52(4).
19. Yu, L., Zhao, R., Huang, J., Zhu, L., and Zhu, X. (2026). A sparse dimension-reduced subspace-based approach for detecting multiple change points in high-dimensional data. *Journal of Multivariate Analysis*.
20. Csiszár, I., and Shields, P. C. (2004). Information theory and statistics: a tutorial. *Foundations and Trends in Communications and Information Theory*, 1(4), 417–528.
21. Rho, M., Tang, H., and Ye, Y. (2010). FragGeneScan: predicting genes in short and error-prone reads. *Nucleic Acids Research*, 38(20), e191.
22. Cancrini, N., Martinelli, F., Roberto, C., and Toninelli, C. (2008). Kinetically constrained spin models. *Probability Theory and Related Fields*, 140, 459–504.

> **Note on references.** Entries 6, 15–19 were added in this version from an
> external prior-art review and have not been independently re-verified against
> the published record. Verify before submission.

---

# Appendix A. Reproducibility

For each configuration, null samples are drawn from a scenario-specific no-change
distribution. The empirical 95th percentile of each detector's null scores is
computed with the higher-quantile convention. Alternative samples are then scored
against those fixed thresholds. The raw score, explicit penalty, population gain,
zero-threshold detection indicator, calibrated detection indicator, optimiser
failure count and selected orbit shift are recorded for every detector.

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

# Appendix B. Reading the empirical slopes

The raw-score regression is the primary test of the theoretical penalty
coefficients. Regressing the mean score on $L G_M$, $\log L$ and effect indicators
should recover a coefficient of one on the first and $-\Delta d/2$ on the second.
The very high $R^2$ values indicate that this affine approximation describes the
simulated mean scores over the tested grid.

Each group-level regression uses only 24 aggregate design points of unequal Monte
Carlo precision, so both ordinary and inverse-variance-weighted fits are reported.
On this grid the two agree to within 0.05 for the constrained detectors; the
largest disagreement is 0.03 for the full detector, and weighting moves the
estimate *toward* prediction in every case.

The calibrated-power crossover slopes answer a different question. Configuration-
specific null calibration adds an empirical threshold that may itself vary with
length, so those crossovers are appropriate for comparing practical sample
requirements at a common false-positive target but their slopes are not expected
to equal the raw MDL coefficient.

# Appendix C. Changes from version 3.1

| Change | Section |
|---|---|
| $g=5$ full-model slope anomaly removed — does not reproduce | 9.1 |
| Higher-mode scenario corrected: the mode-2 flip is the whole change | 8.2, 9.5 |
| Independent-fundamental scenario shown defective; corrected variant added | 9.8 |
| Model D (approximate orbit) added | 3.4, 11 |
| Model selection without an oracle added | 10 |
| $K^\ast$ named; hierarchy restated as a counting statement | 4.4, 13 |
| Caution added: $K^\ast$ is definitional, not a discovery | 13 |
| Block families implemented; dimension separation measured | 12 |
| Penalty-slope decomposition $\Delta d/2 - s$ established | 9.3 |
| Non-comparability of raw scores across models made explicit | 9.6 |
| Split-fraction prediction verified | 9.7 |
| Weighted regressions and bootstrap intervals reported | 8.3, 9.4 |
| Optimiser convergence criterion corrected | 7.5 |
| Related-work positioning added | 1.1 |
| Scope narrowed: known-boundary comparison, simulation study | 1, 14 |

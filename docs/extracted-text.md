# Manuscript text (extracted)

Text extraction of `Geometric_Complexity_Cyclic_Regime_Changes_v3_1.docx`, for
grepping and for reviewing diffs. **The .docx is the source of truth.**

Equations are Office Math (OMML), not images, and are included here inline in
backticks. OMML carries no spacing or delimiters, so a rendered expression like
"log(g - 1)" appears as `logg-1`; read them as symbol sequences, not as
typeset formulas. The five figures are images and do not appear.

---

Geometric Complexity in Cyclic Regime Changes

Full, Fundamental-Subspace, and Shared-Orbit Models under Minimum Description Length

Mac Mayo

July 2026

Abstract

A categorical regime change can be modeled at several levels of structural constraint. An unrestricted detector allows each segment to occupy the full parameter space. A representation-constrained detector restricts each segment independently to a selected invariant subspace. A shared-orbit detector imposes the stronger requirement that both segment distributions arise from one continuous state and differ only by a finite group action. These are different statistical hypotheses and have different minimum-description-length (MDL) complexity increments.

For a known boundary and a regular full model of dimension `dfull`, independently fitting both segments introduces the leading penalty `dfull2logL`. Restricting both segments independently to a `dfund`-dimensional fundamental representation gives `dfund2logL`. If the segments share the same continuous orbit parameter and differ only by a relative element of a fixed cyclic group `Cg`, the continuous-dimension increment is zero; a two-part code pays only the discrete relative-label cost `logg−1`, subject to finite-sample and singular corrections near orbit-collapse points.

We derive these three laws, construct a Fisher-orthonormal Fourier family for direct cyclic categorical models, implement all three detectors, and evaluate them in a 468,000-dataset Monte Carlo study over `C2` through `C6`. The empirical full-model penalty slopes were 1.515, 2.119, and 2.468 for `C4,C5,C6`, compared with predictions 1.5, 2.0, and 2.5. Fundamental-subspace slopes were 0.457 for `C2` and 0.967–1.047 for `C3` through `C6`, compared with predictions 0.5 and 1.0. Shared-orbit score regressions retained small residual slopes, while direct crossover regressions for `C5` and `C6` were 0.037 and 0.062, supporting a near-zero leading logarithmic coefficient rather than establishing exact finite-sample constancy. Under a common 5% null calibration, the shared-orbit detector required approximately 30%, 33%, and 39% fewer observations than the full detector on exact-orbit data for `C4,C5,C6`. Under higher-mode misspecification, the full detector retained high power while the constrained detectors did not. The results confirm that full-space change, invariant-subspace change, and exact-orbit transition are statistically distinct regimes of geometric constraint.

Keywords: minimum description length; cyclic groups; changepoint models; representation theory; information geometry; categorical data; model selection; symmetry.

1. Introduction

A standard two-segment regime comparison asks whether one distribution explains an entire sample or whether separate distributions are needed on either side of a proposed boundary. In an unrestricted categorical model, the alternative assigns one probability vector to the left segment and another to the right. The additional complexity is governed by the dimension of the entire simplex.

Many regime structures are not arbitrary. Physical phases, reading-frame states, rotational sectors, and cyclic operating modes may be related by symmetry. Once symmetry is admitted, however, there are at least two fundamentally different ways to use it.

First, symmetry may identify a low-dimensional invariant subspace in which each regime is allowed to vary independently. This is an ordinary reduction of the admissible parameter space. Second, the regimes may be required to be exact group transforms of one shared state. This is a relational constraint: it removes the independent continuous parameter vector introduced by an ordinary split and replaces it with a discrete group label.

These possibilities lead to three models:

Full independent model. Each segment receives an unrestricted parameter vector.

Independent fundamental-subspace model. Each segment receives its own parameter vector inside a selected invariant subspace.

Shared exact-orbit model. The segments share one continuous parameter vector and differ only by a relative group action.

The models have continuous-dimension increments

`ΔdA=dfull,  ΔdB=dfund,  ΔdC=0.`

The distinction is not semantic. It determines the leading MDL penalty and the sample length required for reliable detection. Earlier formulations of the present project conflated Models B and C: the prose described a shared orbit, while the proposed penalty corresponded to independently fitted subspace coordinates. The corrected framework treats the hypotheses separately and makes their different complexity laws the subject of the analysis.

The contributions of this paper are fourfold.

First, we derive the known-boundary MDL penalties for the three model classes. Second, we construct a correct Fisher-orthonormal Fourier parameterization for direct cyclic categorical families. Third, we implement the three detectors with consistent known-boundary penalties and detector-specific population gains. Fourth, we perform a large Monte Carlo comparison testing both asymptotic score behavior and calibrated practical power.

The empirical study is designed to verify theoretically specified coefficients rather than infer the model dimensions from unconstrained slope fitting. For the full and fundamental models, the observed coefficients closely track their predicted dimensions. The shared-orbit model displays a much smaller logarithmic dependence, along with finite-sample residuals that are consistent with maximization over the relative shift and singular behavior near the uniform state. The misspecification experiments show that the constrained advantage is conditional on the relevant geometry being approximately correct.

The analysis is explicitly an offline known-boundary model comparison. It is not yet a sequential stopping-time procedure. Unknown-boundary scanning and online detection are discussed separately because they add a location or multiplicity cost to every detector.

2. Statistical setup

2.1 Known-boundary two-segment comparison

Let

`X1,…,XL`

be independent categorical observations with a proposed boundary after `L1` observations, where

`L1+L2=L,  L1L→ρ∈0,1.`

The left and right segments are denoted `XL` and `XR`. The boundary is supplied to the detector. No search over candidate locations is performed in the principal analysis.

For model class `M∈{A,B,C}`, write `H0,M` for the one-regime null and `H1,M` for the two-regime alternative. The detector score is

`SM=ℓ1,Mθ1−ℓ0,Mθ0−penML1,L2,`

where `ℓ0,M` and `ℓ1,M` are maximized log likelihoods in nats. The raw MDL rule declares a change when `SM>0`.

2.2 Full parameter family

Let

`P={pθ:θ∈Θ}`

be a regular parametric family of continuous dimension

`dfull=dimΘ.`

For an `m`-category multinomial model,

`P=Δm−1,  dfull=m−1.`

The notation `dfull` is retained because applications with several categorical blocks can have dimensions larger than `m−1`.

2.3 Cyclic group action

Let

`Cg={e,τ,τ2,…,τg−1}`

act on the family through transformations `Tk:P→P`. In the direct categorical model used for the simulations, the alphabet size equals the group order and `Tk` cyclically permutes the category coordinates. In block models, the group may instead permute phase blocks while leaving the within-block alphabet unchanged.

At a symmetric reference distribution `p0`, the tangent space carries an induced real representation of `Cg`. Let

`Vfund⊆Tp0P`

be the selected fundamental invariant component, with real dimension `dfund`. For a direct cyclic categorical model,

`dfund=1,g=2,2,g≥3.`

3. Three model classes

3.1 Model A: full independent change

Under the null, both segments share one unrestricted parameter:

`H0,A: XL,XR∼pθ.`

Under the alternative,

`H1,A: XL∼pθL,  XR∼pθR,`

where `θL` and `θR` are independently fitted. Hence

`d0,A=dfull,  d1,A=2dfull,`

and

`ΔdA=dfull.`

For a direct `g`-category model, `dfull=g−1`.

3.2 Model B: independent fundamental-subspace change

Let `η∈Rdfund` parameterize a smooth exponential-family chart inside the selected invariant subspace. A convenient form is

`qηj=exp{Zηj}r​exp{Zηr},`

where the columns of `Z` span the fundamental component in logit coordinates.

Under the null,

`H0,B: XL,XR∼qη.`

Under the alternative,

`H1,B: XL∼qηL,  XR∼qηR,`

where `ηL` and `ηR` are independent. Therefore,

`ΔdB=dfund.`

Model B states that both regimes lie in the same invariant subspace. It does not require the right regime to be a group transform of the left.

3.3 Model C: shared exact-orbit transition

Let

`qkη=Tkq0η.`

A simultaneous shift of both labels is observationally redundant, so the left label can be fixed as a gauge. Under the null,

`H0,C: XL,XR∼q0η.`

Under the alternative,

`H1,C: XL∼q0η,  XR∼qrη,`

where the relative shift satisfies

`r∈{1,…,g−1}.`

The same continuous vector `η` appears on both sides. Therefore,

`ΔdC=0.`

The alternative introduces only the discrete relative shift. Under a uniform two-part label code, its cost is `logg−1` nats. For fixed `g`, this is constant in `L`.

Table 1. Three levels of geometric constraint for a known boundary.

Model

Alternative relation between segments

`Δd`

Leading incremental penalty

A. Full independent

Arbitrary separate parameters

`dfull`

`dfull2logL+O1`

B. Independent fundamental

Separate parameters in `Vfund`

`dfund`

`dfund2logL+O1`

C. Shared exact orbit

One shared state plus relative `r∈Cg`

`0`

`logg−1+O1`

4. MDL complexity laws

4.1 Regular codelength expansion

For a regular `d`-dimensional family fitted to `N` observations, BIC, regular Laplace marginal likelihood, and standard parametric-complexity expansions share the leading form [1–5]

`−logpxN=−ℓNθ+d2logN+O1.`

The `O1` term depends on the coding convention, prior, Fisher information, and parameter-space geometry.

4.2 Exact split increment

For a known split, the incremental regular complexity produced by replacing one `d`-dimensional fit with two is

`pensplitd;L1,L2=d2logL1+logL2−logL1+L2.`

When `L1/L→ρ`,

`pensplit=d2logL+d2log{ρ1−ρ}+O1.`

Thus the coefficient of `logL` is `d/2`, while the split fraction affects only the bounded term.

For the three models,

`penA=dfull2logL+O1,`

`penB=dfund2logL+O1,`

and

`penC=logg−1+O1.`

4.3 Unknown boundaries

If the boundary is unknown, a detector must encode or search over approximately `L−1` candidate locations. A simple two-part location code adds

`logL−1=logL+O1`

to every model. The corresponding leading coefficients become

`dfull2+1,  dfund2+1,  1.`

A location cost cannot be applied to only one detector in a known-boundary comparison. Doing so confounds changepoint multiplicity with model dimension.

4.4 Units

In bits, a regular penalty is

`Δd2log2L=Δd2ln2lnL.`

The quantity `Δd/2ln2` is therefore a coefficient multiplying `lnL` after conversion to bits. It is not a fixed rate in bits per observation. The per-observation penalty,

`Δd2Llog2L,`

vanishes as `L→∞`.

5. Cyclic Fourier geometry

5.1 Fisher metric at the uniform distribution

For the direct `g`-category model, let

`p0=1g,…,1g.`

The tangent space is

`Tp0Δg−1=v∈Rg:j=0g−1vj=0.`

At `p0`, the Fisher inner product is

`⟨u,v⟩F=j​ujvjp0,j=g uTv.`

5.2 Fundamental Fourier basis

For `g≥3`, define

`cj=2gcos2πjg,  sj=2gsin2πjg.`

Then

`⟨c,c⟩F=⟨s,s⟩F=1,  ⟨c,s⟩F=0.`

Their span is invariant under cyclic permutation. In coefficient space, a one-step cyclic shift acts as a planar rotation by `2π/g`. For `g=2`, the Fisher-unit basis is

`c=121,−1.`

The implementation uses Cartesian coordinates rather than amplitude-angle coordinates:

`vη=ac+bs,  η=a,b.`

This avoids the unidentifiable angular coordinate at zero amplitude.

5.3 Positivity-preserving family

The direct fundamental family is defined through softmax logits:

`qηj=exp{Zηj}r​exp{Zηr},`

where `Z=gc s` for `g≥3` and `Z=gc` for `g=2`. The scaling makes the derivative at the uniform distribution equal to the Fisher-orthonormal tangent basis. If `Rg` is the fundamental rotation matrix, the family satisfies the equivariance identity

`qRgkη=Tkqη.`

The automated tests verified Fisher orthonormality and numerical equivariance through `C10`.

5.4 Local Jensen–Shannon geometry

For small perturbations `p=p0+v` and `q=p0+Rv`,

`JSDp,q=18∥v−Rv∥F2+O∥v∥F4.`

If `R` is the fundamental cyclic rotation and `∥v∥F=ε`,

`JSDp,q=ε241−cos2πg+Oε4.`

Therefore the leading coefficients are `1/2`, `3/8`, and `1/4` for `g=2,3,4`, respectively. No extra factor of `g` appears when `ε` is measured in Fisher norm.

6. Population gains and detection boundaries

6.1 Detector-specific gain

The relevant signal strength is the expected log-likelihood advantage within the model being fitted. Define

`GM=maxθ∈H1,MElogpθX−maxθ∈H0,MElogpθX.`

The expected score has the form

`ESM≈LGM−penML.`

For the full multinomial model with equal segment sizes, the pooled null is the arithmetic mixture and

`GA=JSDpL,pR.`

For a restricted family, the pooled null is generally a KL projection rather than the arithmetic mixture. Consequently, ordinary JSD should not be used automatically for Models B and C. The implementation computes each population gain by optimizing the corresponding population log likelihood.

6.2 Local detection laws

If

`GM=aMε2+oε2,`

then the regular models have boundaries

`εA2≍dfulllogL2aAL,`

and

`εB2≍dfundlogL2aBL.`

For a regular shared-orbit stratum with fixed `g` and an explicit label code,

`εC2≍logg−1+O1aCL.`

The stronger `1/L` scaling of Model C arises not merely from a smaller tangent space, but from sharing the continuous state across the boundary.

6.3 Singular qualification

At `η=0`, all orbit elements coincide. The relative label is then unidentifiable, and the shared-orbit model is singular. The two-part code still adds no independent continuous parameter vector, but exact Bayesian asymptotics near orbit collapse may contain nonregular corrections [6]. The empirical study therefore treats a zero leading coefficient as a theoretical structural prediction, while allowing finite-sample residual length dependence.

7. Detector implementation

7.1 Full detector

The production rerun used a BIC-scored unrestricted multinomial detector so that all three detectors could be compared through maximized likelihood plus explicit complexity increments. For count vectors `cL` and `cR`, the raw gain is

`ℓpL;cL+ℓpR;cR−ℓp;cL+cR.`

The penalty is the exact known-split increment

`g−12logL1+logL2−logL.`

An exact KT/Dirichlet-`1/2` implementation is also available in the accompanying code, but it was not used for the reported production slopes.

7.2 Fundamental detector

The null fits one coordinate vector `η` to the combined counts. The alternative fits `ηL` and `ηR` separately. Optimization uses L-BFGS-B with analytical gradients in Cartesian Fourier coordinates and a smooth softmax map. The penalty is the exact known-split BIC increment with dimension `dfund`.

7.3 Shared-orbit detector

The null fits one `η` to the combined counts. For each nonidentity shift `r`, the right counts are aligned by `−r`, pooled with the left counts, and fitted with one shared `η`. The alternative chooses the shift with the largest shared-state likelihood. Its penalty is

`logg−1,`

with no location cost and no continuous-dimension increment.

For `g=2`, there is only one nonidentity shift and the label cost is `log1=0`.

7.4 Structural validation

Before the production run, six automated tests were executed. They verified:

Fisher orthonormality of the Fourier basis;

exact cyclic equivariance;

the intended continuous-dimension increments;

equality of Models A and B for `C2` and `C3`, where the fundamental component spans the full nontrivial tangent space;

recovery of planted relative shifts in smoke data;

absence of a sample-length-dependent penalty in Model C.

All six tests passed.

8. Monte Carlo design

8.1 Simulation grid

The production experiment used the design in Table 2.

Table 2. Production Monte Carlo design.

Component

Values

Cyclic groups

`C2,C3,C4,C5,C6`

Effect coordinates

0.08, 0.12, 0.18, 0.25

Segment length per side

100, 200, 400, 800, 1,600, 3,200

Total length `L`

200 through 6,400

Alternative trials

500 per configuration

Null calibration trials

1,000 per configuration

Calibration target

`α=0.05`

Full detector

Multinomial BIC

Configurations

312

Detector-level result rows

936

Simulated two-segment datasets

468,000

Base random seed

20260713

The 312 configurations were run with deterministic configuration-specific seeds. Independent configurations were parallelized over 16 worker processes. Results were checkpointed after every completed configuration.

8.2 Data-generating scenarios

Three scenarios were used.

Exact orbit

The left state was generated from a fundamental coordinate with norm equal to the specified effect. The right state was its one-step cyclic transform. This scenario matches Model C and is also contained in Models A and B.

Independent fundamental change

Both regimes were generated inside the fundamental family but were not generally related by a cyclic shift. For `g≥3`, the right coordinate had radius `0.72` times the left radius and angle 0.713 radians. For `g=2`, the right coordinate was `−0.55` times the left coordinate. This scenario matches Model B but generally violates Model C.

Full-space higher-mode change

For `g≥4`, a mode-2 Fourier component with amplitude `0.85` times the effect was added with opposite signs on the two sides of the boundary. For `g=4`, this is the one-dimensional sign representation. For `g=5,6`, it is a higher two-dimensional Fourier mode. This scenario contains signal outside the fundamental component and is used as a misspecification test.

8.3 Two distinct analyses

The experiment separates asymptotic score validation from practical power comparison.

Raw MDL score analysis

For each matched detector and scenario, the mean uncalibrated score was regressed as

`S=β0+βGLGM+βLlogL+effect fixed effects+ϵ.`

The expected coefficient is `βG=1`. The empirical penalty slope is defined as

`cM=−βL.`

Each group-level regression used 24 design points: four effects by six lengths.

Calibrated power analysis

At every group, scenario, effect, and length, an additive critical value was set to the empirical 95th percentile of 1,000 null scores. Power was then estimated using 500 alternative samples. These calibrated curves support fair practical comparisons at a common nominal false-positive target, but they do not directly expose the raw MDL penalty coefficient.

For each effect, the 50% power crossover was estimated by monotone stabilization of the empirical power curve followed by interpolation on log total length. Median crossover ratios summarize relative sample requirements.

9. Results

9.1 Full-model penalty

For the full-space matched scenario, the observed coefficients were close to the theoretical values `g−1/2`.

Table 3. Raw-score regression estimates for the full model.

`g`

`βG` (SE)

Empirical penalty slope (SE)

Predicted slope

`R2`

4

1.003 (0.003)

1.515 (0.066)

1.500

0.9999

5

1.012 (0.004)

2.119 (0.058)

2.000

0.9999

6

0.990 (0.004)

2.468 (0.050)

2.500

0.9999

The gain coefficients were within approximately 1.2% of one. The observed penalty slopes followed the predicted increase with group order. The `C5` estimate exceeded 2.0 by about 0.119, a finite-sample deviation of roughly two standard errors; the overall pattern nevertheless tracks the full dimension rather than a constant fundamental dimension.


Figure 1. Predicted and observed full-model coefficients. Error bars show one regression standard error.

9.2 Fundamental-subspace penalty

Model B predicts a coefficient of `1/2` for `C2` and 1 for every `Cg` with `g≥3`.

Table 4. Raw-score regression estimates for the fundamental-subspace model.

`g`

`βG` (SE)

Empirical penalty slope (SE)

Predicted slope

`R2`

2

0.987 (0.003)

0.457 (0.050)

0.500

1.0000

3

1.019 (0.009)

1.047 (0.031)

1.000

0.9994

4

0.996 (0.004)

1.001 (0.015)

1.000

0.9998

5

1.015 (0.009)

1.019 (0.032)

1.000

0.9993

6

1.002 (0.008)

0.967 (0.027)

1.000

0.9995

The empirical coefficients remain approximately constant from `C3` through `C6`, even though the full-simplex dimension increases from two to five. This is the clearest direct evidence that an independently fitted fundamental family has a different complexity law from the unrestricted multinomial family.

For `C2` and `C3`, the full and fundamental model spaces coincide because

`g−1=dfund.`

The informative separation begins at `C4`.


Figure 2. Predicted and observed fundamental-subspace coefficients. Error bars show one regression standard error.

9.3 Shared exact-orbit penalty

Model C predicts no leading continuous-dimension term. The raw-score regressions yielded the residual estimates in Table 5.

Table 5. Residual raw-score log-length slopes for the shared-orbit model.

`g`

`βG` (SE)

Residual penalty slope (SE)

Structural prediction

`R2`

2

0.994 (0.003)

-0.093 (0.110)

0

0.9999

3

1.005 (0.002)

0.144 (0.064)

0

1.0000

4

1.009 (0.004)

0.160 (0.062)

0

0.9999

5

0.996 (0.004)

0.127 (0.050)

0

0.9999

6

1.008 (0.006)

0.206 (0.047)

0

0.9998

These residual slopes are far below the fundamental and full coefficients, but they are not uniformly zero. Direct raw-threshold crossover regressions, available where at least three effects crossed 50% power inside the grid, gave slopes 0.037 for `C5` and 0.062 for `C6`. Their `R2` values were modest because only three or four crossover points were available.

The evidence therefore supports the qualified statement:

The shared exact-orbit detector has a near-zero leading logarithmic coefficient relative to the regular split models, while finite-sample score behavior retains small group-dependent drift.

The present simulation does not establish exact finite-sample constancy. Plausible sources of the residual include maximization over nonidentity shifts, finite-sample MLE bias, and the singular orbit-collapse point at the uniform distribution.


Figure 3. Residual log-length coefficients for the shared exact-orbit detector. The structural prediction is zero; observed values show small finite-sample drift.

9.4 Calibrated crossover advantage

The common 5% calibration permits direct comparison of sample length at equal nominal false-positive control. Table 6 reports median 50%-power crossover ratios. A ratio below one favors the numerator.

Table 6. Median calibrated crossover-length ratios.

`g`

Shared / full

Shared / fundamental

Fundamental / full

2

0.883

0.883

1.000

3

0.943

0.943

1.000

4

0.705

0.812

0.896

5

0.672

0.829

0.771

6

0.610

0.814

0.793

For exact-orbit data, the shared detector required approximately 30%, 33%, and 39% fewer observations than the full detector for `C4,C5,C6`. Relative to the fundamental detector, the reduction was approximately 17%–19% over the same groups.

For independent fundamental changes, Model B and Model A were identical at `C2` and `C3`. At `C4,C5,C6`, the fundamental detector required approximately 10%, 23%, and 21% fewer observations than the full detector among effects with internal crossover estimates.


Figure 4. Calibrated sample-length advantage on exact-orbit data. The dashed line denotes equal sample requirements.

9.5 Misspecification

The constrained detectors should not dominate when their geometric assumptions are false. The higher-mode scenario tests this requirement.

Table 7. Mean calibrated power for higher-mode full-space changes at total length 6,400.

`g`

Full

Fundamental

Shared orbit

4

0.999

0.399

0.082

5

0.975

0.211

0.044

6

0.962

0.182

0.055

The full detector retained high power because it could represent the added mode. The fundamental detector captured only the component lying in its selected subspace, while the exact-orbit detector was more severely misspecified. This result is important: the sample-efficiency advantage of the constrained detectors is not a consequence of generally lower thresholds. It is conditional on structural correctness.


Figure 5. Calibrated power under a higher-mode change outside the fundamental subspace.

9.6 Relative-shift recovery

On exact-orbit data at total length 6,400, the shared detector also recovered the planted one-step shift with high accuracy.

Table 8. Relative-shift recovery at total length 6,400.

`g`

Mean accuracy

Minimum across effects

Maximum across effects

2

1.000

1.000

1.000

3

1.000

1.000

1.000

4

0.997

0.988

1.000

5

0.9945

0.980

1.000

6

0.989

0.958

1.000

The slight decline with `g` is expected because the detector maximizes over more candidate shifts, while adjacent orbit states become geometrically closer for some directions as group order increases.

10. Interpretation

10.1 The three penalties are empirically distinct

The central empirical result is not a single universal cyclic threshold. It is the separation of three model-dependent complexity laws.

The full detector follows the growing simplex dimension. The fundamental detector follows the fixed dimension of the selected representation. The shared-orbit detector is qualitatively different because the alternative does not introduce a second continuous state.

For `g≥4`,

`g−12>dfund2=1>0.`

The Monte Carlo slopes display the same ordering.

10.2 Dimension reduction versus parameter sharing

Model B gains efficiency through dimension reduction. It remains a conventional split model: the left and right coordinates are independently estimated.

Model C gains efficiency through parameter sharing. Its hypothesis is relational. If the relation is exact, the detector can pool aligned observations from both segments to estimate one state. This is stronger than merely projecting both segments into a smaller space.

The distinction suggests a hierarchy of prior structural knowledge:

`arbitrary change⊃subspace change⊃exact orbit transition.`

Greater structural commitment reduces codelength only when the data support that commitment.

10.3 Why the shared result is qualified

The structural dimension increment of Model C is exactly zero, but a zero dimension increment does not imply that every finite-sample score statistic is independent of `L`. The selected shift is unidentifiable at orbit collapse, the maximum is taken over several discrete alternatives, and likelihood estimators carry finite-sample bias. The observed residual coefficients should therefore be treated as empirical corrections rather than reinterpreted as a two-dimensional regular penalty.

A future exact marginal-likelihood analysis should determine whether these corrections approach a constant, a smaller logarithmic term arising from singular learning, or another slowly varying form in local alternatives.

11. Codon-phase application

11.1 Separating group order from alphabet size

A codon-phase model has

`g=3`

reading-frame phases and

`m=4`

nucleotide symbols. These quantities are not interchangeable.

A phase-specific nucleotide model contains

`p0,p1,p2∈Δ3,`

so the full parameter space is

`Δ33`

with dimension

`dfull=34−1=9.`

The group `C3` acts by permuting the three phase blocks.

11.2 Fundamental isotypic dimension

The phase representation decomposes as

`Rphase3=Vtriv⊕Vfund,`

where `dimVfund=2`. Each phase coordinate carries a three-dimensional nucleotide-contrast tangent space. The nontrivial phase isotypic component therefore has dimension

`dphase-fund=2m−1=6.`

The corresponding known-boundary penalties are

`penfull=92logL+O1,`

`penphase-fund=3logL+O1,`

and, for a shared exact phase orbit,

`penorbit=log2+O1.`

A two-dimensional codon-phase model would require an additional biological restriction, such as a single fixed nucleotide-contrast direction shared across phases.

11.3 Proposed experiment

A corrected frameshift experiment would:

estimate phase-specific nucleotide or codon-emission parameters from verified coding sequences;

generate or curate sequences containing known `+1` and `−1` frameshifts;

compare full, phase-fundamental, and shared-phase-orbit detectors at known boundaries;

repeat with hidden boundaries using the same location treatment for all models;

measure false-positive rate, detection probability, localization error, and sample length.

Frameshift methods based on phase-specific emissions and hidden-state structure are established [9]. The contribution proposed here is not frameshift detection itself, but a representation-theoretic MDL decomposition of three different levels of phase constraint.

12. Relation to the broader geometric framework

The broader Mayo manifold program distinguishes three conceptual levels.

Topology. A flat bundle can carry nontrivial global holonomy despite vanishing local curvature.

Dynamics. A twisted or constrained generator can exhibit a spectral response to the global sector.

Information. A statistical observer faces model-dependent description lengths determined by parameter dimension, invariant subspaces, and shared group relations.

The present paper establishes a result only at the third level. It does not prove that topological holonomy, dynamical spectral gaps, and MDL coefficients are numerically identical.

A motivating numerical coincidence remains: the East model contains an inverse-gap asymptotic involving `1/2ln2`, while a one-dimensional regular BIC increment has coefficient `1/2ln2` when expressed in bits per `lnL`. This resemblance is suggestive but not explanatory. A genuine bridge would need to derive the observer’s statistical family from the dynamics and identify the operating regime in which relaxation and detection asymptotics are comparable.

13. Limitations

Several limitations define the scope of the current result.

First, the observations are independent and categorical. Markov, hidden-state, and continuous-observation extensions may change both the population gain and effective sample size.

Second, the boundary is known. A scanning or online procedure requires multiplicity correction, a stopping rule, false-alarm control over time, and detection-delay analysis.

Third, the direct simulation identifies group order with alphabet size. The codon example requires a block-family implementation in which the group acts on phases while each phase contains a separate observation alphabet.

Fourth, the shared-orbit family is singular at the uniform state. The present two-part code establishes the absence of an added continuous vector but does not derive the exact singular marginal-likelihood expansion.

Fifth, the calibrated critical values were estimated separately for every configuration. This is appropriate for practical power comparison but means that calibrated crossover curves should not be used to infer the raw MDL coefficient.

Sixth, only one fundamental Fourier mode and one higher-mode misspecification were examined. Non-Abelian groups, representation multiplicities, stabilizers, and approximate orbit relations remain open.

Finally, the novelty claim should be kept narrow. Symmetry reduction, MDL dimension penalties, and constrained changepoint models are established ideas individually. The contribution here is their explicit three-way organization into full independent, independent invariant-subspace, and shared exact-orbit hypotheses, together with a matched empirical comparison.

14. Extensions

14.1 Approximate orbit model

A natural interpolation between Models B and C is

`ηR=RgrηL+δ,`

with a shrinkage prior or code on `δ`. Setting `δ=0` gives the exact orbit; allowing unrestricted `δ` recovers the independent subspace model. This would quantify how much deviation from exact symmetry can be tolerated before the relational advantage disappears.

14.2 Stabilizer-adaptive labels

If `η` is invariant under a nontrivial subgroup, the orbit contains fewer than `g` distinct states. The effective label cost should depend on

`Cg⋅η=gStabη.`

A stabilizer-adaptive code could improve finite-sample behavior near symmetric strata.

14.3 Exact universal codes

The unrestricted multinomial family admits the KT mixture. Comparable exact or numerically integrated universal codes for the fundamental and shared-orbit families would eliminate reliance on regular BIC constants and permit direct finite-sample regret comparisons.

14.4 Sequential extension

An online extension would replace the fixed split with a stopping time. The correct questions would then include average run length under the null, expected detection delay under each geometric alternative, and the cost of repeatedly searching over both boundary locations and relative group elements.

15. Conclusion

Cyclic structure constrains a regime-change problem in two different ways. It can restrict each regime independently to a low-dimensional invariant subspace, or it can require the regimes to be exact transforms of one shared state.

For a known boundary, these produce three leading MDL structures:

`penAL=dfull2logL+O1,`

`penBL=dfund2logL+O1,`

and

`penCL=logg−1+O1,`

with the final expression understood as a two-part structural code on a regular orbit stratum.

The Monte Carlo results strongly validate the full and independent-fundamental coefficients. They also support a near-zero leading coefficient for the shared exact-orbit detector, while revealing finite-sample residuals that require a singular or exact coding analysis. Under matched geometry, the constrained detectors reduce required sample length. Under higher-mode misspecification, the advantage disappears and the full detector dominates.

The essential conclusion is therefore not that cyclic groups possess one universal detection threshold. It is that different levels of geometric knowledge define different statistical questions:

`Do the distributions differ arbitrarily?`

`Do they differ within a selected representation?`

`Is one an exact symmetry transform of the other?`

Their answers carry different complexity increments. That distinction is the correct foundation for future work on cyclic changepoints, phase transitions, biological reading frames, and symmetry-aware sequential detection.

Code and data availability

The implementation and complete rerun outputs accompanying this manuscript are:

mayo_cyclic_detectors_v3.py — detector families, fitting routines, population gains, and structural tests;

v3_simulation.py — scenario construction, calibration, and Monte Carlo evaluation;

run_v3_full_parallel.py — deterministic parallel production runner;

test_mayo_cyclic_detectors_v3.py — automated structural tests;

v3_full_results.csv — all 936 detector-level result rows;

v3_score_regression_summary.csv — raw-score regression coefficients;

v3_crossover_estimates.csv — estimated 50%-power crossover lengths;

v3_crossover_ratio_summary.csv — median practical sample-length ratios.

The production run used base seed 20260713.

References

Schwarz, G. (1978). Estimating the dimension of a model. The Annals of Statistics, 6(2), 461–464.

Krichevsky, R. E., and Trofimov, V. K. (1981). The performance of universal encoding. IEEE Transactions on Information Theory, 27(2), 199–207.

Barron, A. R., Rissanen, J., and Yu, B. (1998). The minimum description length principle in coding and modeling. IEEE Transactions on Information Theory, 44(6), 2743–2760.

Grünwald, P. D. (2007). The Minimum Description Length Principle. MIT Press.

Rissanen, J. (2007). Information and Complexity in Statistical Modeling. Springer.

Watanabe, S. (2009). Algebraic Geometry and Statistical Learning Theory. Cambridge University Press.

Amari, S., and Nagaoka, H. (2000). Methods of Information Geometry. American Mathematical Society and Oxford University Press.

Serre, J.-P. (1977). Linear Representations of Finite Groups. Springer.

Rho, M., Tang, H., and Ye, Y. (2010). FragGeneScan: predicting genes in short and error-prone reads. Nucleic Acids Research, 38(20), e191.

Eaton, M. L. (1989). Group Invariance Applications in Statistics. Institute of Mathematical Statistics.

Csiszár, I., and Shields, P. C. (2004). Information theory and statistics: a tutorial. Foundations and Trends in Communications and Information Theory, 1(4), 417–528.

Cancrini, N., Martinelli, F., Roberto, C., and Toninelli, C. (2008). Kinetically constrained spin models. Probability Theory and Related Fields, 140, 459–504.

Appendix A. Reproducibility details

For each configuration, null samples were drawn from a scenario-specific no-change distribution. The empirical 95th percentile of the three detector scores was computed using the higher quantile convention. Alternative samples were then scored against those fixed thresholds. The raw score, explicit penalty, population gain, zero-threshold detection indicator, calibrated detection indicator, and selected orbit shift were recorded for every detector.

Population gains were computed at the distribution level. The full gain was weighted JSD. The fundamental gain was the difference between the best two-segment and pooled expected log likelihoods inside the fundamental family. The shared-orbit gain maximized the aligned shared-state expected likelihood over nonidentity shifts.

The crossover estimator first replaced the empirical power sequence by its cumulative maximum to suppress Monte Carlo reversals. When power crossed 0.5 within the tested grid, interpolation was performed linearly in log total length. Crossovers below or above the grid were flagged and excluded from slope regressions requiring internal estimates.

Appendix B. Reading the empirical slopes

The raw-score regression is the primary test of the theoretical penalty coefficients. Suppose

`S=LG−clogL+b+ϵ.`

Regressing the mean score on `LG`, `logL`, and effect indicators should recover coefficient one on `LG` and coefficient `−c` on `logL`. The extremely high `R2` values indicate that this affine approximation describes the simulated mean scores over the tested grid.

The calibrated-power crossover slopes answer a different question. Configuration-specific null calibration adds an empirical threshold that may itself vary with length. Those crossovers are appropriate for comparing practical sample requirements at a common false-positive target, but their slopes are not expected to equal the raw MDL coefficient.

Acknowledgment

This version was developed after methodological review identified that an unrestricted simplex model, an independently fitted invariant-subspace model, and a shared exact-orbit model had previously been treated as though they implied the same complexity law. Distinguishing those hypotheses led to the present framework and simulation design.


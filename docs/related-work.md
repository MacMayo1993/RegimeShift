# Related work

A map of the literatures this method sits between, and a precise statement of
where it differs.

> **Provenance.** This list was compiled from the reference list of the
> manuscript itself and from two external prior-art reviews commissioned on it.
> The entries have **not** been independently re-verified against the published
> record by this repository. Check each citation before quoting it in a
> submission — treat this file as a map, not as a bibliography.

## Where this method differs

Stated as narrowly as the evidence supports:

> The individual ingredients are standard. The contribution is combining them
> into a known-boundary MDL comparison that separates *independent
> invariant-subspace change* from *shared-orbit transition*. In that setting the
> shared-orbit alternative introduces no new continuous parameter across the
> split, and therefore carries no leading continuous-dimension penalty — it pays
> only a discrete shift-label cost.

Neither external review located a paper combining all four of: categorical
changepoints, a finite cyclic group acting directly on the statistical family,
an explicit contrast between independently fitted invariant-subspace models and
a shared exact orbit, and a matched MDL penalty comparison. That is a
search-based inference about absence, not proof of it, and the novelty claim
should be worded accordingly.

## 1. MDL and BIC as model-selection foundations

The leading `(d/2) log n` term that Models A and B rely on, and the framing of
MDL as a coding principle rather than a fixed penalty formula.

| Work | Relevance |
|---|---|
| Schwarz (1978), *Estimating the Dimension of a Model*, Annals of Statistics | The regular-model dimension penalty. Direct foundation for Models A and B; contains no changepoint or group structure. |
| Barron, Rissanen & Yu (1998), IEEE Trans. Inf. Theory | Universal coding, stochastic complexity, NML, mixture and predictive codes. Relevant because this work moves between BIC, Laplace and two-part reasoning without fixing one code. |
| Grünwald (2007), *The Minimum Description Length Principle*, MIT Press | Book-length treatment. |
| Rissanen (2007), *Information and Complexity in Statistical Modeling*, Springer | Book-length treatment. |
| Grünwald & Roos (2019), *Minimum Description Length Revisited*, arXiv:1908.08484 | Modern survey; useful for stating what "MDL" does and does not commit you to. |
| Krichevsky & Trofimov (1981), IEEE Trans. Inf. Theory | The KT mixture for multinomials — the exact universal code this implementation does *not* use. |

## 2. Changepoint detection with information criteria

Model selection for segment counts and locations. The closest antecedents for
the scoring framework, all in the unknown-boundary setting this work brackets.

| Work | Relevance |
|---|---|
| Yao (1988), Statistics & Probability Letters | Schwarz criterion for the number of changes; consistency. Establishes BIC-type changepoint selection decades earlier, without categorical or symmetry constraints. |
| Davis, Lee & Rodriguez-Yam (2006), JASA | MDL selection of segment counts, orders and locations for piecewise autoregressions. Close on MDL segmentation; different observation model, unknown boundaries, no orbit relation. |
| Niu, Hao & Zhang (2016), *Multiple Change-Point Detection: A Selective Overview*, Statistical Science | General survey. |
| Recent selective review of information criteria in multiple changepoint detection (PMC10813938) | Survey of the IC-based branch specifically. |

## 3. Categorical and compositional changepoints

The observation model closest to this one.

| Work | Relevance |
|---|---|
| Wang, Zou & Yin (2018), Annals of Statistics | High-dimensional multinomial homogeneity, penalised segmentation, consistency. The closest categorical antecedent; tests unrestricted probability change rather than invariant-subspace or orbit alternatives. |
| Truong & Runge (2024), *Stat* | Exact segmentation for simplex-valued and one-hot categorical observations. Close in geometry; its penalty counts changes and encodes no cyclic relational structure. |
| Li, Lu & Wang (2026), arXiv | Marginalised transition model for periodic, serially correlated categorical data. More realistic than the i.i.d. model here, but periodicity is temporal rather than an exact finite-group orbit shared across two segments. |

## 4. Group invariance and equivariant testing

The literature nearest to Model C's relational hypothesis.

| Work | Relevance |
|---|---|
| Pérez-Ortiz, Lardy, de Heide & Grünwald (2024), Annals of Statistics | Testing between group models via maximal invariants, e-statistics, Haar-prior Bayes factors. The closest modern group-testing work; does not formulate a changepoint as independent-subspace versus shared-orbit, nor derive the three-way MDL hierarchy. |
| Eaton (1989), *Group Invariance Applications in Statistics*, IMS | Standard reference for invariance arguments. |
| Serre (1977), *Linear Representations of Finite Groups*, Springer | Standard reference for the representation theory used in Section 5. |
| Amari & Nagaoka (2000), *Methods of Information Geometry* | The Fisher geometry underlying the orthonormal Fourier basis. |

## 5. Dimension-reduced changepoints

The nearest analogue to Model B's efficiency argument.

| Work | Relevance |
|---|---|
| Yu, Zhao, Huang, Zhu & Zhu (2026), J. Multivariate Analysis | Sparse dimension-reduced subspace for changepoint detection. Close to Model B's motivation, but the subspace is *learned* from sparse structure rather than fixed by a representation, and there is no exact-orbit model. |

## 6. Singular learning theory

Directly relevant to the orbit-collapse point, where the relative shift is
unidentifiable and ordinary BIC reasoning does not automatically apply.

| Work | Relevance |
|---|---|
| Watanabe (2009), *Algebraic Geometry and Statistical Learning Theory* | RLCT-based asymptotics for singular models. |
| Watanabe (2013), *A Widely Applicable Bayesian Information Criterion*, JMLR | WBIC. |
| Drton & Plummer (2017), JRSS B | A BIC for singular models. |

Together these are the reason the shared-orbit claim in this repository is
stated as *no leading continuous-dimension penalty under the known-boundary
description-length accounting*, rather than as an exact finite-sample constant.
See `paper-notes.md` for the empirical counterpart: the residual log-length
slope decomposes exactly as `-s`, the raw gain's departure from `n G`.

## 7. Application context

| Work | Relevance |
|---|---|
| Rho, Tang & Ye (2010), *FragGeneScan*, Nucleic Acids Research | Established frameshift/phase methods. The proposed contribution in Section 11 is not frameshift detection but a representation-theoretic MDL decomposition of phase constraint. |
| Csiszár & Shields (2004), *Information Theory and Statistics: A Tutorial* | Background. |
| Cancrini, Martinelli, Roberto & Toninelli (2008), PTRF | Kinetically constrained spin models; the source of the East-model remark in Section 12, which the manuscript correctly presents as a numerical coincidence rather than a derivation. |

## What a reviewer will ask next

Both external reviews converged on the same open questions, none of which this
repository answers:

1. Is the zero continuous-dimension increment justified under an *exact* code
   rather than a BIC-style approximation?
2. What happens at the singular orbit-collapse stratum — is there an RLCT or
   multiplicity correction?
3. Does unknown-boundary scanning change the leading comparison, given the
   location cost applies symmetrically to all three detectors?
4. Would KT, Dirichlet-mixture or NML codes materially change the ordering?
5. Which application makes the cyclic-orbit assumption compelling enough to
   test empirically?

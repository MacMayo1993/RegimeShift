# Reproduction notes

What this repository reproduces from
`Geometric_Complexity_Cyclic_Regime_Changes_v3_1.docx`, and where it had to make
its own choices.

## Recovered exactly from the manuscript

* The three model classes, their continuous-dimension increments, and the
  known-split penalty `(d/2) * log(n_L n_R / n)` (Sections 3–4).
* `d_fund = 1` for `m = 2` and `2` for `m >= 3`; `d_full = m - 1`.
* The Fisher-orthonormal Cartesian Fourier parameterisation, chosen over
  amplitude–angle coordinates to avoid the unidentifiable angle at zero
  amplitude (Section 5.2).
* The production design of Table 2: groups 2–6, effects
  `0.08, 0.12, 0.18, 0.25`, per-side lengths `100 … 3200`, 500 alternative and
  1,000 null trials per configuration, 5% calibration, 312 configurations,
  936 detector rows, 468,000 datasets, base seed 20260713. The test
  `test_production_grid_matches_the_manuscript_design` pins all of these.
* The two-analysis split of Section 8.3 and the crossover estimator of
  Appendix A (cumulative-maximum stabilisation, then linear interpolation in
  log total length, with out-of-grid crossovers flagged and excluded).

## Constants taken from the manuscript

| Quantity | Value | Source |
|---|---|---|
| Model C label cost | `log(m - 1)` nats | §3.3: "Under a uniform two-part label code, its cost is log(g − 1) nats"; §7.3 confirms `log 1 = 0` at `g = 2`. |
| Independent-fundamental radius ratio | `0.72` | §8.2 |
| Independent-fundamental angle | `0.713` rad | §8.2 |
| Independent-fundamental ratio at `m = 2` | `-0.55` | §8.2 |
| Higher-mode amplitude | `0.85 x` effect | §8.2 |

All five are quoted directly from the document, and their provenance is also
exported machine-readably as `regimeshift.MANUSCRIPT_CONSTANTS`. A test asserts
every entry matches the live module constant, so the table cannot drift.

### A correction worth recording

Earlier versions of this file claimed the manuscript "renders its equations as
images", and three of these constants were *guessed* on that basis — `0.85`,
`-0.6` and `0.8` in place of the true `0.72`, `-0.55` and `0.85`.

That premise was wrong. The document contains no equation images at all: its 400
displayed and inline expressions are Office Math (OMML), and the five embedded
PNGs are Figures 1–5. The first extraction pass read only `w:t` elements, which
silently drops every `m:t` inside an `<m:oMath>` node — so the equations
vanished from the extracted text and their absence was misread as evidence they
were pictures.

The lesson generalises beyond this repository: an extractor that drops content
silently is worse than one that fails, because the gap gets rationalised. The
extraction now includes math inline, and `docs/extracted-text.md` says which
parts of the document it covers.

The formulas recovered alongside the constants all confirm the implementation
that had been written without them — the known-split penalty
`(d/2)[log L1 + log L2 − log(L1+L2)]`, its `(d/2) log(ρ(1−ρ))` split-fraction
term, and the `log(g − 1)` label cost.

## Deliberate deviations

* **`higher_mode` requires `m >= 4`.** At `m = 2` the mode-2 component is
  trivial, and at `m = 3` mode 2 is the conjugate of mode 1 — it lies *inside*
  the fundamental component, so it would not be a misspecification at all.
  The manuscript only reports this scenario for `m = 4, 5, 6`. Requesting it
  below `m = 4` raises rather than silently generating in-subspace data.
* **The higher-mode change is the mode-2 flip alone.** Section 8.2 says a
  mode-2 component "was added with opposite signs on the two sides of the
  boundary". This implementation initially read that as *rotate the fundamental
  coordinate and also add the mode*, which left the change carrying a
  full-strength exact-orbit component. Model C then retained 51–56% of the
  population gain and stayed competitive under "misspecification", flatly
  contradicting Table 7.

  Confining the change to the higher mode — both segments share one fundamental
  coordinate, which only sets the operating point away from uniform —
  reproduces Table 7's structure:

  | | full | fundamental | shared orbit |
  |---|---|---|---|
  | population gain, `m = 6`, effect 0.25 | 0.01005 | 0.00030 | **-0.00768** |
  | share of full gain | 100% | 3.0% | negative |

  The fundamental family keeps a few percent of the signal and the shared-orbit
  gain goes *negative* — aligned pooling is worse than not aligning — which is
  what produces Table 7's below-nominal shared-orbit power of 0.044–0.082. The
  worked example shows the same reversal at calibrated 5%. Pinned by
  `test_higher_mode_reproduces_the_manuscript_misspecification_pattern`.
* **The full detector is BIC-scored**, matching Section 7.1's production rerun.
  The exact KT/Dirichlet mixture code is not implemented; Section 14.3 lists
  exact universal codes for all three families as future work.

## Properties the code has that the manuscript does not state

Two consequences of the Model C construction are worth naming, because they
look like bugs until you see why they are not:

1. **Model C's raw gain can be negative.** Its alternative maximises over
   *nonidentity* shifts only, so unlike Models A and B it does not nest its own
   null. On no-change data the best aligned pooling is typically worse than the
   unaligned pooling. This is exactly why the model pays only a discrete label
   cost.
2. **Model C's raw zero-threshold rule is not conservative.** A penalty of
   `log(m - 1)` does not grow with `n`, so unlike Models A and B it provides no
   increasing protection under the null — and at `m = 2` it is exactly zero,
   leaving no protection at all. In the validation grid the shared-orbit
   detector's zero-threshold null rate reaches ~0.29 at `m = 2` with short
   segments and weak effects, against ~0.02 for the regular detectors. This is
   the same singular qualification as Section 6.3 (the relative label is
   unidentifiable at orbit collapse) and it is precisely why the study compares
   detectors at a *common calibrated 5%* rather than by raw rule.
3. **Model C can still win under a mild departure from exact symmetry.** In the
   independent-fundamental scenario its hypothesis is strictly false and its
   population gain is strictly smaller than Model B's — but at finite samples
   its constant label cost is cheap enough that it still detects at least as
   well at `m = 6`. That is the approximate-orbit regime of Section 14.1: the
   relational advantage degrades gradually. The honest counterweight is the
   higher-mode scenario, where the deviation leaves the fundamental subspace
   entirely and Model C collapses to near-chance power, exactly as in Table 7.
4. **Model A's raw gain does not dominate Model B's.** The unrestricted
   alternative is larger, but so is the unrestricted *null*, and the gain is a
   difference of the two. What is guaranteed is the likelihood nesting
   (`alt_A >= alt_B >= alt_C`, `null_A >= null_B = null_C`) and, because B and C
   share a null, `gain_B >= gain_C`. Both are asserted in
   `tests/test_detectors.py`.

## K* = 1/(2 ln 2) as the penalty quantum

Section 4.4 converts the regular penalty to bits: `(Δd/2) log2 L = Δd/(2 ln 2) ln L`,
and Section 12 notes that `1/(2 ln 2)` also appears in an East-model inverse-gap
asymptotic, calling the resemblance suggestive rather than explanatory.

Naming that constant `K*` makes the paper's hierarchy a counting statement.
Every leading coefficient in the framework is an integer multiple of it:

| model | d | nats per ln n | bits per ln n |
|---|---|---|---|
| A. full | `m - 1` | `(m-1)/2` | `(m-1) K*` |
| B. fundamental | `d_fund` | `d_fund/2` | `d_fund K*` |
| C. shared orbit | 0 | 0 | **0** |

So "Model C has no leading continuous-dimension penalty" becomes "Model C pays
zero `K*` per e-fold". `regimeshift.K_STAR` exports the constant, `predicted_slope`
takes a `units` argument, and the reports accept `--units bits`, which adds a
`k_star_multiple` column.

**A caution about Section 12.** `K*` arising as the per-dimension BIC rate in
bits is *definitional*, not a discovery: it is Schwarz's one-half expressed in
base 2, and any quantity counting half a parameter per e-fold in bits produces
it. The East-model resemblance is therefore only meaningful if the dynamical
occurrence is not likewise a units artefact — if `1/(2 ln 2)` enters there from
the structure of the constrained generator rather than from a choice of base.
That is a sharper form of the manuscript's own caveat, and it is the question a
bridge between the levels would have to answer first.

## Analysis-integrity work beyond the manuscript

An external methodological review raised five concerns about the *analysis*
rather than the models. All five are now addressed in code, and the fixes are
pinned by `tests/test_analysis_integrity.py` plus Monte Carlo counterparts in
the slow suite.

**Weighted regression.** Each group-level regression has only
`n_effects x n_lengths` aggregate design points, and those means have unequal
Monte Carlo precision, so an unweighted fit can misstate both the slope and its
standard error. `score_regression(..., weighted=True)` uses inverse
`sd_score^2 / n_alt` weights, and the report now carries OLS and WLS side by
side plus their disagreement (`weighting_shift`). On the validation grid the
constrained detectors move by 0.02–0.05 — the reported slopes are not artifacts
of the specification. The full detector at `m = 4` moves more, from 1.23 to
1.46 against a prediction of 1.5, i.e. weighting *improves* agreement with
theory. Design condition numbers (130–700) are now reported too.

**The penalty slope is decomposable, exactly.** Because this implementation
subtracts a penalty it computes exactly, the gain and penalty coefficients need
not be estimated jointly. Writing the mean raw gain as `n G + a + s log n`,

    penalty_slope = d/2 - s

holds as an algebraic identity, where `s` (`gain_residual_regression`) is the
only empirical quantity in it. This is worth stating plainly: **for Model C,
`d = 0`, so the residual slopes of Table 5 are `-s` and nothing else.** They are
a property of the maximised likelihood gain — shift maximisation and
finite-sample MLE bias, exactly the causes Section 9.3 hypothesises — and not
evidence of a hidden continuous-dimension penalty. The identity is asserted to
1e-8 on both synthetic and simulated output.

**Crossover uncertainty.** `crossover_bootstrap` and
`crossover_ratio_bootstrap` resample the whole pipeline — power draw,
cumulative-maximum stabilisation, interpolation, median across effects. Two
limits are documented rather than hidden: the resampling covers binomial power
noise but not the variability of the empirical 95th-percentile critical value,
and detectors are resampled independently although they score the same datasets
and are positively correlated, which makes the ratio intervals conservative.

**Optimiser audit.** `fit_fundamental` now checks convergence and
`run_config` reports `optimizer_failures` per configuration, asserted to be zero
across the grid. Note the subtlety this exposed: judging convergence by
`OptimizeResult.success` produces ~2.5% false alarms, because under our tight
`ftol` L-BFGS-B reports `ABNORMAL` whenever its line search cannot improve at
machine precision — which happens *at* the optimum, with observed gradient norms
around 5e-8. Convergence is therefore judged on the first-order condition,
relative to `n` since the gradient is `B^T (counts - n p)`.

The audit also surfaced a second, separate property. When a category has zero
count — routine on short segments — the fundamental MLE does not exist: the
likelihood rises toward the simplex boundary and is asymptotically flat along
that direction. Different optimiser starts then halt at very different
coordinates (observed `|theta|` of 8 versus 15) **while agreeing on the
log-likelihood to six decimals**. Since every detector consumes only
likelihoods, scores are unaffected and start-independent, and both facts are
pinned by tests. It does mean a fitted *coordinate* from a short segment should
not be interpreted, and it is the reason the fit retains two starts even though
the log-likelihood is concave in theta.

**Split fraction.** `Config.split_fraction` makes rho a design dimension
(`segment_length` is the left segment; the default 0.5 reproduces the balanced
design and leaves every existing seed and checkpoint unchanged). This tests a
prediction the manuscript makes but never checked: rho shifts only the bounded
term `(d/2) log(rho (1 - rho))`, leaving the `log n` coefficient at `d/2`. Both
the closed-form and Monte Carlo versions are asserted.

## Model D: the approximate-orbit model (Section 14.1)

Implemented beyond the manuscript's production comparison, because Section 14.1
poses a question the three-model study cannot answer: how far from an exact
orbit a change can drift before the relational code stops paying.

The manuscript gives the form `eta_R = R_g^r eta_L + delta` "with a shrinkage
prior or code on `delta`" but does not fix the code. This implementation takes
a Gaussian prior `delta ~ N(0, tau^2 I)`; the Laplace approximation then gives
a deviation cost of `(d/2) log(1 + n_R tau^2)`, since the Fisher-orthonormal
coordinates make the right segment's information `n_R` per unit per direction.

Two properties make the nesting exact rather than approximate:

* `tau = 0` pins `delta` out and maximises over the same nonidentity shifts, so
  Model D **is** Model C — same gain, same penalty, same selected shift, to
  1e-9. Asserted directly.
* Large `tau` leaves `delta` free, and the maximised gain matches Model B's to
  1e-4, since the alternative can then reach any pair of coordinates.

**The caveat that belongs with it.** For any *fixed* `tau > 0` the leading
coefficient is `d/2` — Model B's rate, not an intermediate one. The
interpolation lives entirely in the bounded term, i.e. at finite `n`. That is
not a defect for the question being asked, since tolerance to imperfect symmetry
is a finite-sample question, but a genuine interpolation of the *leading*
coefficient would require `tau` shrinking with `n`, and this implementation does
not do that.

The `approximate_orbit` scenario supplies matching data, displacing the right
state perpendicular to the rotated state so the deviation is a departure from
the orbit rather than a rescaling along it. Sweeping it at `m = 6` shows three
regimes: the rigid code wins out to a deviation of about 0.25 effects, Model D
wins from roughly 0.5 to 1.0 by beating *both* endpoints, and beyond about 1.5
the relation is not worth encoding at all. See the README table.

Model D is deliberately **not** wired into `run_all_detectors` or the production
grid, which reproduce the manuscript's three-model comparison. It is available
through the API and covered by `tests/test_approximate_orbit.py`.

## The committed production run

`results/v3-production/` holds one complete run of the Table 2 design — 312
configurations, 936 detector rows, 468,000 datasets, base seed 20260713, 17.4
minutes on 4 workers — with a `run_manifest.json` carrying the commit,
environment and per-file checksums.

Its numbers reproduce the manuscript's, and in most cells land closer to theory:

| detector | m | this run | (WLS) | manuscript | predicted |
|---|---:|---:|---:|---:|---:|
| full | 4 | 1.488 | 1.504 | 1.515 | 1.5 |
| full | 5 | **1.967** | 1.994 | *2.119* | 2.0 |
| full | 6 | 2.488 | 2.518 | 2.468 | 2.5 |
| fundamental | 2 | 0.522 | 0.474 | 0.457 | 0.5 |
| fundamental | 3 | 0.996 | 0.979 | 1.047 | 1.0 |
| fundamental | 4 | 0.971 | 0.961 | 1.001 | 1.0 |
| fundamental | 5 | 0.964 | 0.988 | 1.019 | 1.0 |
| fundamental | 6 | 1.036 | 1.017 | 0.967 | 1.0 |
| shared orbit | 2 | -0.228 | -0.074 | -0.093 | 0 |
| shared orbit | 3 | 0.041 | 0.090 | 0.144 | 0 |
| shared orbit | 4 | 0.083 | 0.165 | 0.160 | 0 |
| shared orbit | 5 | 0.139 | 0.142 | 0.127 | 0 |
| shared orbit | 6 | 0.122 | 0.179 | 0.206 | 0 |

**One result worth the manuscript's attention.** Section 9.1 singles out the
`m = 5` full-model slope, which the manuscript reports as 2.119 — "exceeded 2.0
by about 0.119, a finite-sample deviation of roughly two standard errors" — and
discusses it as a real feature. This run gets **1.967** (WLS 1.994) at the same
design and seed convention, i.e. right on 2.0, and the gain-residual slope there
is 0.033, so the raw gain tracks `n G` with no drift to explain. The anomaly
does not reproduce. Section 9.1 appears to be interpreting that run's Monte
Carlo noise, and the sentence should probably go.

The calibrated crossover ratios also track Table 6 closely at the group orders
where the models differ — shared/full of 0.687, 0.630, 0.608 for `m = 4, 5, 6`
against the manuscript's 0.705, 0.672, 0.610 — while coming out lower at
`m = 2, 3` (0.758, 0.788 against 0.883, 0.943).

Two caveats about the shipped files. The manifest records `dirty: true` because
the run was launched with the results README uncommitted; the commit it names is
the infrastructure commit, and no code differed. And the `higher_mode` rows of
`crossover_ratio_summary.csv` are meaningless — under the corrected
misspecification scenario the constrained detectors barely reach 50% power
inside the grid, so only one effect yields an internal crossover and the ratios
(18.95, 25.12, NaN) are dividing near-noise by near-noise. The `_n` columns show
this; the exact-orbit and independent-fundamental rows are the usable ones.

## Model selection, and a scenario defect it exposed

Every efficiency number in the manuscript is an *oracle* number: the detector
matching the generating geometry is chosen in advance. `regimeshift.selection`
chooses it from the data instead.

Selection cannot use the detector scores. A score is `gain - penalty` against
*its own* null, and those nulls differ — Model A pools an unrestricted
multinomial, Models B/C/D pool a fundamental coordinate — so the scores have
different origins. `code_lengths` returns the total description length of the
same data under each hypothesis, which is comparable, and the detector scores
fall out exactly as differences (asserted to 1e-9):

    score_A = L(null_full) - L(full)
    score_B = L(null_fundamental) - L(fundamental)
    score_C = L(null_fundamental) - L(shared_orbit)

Selecting the shortest code answers both questions at once — whether a change
occurred and what kind — and the six candidates include the two nulls, so no
separate detection step is needed.

At `m = 6`, effect 0.25, 200 trials:

| generated from | n/side | picks the generating family | false-change rate |
|---|---:|---:|---:|
| exact orbit | 200 / 800 / 3200 | 66% / 90% / 94% | — |
| higher mode | 200 / 800 / 3200 | 1% / 25% / 100% | — |
| no change | 200 / 800 / 3200 | — | 4.5% / 0% / 0% |

### The defect

The `independent_fundamental` scenario selects `fundamental` almost never at
`m = 6` — 3% — choosing `approximate_orbit` 93% of the time instead. That is
the selector being right, not wrong.

Section 8.2 fixes the angular offset at **0.713 rad** while the one-step
rotation is `2 pi / m`, which *shrinks* as `m` grows. The scenario therefore
slides toward being an orbit:

| m | one-step rotation | distance from nearest orbit | selector picks |
|---:|---:|---:|---|
| 3 | 2.094 rad | 1.12 effects | `fundamental` 100% |
| 4 | 1.571 rad | 0.76 effects | `fundamental` 91% |
| 5 | 1.257 rad | 0.53 effects | `approximate_orbit` 73% |
| 6 | 1.047 rad | 0.40 effects | `approximate_orbit` 93% |

The crossover near 0.6 effects matches Model D's winning band independently
measured in the README sweep.

**Why this matters for the manuscript.** The scenario meant to represent "Model
B territory" is not holding its distance from Model C's hypothesis constant
across `m`. Any `m`-dependence in the independent-fundamental results — Table 6's
row included — is therefore partly an artifact of the scenario drifting toward
an orbit as `m` increases, not a property of the detectors. It also explains, and
quantifies, the earlier observation that Model C stays competitive there at
`m = 6`: the data is only 0.40 effects from its hypothesis.

The fix is to hold the angular offset at a fixed *fraction* of `2 pi / m`, or to
hold the orbit distance itself constant, rather than fixing it in radians. This
implementation keeps the manuscript's constant so the reproduction stays
faithful; `tests/test_selection.py` pins the drift so it cannot go unnoticed.

## Not implemented

Scope limits carried over from Sections 11 and 13 of the manuscript:

* Unknown-boundary scanning and the `log n` location cost of Section 4.3, and
  any sequential or online stopping rule (Section 14.4).
* The block/codon-phase family in which `C_3` permutes phase blocks while each
  phase carries its own nucleotide alphabet (Section 11). The direct simulation
  here identifies group order with alphabet size.
* The approximate-orbit interpolation and stabilizer-adaptive label codes of
  Sections 14.1–14.2.
* Markov, hidden-state and continuous-observation extensions.

## Figures

The manuscript's five figures are plots of the four report CSVs. They are not
regenerated here; `python -m regimeshift analyse` rebuilds the underlying tables
from any results file.

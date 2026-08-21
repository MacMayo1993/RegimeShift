# Independent audit of the v4 framework

An accuracy audit of the idea in `docs/paper/geometric-complexity-v4.md` and its
reference implementation, followed by a plain-English explanation and a map of
what the work builds on.

The audit was done against the repository at commit `bc72290`, by:

* re-deriving every closed-form penalty in the framework from scratch, rather
  than reading it off the manuscript;
* running the full test suite — **625 fast tests and all 52 slow statistical
  tests pass** (the slow layer takes 1h37m on this four-core machine, not the
  advertised five minutes; see F12);
* recomputing **all** of the manuscript's headline tables directly from
  `results/v3-production/full_results.csv`, without using
  `regimeshift.analysis`;
* running three simulations that are **not** in the repository, chosen to test
  the central claim in a way the repository's own regressions structurally
  cannot;
* spot-checking the citations that carry weight against the published record.

**Verdict: the mathematics is correct, the code implements it faithfully, and
the reported numbers reproduce.** The framework is also unusually candid about
its own limits — several corrections that an auditor would normally have to
find are already recorded in the manuscript and in `docs/paper-notes.md`. The
findings below are therefore mostly about *emphasis and framing*, plus one small
code defect, which is fixed in this branch.

---

## 1. What was checked, and what it showed

### 1.1 The mathematics — re-derived independently

| Claim | Status | How it was checked |
|---|---|---|
| `d_full = m − 1`, `d_fund = 1` (m=2), `2` (m≥3) | **Correct** | Real irreducible decomposition of the permutation representation of `C_m` on the sum-zero tangent space. Modes `k` and `m−k` pair into 2-D real components; `k = 1` is 1-D only when `m = 2`. |
| Known-split penalty `(d/2)·log(n_L n_R / n)` | **Correct** | Two `d`-dimensional fits at `(d/2)log n_L` and `(d/2)log n_R` minus one at `(d/2)log n`. With `n_L = ρn` this is `(d/2)log n + (d/2)log ρ(1−ρ)`, so the `log n` coefficient is `d/2` and `ρ` only moves the bounded term — exactly as claimed. |
| Model C's alternative = fitting **one** coordinate to `c_L + roll(c_R, −s)` | **Correct, and exact** | Because the chart is equivariant (`p(R^s η) = roll(p(η), s)`), scoring `c_R` against `roll(p, s)` equals scoring `roll(c_R, −s)` against `p`. Summing the counts is therefore the exact profile likelihood, not an approximation. |
| Model C's continuous increment is 0 | **Correct** | The alternative and its null both carry one `d_fund`-vector; only a discrete label is added. |
| `K* = 1/(2 ln 2) = 0.7213475…` | **Correct, and definitional** | It is Schwarz's one-half in base 2. The repo says so itself. See finding F7. |
| Local JSD coefficient `(1 − cos 2π/m)/4` | **Correct** | Analytically: `JS ≈ ⅛‖η − Rη‖²_F = ¼(1 − cos 2π/m)‖η‖²`. Numerically at `ε = 10⁻³`: 0.50000, 0.37513, 0.25000, 0.12500 for `m = 2, 3, 4, 6` against 1/2, 3/8, 1/4, 1/8. |
| Model D's Schur complement `J_eff = L₁L₂/(L₁+L₂)` | **Correct** | Joint information `[[L₁+L₂, L₂],[L₂, L₂]]`; Schur complement for `δ` is `L₂ − L₂²/(L₁+L₂) = L₁L₂/(L₁+L₂)`. |
| Model D's penalty `(d/2)log(1 + τ²J_eff)` | **Correct, with a caveat worth writing down** | See below. |

The Model D correction deserves a note, because it is the one place where a
reader can talk themselves into the superseded formula. Laplace on the *joint*
`(η, δ)` gives

```
½·logdet(H + diag(0, τ⁻²)) = ½·log(L₁+L₂) + ½·log(1 + τ²·L₁L₂/(L₁+L₂))
```

per direction — that is, exactly `(d/2)log n` for the shared state **plus**
`pen_δ`. The corrected formula is right *provided* the shared state is charged
the full `(d/2)log n` with `n = L₁ + L₂`, which is what `selection.code_lengths`
does. If instead you marginalise `δ` first at fixed `η` and then profile `η`,
the naive `(d/2)log(1 + τ²L₂)` appears correct — I reproduced that decomposition
numerically and it matches the superseded formula to 0.001 nats. The two
decompositions must give the same total, and they do; they simply split it
differently. The repository's own brute-force test
(`test_deviation_penalty_matches_brute_force_marginalisation`) isolates the right
quantity — the ratio of the D and C marginals, where `η`'s measure cancels — and
so lands on the corrected formula. **Suggested doc change:** say this explicitly
in §3.4 or `paper-notes.md`, because "which decomposition" is precisely the trap
v3.1 fell into.

### 1.2 The reported numbers — recomputed from the committed run

Every headline table reproduces exactly from `full_results.csv` using an
independent regression written for this audit:

* **Penalty slopes** (Tables 3–5): 1.4882, 1.9670, 2.4876 (full, `m=4,5,6`);
  0.5217, 0.9962, 0.9707, 0.9636, 1.0357 (fundamental, `m=2…6`); −0.2280,
  0.0406, 0.0832, 0.1392, 0.1219 (shared orbit). Identical to the published
  values.
* **The identity `penalty slope = Δd/2 − s`** holds to floating point. Note it is
  algebraically forced — the score *is* the gain minus an exactly computed
  penalty — so it is a consistency check, not evidence. The manuscript says this
  (§8.3, Appendix B); the README is less careful (finding F1).
* **Zero-threshold null rates** (§9.6): mean/worst 0.0053/0.063 (full),
  0.0090/0.063 (fundamental), 0.0531/0.281 (shared orbit). Reproduces.
* **Misspecification power** (Table 7) and **shift recovery** (Table 8):
  reproduce to the last digit.
* **Crossover ratios** (Table 6): reproduce. See finding F3 about how thin they
  are.
* The README's worked example reproduces bit-for-bit
  (`gain=23.629 / 23.481 / 22.832`, `score=6.917 / 16.796 / 21.223`).
* **The two selection tables** of §10 and §11 — which are *not* produced by the
  production grid — regenerate exactly from `scripts/regenerate_selection_tables.py`,
  every cell, including the paired-advantage table that withdraws Model D's
  "middle band" claim (58/77/90%, 1/21/100%, 4/2/0%; and −0.31, +0.28, −0.18,
  −2.55, −7.38 nats). This is the strongest provenance in the repository: the
  script is deterministic and the manuscript quotes its output verbatim.
* 936 detector rows = 312 configurations × 3; 312 × (500 + 1000) = 468,000
  datasets. The arithmetic in Table 2 is right.

### 1.3 Three checks the repository does not run

The score regressions cannot independently establish the complexity law, because
the penalty is inserted into the score before regressing. So I measured the
dimensions a different way: under the null, `2 × raw gain` for a regular
`d`-dimensional split is `χ²_d`, so the **mean null raw gain is `d/2`** — no
penalty involved, nothing inserted. At `m = 4`, `η = 0`, 200 trials per point:

| n per side | full (predict 1.5) | fundamental (predict 1.0) | shared orbit (predict 0) |
|---:|---:|---:|---:|
| 200 | 1.519 | 0.995 | 0.626 |
| 800 | 1.510 | 0.978 | 0.645 |
| 3,200 | 1.616 | 1.062 | 0.654 |
| 12,800 | 1.554 | 1.023 | 0.668 |

Over a 64-fold range of sample size the full and fundamental gains sit at their
predicted `d/2` and the shared-orbit gain is **flat** — it does not grow with
`log n` at any rate. This is the cleanest evidence in or out of the repository
that Model C's continuous increment really is zero, and it is worth adding to
the test suite.

The second check explains §9.6 mechanically. Repeat the same experiment at a
*regular* orbit point (`η = 0.25`, still no change between segments):

| n per side | shared-orbit mean null gain |
|---:|---:|
| 200 | −4.25 |
| 800 | −20.4 |
| 3,200 | −90.5 |
| 12,800 | −379.0 |

Model C's protection under the null is **linear in `n` and comes from the
likelihood, not from the penalty** — mis-aligned pooling is a strictly wrong
model and is punished at rate `n`. At the collapse point `η = 0` that protection
vanishes entirely (the table above), leaving only `log(m−1)`, which is *zero* at
`m = 2`. That is the whole story behind the 0.281 worst-case null rate, and it
is a sharper statement than "the raw rule is not conservative".

The third check is the Model D marginalisation described in §1.1.

### 1.4 Citations

Spot-checked against the published record: Wang, Zou & Yin, *Change-point
detection in multinomial data with a large number of categories*, **Annals of
Statistics** 46(5), 2018 — confirmed. Pérez-Ortiz, Lardy, de Heide & Grünwald,
*E-statistics, group invariance and anytime-valid testing*, **Annals of
Statistics** 52(4), 1410–1432, 2024 — confirmed, including volume and pagination.
The East-model constant is real: `ln(1/gap) ∼ ln²(1/q)/(2 ln 2)` as `q → 0`
(Aldous–Diaconis; sharpened in the kinetically-constrained-models literature).
See F7 for what that does and does not license.

---

## 2. Findings

Twelve findings, ordered by how much they should change what the work says about
itself. None of them is a mathematical error.

**F1 — The README oversells the regression evidence.** `README.md` says the full
and fundamental slopes are "the clearest evidence that an independently fitted
fundamental family obeys a different complexity law". The manuscript's own §8.3
explicitly denies this: the penalty is inserted, so the regression measures only
the gain residual `s`. The clearest evidence is elsewhere — the `χ²_d` null-gain
reading in `blocks.py` §12.3, and the table in §1.3 above. *Fix: soften the
README bullet and point it at those.*

**F2 — "Near-zero" residual slopes are, at some group orders, statistically
non-zero.** `m = 2`: −0.228 ± 0.084 (2.7 se). `m = 6`: 0.122 ± 0.038 (3.2 se).
The manuscript's qualified wording is defensible, and the `Δd/2 − s`
decomposition correctly localises the drift in the likelihood gain rather than
in a hidden penalty. But Table 5 prints the slopes without their standard
errors, which are already in `score_regression_summary.csv`. *Fix: print them.*

**F3 — The crossover medians are thinner than they look.** Table 6's headline —
"31%, 37% and 39% fewer observations" — is a median over the *interior* crossover
estimates in each cell, and `crossover_ratio_summary.csv` records how many that
is: **2 to 4 points per cell**. §9.4 discloses the aggregate (44 of 156 estimates
out of grid) but not the per-cell count. A median of three numbers should be
labelled as one. *Fix: add the `_n` column to Table 6.*

**F4 — Model D's `τ = 0.05` is an unpriced free parameter.** It is not estimated,
not encoded, and not swept. Any selection result that returns `approximate_orbit`
— including the statement that exact-orbit recovery is "58% / 77% / 90%, with
essentially all the remainder going to `approximate_orbit`" — is conditional on
that number. §14.7 lists a data-chosen `τ` as future work, which is right, but
the recovery percentages are quoted in the meantime as if they were properties of
the selector. *Fix: sweep `τ ∈ {0.02, 0.05, 0.1, 0.2}` and report the sensitivity
next to the recovery table; it is a cheap run.*

**F5 — "Distance from the nearest orbit" means nearest *non-identity* orbit
point.** The tabulated 1.12 / 0.76 / 0.53 / 0.40 are distances to the nearest
`R^s η_L` with `s ≠ 0`, which is the right quantity (Model C never considers the
identity). But at `m = 3` and `m = 4` the *un-shifted* left state is closer
(0.655 in both cases) than any shifted one, so the phrase is loose. *Fix: one
word — "nearest non-identity orbit point".*

**F6 — One comparison is more robust than the paper admits, and it is the
important one.** §10 says cross-model lengths omit `O(1)` terms and that
"comparisons among C and D at fixed `L` … should be read as provisional". True
for C-vs-D. But **C against its own null is not provisional**: both hypotheses use
the same `η` chart, and both carry total information `(L₁+L₂)·I`, so the
Jeffreys/Fisher-volume terms are identical and cancel exactly. The entire
difference is `gain − log(m−1)`, with no convention-dependent residue. The
framework's central claim therefore rests on its *most* robust comparison, not
its most fragile one. *Fix: say so — it is a free strengthening.*

**F7 — The `K*` / East-model resemblance should be demoted further.** The
manuscript already says `K*` is definitional and that the coincidence only
carries weight if the East-model occurrence is not likewise a units artefact.
Two observations settle it in the direction the paper suspects. First, the
constants have different **dimensions**: `K*` is bits per e-fold of *sample
size*; the East constant is the coefficient of `ln²(1/q)` in a *relaxation time*.
Second, the `2 ln 2` there arises from the binary hierarchy of the energy barrier
(`log₂` of a length), i.e. from counting in base 2 — the same kind of artefact.
*Fix: move Appendix D's bridge to a footnote and state the dimensional mismatch;
it costs nothing and removes the one place a hostile reader would start.*

**F8 — Calibration is assumed available.** Every honest comparison in the paper
is made at a calibrated 5%, which requires 1,000 null draws from the correct
no-change distribution at every configuration. In an application you do not have
that; you have one dataset and one boundary. The raw rule is the thing a
practitioner would actually run, and §9.6 shows it is the thing that
misbehaves. *Fix: state the practical consequence — Model C needs either a
calibration source or a genuinely coded null, and at `m = 2` it needs one badly.*

**F9 — Most of the measured advantage comes from Model B, not Model C.** From
Table 6 at `m = 6`: fundamental/full = 0.708, shared/fundamental = 0.852. Knowing
the *subspace* buys ~29% fewer observations; adding the much stronger claim that
the change is an *exact rotation of a shared state* buys a further ~15%. The
paper's headline ratio (0.608) is the product of the two, and the framing
naturally reads as if the orbit constraint is doing the work. It is doing the
smaller half. *Fix: report the decomposition; it is more interesting than the
product, and it is what tells a practitioner whether Model C is worth the risk.*

**F10 — And the risk is large and quantified.** The same file gives the
misspecification crossover ratios: under a higher-mode change the constrained
detectors need roughly **19× and 25×** the sample length of the full detector at
`m = 4, 5` (single interior estimates, so order-of-magnitude only). The asymmetry
— correct structure buys ~1.6×, wrong structure costs ~20× — is the single most
decision-relevant number in the study and it does not appear in the manuscript.
*Fix: put it in §9.5.*

**F11 — A code defect: the block detectors did not validate their inputs.**
`validate_pair` in `detectors.py` exists precisely because a mismatched pair
"yields a score built from two different alphabet sizes — a wrong number rather
than an error". The block module had no equivalent:
`block_shared_orbit_detector` and `block_fundamental_detector` happened to raise
from inside their fits, but `block_full_detector` computed its likelihood from
whatever array it was handed and its penalty from `geometry`, so a `(3, 4)` count
array scored against `BlockGeometry(6, 4)` returned a plausible score with
`dimension_increment = 18` instead of 9. **Fixed in this branch**:
`validate_block_pair` now guards all three block detectors, with tests.

**F12 — The advertised slow-suite runtime is off by an order of magnitude.**
`README.md` says `pytest -m slow` takes "~5 minutes on 4 cores". On an idle
four-core machine it took **1 hour 36 minutes**, essentially all of it in the
module-scoped grid fixture (4,532 s) and the split-fraction test (1,257 s). The
arithmetic says the longer figure is the honest one: 160 configurations × 600
datasets, each costing up to ten L-BFGS fits at `m = 6` — on the order of a
million optimisations. Everything passes; only the estimate is wrong, and a
contributor who trusts it will assume the run has hung. *Fix: quote a measured
range, or note the hardware the five minutes was measured on.*

---

## 3. The idea in plain English

### The situation

You are watching something that cycles through a fixed set of labelled states —
the three reading frames of a DNA sequence, the six phases of a motor, the
operating modes of a sensor. All you record is *how often each state occurs*.
At some known moment, the pattern of frequencies changes. You want to know
whether it really changed, and what kind of change it was.

### The three (now four) stories

1. **"Something changed."** The frequencies are just different numbers now.
   Nothing links the before to the after.
2. **"It moved within its natural family."** The system stayed the kind of thing
   it was — the same cyclic-looking shape — but moved to a different one of those
   shapes.
3. **"Nothing changed except the phase."** The system is doing *exactly* what it
   was doing before, one step around the cycle. A motor slipped a phase; a DNA
   sequence lost a base and shifted its reading frame.
4. **"Almost the phase."** Story 3, allowing a small wobble.

### Why the distinction is worth money

Think of it as what you have to *say* to describe the change to someone else.

* Story 1: you must transmit a whole new set of frequencies. Expensive, and the
  cost grows as you collect more data, because more data means you must specify
  those numbers more precisely.
* Story 2: you must transmit a new pair of coordinates. Cheaper — two numbers
  instead of five — but still growing with the amount of data.
* Story 3: you transmit **nothing new about the shape at all**. The shape is
  already known from the first segment. All you say is *which way it rotated* —
  one of `m − 1` choices, a single word. That cost does not grow with the amount
  of data. Ever.

That last sentence is the whole idea. The usual way to make a model cheap is to
give it fewer parameters (**dimension reduction**, story 2). Story 3 does
something categorically different: it gives the second segment **no new
parameters at all**, only a shared one plus a label (**parameter sharing**). In
the standard accounting — half a "unit of cost" per parameter per e-fold of data
— story 1 pays `m−1` units, story 2 pays 2 units, and story 3 pays **zero**.
These two ways of exploiting symmetry are usually lumped together; the work's
contribution is separating them and showing they sit at different levels of a
hierarchy.

### The catch, which the work is honest about

A sharper claim is cheaper *when it is true* and disastrous when it is not.
The measured trade at `m = 6`:

* Change really is a phase slip → you need about **40% less data** than the
  unrestricted method to detect it. (Of which roughly two-thirds comes from
  knowing the subspace, one-third from the sharing.)
* Change is of a kind the cyclic story cannot represent at all → you need about
  **twenty times more** data. The constrained detector isn't merely weaker; its
  "evidence" can go negative, because forcing an alignment that isn't there is
  worse than not aligning.

There is also a specific failure mode with a clean geometric cause. When the
system sits exactly at the perfectly-uniform state, *every* rotation does
nothing — the "which way did it rotate" question has no answer. The detector
still asks it, picks whichever rotation noise favours, and calls the result
evidence. With only two states there is no label cost at all to offset that, and
the false-alarm rate blows out to 28%. This is why every comparison in the paper
is made at a calibrated threshold rather than by the raw rule.

### What is actually established

That the accounting is right, that the code implements it, and that simulated
data behaves the way the accounting says. Not that any real system is a cyclic
orbit — **no real dataset has been analysed**, which the paper states as its
largest gap. The proposed test case (detecting DNA frameshifts as a phase
rotation) is described but has not been run.

---

## 4. What this builds on

Nothing in the toolkit is new. The claim — stated this way in the paper itself —
is that the *combination* is.

| Layer | Source it comes from | What is taken |
|---|---|---|
| "Half a `log n` per parameter" | **Schwarz (1978)**, BIC | The entire penalty arithmetic for Models A and B. `K*` is this constant in base 2. |
| Description length as a *coding* principle | **Rissanen**; **Barron, Rissanen & Yu (1998)**; **Grünwald (2007)** | The framing that lets a *discrete label* (`log(m−1)`) and a *continuous parameter* (`(d/2)log n`) be priced in the same units. Model C's whole argument lives here. |
| Two-part codes | Same lineage; **Krichevsky–Trofimov (1981)** for the exact multinomial code *not* used | The `log(m−1)` label cost is a two-part code, one legitimate choice among several. |
| Information criteria for changepoints | **Yao (1988)**; **Davis, Lee & Rodriguez-Yam (2006)** | The idea of scoring a segmentation by likelihood plus a dimension penalty. |
| Categorical changepoints | **Wang, Zou & Yin (2018)**; **Truong & Runge (2024)** | The observation model — multinomial counts on either side of a boundary. Both test *unrestricted* change; neither has a subspace or orbit alternative. |
| Dimension-reduced changepoints | **Yu et al. (2026)** | The nearest analogue of Model B — but their subspace is *learned* from sparsity, where here it is *fixed by the group representation*. |
| Which subspace, and why | **Serre (1977)**, representation theory of finite groups | The decomposition `R^m = V_triv ⊕ V_fund ⊕ …` and hence `d_fund = 1` or `2`. This is what makes the subspace canonical rather than chosen. |
| The metric that makes `‖η‖` meaningful | **Amari & Nagaoka (2000)**, information geometry | Fisher-orthonormal coordinates, so "effect size" is a distance and the rotation is a literal rotation. |
| Group invariance in testing | **Eaton (1989)**; **Pérez-Ortiz, Lardy, de Heide & Grünwald (2024)** | The nearest living relative. Their maximal invariant for `C_g` on the fundamental subspace is exactly what Model C conditions on — but they ask a sequential testing question with no coding component. This is also the most credible route to a sequential version. |
| What happens when the orbit collapses | **Watanabe (2009, 2013)**; **Drton & Plummer (2017)** | Singular learning theory. Cited to *fence off* the `η = 0` stratum, where BIC dimension counting is not the right theory — not to solve it. |
| The application target | **Rho, Tang & Ye (2010)**, FragGeneScan | Frameshift detection via phase-specific emissions already exists. The proposed contribution is the MDL decomposition of phase constraint, not the detection. |
| The speculative bridge | **Cancrini et al. (2008)**, East model | The `1/(2 ln 2)` coincidence. See F7 — treat as a footnote. |

**The genuinely new part**, stated as narrowly as the evidence supports: the
explicit four-way MDL separation of *unrestricted change*, *independent
invariant-subspace change*, *shared exact-orbit change*, and *shrinkage toward
an orbit*, with the observation that the third carries **no continuous-dimension
increment at all** on a regular orbit stratum while the second carries
`d_fund`. Two external prior-art reviews failed to find that combination
elsewhere, which is evidence of absence rather than proof of it — and the
repository says so.

---

## 5. Recommended next steps, in priority order

1. **Analyse one real dataset.** This is the paper's own limitation 7 and it is
   correct. Every number here comes from data generated by the models being
   scored. A controlled categorical phase dataset — it need not be the frameshift
   application — would do more for the argument than any further simulation.
2. **Report the B-vs-C decomposition and the misspecification cost** (F9, F10).
   Two numbers already in the committed results; they change how a reader weighs
   the method.
3. **Sweep `τ`** (F4) before quoting selector recovery rates again.
4. **Add the null-gain dimension check** of §1.3 to the test suite: it is the
   only non-circular measurement of the complexity law in the framework, and it
   costs a few seconds.
5. **Claim the robustness of the C-vs-null comparison** (F6), and demote the
   East-model bridge (F7).

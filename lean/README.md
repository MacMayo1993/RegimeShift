# Lean verification

Machine-checked companion to `docs/paper/geometric-complexity-v4.md`. The scope
is deliberately narrow: Lean verifies the **geometry and the dimension
counting**, which is where the paper's novelty lies. The asymptotic coding law
is assumed, as a classical theorem, exactly where the manuscript assumes it.

> **Status: compiling.** `.github/workflows/lean.yml` builds this library
> against mathlib on every push touching `lean/`, and fails the job if any
> `sorry` appears. Every entry marked *proved* below is checked by the
> compiler.

## Claim-by-claim status

| Manuscript claim | Section | Lean | Status |
|---|---|---|---|
| `p(R_g^r η) = T_r p(η)` (equivariance of the chart) | 5.3 | `RegimeShift.equivariance`, `equivariance_std` | proved (for the chart as scaled here; see the `g = 2` note below) |
| Softmax chart lands in the simplex, is positive | 5.3 | `softmax_sum_one`, `softmax_pos` | proved |
| Shift is a `C_g` action on categories | 2.4 | `shift_shift`, `shift_zero` | proved |
| Standard basis `e j = exp(2πij/g)` is a unitary character | 5.2 | `zeta_map_add`, `zeta_norm_one` | proved |
| Regular orbit ⇒ trivial stabiliser, `g` distinct points, `r` identified | 3.3, Prop 1 asm 2 | `act_injective`, `stabilizer_trivial` | proved |
| Orbit collapse at `η = 0` | 3.3, 6.3 | `act_zero_state` | proved |
| `pen_split(d;L₁,L₂)` and its `ρ(1-ρ)` form | 4.2 | `penSplit_rho` | proved |
| Leading coefficients are integer multiples of `K*` | 4.4 | `leading_coeff_multiple_of_Kstar` | **proved, but near-vacuous** — reduces to `n / (2 log 2) = n * K*` for `n : ℕ`, so it verifies only that the increments are whole numbers, which is a definition and not a theorem |
| `Δd_A = g-1`, `Δd_B = Δd_D = d_fund`, `Δd_C = 0` | 3.1–3.4 | `Model.increment` | **definition, not theorem** |
| `d_fund = 1` for `g=2`, `2` for `g≥3` | 2.4 | `dFund`; `finrank_fund_two`, `finrank_fund_ge_three` | **partly proved** — the span of the two logit directions has dimension `dFund g`; that this span *is* the fundamental isotypic component is not formalised |
| `cos φⱼ`, `sin φⱼ` independent iff `g ≥ 3`; sine vanishes at `g=2` | 5.2 | `cosSin_linearIndependent`, `sinVec_two` | proved |
| `(d/2)·log L` leading penalty from a `d`-dimensional fit | 4.1 | — | **assumed** (classical BIC/Laplace; not in mathlib) |
| Fisher orthonormality of the fundamental basis | 5.1–5.2 | — | not formalised |
| Local JSD coefficient `(1-cos(2π/g))/4` | 5.4 | — | not formalised |
| Singular asymptotics / RLCT at orbit collapse | 4.4, 6.3 | — | out of scope |
| Monte Carlo slopes, calibrated power, misspecification | 9 | — | empirical, out of scope |

**One scaling caveat.** `Fourier.logit` applies `√2` uniformly. That is the
Fisher-orthonormal scaling for `g ≥ 3` only; the `g = 2` sign representation needs
a factor of `1`, which is what `fourier_design_matrix(2)` returns in the Python
implementation. So at `g = 2` the formalised chart is `√2` times the code's. No
proof here depends on it — equivariance is invariant under positive rescaling of
the logits — but the two are not the same chart, and Fisher orthonormality, the
one property that would distinguish them, is unformalised.

The three most valuable next targets, in order: identifying the span of
`cosVec`/`sinVec` with the fundamental isotypic component of the permutation
representation of `C_g` on the sum-zero tangent space, which would finish the
`dFund` story; Fisher orthonormality (root-of-unity sums); the local JSD limit.

## Building

```
cd lean && lake exe cache get && lake build
```

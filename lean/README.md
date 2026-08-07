# Lean verification

Machine-checked companion to `docs/paper/geometric-complexity-v4.md`. The scope
is deliberately narrow: Lean verifies the **geometry and the dimension
counting**, which is where the paper's novelty lies. The asymptotic coding law
is assumed, as a classical theorem, exactly where the manuscript assumes it.

> **Status: not yet compiled.** These files were written in an environment with
> no Lean toolchain and no network access to `lean-lang.org`, so nothing here
> has been checked by the compiler. The CI workflow
> (`.github/workflows/lean.yml`) is the first thing that will actually build
> them; treat any claim below as provisional until that workflow is green.
> There are no `sorry`s — but an uncompiled proof is not a proof.

## Claim-by-claim status

| Manuscript claim | Section | Lean | Status |
|---|---|---|---|
| `p(R_g^r η) = T_r p(η)` (equivariance of the chart) | 5.3 | `RegimeShift.equivariance`, `equivariance_std` | proved |
| Softmax chart lands in the simplex, is positive | 5.3 | `softmax_sum_one`, `softmax_pos` | proved |
| Shift is a `C_g` action on categories | 2.4 | `shift_shift`, `shift_zero` | proved |
| Standard basis `e j = exp(2πij/g)` is a unitary character | 5.2 | `zeta_map_add`, `zeta_norm_one` | proved |
| Regular orbit ⇒ trivial stabiliser, `g` distinct points, `r` identified | 3.3, Prop 1 asm 2 | `act_injective`, `stabilizer_trivial` | proved |
| Orbit collapse at `η = 0` | 3.3, 6.3 | `act_zero_state` | proved |
| `pen_split(d;L₁,L₂)` and its `ρ(1-ρ)` form | 4.2 | `penSplit_rho` | proved |
| Leading coefficients are integer multiples of `K*` | 4.4 | `leading_coeff_multiple_of_Kstar` | proved (given the increments of §3) |
| `Δd_A = g-1`, `Δd_B = Δd_D = d_fund`, `Δd_C = 0` | 3.1–3.4 | `Model.increment` | **definition, not theorem** |
| `d_fund = 1` for `g=2`, `2` for `g≥3` | 2.4 | `dFund` | **definition, not theorem** — the isotypic decomposition is not formalised |
| `(d/2)·log L` leading penalty from a `d`-dimensional fit | 4.1 | — | **assumed** (classical BIC/Laplace; not in mathlib) |
| Fisher orthonormality of the fundamental basis | 5.1–5.2 | — | not formalised |
| Local JSD coefficient `(1-cos(2π/g))/4` | 5.4 | — | not formalised |
| Singular asymptotics / RLCT at orbit collapse | 4.4, 6.3 | — | out of scope |
| Monte Carlo slopes, calibrated power, misspecification | 9 | — | empirical, out of scope |

The three most valuable next targets, in order: `dFund` as a theorem about the
isotypic decomposition of the permutation representation of `C_g` on the
sum-zero tangent space; Fisher orthonormality (root-of-unity sums); the local
JSD limit.

## Building

```
cd lean && lake exe cache get && lake build
```

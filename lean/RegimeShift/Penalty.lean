/-
# Split penalties and the constant `K*`

Sections 4.2 and 4.4. Everything here is *definitional* on the statistics side:
`penSplit` is the BIC/Laplace split increment as the manuscript defines it, not
something derived from a coding theorem. Mathlib has no Laplace expansion of the
marginal likelihood, so the step from "the increment has continuous dimension
`d`" to "the penalty is `(d/2) log L`" is assumed here exactly as the manuscript
assumes it (Section 4.1). What is proved is the algebra downstream of that
assumption, including the `ρ(1-ρ)` reformulation that Section 9.7 checks
empirically.
-/
import Mathlib

namespace RegimeShift

open Real

/-- Incremental complexity of replacing one `d`-dimensional fit by two, at a
known boundary (Section 4.2). -/
noncomputable def penSplit (d L₁ L₂ : ℝ) : ℝ :=
  d / 2 * (Real.log L₁ + Real.log L₂ - Real.log (L₁ + L₂))

/-- With `L₁ = ρL` and `L₂ = (1-ρ)L`, the penalty is `(d/2) log L` plus a term
that depends on the split fraction but not on `L`: the split fraction affects
only the bounded term. -/
theorem penSplit_rho (d L ρ : ℝ) (hL : 0 < L) (hρ : 0 < ρ) (hρ' : ρ < 1) :
    penSplit d (ρ * L) ((1 - ρ) * L)
      = d / 2 * Real.log L + d / 2 * Real.log (ρ * (1 - ρ)) := by
  have h1 : (0:ℝ) < 1 - ρ := by linarith
  have hsum : ρ * L + (1 - ρ) * L = L := by ring
  rw [penSplit, hsum, Real.log_mul (ne_of_gt hρ) (ne_of_gt hL),
    Real.log_mul (ne_of_gt h1) (ne_of_gt hL),
    Real.log_mul (ne_of_gt hρ) (ne_of_gt h1)]
  ring

/-- Per-dimension penalty rate in bits per e-fold: Schwarz's one-half expressed
in base 2 (Section 4.4). -/
noncomputable def Kstar : ℝ := 1 / (2 * Real.log 2)

/-- Continuous-dimension increments of the four model classes on a regular
stratum (Sections 3.1–3.4). `d_fund` is `1` for `g = 2` and `2` for `g ≥ 3`
(Section 2.4). Model D's increment is the asymptotic one. -/
def dFund (g : ℕ) : ℕ := if g ≤ 2 then 1 else 2

inductive Model | A | B | C | D

/-- `Δd` for each model class. Read only at `g ≥ 2`: the truncated subtraction
below returns `0` at `g ∈ {0, 1}`, where the cyclic model is not defined and
`dFund` is likewise a placeholder. -/
def Model.increment : Model → ℕ → ℕ
  | .A, g => g - 1
  | .B, g => dFund g
  | .C, _ => 0
  | .D, g => dFund g

/-- **The counting statement of Section 4.4.** On regular strata and under the
BIC/Laplace convention, every leading coefficient in bits per e-fold is an
integer multiple of `K*`.

Be clear about how little this verifies. It reduces to `∀ n : ℕ, n / (2 log 2) =
n * K*`, which is arithmetic, and would hold for *any* `ℕ`-valued `increment`
whatsoever — the mathematical content is entirely in `Model.increment` being
`ℕ`-valued, and that is a definition transcribed from §3, not a theorem. What the
statement rules out is only a non-integer count; it says nothing about whether
the increments are the right ones. On singular strata the premise fails outright,
since the leading coefficient is a real log canonical threshold that need not be
a half-integer. -/
theorem leading_coeff_multiple_of_Kstar (M : Model) (g : ℕ) :
    ∃ n : ℕ, (M.increment g : ℝ) / (2 * Real.log 2) = n * Kstar :=
  ⟨M.increment g, by rw [Kstar]; ring⟩

end RegimeShift

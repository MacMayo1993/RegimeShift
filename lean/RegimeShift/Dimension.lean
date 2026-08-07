/-
# The dimension of the fundamental component

Section 2.4: `d_fund = 1` for `g = 2` (the sign representation) and `2` for
`g ≥ 3`. What is proved here is the dimension of the *selected* fundamental
component — the span of the two logit directions `cos φⱼ` and `sin φⱼ` of
Section 5.2 — namely that they are linearly independent exactly when `g ≥ 3`,
and that at `g = 2` the sine direction vanishes identically, leaving the one
sign direction.

This is not yet the full representation-theoretic statement. It does not prove
that this span *is* the fundamental isotypic component of the permutation
representation of `C_g` on the sum-zero tangent space; that decomposition
remains unformalised, and `dFund` remains the definition the rest of the
development uses. What it does establish is that the definition is not an
arbitrary numeral: it is the real dimension of the space the implementation
actually fits in, with the `g = 2` collapse to one dimension proved rather than
asserted.
-/
import RegimeShift.Basic
import RegimeShift.Penalty

namespace RegimeShift

open Real

variable {g : ℕ} [NeZero g]

/-- The cosine logit direction of Section 5.2, `φⱼ = 2πj/g`. -/
noncomputable def cosVec (g : ℕ) : ZMod g → ℝ :=
  fun j => Real.cos (2 * π * j.val / g)

/-- The sine logit direction. At `g = 2` this vanishes identically, which is
why the fundamental component drops to the one-dimensional sign
representation. -/
noncomputable def sinVec (g : ℕ) : ZMod g → ℝ :=
  fun j => Real.sin (2 * π * j.val / g)

@[simp] lemma cosVec_zero : cosVec g 0 = 1 := by
  simp [cosVec]

@[simp] lemma sinVec_zero : sinVec g 0 = 0 := by
  simp [sinVec]

/-- For `g ≥ 3` the angle `2π/g` lies strictly between `0` and `π`, so the sine
direction is genuinely nonzero — this is the whole reason the fundamental
component is two-dimensional there. -/
lemma sin_two_pi_div_ne_zero (hg : 3 ≤ g) : Real.sin (2 * π / g) ≠ 0 := by
  have hg3 : (3 : ℝ) ≤ (g : ℝ) := by exact_mod_cast hg
  have hgpos : (0 : ℝ) < (g : ℝ) := by linarith
  refine ne_of_gt (Real.sin_pos_of_pos_of_lt_pi ?_ ?_)
  · positivity
  · rw [div_lt_iff hgpos]
    nlinarith [Real.pi_pos]

/-- **The two directions are independent exactly when `g ≥ 3`.** -/
theorem cosSin_linearIndependent (hg : 3 ≤ g) :
    LinearIndependent ℝ ![cosVec g, sinVec g] := by
  haveI : Fact (1 < g) := ⟨by omega⟩
  rw [LinearIndependent.pair_iff]
  intro a b hab
  have h0 := congrFun hab 0
  have h1 := congrFun hab 1
  simp at h0
  have ha : a = 0 := by simpa using h0
  rw [ha] at h1
  simp [cosVec, sinVec, ZMod.val_one] at h1
  rcases h1 with h1 | h1
  · exact ⟨ha, h1⟩
  · exact absurd h1 (sin_two_pi_div_ne_zero hg)

/-- At `g = 2` the sine direction vanishes identically: `sin(πj) = 0`. -/
theorem sinVec_two : sinVec 2 = 0 := by
  funext j
  have h : 2 * π * (j.val : ℝ) / 2 = (j.val : ℝ) * π := by ring
  simp [sinVec, h, Real.sin_nat_mul_pi]

lemma cosVec_ne_zero : cosVec g ≠ 0 := by
  intro h
  have h0 := congrFun h 0
  rw [cosVec_zero] at h0
  exact one_ne_zero h0

/-- **`d_fund = 2` for `g ≥ 3`.** -/
theorem finrank_fund_ge_three (hg : 3 ≤ g) :
    Module.finrank ℝ (Submodule.span ℝ (Set.range ![cosVec g, sinVec g]))
      = dFund g := by
  have hg2 : ¬ (g ≤ 2) := by omega
  rw [finrank_span_eq_card (cosSin_linearIndependent hg)]
  simp [dFund, hg2]

/-- **`d_fund = 1` for `g = 2`**: the sine direction is identically zero, so the
span collapses to the one-dimensional sign representation. -/
theorem finrank_fund_two :
    Module.finrank ℝ (Submodule.span ℝ (Set.range ![cosVec 2, sinVec 2]))
      = dFund 2 := by
  have hrange : Set.range ![cosVec 2, sinVec 2] = {cosVec 2, (0 : ZMod 2 → ℝ)} := by
    rw [sinVec_two]
    simp
  rw [hrange, Set.pair_comm, Submodule.span_insert_zero,
    finrank_span_singleton (cosVec_ne_zero (g := 2))]
  simp [dFund]

end RegimeShift

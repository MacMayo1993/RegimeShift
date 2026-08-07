/-
# Cyclic shift action and the softmax chart

Formalises the objects of Sections 2.4 and 5.3 of the manuscript: the action of
`C_g` on categorical distributions by cyclic relabelling, and the
positivity-preserving softmax chart `p(η) = softmax(B η)`.

The index type is `ZMod g`, so the group and the alphabet coincide — this is the
*direct* cyclic categorical model of Sections 3–11, where the group acts on the
categories themselves and `a = g`. Block families (Section 12), where the group
permutes phase blocks and each block carries its own alphabet, are not modelled
here.
-/
import Mathlib

open Finset

namespace RegimeShift

variable {g : ℕ} [NeZero g]

/-- Category index of the direct cyclic model: the group order equals the
alphabet size. -/
abbrev Idx (g : ℕ) := ZMod g

instance : Nonempty (Idx g) := ⟨0⟩

/-- The action `T_r` of the group element `r` on a function of the categories,
by cyclic relabelling. Applied to a probability vector this is the `T_r` of
Section 2.4; applied to a logit vector it is the same permutation acting one
level down. -/
def shift (r : Idx g) (p : Idx g → ℝ) : Idx g → ℝ := fun j => p (j - r)

omit [NeZero g] in
@[simp] lemma shift_apply (r : Idx g) (p : Idx g → ℝ) (j : Idx g) :
    shift r p j = p (j - r) := rfl

omit [NeZero g] in
@[simp] lemma shift_zero (p : Idx g → ℝ) : shift 0 p = p := by
  funext j; simp [shift]

omit [NeZero g] in
/-- `shift` is an action of the additive group `ZMod g`. -/
lemma shift_shift (r s : Idx g) (p : Idx g → ℝ) :
    shift r (shift s p) = shift (r + s) p := by
  funext j
  simp [shift, sub_sub]

omit [NeZero g] in
/-- `shift r` is precomposition with the equivalence `j ↦ j - r`. -/
lemma shift_eq_comp (r : Idx g) (p : Idx g → ℝ) :
    shift r p = p ∘ (Equiv.subRight r) := rfl

/-- The softmax chart. With `ℓ = B η` this is `p(η)` of Section 5.3. -/
noncomputable def softmax (ℓ : Idx g → ℝ) : Idx g → ℝ :=
  fun j => Real.exp (ℓ j) / ∑ k, Real.exp (ℓ k)

lemma sum_exp_pos (ℓ : Idx g → ℝ) : 0 < ∑ k, Real.exp (ℓ k) :=
  Finset.sum_pos (fun k _ => Real.exp_pos (ℓ k)) univ_nonempty

lemma softmax_pos (ℓ : Idx g → ℝ) (j : Idx g) : 0 < softmax ℓ j :=
  div_pos (Real.exp_pos _) (sum_exp_pos ℓ)

/-- The chart lands in the simplex. -/
lemma softmax_sum_one (ℓ : Idx g → ℝ) : ∑ j, softmax ℓ j = 1 := by
  simp only [softmax]
  rw [← Finset.sum_div]
  exact div_self (ne_of_gt (sum_exp_pos ℓ))

/-- Softmax commutes with any relabelling of the categories: the normalising
constant is invariant under a bijection of the index set. -/
theorem softmax_comp (e : Idx g ≃ Idx g) (ℓ : Idx g → ℝ) :
    softmax (ℓ ∘ e) = (softmax ℓ) ∘ e := by
  funext j
  simp only [softmax, Function.comp_apply]
  rw [Equiv.sum_comp e (fun k => Real.exp (ℓ k))]

/-- Softmax intertwines the shift action on logits with the shift action on
probabilities. This is the level-crossing step of the equivariance identity
`p(R_g^r η) = T_r p(η)`: it reduces that identity to the purely linear claim
that `B` is equivariant, which is `Fourier.logit_rotate`. -/
theorem softmax_shift (r : Idx g) (ℓ : Idx g → ℝ) :
    softmax (shift r ℓ) = shift r (softmax ℓ) := by
  rw [shift_eq_comp, shift_eq_comp, softmax_comp]

end RegimeShift

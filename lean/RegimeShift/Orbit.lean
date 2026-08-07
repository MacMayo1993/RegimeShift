/-
# Regular orbits (Proposition 1, assumption 2)

Proposition 1 has two halves. The statistical half — that a zero-dimensional
increment yields a penalty constant in `L` — rests on the BIC/Laplace expansion
and is *assumed*, not proved (see `Penalty.lean`). The half that is provable
here is the identifiability content of assumption 2: on a regular orbit the
stabiliser is trivial, the `g` orbit points are distinct, and the relative shift
`r` is therefore recoverable. That is what makes the two-part code over the
`g - 1` nonidentity elements a valid code and the continuous increment zero.

At `η = 0` every group element acts trivially, the orbit collapses to a point,
and `r` is not identifiable — the singular stratum where the regular dimension
count does not apply at all.
-/
import RegimeShift.Fourier

namespace RegimeShift

variable {g : ℕ} [NeZero g]

/-- The action of `C_g` on the fundamental coordinate. -/
noncomputable def act (χ : FundChar g) (r : ZMod g) (η : ℂ) : ℂ := χ.e r * η

/-- The orbit is a genuine group action. -/
lemma act_add (χ : FundChar g) (r s : ZMod g) (η : ℂ) :
    act χ r (act χ s η) = act χ (r + s) η := by
  simp [act, χ.map_add]; ring

lemma act_zero (χ : FundChar g) (η : ℂ) : act χ 0 η = η := by
  have hne : χ.e 0 ≠ 0 := by
    intro h
    have h1 := χ.norm_one 0
    rw [h] at h1
    simp at h1
  have h := χ.map_add 0 0
  rw [add_zero] at h
  have h1 : (1 : ℂ) * χ.e 0 = χ.e 0 * χ.e 0 := by rw [one_mul]; exact h
  simp [act, (mul_right_cancel₀ hne h1).symm]

/-- **Regular orbits are free.** If the character is faithful and `η ≠ 0`, the
`g` orbit points `{R^s η}` are pairwise distinct, so the relative shift is
identified by the pair of segment states. -/
theorem act_injective (χ : FundChar g) (hχ : Function.Injective χ.e)
    (η : ℂ) (hη : η ≠ 0) : Function.Injective (fun r => act χ r η) := by
  intro r s hrs
  exact hχ (mul_right_cancel₀ hη hrs)

/-- **Trivial stabiliser at a regular point** — assumption 2 of Proposition 1,
in the form it is used: no nonidentity element fixes `η`. -/
theorem stabilizer_trivial (χ : FundChar g) (hχ : Function.Injective χ.e)
    (η : ℂ) (hη : η ≠ 0) (r : ZMod g) (hr : act χ r η = η) : r = 0 := by
  refine act_injective χ hχ η hη ?_
  show act χ r η = act χ 0 η
  rw [hr, act_zero]

/-- **Orbit collapse at the singular point.** At `η = 0` every group element
acts trivially, so the relative shift carries no information and the two-part
code of Proposition 1 is decoding an unidentifiable label. -/
theorem act_zero_state (χ : FundChar g) (r : ZMod g) : act χ r 0 = 0 := by
  simp [act]

end RegimeShift

/-
# The standard character `e j = exp(2πij/g)`

Instantiates `FundChar` with the concrete Fourier basis of Section 5.2, so that
`equivariance` is a statement about the actual family the implementation uses
rather than about an abstract character. The only content is that the map is a
homomorphism despite `j ↦ j.val` not being one: the discrepancy is a whole
number of turns, which the exponential kills.
-/
import RegimeShift.Fourier

open Real Complex

namespace RegimeShift

variable {g : ℕ} [NeZero g]

/-- `ζ_g^j`, with the exponent read through `ZMod.val`. -/
noncomputable def zeta (g : ℕ) (j : ZMod g) : ℂ :=
  Complex.exp (2 * Real.pi * Complex.I * (j.val : ℂ) / (g : ℂ))

lemma zeta_map_add (x y : ZMod g) : zeta g (x + y) = zeta g x * zeta g y := by
  have hg : (g : ℂ) ≠ 0 := Nat.cast_ne_zero.mpr (NeZero.ne g)
  set q : ℕ := (x.val + y.val) / g with hq
  have hsplit : (x.val + y.val : ℕ) = g * q + (x + y).val := by
    rw [ZMod.val_add, hq]
    exact (Nat.div_add_mod (x.val + y.val) g).symm
  rw [zeta, zeta, zeta, ← Complex.exp_add]
  have : 2 * Real.pi * Complex.I * (x.val : ℂ) / (g : ℂ)
        + 2 * Real.pi * Complex.I * (y.val : ℂ) / (g : ℂ)
      = (q : ℤ) * (2 * Real.pi * Complex.I)
        + 2 * Real.pi * Complex.I * ((x + y).val : ℂ) / (g : ℂ) := by
    have hcast : ((x.val + y.val : ℕ) : ℂ) = (g : ℂ) * (q : ℂ) + ((x + y).val : ℂ) := by
      exact_mod_cast hsplit
    field_simp at hcast ⊢
    push_cast
    linear_combination (2 * Real.pi * Complex.I) * hcast
  rw [this, Complex.exp_add, Complex.exp_int_mul_two_pi_mul_I, one_mul]

lemma zeta_norm_one (x : ZMod g) : ‖zeta g x‖ = 1 := by
  rw [zeta, Complex.norm_exp]
  have : (2 * Real.pi * Complex.I * (x.val : ℂ) / (g : ℂ)).re = 0 := by
    simp [Complex.div_re, Complex.mul_re, Complex.mul_im]
  rw [this, Real.exp_zero]

/-- The fundamental Fourier basis of Section 5.2 as a `FundChar`. -/
noncomputable def stdChar (g : ℕ) [NeZero g] : FundChar g where
  e := zeta g
  map_add := zeta_map_add
  norm_one := zeta_norm_one

/-- Equivariance for the concrete family: what
`tests/test_fourier_geometry.py::test_equivariance` checks numerically. -/
theorem equivariance_std (η : ℂ) (r : ZMod g) :
    softmax (logit (stdChar g) (zeta g r * η))
      = shift r (softmax (logit (stdChar g) η)) :=
  equivariance (stdChar g) η r

end RegimeShift

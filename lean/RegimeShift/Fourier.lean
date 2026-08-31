/-
# Fundamental Fourier family and the equivariance identity

Section 5.2–5.3 of the manuscript. The fundamental invariant component of the
tangent space at the uniform distribution is two-dimensional for `g ≥ 3` (one
cosine and one sine direction) and one-dimensional for `g = 2` (the sign
direction). Both cases are captured uniformly here by taking the coordinate
`η` to be a *complex* number and the basis to be a character `e : ZMod g → ℂ`:

* for `g ≥ 3`, `e j = exp(2πij/g)` and multiplication by `e r` is exactly the
  planar rotation `R_g^r` by `2πr/g` of Section 5.2;
* for `g = 2`, `e j = (-1)^j` is real and `η` may be taken real, recovering the
  sign representation.

The main theorem `equivariance` is the identity `p(R_g^r η) = T_r p(η)` that
`tests/test_fourier_geometry.py` checks numerically to `1e-13` through `g = 8`.

Scope: this file proves equivariance, which is an exact algebraic identity. It
says nothing about Fisher orthonormality of the basis (Section 5.1–5.2) or the
local Jensen–Shannon coefficient (Section 5.4); see `lean/README.md` for the
claim-by-claim status table.
-/
import RegimeShift.Basic

namespace RegimeShift

variable {g : ℕ} [NeZero g]

/-- A unitary additive character of `ZMod g`, playing the role of the
fundamental Fourier basis. Multiplication by `e r` is the action of the group
element `r` on the fundamental coordinate — the matrix `R_g^r` of Section 5.3. -/
structure FundChar (g : ℕ) where
  e : ZMod g → ℂ
  map_add : ∀ x y, e (x + y) = e x * e y
  norm_one : ∀ x, ‖e x‖ = 1

namespace FundChar

variable (χ : FundChar g)

omit [NeZero g] in
lemma mul_star_self (x : ZMod g) : χ.e x * (starRingEnd ℂ) (χ.e x) = 1 := by
  have habs : Complex.abs (χ.e x) = 1 := by
    simpa [Complex.norm_eq_abs] using χ.norm_one x
  rw [Complex.mul_conj, Complex.normSq_eq_abs, habs]
  norm_num

omit [NeZero g] in
/-- The character carries subtraction in the index to multiplication by the
conjugate — the only fact about the basis that equivariance needs. -/
lemma e_sub (j r : ZMod g) : χ.e (j - r) = χ.e j * (starRingEnd ℂ) (χ.e r) := by
  have h : χ.e j = χ.e (j - r) * χ.e r := by
    rw [← χ.map_add]; ring_nf
  calc χ.e (j - r) = χ.e (j - r) * (χ.e r * (starRingEnd ℂ) (χ.e r)) := by
        rw [χ.mul_star_self]; ring
    _ = χ.e j * (starRingEnd ℂ) (χ.e r) := by rw [h]; ring

end FundChar

/-- The fundamental logit vector `B η` of Section 5.3. The real and imaginary
parts of `η` are the cosine and sine coordinates.

The `√2` is the Fisher-orthonormal scaling **for `g ≥ 3`**, where it makes `‖η‖₂`
the Fisher norm of the induced perturbation. It is *not* the `g = 2` scaling: the
sign representation is one-dimensional and needs a factor of `1`, which is what
`fourier_design_matrix(2)` returns in the reference implementation. So at `g = 2`
this chart is `√2` times the implementation's, and `‖η‖₂` is not its Fisher norm.

Nothing below depends on the choice. Equivariance is invariant under a positive
rescaling of the logits, so `logit_rotate` and `equivariance` hold for either;
Fisher orthonormality, which is where the two differ, is not formalised here. -/
noncomputable def logit (χ : FundChar g) (η : ℂ) : Idx g → ℝ :=
  fun j => Real.sqrt 2 * (η * (starRingEnd ℂ) (χ.e j)).re

/-- **Linear equivariance.** Rotating the fundamental coordinate by the group
element `r` shifts the logit vector by `r`. This is `B (R_g^r η) = T_r (B η)`. -/
theorem logit_rotate (χ : FundChar g) (η : ℂ) (r : Idx g) :
    logit χ (χ.e r * η) = shift r (logit χ η) := by
  funext j
  simp only [logit, shift_apply, χ.e_sub, map_mul, starRingEnd_self_apply]
  congr 2
  ring

/-- **Equivariance of the chart (Section 5.3): `p(R_g^r η) = T_r p(η)`.**

Together with `Proposition 1`'s regularity hypothesis this is what makes the
shared-orbit model of Section 3.3 well posed: the right segment's distribution
is determined by the left segment's continuous coordinate and the discrete
relative label `r`, with no further continuous parameter. -/
theorem equivariance (χ : FundChar g) (η : ℂ) (r : Idx g) :
    softmax (logit χ (χ.e r * η)) = shift r (softmax (logit χ η)) := by
  rw [logit_rotate, softmax_shift]

end RegimeShift

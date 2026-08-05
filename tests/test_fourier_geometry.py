"""Structural tests for the cyclic Fourier geometry (Section 5).

These cover two of the six structural validations listed in Section 7.4:
Fisher orthonormality of the Fourier basis, and exact cyclic equivariance.
"""

from __future__ import annotations

import numpy as np
import pytest

from regimeshift.fourier import (
    fisher_inner_product,
    fourier_design_matrix,
    full_dimension,
    fundamental_dimension,
    fundamental_tangent_basis,
    higher_mode_logits,
    probabilities,
    rotation_matrix,
)

GROUPS = [2, 3, 4, 5, 6, 7, 8]


@pytest.mark.parametrize("m", GROUPS)
def test_fundamental_dimension_matches_representation_theory(m):
    # The sign representation at m = 2 is one-dimensional; every other
    # fundamental Fourier mode is a two-dimensional real representation.
    assert fundamental_dimension(m) == (1 if m == 2 else 2)
    assert full_dimension(m) == m - 1


@pytest.mark.parametrize("m", GROUPS)
def test_fisher_orthonormal_basis(m):
    basis = fundamental_tangent_basis(m)
    d = fundamental_dimension(m)
    assert basis.shape == (d, m)
    gram = np.array([[fisher_inner_product(u, v, m) for v in basis] for u in basis])
    np.testing.assert_allclose(gram, np.eye(d), atol=1e-12)


@pytest.mark.parametrize("m", GROUPS)
def test_basis_vectors_are_tangent_to_the_simplex(m):
    basis = fundamental_tangent_basis(m)
    np.testing.assert_allclose(basis.sum(axis=1), 0.0, atol=1e-14)


@pytest.mark.parametrize("m", GROUPS)
def test_softmax_derivative_equals_the_tangent_basis(m):
    """The logit scaling is chosen so d p / d theta at zero is the Fisher basis."""
    d = fundamental_dimension(m)
    eps = 1e-6
    for a in range(d):
        step = np.zeros(d)
        step[a] = eps
        numeric = (probabilities(step, m) - probabilities(-step, m)) / (2 * eps)
        np.testing.assert_allclose(numeric, fundamental_tangent_basis(m)[a], atol=1e-8)


@pytest.mark.parametrize("m", GROUPS)
@pytest.mark.parametrize("steps", [1, 2, 3])
def test_exact_cyclic_equivariance(m, steps):
    """p(R^s theta) = g^s p(theta), i.e. rotating coordinates permutes categories."""
    rng = np.random.default_rng(11 + m)
    for _ in range(20):
        theta = rng.normal(scale=0.6, size=fundamental_dimension(m))
        rotated = probabilities(rotation_matrix(m, steps) @ theta, m)
        shifted = np.roll(probabilities(theta, m), steps)
        np.testing.assert_allclose(rotated, shifted, atol=1e-13)


@pytest.mark.parametrize("m", GROUPS)
def test_rotation_by_full_period_is_the_identity(m):
    np.testing.assert_allclose(rotation_matrix(m, m), np.eye(fundamental_dimension(m)), atol=1e-13)


@pytest.mark.parametrize("m", GROUPS)
def test_fisher_norm_of_theta_is_preserved_by_the_group_action(m):
    rng = np.random.default_rng(5)
    theta = rng.normal(size=fundamental_dimension(m))
    for s in range(m):
        np.testing.assert_allclose(
            np.linalg.norm(rotation_matrix(m, s) @ theta), np.linalg.norm(theta), atol=1e-13
        )


@pytest.mark.parametrize("m", GROUPS)
def test_probabilities_are_a_valid_distribution(m):
    rng = np.random.default_rng(3)
    for _ in range(10):
        p = probabilities(rng.normal(scale=2.0, size=fundamental_dimension(m)), m)
        assert p.shape == (m,)
        assert np.all(p > 0)
        np.testing.assert_allclose(p.sum(), 1.0, atol=1e-14)


def test_zero_coordinate_gives_the_uniform_distribution():
    for m in GROUPS:
        np.testing.assert_allclose(probabilities(np.zeros(fundamental_dimension(m)), m), np.full(m, 1 / m))


@pytest.mark.parametrize("m", [4, 5, 6, 8])
def test_higher_mode_is_fisher_orthogonal_to_the_fundamental(m):
    """The misspecification direction must genuinely leave the fundamental subspace."""
    extra = higher_mode_logits(m, 1.0, mode=2)
    # Convert the logit perturbation to a tangent vector at the uniform state.
    tangent = (extra - extra.mean()) / m
    for basis_vector in fundamental_tangent_basis(m):
        assert abs(fisher_inner_product(tangent, basis_vector, m)) < 1e-12
    assert np.linalg.norm(tangent) > 0


def test_higher_mode_rejects_out_of_range_modes():
    with pytest.raises(ValueError):
        higher_mode_logits(4, 1.0, mode=4)
    with pytest.raises(ValueError):
        higher_mode_logits(4, 1.0, mode=0)


def test_invalid_group_orders_are_rejected():
    for bad in (1, 0, -3):
        with pytest.raises(ValueError):
            fourier_design_matrix(bad)


@pytest.mark.parametrize("m", GROUPS)
def test_local_jensen_shannon_coefficient(m):
    """Section 5.4: for small Fisher-norm perturbations the weighted JSD between
    p(theta) and p(R theta) is quadratic with the coefficient set by the rotation
    angle, with no extra factor of m."""
    from regimeshift.gains import weighted_jensen_shannon

    R = rotation_matrix(m, 1)
    d = fundamental_dimension(m)
    unit = np.eye(d)[0]
    # JSD(p, q) -> (1/8) ||delta||_F^2, and theta carries the Fisher norm, so the
    # coefficient is (1/8) ||R e - e||^2 = (1 - cos(2 pi / m)) / 4.
    predicted = 0.125 * np.linalg.norm(R @ unit - unit) ** 2
    np.testing.assert_allclose(predicted, (1 - np.cos(2 * np.pi / m)) / 4, atol=1e-14)

    for eps in (5e-3, 2.5e-3):
        theta = eps * unit
        jsd = weighted_jensen_shannon(probabilities(theta, m), probabilities(R @ theta, m))
        np.testing.assert_allclose(jsd / eps**2, predicted, rtol=3e-3)

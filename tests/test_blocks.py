"""Block (codon-phase) families, where the group acts on phases not the alphabet.

Section 11. The direct model identifies group order with alphabet size, so
raising ``m`` changes the number of candidate shifts, the simplex dimension, the
category sparsity and the orbit spacing all at once. Both external reviews named
that confound. These tests check that the block family separates them, and that
its dimension arithmetic reproduces Section 11's codon numbers exactly.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.optimize import approx_fprime

from regimeshift.blocks import (
    CODON,
    BlockGeometry,
    _contrast_basis,
    _neg_block_loglik_and_grad,
    block_full_detector,
    block_fundamental_detector,
    block_multinomial_loglik,
    block_probabilities,
    block_rotation,
    block_shared_orbit_detector,
    fit_block_fundamental,
    run_block_detectors,
)
from regimeshift.detectors import fit_failure_count, label_cost, reset_fit_failures
from regimeshift.fourier import fourier_design_matrix

GEOMETRIES = [(2, 4), (3, 4), (4, 4), (6, 4), (3, 8), (5, 3)]


def uniform_counts(geometry, n, rng):
    p = block_probabilities(
        np.zeros((geometry.mode_dimension, geometry.contrast_dimension)), geometry
    )
    return np.array([rng.multinomial(n, p[j]) for j in range(geometry.g)])


# --------------------------------------------------------------------------
# Section 11's arithmetic
# --------------------------------------------------------------------------


def test_codon_case_reproduces_section_eleven():
    """``g = 3`` reading frames, ``a = 4`` nucleotides: d_full = 3(4-1) = 9,
    dim V_fund = 2, d_phase_fund = 2(4-1) = 6."""
    assert (CODON.g, CODON.a) == (3, 4)
    assert CODON.full_dimension == 9
    assert CODON.mode_dimension == 2
    assert CODON.fundamental_dimension == 6
    assert CODON.contrast_dimension == 3


@pytest.mark.parametrize("g,a", GEOMETRIES)
def test_dimensions_are_the_product_of_phase_and_contrast_parts(g, a):
    geometry = BlockGeometry(g, a)
    assert geometry.full_dimension == g * (a - 1)
    assert geometry.fundamental_dimension == geometry.mode_dimension * (a - 1)
    assert geometry.mode_dimension == (1 if g == 2 else 2)


def test_group_order_and_alphabet_size_move_independently():
    """The point of the block family. At fixed alphabet, raising the group order
    raises the full dimension but leaves the fundamental one alone; raising the
    alphabet raises both. In the direct model these could never be separated."""
    fixed_alphabet = [BlockGeometry(g, 4) for g in (3, 4, 5, 6)]
    full = [geometry.full_dimension for geometry in fixed_alphabet]
    fundamental = [geometry.fundamental_dimension for geometry in fixed_alphabet]
    assert full == [9, 12, 15, 18]
    assert fundamental == [6, 6, 6, 6]

    fixed_group = [BlockGeometry(3, a) for a in (4, 6, 8)]
    assert [geometry.full_dimension for geometry in fixed_group] == [9, 15, 21]
    assert [geometry.fundamental_dimension for geometry in fixed_group] == [6, 10, 14]


def test_invalid_geometries_are_rejected():
    for g, a in ((1, 4), (0, 4), (3, 1), (3, 0)):
        with pytest.raises(ValueError):
            BlockGeometry(g, a)


# --------------------------------------------------------------------------
# geometry
# --------------------------------------------------------------------------


@pytest.mark.parametrize("a", [2, 3, 4, 8])
def test_contrast_basis_is_orthonormal_and_mean_zero(a):
    basis = _contrast_basis(a)
    assert basis.shape == (a - 1, a)
    np.testing.assert_allclose(basis @ basis.T, np.eye(a - 1), atol=1e-13)
    np.testing.assert_allclose(basis.sum(axis=1), 0.0, atol=1e-13)


@pytest.mark.parametrize("g,a", GEOMETRIES)
def test_probabilities_are_valid_per_block(g, a):
    geometry = BlockGeometry(g, a)
    rng = np.random.default_rng(g * a)
    theta = rng.normal(scale=0.5, size=(geometry.mode_dimension, geometry.contrast_dimension))
    p = block_probabilities(theta, geometry)
    assert p.shape == (g, a)
    assert np.all(p > 0)
    np.testing.assert_allclose(p.sum(axis=1), 1.0, atol=1e-13)


@pytest.mark.parametrize("g,a", GEOMETRIES)
def test_zero_coordinates_give_uniform_blocks(g, a):
    geometry = BlockGeometry(g, a)
    p = block_probabilities(np.zeros((geometry.mode_dimension, geometry.contrast_dimension)), geometry)
    np.testing.assert_allclose(p, 1.0 / a, atol=1e-14)


@pytest.mark.parametrize("g,a", GEOMETRIES)
def test_phase_shift_equals_coordinate_rotation(g, a):
    """The defining equivariance: permuting phase blocks acts on the fundamental
    coordinates as ``R (x) I`` -- a rotation of the phase-mode index, identical in
    every alphabet-contrast direction. The alphabet is untouched."""
    geometry = BlockGeometry(g, a)
    rng = np.random.default_rng(g + a)
    for _ in range(10):
        theta = rng.normal(scale=0.4, size=(geometry.mode_dimension, geometry.contrast_dimension))
        for shift in range(1, g):
            rotated = block_probabilities(block_rotation(theta, geometry, shift), geometry)
            shifted = np.roll(block_probabilities(theta, geometry), shift, axis=0)
            np.testing.assert_allclose(rotated, shifted, atol=1e-13)


@pytest.mark.parametrize("g,a", [(3, 4), (6, 4), (3, 8)])
def test_block_loglik_gradient_is_correct(g, a):
    geometry = BlockGeometry(g, a)
    rng = np.random.default_rng(g * a + 1)
    counts = uniform_counts(geometry, 300, rng).astype(float)
    design = fourier_design_matrix(g)
    contrasts = _contrast_basis(a)
    shape = (geometry.mode_dimension, geometry.contrast_dimension)
    for _ in range(6):
        x = rng.normal(scale=0.3, size=int(np.prod(shape)))
        _, analytic = _neg_block_loglik_and_grad(x, counts, design, contrasts, shape)
        numeric = approx_fprime(
            x, lambda z: _neg_block_loglik_and_grad(z, counts, design, contrasts, shape)[0], 1e-7
        )
        np.testing.assert_allclose(analytic, numeric, atol=1e-3, rtol=1e-4)


@pytest.mark.parametrize("g,a", GEOMETRIES)
def test_fit_recovers_a_planted_coordinate(g, a):
    geometry = BlockGeometry(g, a)
    rng = np.random.default_rng(g * a)
    theta = rng.normal(scale=0.25, size=(geometry.mode_dimension, geometry.contrast_dimension))
    p = block_probabilities(theta, geometry)
    counts = np.array([rng.multinomial(200_000, p[j]) for j in range(g)])
    fitted, _ = fit_block_fundamental(counts, geometry)
    np.testing.assert_allclose(fitted, theta, atol=0.02)


# --------------------------------------------------------------------------
# the detectors, and the dimensions they actually pay
# --------------------------------------------------------------------------


@pytest.mark.parametrize("g,a", GEOMETRIES)
def test_detectors_report_the_intended_increments(g, a):
    geometry = BlockGeometry(g, a)
    rng = np.random.default_rng(g + a)
    cL = uniform_counts(geometry, 300, rng)
    cR = uniform_counts(geometry, 300, rng)
    results = run_block_detectors(cL, cR, geometry)
    assert results["block_full"].dimension_increment == geometry.full_dimension
    assert results["block_fundamental"].dimension_increment == geometry.fundamental_dimension
    assert results["block_shared_orbit"].dimension_increment == 0.0
    assert results["block_shared_orbit"].penalty == pytest.approx(label_cost(g))


@pytest.mark.parametrize("g,a,expected_full,expected_fund", [
    (3, 4, 4.5, 3.0),
    (6, 4, 9.0, 3.0),
    (3, 8, 10.5, 7.0),
])
def test_null_gains_recover_the_dimensions_empirically(g, a, expected_full, expected_fund):
    """Under the null a regular ``d``-dimensional split has ``2 x gain ~ chi^2_d``,
    so the mean raw gain is ``d/2``. This measures the dimension arithmetic
    rather than restating it -- and shows the separation directly: at ``a = 4``,
    going from ``g = 3`` to ``g = 6`` doubles the full model's mean gain from 4.5
    to 9.0 while the fundamental model's stays at 3.0.
    """
    geometry = BlockGeometry(g, a)
    rng = np.random.default_rng(20260713 + g * 100 + a)
    full, fundamental = [], []
    for _ in range(500):
        cL = uniform_counts(geometry, 400, rng)
        cR = uniform_counts(geometry, 400, rng)
        full.append(block_full_detector(cL, cR, geometry).raw_gain)
        fundamental.append(block_fundamental_detector(cL, cR, geometry).raw_gain)
    assert np.mean(full) == pytest.approx(expected_full, rel=0.12)
    assert np.mean(fundamental) == pytest.approx(expected_fund, rel=0.12)


@pytest.mark.parametrize("g,a", [(3, 4), (4, 5), (6, 4)])
def test_shared_orbit_recovers_a_planted_phase_shift(g, a):
    """A frameshift, in the codon reading."""
    geometry = BlockGeometry(g, a)
    rng = np.random.default_rng(g * a + 7)
    theta = rng.normal(scale=0.35, size=(geometry.mode_dimension, geometry.contrast_dimension))
    for shift in range(1, g):
        p_left = block_probabilities(theta, geometry)
        p_right = np.roll(p_left, shift, axis=0)
        hits = 0
        for _ in range(15):
            cL = np.array([rng.multinomial(3000, p_left[j]) for j in range(g)])
            cR = np.array([rng.multinomial(3000, p_right[j]) for j in range(g)])
            hits += block_shared_orbit_detector(cL, cR, geometry).selected_shift == shift
        assert hits >= 14, f"g={g} a={a} shift={shift}: {hits}/15"


@pytest.mark.parametrize("g,a", [(3, 4), (6, 4)])
def test_shared_orbit_beats_the_alternatives_on_a_phase_shift(g, a):
    """The block analogue of the paper's central comparison: on an exact phase
    orbit the relational code wins, and by more as the group order grows."""
    geometry = BlockGeometry(g, a)
    rng = np.random.default_rng(g)
    theta = np.zeros((geometry.mode_dimension, geometry.contrast_dimension))
    theta[0, 0] = 0.3
    p_left = block_probabilities(theta, geometry)
    p_right = np.roll(p_left, 1, axis=0)
    wins = 0
    for _ in range(25):
        cL = np.array([rng.multinomial(800, p_left[j]) for j in range(g)])
        cR = np.array([rng.multinomial(800, p_right[j]) for j in range(g)])
        results = run_block_detectors(cL, cR, geometry)
        wins += results["block_shared_orbit"].score > max(
            results["block_full"].score, results["block_fundamental"].score
        )
    assert wins >= 23, f"{wins}/25"


@pytest.mark.parametrize("g,a", GEOMETRIES)
def test_block_multinomial_loglik_matches_per_block_frequencies(g, a):
    geometry = BlockGeometry(g, a)
    rng = np.random.default_rng(g * a + 3)
    counts = uniform_counts(geometry, 500, rng).astype(float)
    expected = 0.0
    for row in counts:
        p = row / row.sum()
        nz = row > 0
        expected += float(np.sum(row[nz] * np.log(p[nz])))
    assert block_multinomial_loglik(counts) == pytest.approx(expected)


def test_block_loglik_handles_empty_blocks():
    counts = np.array([[0.0, 0, 0, 0], [5, 5, 0, 0], [1, 2, 3, 4]])
    assert np.isfinite(block_multinomial_loglik(counts))


@pytest.mark.parametrize("g,a", GEOMETRIES)
def test_fits_converge_across_the_geometries(g, a):
    geometry = BlockGeometry(g, a)
    reset_fit_failures()
    rng = np.random.default_rng(g * a)
    for n in (20, 200, 5000):
        for _ in range(10):
            fit_block_fundamental(uniform_counts(geometry, n, rng), geometry)
    assert fit_failure_count() == 0


def test_wrong_shaped_counts_are_rejected():
    with pytest.raises(ValueError):
        fit_block_fundamental(np.ones((2, 4)), CODON)
    with pytest.raises(ValueError):
        block_probabilities(np.zeros((1, 3)), CODON)


@pytest.mark.parametrize(
    "detector",
    [block_full_detector, block_fundamental_detector, block_shared_orbit_detector],
)
def test_block_detectors_reject_counts_of_the_wrong_geometry(detector):
    """A mismatched geometry must raise, not return a score built from two of them.

    ``block_full_detector`` reads its likelihood off the array it is handed and
    its penalty off ``geometry``, so before ``validate_block_pair`` a (3, 4)
    array scored against ``BlockGeometry(6, 4)`` returned a plausible-looking
    number with dimension increment 18 instead of 9.
    """
    counts = np.full((3, 4), 25.0)
    with pytest.raises(ValueError):
        detector(counts, counts, BlockGeometry(g=6, a=4))


@pytest.mark.parametrize(
    "detector",
    [block_full_detector, block_fundamental_detector, block_shared_orbit_detector],
)
def test_block_detectors_reject_empty_or_negative_segments(detector):
    good = np.full(CODON.shape, 25.0)
    with pytest.raises(ValueError):
        detector(good, np.zeros(CODON.shape), CODON)
    bad = good.copy()
    bad[0, 0] = -1.0
    with pytest.raises(ValueError):
        detector(good, bad, CODON)

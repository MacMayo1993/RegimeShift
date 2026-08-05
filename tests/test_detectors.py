"""Structural tests for the three detectors (Sections 3, 4 and 7).

Together with ``test_fourier_geometry.py`` these cover the six structural
validations of Section 7.4:

1. Fisher orthonormality of the Fourier basis        (test_fourier_geometry)
2. exact cyclic equivariance                         (test_fourier_geometry)
3. the intended continuous-dimension increments      (here)
4. equality of Models A and B for m = 2 and m = 3    (here)
5. recovery of planted relative shifts in smoke data (here)
6. absence of a sample-length-dependent penalty in C  (here)
"""

from __future__ import annotations

import numpy as np
import pytest

from regimeshift.detectors import (
    fit_fundamental,
    full_detector,
    fundamental_detector,
    fundamental_loglik,
    label_cost,
    multinomial_loglik,
    run_all_detectors,
    shared_orbit_detector,
    split_penalty,
)
from regimeshift.fourier import fundamental_dimension, probabilities, rotation_matrix
from regimeshift.scenarios import build_segments

GROUPS = [2, 3, 4, 5, 6]


def counts_from(p, n, rng):
    return rng.multinomial(n, p)


# --------------------------------------------------------------------------
# 3. continuous-dimension increments
# --------------------------------------------------------------------------


@pytest.mark.parametrize("m", GROUPS)
def test_dimension_increments_match_the_three_model_classes(m):
    rng = np.random.default_rng(m)
    cL = counts_from(np.full(m, 1 / m), 500, rng)
    cR = counts_from(np.full(m, 1 / m), 500, rng)
    results = run_all_detectors(cL, cR, m)
    assert results["full"].dimension_increment == m - 1
    assert results["fundamental"].dimension_increment == fundamental_dimension(m)
    assert results["shared_orbit"].dimension_increment == 0.0


@pytest.mark.parametrize("m", GROUPS)
def test_penalty_ordering_is_a_strict_hierarchy(m):
    """A > B >= C in complexity, with equality of A and B exactly when the
    fundamental component spans the whole nontrivial tangent space."""
    rng = np.random.default_rng(m + 100)
    cL = counts_from(np.full(m, 1 / m), 1000, rng)
    cR = counts_from(np.full(m, 1 / m), 1000, rng)
    results = run_all_detectors(cL, cR, m)
    if m <= 3:
        assert results["full"].penalty == pytest.approx(results["fundamental"].penalty)
    else:
        assert results["full"].penalty > results["fundamental"].penalty
    assert results["fundamental"].penalty > results["shared_orbit"].penalty


def test_split_penalty_leading_coefficient_is_half_the_dimension():
    """(d/2) log n is the leading term; the split fraction only shifts the
    bounded remainder."""
    for dim in (1, 2, 5):
        for fraction in (0.5, 0.25, 0.1):
            slopes = []
            for n in (10**5, 10**6, 10**7):
                nL = int(n * fraction)
                slopes.append(split_penalty(dim, nL, n - nL))
            # Successive increments over a factor-of-10 length increase.
            diffs = np.diff(slopes) / np.log(10)
            np.testing.assert_allclose(diffs, dim / 2, rtol=1e-6)


def test_split_penalty_rejects_empty_segments():
    for bad in ((0, 10), (10, 0)):
        with pytest.raises(ValueError):
            split_penalty(2, *bad)


@pytest.mark.parametrize("m", GROUPS)
def test_label_cost_matches_the_number_of_nonidentity_shifts(m):
    assert label_cost(m) == pytest.approx(np.log(m - 1))
    if m == 2:
        assert label_cost(m) == 0.0


# --------------------------------------------------------------------------
# 4. Models A and B coincide when the fundamental spans the tangent space
# --------------------------------------------------------------------------


@pytest.mark.parametrize("m", [2, 3])
def test_full_and_fundamental_coincide_for_small_groups(m):
    """For m = 2 and m = 3, d_fund = m - 1, so the two model spaces are equal
    and the detectors must agree in gain, penalty and score."""
    rng = np.random.default_rng(7)
    for _ in range(10):
        p_left = rng.dirichlet(np.ones(m))
        p_right = rng.dirichlet(np.ones(m))
        cL = counts_from(p_left, 400, rng)
        cR = counts_from(p_right, 400, rng)
        a = full_detector(cL, cR, m)
        b = fundamental_detector(cL, cR, m)
        assert a.raw_gain == pytest.approx(b.raw_gain, abs=1e-8)
        assert a.penalty == pytest.approx(b.penalty)
        assert a.score == pytest.approx(b.score, abs=1e-8)


@pytest.mark.parametrize("m", [4, 5, 6])
def test_maximised_likelihoods_are_nested(m):
    """A's alternative contains B's, which contains C's; and A's null contains
    the common null of B and C. Note that this does *not* order the raw gains of
    A and B, because A also raises its own null."""
    rng = np.random.default_rng(m * 13)
    for _ in range(10):
        cL = counts_from(rng.dirichlet(np.ones(m)), 500, rng)
        cR = counts_from(rng.dirichlet(np.ones(m)), 500, rng)

        null_full = multinomial_loglik(cL + cR)
        _, null_constrained = fit_fundamental(cL + cR, m)
        assert null_full >= null_constrained - 1e-8

        alt_full = multinomial_loglik(cL) + multinomial_loglik(cR)
        _, ll_left = fit_fundamental(cL, m)
        _, ll_right = fit_fundamental(cR, m)
        alt_fundamental = ll_left + ll_right
        alt_shared = shared_orbit_detector(cL, cR, m).raw_gain + null_constrained

        assert alt_full >= alt_fundamental - 1e-8
        assert alt_fundamental >= alt_shared - 1e-8


@pytest.mark.parametrize("m", GROUPS)
def test_regular_gains_are_non_negative(m):
    """Models A and B each nest their own null, so their maximised gains cannot
    be negative."""
    rng = np.random.default_rng(m + 71)
    for _ in range(15):
        cL = counts_from(np.full(m, 1 / m), 300, rng)
        cR = counts_from(np.full(m, 1 / m), 300, rng)
        results = run_all_detectors(cL, cR, m)
        assert results["full"].raw_gain >= -1e-7
        assert results["fundamental"].raw_gain >= -1e-7


@pytest.mark.parametrize("m", GROUPS)
def test_shared_orbit_gain_may_be_negative(m):
    """Model C's alternative ranges over *nonidentity* shifts only, so it does
    not nest its own null: on no-change data the best aligned pooling is
    typically worse than the unaligned pooling. This is a property of the
    hypothesis, not a defect -- it is why Model C pays only a discrete label
    cost and no continuous-dimension increment."""
    rng = np.random.default_rng(m + 71)
    gains = []
    for _ in range(30):
        cL = counts_from(np.full(m, 1 / m), 300, rng)
        cR = counts_from(np.full(m, 1 / m), 300, rng)
        gains.append(shared_orbit_detector(cL, cR, m).raw_gain)
    assert min(gains) < 0
    assert np.isfinite(gains).all()


@pytest.mark.parametrize("m", [4, 5, 6])
def test_shared_orbit_gain_never_exceeds_the_fundamental_gain(m):
    """B and C share a null and C's alternative is contained in B's, so this
    ordering *is* guaranteed."""
    rng = np.random.default_rng(m * 17)
    for _ in range(10):
        cL = counts_from(rng.dirichlet(np.ones(m)), 400, rng)
        cR = counts_from(rng.dirichlet(np.ones(m)), 400, rng)
        results = run_all_detectors(cL, cR, m)
        assert results["fundamental"].raw_gain >= results["shared_orbit"].raw_gain - 1e-8


# --------------------------------------------------------------------------
# fitting routines
# --------------------------------------------------------------------------


@pytest.mark.parametrize("m", GROUPS)
def test_fundamental_fit_recovers_the_generating_coordinate(m):
    """With a large sample the MLE converges to the planted coordinate."""
    rng = np.random.default_rng(m)
    theta = rng.normal(scale=0.3, size=fundamental_dimension(m))
    counts = counts_from(probabilities(theta, m), 400_000, rng)
    theta_hat, _ = fit_fundamental(counts, m)
    np.testing.assert_allclose(theta_hat, theta, atol=0.03)


@pytest.mark.parametrize("m", GROUPS)
def test_fundamental_fit_beats_perturbations_of_itself(m):
    """A first-order optimality check on the fitted coordinate."""
    rng = np.random.default_rng(m + 4)
    counts = counts_from(rng.dirichlet(np.ones(m)), 2000, rng)
    theta_hat, ll_hat = fit_fundamental(counts, m)
    for _ in range(20):
        step = rng.normal(scale=0.05, size=theta_hat.shape)
        assert fundamental_loglik(theta_hat + step, counts, m) <= ll_hat + 1e-7


@pytest.mark.parametrize("m", GROUPS)
def test_fundamental_fit_is_equivariant(m):
    """Rolling the counts rotates the fitted coordinate and leaves the
    likelihood unchanged."""
    rng = np.random.default_rng(m + 30)
    counts = counts_from(rng.dirichlet(np.ones(m)), 3000, rng)
    theta_hat, ll = fit_fundamental(counts, m)
    for s in range(1, m):
        theta_s, ll_s = fit_fundamental(np.roll(counts, s), m)
        assert ll_s == pytest.approx(ll, abs=1e-6)
        np.testing.assert_allclose(theta_s, rotation_matrix(m, s) @ theta_hat, atol=1e-4)


@pytest.mark.parametrize("m", GROUPS)
def test_multinomial_loglik_matches_the_empirical_frequencies(m):
    rng = np.random.default_rng(2)
    counts = counts_from(rng.dirichlet(np.ones(m)), 1000, rng)
    p_hat = counts / counts.sum()
    expected = float(np.sum(counts[counts > 0] * np.log(p_hat[counts > 0])))
    assert multinomial_loglik(counts) == pytest.approx(expected)


def test_multinomial_loglik_handles_zero_counts():
    assert multinomial_loglik(np.array([0.0, 0.0])) == 0.0
    assert np.isfinite(multinomial_loglik(np.array([10.0, 0.0, 5.0])))


# --------------------------------------------------------------------------
# 5. planted relative-shift recovery
# --------------------------------------------------------------------------


@pytest.mark.parametrize("m", [3, 4, 5, 6])
@pytest.mark.parametrize("shift", [1, 2])
def test_shared_orbit_recovers_a_planted_shift(m, shift):
    rng = np.random.default_rng(m * 100 + shift)
    theta = np.array([0.25, 0.0])[: fundamental_dimension(m)]
    p_left = probabilities(theta, m)
    p_right = probabilities(rotation_matrix(m, shift) @ theta, m)
    hits = 0
    trials = 25
    for _ in range(trials):
        cL = counts_from(p_left, 4000, rng)
        cR = counts_from(p_right, 4000, rng)
        hits += shared_orbit_detector(cL, cR, m).selected_shift == shift
    assert hits >= trials - 1


@pytest.mark.parametrize("m", [4, 5, 6])
def test_shared_orbit_recovers_the_scenario_shift(m):
    """Smoke data straight from the exact-orbit scenario."""
    segments = build_segments(m, "exact_orbit", 0.25)
    rng = np.random.default_rng(m)
    hits = [
        shared_orbit_detector(
            counts_from(segments.p_left, 3200, rng), counts_from(segments.p_right, 3200, rng), m
        ).selected_shift
        == segments.planted_shift
        for _ in range(30)
    ]
    assert np.mean(hits) >= 0.9


@pytest.mark.parametrize("m", [3, 4, 5])
def test_shared_orbit_never_selects_the_identity(m):
    rng = np.random.default_rng(9)
    result = shared_orbit_detector(
        counts_from(np.full(m, 1 / m), 200, rng), counts_from(np.full(m, 1 / m), 200, rng), m
    )
    assert result.selected_shift in range(1, m)


# --------------------------------------------------------------------------
# 6. the Model C penalty does not depend on sample length
# --------------------------------------------------------------------------


@pytest.mark.parametrize("m", GROUPS)
def test_shared_orbit_penalty_is_constant_in_sample_length(m):
    rng = np.random.default_rng(m + 500)
    penalties = set()
    for n in (50, 500, 5_000, 50_000):
        cL = counts_from(np.full(m, 1 / m), n, rng)
        cR = counts_from(np.full(m, 1 / m), n, rng)
        penalties.add(round(shared_orbit_detector(cL, cR, m).penalty, 12))
    assert len(penalties) == 1
    assert penalties.pop() == pytest.approx(label_cost(m))


@pytest.mark.parametrize("m", [4, 5, 6])
def test_regular_penalties_do_grow_with_sample_length(m):
    """The contrast that makes the previous test meaningful."""
    rng = np.random.default_rng(m)
    for detector, dim in ((full_detector, m - 1), (fundamental_detector, fundamental_dimension(m))):
        penalties = []
        for n in (1_000, 10_000, 100_000):
            cL = counts_from(np.full(m, 1 / m), n, rng)
            cR = counts_from(np.full(m, 1 / m), n, rng)
            penalties.append(detector(cL, cR, m).penalty)
        np.testing.assert_allclose(np.diff(penalties) / np.log(10), dim / 2, rtol=1e-6)


@pytest.mark.parametrize("m", GROUPS)
def test_scores_are_finite_at_extreme_counts(m):
    """Degenerate segments (all mass in one category) must not produce NaNs."""
    cL = np.zeros(m)
    cL[0] = 100
    cR = np.zeros(m)
    cR[-1] = 100
    for result in run_all_detectors(cL, cR, m).values():
        assert np.isfinite(result.score)
        assert np.isfinite(result.raw_gain)

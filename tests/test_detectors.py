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
    approximate_orbit_detector,
    fit_fundamental,
    full_detector,
    fundamental_detector,
    fundamental_loglik,
    fundamental_mle_exists,
    label_cost,
    multinomial_loglik,
    run_all_detectors,
    shared_orbit_detector,
    split_penalty,
    validate_pair,
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


# --------------------------------------------------------------------------
# Uniform input validation
# --------------------------------------------------------------------------

ALL_DETECTORS = (
    full_detector,
    fundamental_detector,
    shared_orbit_detector,
    approximate_orbit_detector,
)


@pytest.mark.parametrize("detector", ALL_DETECTORS)
def test_detectors_reject_alphabet_mismatch(detector):
    """A vector whose length disagrees with ``m`` must raise, not score.

    Without the check the likelihood is computed from the vector length and the
    penalty from ``m``, so the detector returns a score assembled from two
    different alphabet sizes.
    """
    m = 4
    cL = np.full(m + 1, 10.0)
    cR = np.full(m + 1, 10.0)
    with pytest.raises(ValueError, match="shape"):
        detector(cL, cR, m)


@pytest.mark.parametrize("detector", ALL_DETECTORS)
def test_detectors_reject_empty_segment(detector):
    m = 4
    cL = np.full(m, 10.0)
    cR = np.zeros(m)
    with pytest.raises(ValueError, match="nonempty"):
        detector(cL, cR, m)


@pytest.mark.parametrize("detector", ALL_DETECTORS)
def test_detectors_reject_negative_and_nonfinite(detector):
    m = 4
    good = np.full(m, 10.0)
    bad_negative = good.copy()
    bad_negative[0] = -1.0
    with pytest.raises(ValueError, match="nonnegative"):
        detector(bad_negative, good, m)

    bad_nan = good.copy()
    bad_nan[0] = np.nan
    with pytest.raises(ValueError, match="finite"):
        detector(good, bad_nan, m)


def test_validate_pair_returns_float_arrays():
    m = 3
    cL, cR = validate_pair([1, 2, 3], (4, 5, 6), m)
    assert cL.dtype == float and cR.dtype == float
    np.testing.assert_allclose(cL, [1.0, 2.0, 3.0])


def test_split_penalty_is_negative_at_tiny_sizes():
    """The asymptotic formula is not a codelength at very small ``n``.

    Documented behaviour, asserted so it cannot change silently: with one
    observation per segment the increment is ``-(dim/2) log 2``.
    """
    assert split_penalty(2, 1, 1) == pytest.approx(-np.log(2.0))


# ---------------------------------------------------------------------------
# Existence of the fundamental MLE (Section 7.5)
# ---------------------------------------------------------------------------


def _t_bar_is_interior(counts: np.ndarray, m: int, tol: float = 1e-12) -> bool:
    """Barndorff-Nielsen condition, computed from the convex hull directly.

    The reference implementation of the criterion :func:`fundamental_mle_exists`
    shortcuts. It builds ``t_bar = sum_j f_j B_j`` and asks whether it is
    strictly inside ``conv{B_j}``, with no appeal to the polygon's face
    structure -- so agreeing with it is a real check on the shortcut.
    """
    from regimeshift.fourier import fourier_design_matrix

    B = fourier_design_matrix(m)
    f = np.asarray(counts, dtype=float)
    f = f / f.sum()
    t = B.T @ f
    if m == 2:
        return abs(t[0]) < np.abs(B[:, 0]).max() - tol
    for j in range(m):
        a, b = B[j], B[(j + 1) % m]
        edge = b - a
        normal = np.array([-edge[1], edge[0]])
        if normal @ (-a) > 0:  # orient outward, away from the centre
            normal = -normal
        if normal @ (t - a) >= -tol * np.linalg.norm(normal):
            return False
    return True


@pytest.mark.parametrize("m", [2, 3, 4, 5, 6, 7])
def test_mle_existence_rule_matches_the_convex_support_condition(m):
    """The O(m) adjacency rule agrees with the convex-hull condition exactly.

    Checked over every support pattern at three weightings, which is the whole
    space the rule can be wrong on: the criterion depends on the support and,
    through ``t_bar``, on the weights.
    """
    import itertools

    for k in range(1, m + 1):
        for support in itertools.combinations(range(m), k):
            for weights in ([1] * k, list(range(1, k + 1)), [1] + [7] * (k - 1)):
                counts = np.zeros(m)
                counts[list(support)] = np.array(weights, dtype=float) * 10.0
                assert fundamental_mle_exists(counts, m) == _t_bar_is_interior(counts, m), (
                    f"m={m} counts={counts}"
                )


@pytest.mark.parametrize("m", [3, 4, 5, 6, 7])
def test_mle_existence_rule_matches_on_random_counts(m):
    rng = np.random.default_rng(20260713 + m)
    for _ in range(400):
        counts = rng.integers(0, 6, size=m).astype(float)
        if counts.sum() == 0:
            continue
        assert fundamental_mle_exists(counts, m) == _t_bar_is_interior(counts, m)


def test_a_zero_cell_does_not_by_itself_break_the_mle():
    """The manuscript's Section 7.5 used to say a zero count implies the MLE
    does not exist. It does not, for ``m >= 3``.

    ``[0, 10, 10, 10]`` has an empty cell and an ordinary interior optimum:
    driving the first cell's probability to zero would drive its *opposite*
    cell up, and that cell has positive count, so the likelihood turns over.
    """
    counts = np.array([0.0, 10.0, 10.0, 10.0])
    assert fundamental_mle_exists(counts, 4)

    theta_hat, ll_hat = fit_fundamental(counts, 4)
    assert np.linalg.norm(theta_hat) < 1.0
    assert np.isfinite(ll_hat)
    for direction in np.eye(2):
        assert fundamental_loglik(theta_hat + direction, counts, 4) < ll_hat
        assert fundamental_loglik(theta_hat - direction, counts, 4) < ll_hat


def test_two_empty_opposite_cells_still_leave_a_finite_optimum():
    """Two zero cells, and the MLE is not merely finite but exactly uniform:
    opposite vertices of the square span no face of it."""
    counts = np.array([10.0, 0.0, 10.0, 0.0])
    assert fundamental_mle_exists(counts, 4)
    theta_hat, ll_hat = fit_fundamental(counts, 4)
    np.testing.assert_allclose(theta_hat, np.zeros(2), atol=1e-8)
    assert np.isfinite(ll_hat)


def test_two_adjacent_cells_are_the_case_that_genuinely_fails():
    """All mass on a cyclically adjacent pair puts ``t_bar`` on an edge.

    The supremum is finite but unattained, so the likelihood is *flat* along
    the escape direction -- which is why an optimiser stops at a large, and
    start-dependent, coordinate while reporting a sensible log-likelihood.
    """
    counts = np.array([10.0, 10.0, 0.0, 0.0])
    assert not fundamental_mle_exists(counts, 4)

    theta_hat, ll_hat = fit_fundamental(counts, 4)
    assert np.linalg.norm(theta_hat) > 5.0
    # flat, not decreasing: going four times further costs nothing
    assert fundamental_loglik(4.0 * theta_hat, counts, 4) == pytest.approx(ll_hat, abs=1e-9)


def test_at_m_two_a_zero_cell_is_exactly_the_failure_condition():
    """At ``m = 2`` the hull is a segment and both endpoints are faces, so the
    blanket claim happens to be right -- which is presumably where it came
    from. It does not generalise."""
    assert not fundamental_mle_exists(np.array([0.0, 10.0]), 2)
    assert not fundamental_mle_exists(np.array([10.0, 0.0]), 2)
    assert fundamental_mle_exists(np.array([1.0, 10.0]), 2)


def test_zero_cells_are_not_routine_on_the_production_grid():
    """Section 7.5 called zero counts "routine on short segments". They are not.

    Computed over the actual design rather than a worst-case corner: every
    configuration's segment probabilities and its own segment lengths, weighted
    by how many datasets it draws. A union bound on P(some cell empty) is an
    overestimate and still lands at a fraction of one segment across the whole
    run.

    The stronger statement is about the criterion rather than the counts. Even
    at the worst corner a *genuine* failure needs the whole segment to land in
    at most two cyclically adjacent cells, and no adjacent pair there carries
    enough mass for that to be reachable at any length on the grid.
    """
    from regimeshift.runner import PRODUCTION_GRID
    from regimeshift.scenarios import build_segments
    from regimeshift.simulation import build_grid

    expected = 0.0
    segments_scored = 0
    worst = 0.0
    smallest_cell = 1.0
    for config in build_grid(**PRODUCTION_GRID):
        segments = build_segments(config.m, config.scenario, config.effect)
        drawn = config.n_alt + config.n_null
        for p, n in ((segments.p_left, config.n_left), (segments.p_right, config.n_right)):
            smallest_cell = min(smallest_cell, float(p.min()))
            p_empty = float(np.sum((1.0 - p) ** n))
            worst = max(worst, p_empty)
            expected += p_empty * drawn
            segments_scored += drawn

    assert segments_scored == 936_000
    assert smallest_cell == pytest.approx(0.091, abs=0.001)  # quoted in Section 7.5
    assert worst == pytest.approx(7.2e-5, rel=0.05)          # worst corner, n = 100
    assert expected < 1.0                                    # about 0.2 segments in 936,000

    # and the worst corner cannot break the criterion at any length on the grid
    p = build_segments(6, "higher_mode", 0.25).p_right
    heaviest_adjacent_pair = max(p[j] + p[(j + 1) % 6] for j in range(6))
    assert heaviest_adjacent_pair**100 < 1e-30

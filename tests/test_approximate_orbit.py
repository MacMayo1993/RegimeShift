"""Model D: the approximate-orbit interpolation of Section 14.1.

    eta_R = R^r eta_L + delta,   delta ~ N(0, tau^2 I)

with ``tau = 0`` pinning the deviation at zero (Model C, exact orbit) and large
``tau`` leaving it free (Model B, independent subspace). The point of the model
is to answer "how much deviation from exact symmetry can be tolerated before the
relational advantage disappears" rather than treating symmetry as all-or-nothing.

These tests pin the two nesting limits, the shape of the deviation code, and the
behaviour across a deviation sweep.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.optimize import approx_fprime

from regimeshift.detectors import (
    _neg_joint_and_grad,
    approximate_orbit_detector,
    deviation_penalty,
    fit_approximate_orbit,
    fit_failure_count,
    fundamental_detector,
    label_cost,
    reset_fit_failures,
    shared_orbit_detector,
    split_penalty,
)
from regimeshift.fourier import fourier_design_matrix, fundamental_dimension, rotation_matrix
from regimeshift.scenarios import build_segments

GROUPS = [2, 3, 4, 5, 6]


def counts(p, n, rng):
    return rng.multinomial(n, p)


# --------------------------------------------------------------------------
# the deviation code
# --------------------------------------------------------------------------


def test_deviation_penalty_vanishes_at_zero_scale():
    for dim in (1, 2):
        for n in (10, 1000, 10**6):
            assert deviation_penalty(dim, n, n, 0.0) == 0.0


@pytest.mark.parametrize("dim", [1, 2])
def test_deviation_penalty_is_increasing_in_scale_and_length(dim):
    scales = [0.0, 0.01, 0.1, 1.0, 10.0]
    values = [deviation_penalty(dim, 1000, 1000, s) for s in scales]
    assert values == sorted(values)
    assert all(a < b for a, b in zip(values, values[1:]))

    lengths = [10, 100, 1000, 10_000]
    by_length = [deviation_penalty(dim, n, n, 0.1) for n in lengths]
    assert by_length == sorted(by_length)


@pytest.mark.parametrize("dim", [1, 2])
def test_deviation_penalty_approaches_the_regular_rate_for_a_wide_prior(dim):
    """With a fixed positive scale the leading coefficient is ``dim/2`` -- Model
    B's rate, not something in between. The interpolation is in the bounded
    term. This is the honest asymptotic statement and the tests should say so."""
    slopes = []
    for n in (10**5, 10**6, 10**7):
        slopes.append(deviation_penalty(dim, n, n, 1.0))
    increments = np.diff(slopes) / np.log(10)
    np.testing.assert_allclose(increments, dim / 2, rtol=1e-4)


def test_deviation_penalty_rejects_invalid_arguments():
    with pytest.raises(ValueError):
        deviation_penalty(2, 100, 100, -0.1)
    with pytest.raises(ValueError):
        deviation_penalty(2, 100, 0, 0.1)


@pytest.mark.parametrize("dim", [1, 2])
def test_deviation_penalty_uses_the_profiled_information(dim):
    """The shared state is estimated, not known, so the information that
    constrains ``delta`` is the Schur complement ``L1 L2 / (L1 + L2)`` -- not
    ``L2``. On a balanced split that halves the effective information, and the
    gap to the known-shared-state value rises to ``(dim/2) log 2`` as the prior
    stops doing the constraining."""
    tau = 0.05
    for half in (500, 5000, 50_000):
        profiled = deviation_penalty(dim, half, half, tau)
        np.testing.assert_allclose(
            profiled, 0.5 * dim * np.log1p(0.5 * half * tau**2), rtol=1e-12
        )
        gap = 0.5 * dim * np.log1p(half * tau**2) - profiled
        assert 0 < gap < 0.5 * dim * np.log(2)

    # the gap approaches (dim/2) log 2 from below as half * tau^2 -> infinity
    gaps = [
        0.5 * dim * np.log1p(h * tau**2) - deviation_penalty(dim, h, h, tau)
        for h in (10**5, 10**7, 10**9)
    ]
    assert gaps == sorted(gaps)
    np.testing.assert_allclose(gaps[-1], 0.5 * dim * np.log(2), rtol=1e-6)


def test_deviation_penalty_matches_brute_force_marginalisation():
    """Ground the closed form against the thing it approximates.

    Integrate the joint likelihood over ``(eta, delta)`` under the Gaussian code
    on ``delta``, and over ``eta`` alone under the exact-orbit model, both by
    quadrature. The penalty is exactly ``-log`` of the ratio of the two
    marginals, normalised at each model's own maximum. At the Fisher reference
    point the closed form should land within a few thousandths of a nat -- and
    the version-3.1 formula, which used ``L2`` rather than the profiled
    information, should be visibly outside that.
    """
    from scipy import integrate

    from regimeshift.detectors import _rotation, fit_fundamental
    from regimeshift.fourier import probabilities

    m, tau = 2, 0.15
    d = fundamental_dimension(m)
    R = _rotation(m, 1)
    theta_true = np.zeros(d)  # the Fisher-orthonormal reference point

    def loglik(theta, c):
        return float(c @ np.log(probabilities(theta, m)))

    for n_left, n_right in [(200, 200), (400, 100)]:
        rng = np.random.default_rng(4)
        empirical = []
        for _ in range(5):
            cL = rng.multinomial(n_left, probabilities(theta_true, m)).astype(float)
            cR = rng.multinomial(n_right, probabilities(R @ theta_true, m)).astype(float)
            _, _, ll_pen = fit_approximate_orbit(cL, cR, m, 1, tau)
            _, ll_orbit = fit_fundamental(cL + np.roll(cR, -1), m)

            def joint(dl, et):
                return np.exp(
                    loglik(np.array([et]), cL)
                    + loglik(R @ np.array([et]) + np.array([dl]), cR)
                    - 0.5 * dl**2 / tau**2
                    - 0.5 * np.log(2 * np.pi * tau**2)
                    - ll_pen
                )

            def orbit(et):
                return np.exp(
                    loglik(np.array([et]), cL)
                    + loglik(R @ np.array([et]), cR)
                    - ll_orbit
                )

            v_d, _ = integrate.dblquad(joint, -6, 6, -6, 6, epsabs=1e-13, epsrel=1e-10)
            v_c, _ = integrate.quad(orbit, -6, 6, epsabs=1e-13, epsrel=1e-10)
            empirical.append(-np.log(v_d / v_c))

        measured = float(np.mean(empirical))
        assert measured == pytest.approx(
            deviation_penalty(d, n_left, n_right, tau), abs=0.01
        )
        superseded = 0.5 * d * np.log1p(n_right * tau**2)
        if n_left == n_right:  # where the two formulas differ most
            assert abs(measured - superseded) > 0.2


@pytest.mark.parametrize("dim", [1, 2])
def test_deviation_penalty_is_symmetric_and_bounded_by_each_segment(dim):
    """The profiled information is symmetric in the two segments -- neither one
    alone identifies the deviation -- and a very long left segment recovers the
    known-shared-state limit, where ``L2`` alone is correct after all."""
    tau = 0.1
    np.testing.assert_allclose(
        deviation_penalty(dim, 300, 900, tau), deviation_penalty(dim, 900, 300, tau)
    )
    for n_left in (10**5, 10**7, 10**9):
        assert deviation_penalty(dim, n_left, 400, tau) < 0.5 * dim * np.log1p(400 * tau**2)
    np.testing.assert_allclose(
        deviation_penalty(dim, 10**9, 400, tau),
        0.5 * dim * np.log1p(400 * tau**2),
        rtol=1e-5,
    )


# --------------------------------------------------------------------------
# the fit
# --------------------------------------------------------------------------


@pytest.mark.parametrize("m", GROUPS)
def test_joint_objective_gradient_is_correct(m):
    rng = np.random.default_rng(m)
    B = fourier_design_matrix(m)
    R = rotation_matrix(m, 1)
    d = fundamental_dimension(m)
    cL = counts(np.full(m, 1 / m), 400, rng).astype(float)
    cR = counts(np.full(m, 1 / m), 400, rng).astype(float)
    for _ in range(10):
        x = rng.normal(scale=0.3, size=2 * d)
        _, analytic = _neg_joint_and_grad(x, cL, cR, B, R, 0.04)
        numeric = approx_fprime(x, lambda z: _neg_joint_and_grad(z, cL, cR, B, R, 0.04)[0], 1e-7)
        np.testing.assert_allclose(analytic, numeric, atol=1e-3, rtol=1e-4)


@pytest.mark.parametrize("m", GROUPS)
def test_zero_scale_pins_the_deviation_at_zero(m):
    rng = np.random.default_rng(m + 3)
    cL = counts(np.full(m, 1 / m), 300, rng)
    cR = counts(np.full(m, 1 / m), 300, rng)
    _, delta, _ = fit_approximate_orbit(cL, cR, m, 1, 0.0)
    np.testing.assert_allclose(delta, 0.0, atol=0)


@pytest.mark.parametrize("m", GROUPS)
def test_a_wider_prior_admits_a_larger_deviation(m):
    """The shrinkage is doing what shrinkage does."""
    segments = build_segments(m, "approximate_orbit", 0.25, deviation=1.0)
    rng = np.random.default_rng(m)
    cL = counts(segments.p_left, 4000, rng)
    cR = counts(segments.p_right, 4000, rng)
    norms = [
        np.linalg.norm(fit_approximate_orbit(cL, cR, m, 1, tau)[1])
        for tau in (0.0, 0.01, 0.05, 0.5)
    ]
    assert norms == sorted(norms)
    assert norms[-1] > norms[0]


# --------------------------------------------------------------------------
# the two nesting limits
# --------------------------------------------------------------------------


@pytest.mark.parametrize("m", GROUPS)
def test_zero_scale_reproduces_the_shared_orbit_detector_exactly(m):
    """Not approximately: Model D at ``tau = 0`` *is* Model C, because both
    maximise over the same nonidentity shifts with the deviation pinned out."""
    segments = build_segments(m, "exact_orbit", 0.25)
    rng = np.random.default_rng(m * 7)
    for _ in range(8):
        cL = counts(segments.p_left, 600, rng)
        cR = counts(segments.p_right, 600, rng)
        c = shared_orbit_detector(cL, cR, m)
        d = approximate_orbit_detector(cL, cR, m, deviation_scale=0.0)
        assert d.raw_gain == pytest.approx(c.raw_gain, abs=1e-9)
        assert d.penalty == pytest.approx(c.penalty, abs=1e-12)
        assert d.score == pytest.approx(c.score, abs=1e-9)
        assert d.selected_shift == c.selected_shift
        assert d.dimension_increment == 0.0


@pytest.mark.parametrize("m", [3, 4, 5, 6])
def test_a_wide_prior_recovers_the_fundamental_gain(m):
    """With the deviation effectively unconstrained, the alternative can reach
    any pair of coordinates, so the maximised gain matches Model B's."""
    segments = build_segments(m, "independent_fundamental", 0.25)
    rng = np.random.default_rng(m * 11)
    cL = counts(segments.p_left, 3000, rng)
    cR = counts(segments.p_right, 3000, rng)
    wide = approximate_orbit_detector(cL, cR, m, deviation_scale=50.0)
    assert wide.raw_gain == pytest.approx(fundamental_detector(cL, cR, m).raw_gain, abs=1e-4)


@pytest.mark.parametrize("m", [4, 6])
def test_gain_increases_monotonically_with_the_prior_width(m):
    """A wider prior can only enlarge the feasible set, so the maximised
    penalised likelihood gain cannot decrease."""
    segments = build_segments(m, "approximate_orbit", 0.25, deviation=0.5)
    rng = np.random.default_rng(m)
    cL = counts(segments.p_left, 2000, rng)
    cR = counts(segments.p_right, 2000, rng)
    gains = [
        approximate_orbit_detector(cL, cR, m, deviation_scale=t).raw_gain
        for t in (0.0, 0.02, 0.1, 0.5, 5.0)
    ]
    assert all(a <= b + 1e-7 for a, b in zip(gains, gains[1:]))


@pytest.mark.parametrize("m", [4, 6])
def test_penalty_sits_between_the_two_models_it_interpolates(m):
    """At a modest prior width Model D's penalty is above Model C's constant and
    below Model B's split increment -- the interpolation, in the only place it
    can live at finite n."""
    n = 2000
    d = fundamental_dimension(m)
    c_penalty = label_cost(m)
    b_penalty = split_penalty(d, n, n)
    d_penalty = label_cost(m) + deviation_penalty(d, n, n, 0.02)
    assert c_penalty < d_penalty < b_penalty + label_cost(m)


# --------------------------------------------------------------------------
# the scenario, and what the model is for
# --------------------------------------------------------------------------


@pytest.mark.parametrize("m", GROUPS)
def test_zero_deviation_scenario_is_the_exact_orbit(m):
    approx = build_segments(m, "approximate_orbit", 0.25, deviation=0.0)
    exact = build_segments(m, "exact_orbit", 0.25)
    np.testing.assert_allclose(approx.p_left, exact.p_left, atol=1e-14)
    np.testing.assert_allclose(approx.p_right, exact.p_right, atol=1e-14)
    assert approx.planted_shift == 1
    assert approx.deviation == 0.0


@pytest.mark.parametrize("m", [3, 4, 5, 6])
def test_deviation_moves_the_right_state_off_the_orbit(m):
    exact = build_segments(m, "exact_orbit", 0.25)
    previous = 0.0
    for dev in (0.1, 0.5, 1.0, 2.0):
        segments = build_segments(m, "approximate_orbit", 0.25, deviation=dev)
        distance = np.linalg.norm(segments.theta_right - exact.theta_right)
        assert distance > previous
        previous = distance
        assert segments.planted_shift is None
        # The displacement is perpendicular to the rotated state, so it is a
        # departure from the orbit rather than a rescaling along it.
        assert abs(float(exact.theta_right @ (segments.theta_right - exact.theta_right))) < 1e-12


def test_deviation_is_rejected_for_other_scenarios():
    for scenario in ("exact_orbit", "independent_fundamental", "higher_mode"):
        with pytest.raises(ValueError, match="approximate_orbit"):
            build_segments(4, scenario, 0.2, deviation=0.5)
    with pytest.raises(ValueError):
        build_segments(4, "approximate_orbit", 0.2, deviation=-1.0)


@pytest.mark.parametrize("m", [4, 6])
def test_relational_advantage_decays_as_the_orbit_relation_breaks(m):
    """The question Section 14.1 poses. Model C's edge over Model B is largest
    on an exact orbit and shrinks as the deviation grows -- gradually, not as a
    cliff, which is what makes an approximate-orbit code worth having."""
    rng = np.random.default_rng(m)
    margins = []
    for dev in (0.0, 0.5, 1.5):
        segments = build_segments(m, "approximate_orbit", 0.25, deviation=dev)
        diffs = []
        for _ in range(40):
            cL = counts(segments.p_left, 800, rng)
            cR = counts(segments.p_right, 800, rng)
            diffs.append(
                shared_orbit_detector(cL, cR, m).score - fundamental_detector(cL, cR, m).score
            )
        margins.append(float(np.mean(diffs)))
    assert margins[0] > margins[1] > margins[2]


@pytest.mark.parametrize("m", [4, 6])
def test_approximate_orbit_is_the_better_code_in_the_middle(m):
    """Where the model earns its place: on a mildly perturbed orbit it should
    beat *both* endpoints -- Model C is too rigid to fit the deviation, Model B
    pays full price for it."""
    segments = build_segments(m, "approximate_orbit", 0.25, deviation=0.6)
    rng = np.random.default_rng(m + 20)
    wins = 0
    trials = 40
    for _ in range(trials):
        cL = counts(segments.p_left, 1500, rng)
        cR = counts(segments.p_right, 1500, rng)
        approx = approximate_orbit_detector(cL, cR, m, deviation_scale=0.05).score
        rigid = shared_orbit_detector(cL, cR, m).score
        loose = fundamental_detector(cL, cR, m).score
        wins += approx > max(rigid, loose)
    assert wins > trials * 0.6, f"approximate orbit won only {wins}/{trials}"


def test_the_fit_converges_across_the_sweep():
    reset_fit_failures()
    rng = np.random.default_rng(0)
    for m in (2, 4, 6):
        for dev in (0.0, 0.5, 2.0):
            segments = build_segments(m, "approximate_orbit", 0.18, deviation=dev)
            for tau in (0.0, 0.05, 1.0):
                cL = counts(segments.p_left, 500, rng)
                cR = counts(segments.p_right, 500, rng)
                result = approximate_orbit_detector(cL, cR, m, deviation_scale=tau)
                assert np.isfinite(result.score)
                assert result.selected_shift in range(1, m)
    assert fit_failure_count() == 0

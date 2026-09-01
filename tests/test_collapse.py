"""The shared-orbit limit law at and near orbit collapse (Section 6.4).

These are the tests behind the one theorem in the paper that is not a dimension
count. They check three things: that the law describes the detector, that exact
collapse is the least favourable local null, and that no constant label cost can
control the raw rule there -- the paper's negative result about naive MDL
thresholding.
"""

from __future__ import annotations

import numpy as np
import pytest

from regimeshift.collapse import (
    collapse_law,
    critical_value,
    false_alarm_rate,
    least_favourable_sweep,
)
from regimeshift.detectors import label_cost, shared_orbit_detector
from regimeshift.fourier import fundamental_dimension, probabilities

GROUPS = [2, 3, 4, 5, 6]


def _simulate(m, h, n_left, n_right, reps, rng):
    """Raw shared-orbit gains with both segments at eta = h / sqrt(L)."""
    d = fundamental_dimension(m)
    eta = np.zeros(d) if h is None else np.asarray(h, float) / np.sqrt(n_left + n_right)
    p = probabilities(eta, m)
    A = rng.multinomial(n_left, p, size=reps).astype(float)
    B = rng.multinomial(n_right, p, size=reps).astype(float)
    return np.array([shared_orbit_detector(A[i], B[i], m).raw_gain for i in range(reps)])


@pytest.mark.slow
@pytest.mark.parametrize("m", GROUPS)
def test_the_law_describes_the_detector_at_exact_collapse(m):
    """Both segments uniform: the detector's raw gain matches the limit law."""
    rng = np.random.default_rng(20260713 + m)
    sim = _simulate(m, None, 2000, 2000, 4000, rng)
    law = collapse_law(m, reps=200_000, rng=np.random.default_rng(11 + m))

    for q in (0.5, 0.9, 0.95, 0.99):
        assert np.quantile(sim, q) == pytest.approx(np.quantile(law, q), abs=0.12), f"quantile {q}"
    assert sim.mean() == pytest.approx(law.mean(), abs=0.05)


@pytest.mark.slow
@pytest.mark.parametrize("m", [4, 6])
@pytest.mark.parametrize("h_norm", [1.5, 3.0])
def test_the_law_describes_the_detector_at_local_alternatives(m, h_norm):
    """And away from collapse, where w_r = R^r h - h stops vanishing."""
    d = fundamental_dimension(m)
    h = np.zeros(d)
    h[0] = h_norm
    rng = np.random.default_rng(4242 + m + int(10 * h_norm))
    sim = _simulate(m, h, 1500, 1500, 4000, rng)
    law = collapse_law(m, h=h, reps=200_000, rng=np.random.default_rng(7 + m))

    assert np.quantile(sim, 0.95) == pytest.approx(np.quantile(law, 0.95), abs=0.12)
    assert sim.mean() == pytest.approx(law.mean(), abs=0.08)


@pytest.mark.slow
@pytest.mark.parametrize("m", [3, 5])
@pytest.mark.parametrize("rho", [0.25, 0.1])
def test_the_law_holds_at_unbalanced_splits(m, rho):
    total, reps = 6000, 4000
    n_left = int(total * rho)
    rng = np.random.default_rng(808 + m + int(100 * rho))
    sim = _simulate(m, None, n_left, total - n_left, reps, rng)
    law = collapse_law(m, rho=n_left / total, reps=200_000, rng=np.random.default_rng(3 + m))

    assert np.quantile(sim, 0.95) == pytest.approx(np.quantile(law, 0.95), abs=0.15)
    assert sim.mean() == pytest.approx(law.mean(), abs=0.06)


@pytest.mark.parametrize("m", GROUPS)
def test_collapse_is_least_favourable_over_the_evaluated_grid(m):
    """A numerical finding, not a corollary.

    Swept over ||h|| <= 6 and a full orbit sector with common random numbers.
    Paired draws matter: estimating a few hundred quantiles independently and
    taking the maximum manufactures an excess of about 0.02 from Monte Carlo
    noise alone, the same size as the effect being looked for. They do not
    remove selection bias entirely, so the sweep also re-checks its own argmax
    against the origin on independent draws -- that difference should straddle
    zero rather than sit above it.

    None of this shows uniform size control. It is evidence at the settings
    tested, and the paper states it that way.
    """
    result = least_favourable_sweep(m, reps=120_000)
    assert result["excess"] < 0.01, result
    assert result["argmax_norm"] < 0.2, result
    assert result["points_exceeding_by_0.01"] == 0, result
    assert abs(result["independent_excess_at_argmax"]) < 0.02, result


@pytest.mark.parametrize("m", GROUPS)
@pytest.mark.parametrize("h_norm", [0.0, 1.5, 3.0])
def test_the_balanced_Y_form_is_the_same_law(m, h_norm):
    """The readable statement of the theorem on a balanced split.

    With ``Y_i = h_n + Z_i`` in the per-segment convention ``eta_n = h_n/sqrt(n)``,

        W = 1/4 max_r [ ||Y_1 + R^-r Y_2||^2 - ||Y_1 + Y_2||^2 ],

    which is the general form evaluated at ``rho = 1/2`` and ``h = sqrt(2) h_n``.
    Driven by the same random numbers the two agree exactly, which is the check
    that the convention bookkeeping is right.
    """
    from regimeshift.fourier import rotation_matrix

    d = fundamental_dimension(m)
    h_n = np.zeros(d)
    h_n[0] = h_norm

    def y_form(rng, reps):
        Y1 = h_n + rng.normal(size=(reps, d))
        Y2 = h_n + rng.normal(size=(reps, d))
        base = np.sum((Y1 + Y2) ** 2, axis=1)
        best = np.full(reps, -np.inf)
        for r in range(1, m):
            best = np.maximum(best, np.sum((Y1 + Y2 @ rotation_matrix(m, -r).T) ** 2, axis=1))
        return (best - base) / 4.0

    reps = 60_000
    y = y_form(np.random.default_rng(4), reps)
    w = collapse_law(m, h=np.sqrt(2) * h_n, rho=0.5, reps=reps, rng=np.random.default_rng(4))
    np.testing.assert_allclose(np.quantile(y, [0.5, 0.9, 0.95, 0.99]),
                               np.quantile(w, [0.5, 0.9, 0.95, 0.99]), atol=1e-9)


def test_the_naive_unbalanced_extension_of_the_Y_form_is_wrong():
    """Guards a trap: weighting ``Y_i`` by ``sqrt(rho_i)`` reproduces the law at
    ``h = 0`` and fails off collapse, so agreement at the collapse point is not
    evidence that a candidate general form is right."""
    from regimeshift.fourier import rotation_matrix

    m, rho, d = 4, 0.25, fundamental_dimension(4)
    a, b = np.sqrt(rho), np.sqrt(1 - rho)

    def naive(h, rng, reps):
        Y1 = h + rng.normal(size=(reps, d)); Y2 = h + rng.normal(size=(reps, d))
        base = np.sum((a * Y1 + b * Y2) ** 2, axis=1)
        best = np.full(reps, -np.inf)
        for r in range(1, m):
            best = np.maximum(best, np.sum((a * Y1 + b * (Y2 @ rotation_matrix(m, -r).T)) ** 2, axis=1))
        return (best - base) / 2.0

    at_zero = np.zeros(d)
    off = np.array([4.0, 0.0])
    q = lambda x: float(np.quantile(x, 0.95))
    # agrees at collapse ...
    assert q(naive(at_zero, np.random.default_rng(1), 80_000)) == pytest.approx(
        q(collapse_law(m, at_zero, rho, 80_000, np.random.default_rng(2))), abs=0.05)
    # ... and is badly wrong away from it
    assert abs(q(naive(off, np.random.default_rng(1), 80_000))
               - q(collapse_law(m, off, rho, 80_000, np.random.default_rng(2)))) > 1.0


@pytest.mark.parametrize("m", GROUPS)
def test_the_limit_does_not_depend_on_the_sample_size(m):
    """Which is why the zero-threshold false-positive rate does not vanish: the
    statistic is O_p(1) at collapse, and no amount of data shrinks it."""
    rng = np.random.default_rng(99 + m)
    small = _simulate(m, None, 400, 400, 2500, rng)
    large = _simulate(m, None, 3200, 3200, 2500, rng)
    assert np.quantile(small, 0.95) == pytest.approx(np.quantile(large, 0.95), abs=0.2)


@pytest.mark.parametrize("m", GROUPS)
def test_no_constant_label_cost_controls_the_raw_rule_at_collapse(m):
    """The paper's negative result.

    The two-part code charges ``log(g-1)``; a code including the identity would
    charge ``log g``. Neither is anywhere near the level a test would need, and
    the gap is not a small-sample effect -- it is the asymptotic law.
    """
    rng = np.random.default_rng(5150 + m)
    draws = collapse_law(m, reps=200_000, rng=rng)

    at_label = float(np.mean(draws > label_cost(m)))
    at_log_g = float(np.mean(draws > np.log(m)))
    assert at_label > 0.13, f"g={m}: {at_label}"       # far above any nominal level
    assert at_log_g < at_label                          # stricter, and still far off
    assert at_log_g > 0.09, f"g={m}: {at_log_g}"

    # the level-0.05 threshold is several nats above either convention
    q95 = float(np.quantile(draws, 0.95))
    assert q95 > label_cost(m) + 0.8
    assert q95 > np.log(m)


def test_the_critical_value_is_a_level_not_a_constant():
    """Guards the distinction the paper insists on: this object moves with alpha,
    so it is a testing threshold and cannot stand in for a codelength constant."""
    rng = np.random.default_rng(31337)
    draws = collapse_law(6, reps=200_000, rng=rng)
    q90, q95, q99 = (float(np.quantile(draws, q)) for q in (0.90, 0.95, 0.99))
    assert q90 < q95 < q99
    assert q99 - q90 > 1.0


def test_collapse_law_rejects_invalid_arguments():
    with pytest.raises(ValueError, match="rho"):
        collapse_law(4, rho=0.0, reps=10)
    with pytest.raises(ValueError, match="rho"):
        collapse_law(4, rho=1.0, reps=10)
    with pytest.raises(ValueError, match="shape"):
        collapse_law(4, h=np.zeros(3), reps=10)


@pytest.mark.parametrize("m", GROUPS)
def test_large_local_states_return_the_detector_to_regular_behaviour(m):
    """Far from collapse the label is identifiable again and the statistic
    collapses toward minus infinity, which is why the failure is local."""
    d = fundamental_dimension(m)
    near = collapse_law(m, reps=60_000, rng=np.random.default_rng(1))
    far_h = np.zeros(d)
    far_h[0] = 8.0
    far = collapse_law(m, h=far_h, reps=60_000, rng=np.random.default_rng(1))
    assert far.mean() < near.mean() - 1.0
    assert np.quantile(far, 0.95) < np.quantile(near, 0.95)


# ---------------------------------------------------------------------------
# Corollary 3: the zero-threshold rate at collapse is exactly (g-1)/g
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("m", [2, 3, 4, 5, 6, 7, 8])
def test_zero_threshold_rate_at_collapse_is_g_minus_one_over_g(m):
    """An exact identity, not a simulated quantity.

    At collapse the aligned energies differ only in their cross term, so the
    argmax over all ``g`` shifts is cyclically exchangeable and each is equally
    likely to win. The raw gain is positive exactly when a nonidentity shift
    wins, which is ``g-1`` times out of ``g``.
    """
    draws = collapse_law(m, reps=400_000, rng=np.random.default_rng(1000 + m))
    assert float(np.mean(draws > 0)) == pytest.approx((m - 1) / m, abs=0.003)


@pytest.mark.parametrize("m", [3, 4, 6])
@pytest.mark.parametrize("rho", [0.5, 0.25, 0.1, 0.75])
def test_that_rate_does_not_depend_on_the_split_fraction(m, rho):
    """The corollary holds for every ``rho``, because
    ``A_r = rho|U1|^2 + (1-rho)|U2|^2 + 2 sqrt(rho(1-rho)) U1'R^-r U2`` and only
    the cross term carries ``r`` -- so the *event* ``{W > 0}`` is free of rho."""
    draws = collapse_law(m, rho=rho, reps=400_000, rng=np.random.default_rng(1000 + m))
    assert float(np.mean(draws > 0)) == pytest.approx((m - 1) / m, abs=0.003)


@pytest.mark.parametrize("m", [3, 5])
def test_the_argmax_over_shifts_is_uniform_at_collapse(m):
    """The mechanism behind Corollary 3, checked directly: every shift including
    the identity is equally likely to maximise the aligned energy."""
    from regimeshift.fourier import rotation_matrix

    d = fundamental_dimension(m)
    reps = 400_000
    rng = np.random.default_rng(7)
    U1 = rng.normal(size=(reps, d))
    U2 = rng.normal(size=(reps, d))
    energies = np.stack(
        [np.sum((U1 + U2 @ rotation_matrix(m, -r).T) ** 2, axis=1) for r in range(m)], axis=1
    )
    shares = np.bincount(energies.argmax(axis=1), minlength=m) / reps
    np.testing.assert_allclose(shares, np.full(m, 1.0 / m), atol=0.004)


@pytest.mark.parametrize("m", [2, 3, 4, 5, 6])
def test_the_level_threshold_sits_well_above_the_codelength(m):
    """Quantifies the gap the paper quotes: 1.0 to 1.8 nats between the
    level-0.05 critical value and the two-part label cost."""
    q95 = critical_value(m, reps=400_000, rng=np.random.default_rng(20260713 + m))
    gap = q95 - label_cost(m)
    assert 1.0 <= gap <= 1.8, f"g={m}: gap {gap}"

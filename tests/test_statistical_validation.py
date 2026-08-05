"""Seeded Monte Carlo validation of the manuscript's empirical claims.

These are the claims that can only fail probabilistically. Everything here runs
a real simulation grid, so the whole module is marked ``slow``:

    pytest -m slow

Tolerances are deliberately looser than the manuscript's own standard errors:
the grids below are a fraction of the 468,000-dataset production run, so they
test the *shape* of each law (which coefficient, which ordering, which
direction) rather than the third decimal place. Every configuration is seeded,
so a failure here is a real regression, not a flaky draw.

The reference values in the docstrings are the manuscript's production numbers
(Tables 3-8), quoted so a drift in this implementation is visible.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from regimeshift.analysis import (
    crossover_bootstrap,
    crossover_estimates,
    crossover_ratio_bootstrap,
    crossover_ratio_summary,
    gain_residual_regression,
    predicted_slope,
    score_regression,
    score_regression_summary,
)
from regimeshift.runner import run_grid
from regimeshift.simulation import build_grid

pytestmark = pytest.mark.slow

EFFECTS = (0.08, 0.12, 0.18, 0.25)
LENGTHS = (100, 200, 400, 800, 1600)
N_ALT = 300
N_NULL = 300


@pytest.fixture(scope="module")
def results() -> pd.DataFrame:
    """One shared grid across groups 2, 4 and 6 and all three scenarios."""
    configs = build_grid(
        groups=(2, 4, 6),
        scenarios=("exact_orbit", "independent_fundamental", "higher_mode"),
        effects=EFFECTS,
        segment_lengths=LENGTHS,
        n_alt=N_ALT,
        n_null=N_NULL,
    )
    return run_grid(configs, workers=4)


def _fit(results: pd.DataFrame, detector: str, scenario: str, m: int) -> dict:
    subset = results[
        (results["detector"] == detector) & (results["scenario"] == scenario) & (results["m"] == m)
    ]
    assert len(subset) == len(EFFECTS) * len(LENGTHS)
    return score_regression(subset)


# --------------------------------------------------------------------------
# the affine score model itself (Appendix B)
# --------------------------------------------------------------------------


@pytest.mark.parametrize("detector,scenario", [
    ("full", "higher_mode"),
    ("fundamental", "independent_fundamental"),
    ("shared_orbit", "exact_orbit"),
])
@pytest.mark.parametrize("m", [4, 6])
def test_affine_score_model_fits_and_recovers_a_unit_gain_coefficient(results, detector, scenario, m):
    """The mean score should be gain*n minus a penalty term: coefficient one on
    the expected total gain, and an extremely high R^2. The manuscript reports
    gain coefficients within about 1.2% of one and R^2 >= 0.999."""
    fit = _fit(results, detector, scenario, m)
    assert fit["r_squared"] > 0.99
    assert fit["beta_gain"] == pytest.approx(1.0, abs=0.08)


# --------------------------------------------------------------------------
# Table 3: the full model tracks the simplex dimension
# --------------------------------------------------------------------------


@pytest.mark.parametrize("m", [2, 4, 6])
def test_full_penalty_slope_tracks_the_full_dimension(results, m):
    """Manuscript Table 3: 1.515, 2.119, 2.468 for m = 4, 5, 6 against
    predictions 1.5, 2.0, 2.5."""
    fit = _fit(results, "full", "higher_mode" if m >= 4 else "exact_orbit", m)
    predicted = predicted_slope("full", m)
    assert fit["penalty_slope"] == pytest.approx(predicted, abs=0.35)


def test_full_penalty_slope_increases_with_group_order(results):
    """The signature of Model A: the coefficient grows with the simplex
    dimension rather than staying fixed."""
    slopes = [_fit(results, "full", "higher_mode", m)["penalty_slope"] for m in (4, 6)]
    assert slopes[1] > slopes[0] + 0.4


# --------------------------------------------------------------------------
# Table 4: the fundamental model tracks the representation dimension
# --------------------------------------------------------------------------


@pytest.mark.parametrize("m", [2, 4, 6])
def test_fundamental_penalty_slope_tracks_the_representation_dimension(results, m):
    """Manuscript Table 4: 0.457 at m = 2 and 0.967-1.047 for m = 3..6, against
    predictions 0.5 and 1.0."""
    fit = _fit(results, "fundamental", "independent_fundamental", m)
    assert fit["penalty_slope"] == pytest.approx(predicted_slope("fundamental", m), abs=0.3)


def test_fundamental_slope_is_flat_in_group_order(results):
    """The clearest separation from Model A: the fundamental coefficient stays
    near 1 from m = 4 to m = 6 while the full-simplex dimension grows from 3
    to 5."""
    fundamental = [_fit(results, "fundamental", "independent_fundamental", m)["penalty_slope"] for m in (4, 6)]
    assert abs(fundamental[1] - fundamental[0]) < 0.35
    full = [_fit(results, "full", "higher_mode", m)["penalty_slope"] for m in (4, 6)]
    assert full[1] - full[0] > abs(fundamental[1] - fundamental[0])


# --------------------------------------------------------------------------
# Table 5: the shared-orbit residual slope is near zero, not exactly zero
# --------------------------------------------------------------------------


@pytest.mark.parametrize("m", [2, 4, 6])
def test_shared_orbit_residual_slope_is_near_zero(results, m):
    """Manuscript Table 5 reports residuals of -0.093 to 0.206 against a
    structural prediction of zero. The claim under test is the qualified one of
    Section 9.3: a *near-zero leading logarithmic coefficient*, not exact
    finite-sample constancy."""
    fit = _fit(results, "shared_orbit", "exact_orbit", m)
    assert abs(fit["penalty_slope"]) < 0.5


@pytest.mark.parametrize("m", [4, 6])
def test_shared_orbit_slope_is_far_below_the_regular_slopes(results, m):
    """The separation that carries the paper's central result."""
    shared = abs(_fit(results, "shared_orbit", "exact_orbit", m)["penalty_slope"])
    fundamental = _fit(results, "fundamental", "independent_fundamental", m)["penalty_slope"]
    full = _fit(results, "full", "higher_mode", m)["penalty_slope"]
    assert shared < fundamental / 2
    assert shared < full / 2
    assert fundamental < full


# --------------------------------------------------------------------------
# Table 6: calibrated sample-length advantage on matched data
# --------------------------------------------------------------------------


def test_calibrated_crossover_advantage_on_exact_orbit_data(results):
    """Manuscript Table 6: on exact-orbit data the shared detector needed about
    30%, 33% and 39% fewer observations than the full detector for m = 4, 5, 6,
    and 17-19% fewer than the fundamental detector."""
    ratios = crossover_ratio_summary(crossover_estimates(results), scenario="exact_orbit")
    ratios = ratios.set_index("m")
    for m in (4, 6):
        assert ratios.loc[m, "shared_orbit/full"] < 1.0
        assert ratios.loc[m, "shared_orbit/fundamental"] <= 1.0
        assert ratios.loc[m, "fundamental/full"] <= 1.0
    # The advantage should widen as the full model's dimension grows.
    assert ratios.loc[6, "shared_orbit/full"] <= ratios.loc[4, "shared_orbit/full"] + 0.1


def test_no_detector_advantage_at_m2_where_all_models_coincide(results):
    """At m = 2 the fundamental component spans the whole tangent space, so
    Models A and B are the same model and their crossovers must agree."""
    ratios = crossover_ratio_summary(crossover_estimates(results), scenario="exact_orbit")
    row = ratios[ratios["m"] == 2]
    if len(row) and np.isfinite(row.iloc[0]["fundamental/full"]):
        assert row.iloc[0]["fundamental/full"] == pytest.approx(1.0, abs=1e-6)


# --------------------------------------------------------------------------
# Table 7: the advantage is conditional on structural correctness
# --------------------------------------------------------------------------


@pytest.mark.parametrize("m", [4, 6])
def test_full_detector_dominates_under_higher_mode_misspecification(results, m):
    """Manuscript Table 7 at total length 6,400: full power 0.999/0.975/0.962
    against fundamental 0.399/0.211/0.182 and shared orbit 0.082/0.044/0.055.

    This is the result that keeps the constrained advantage honest: it is not a
    generally lower threshold, it is conditional on the geometry being right."""
    longest = max(LENGTHS) * 2
    subset = results[
        (results["scenario"] == "higher_mode")
        & (results["m"] == m)
        & (results["total_length"] == longest)
    ]
    power = subset.groupby("detector")["power_calibrated"].mean()
    assert power["full"] > power["fundamental"]
    assert power["full"] > power["shared_orbit"]
    assert power["full"] > 0.6


@pytest.mark.parametrize("m", [4, 6])
def test_constrained_detectors_win_when_their_geometry_is_correct(results, m):
    """The other direction of the same claim, on matched exact-orbit data."""
    subset = results[
        (results["scenario"] == "exact_orbit")
        & (results["m"] == m)
        & (results["total_length"] == 2 * LENGTHS[1])
    ]
    power = subset.groupby("detector")["power_calibrated"].mean()
    assert power["shared_orbit"] >= power["full"]
    assert power["fundamental"] >= power["full"]


def test_shared_orbit_survives_a_mild_departure_from_exact_symmetry(results):
    """A finite-sample effect worth pinning down, and *not* a contradiction.

    In the independent-fundamental scenario the right coordinate is a shrunk,
    rotated version of the left, so Model C's hypothesis is strictly false and
    its population gain is strictly smaller than Model B's (asserted exactly in
    ``test_shared_orbit_gain_falls_short_on_independent_changes``). At finite
    samples, however, Model C's constant label cost is so much cheaper than a
    growing regular penalty that it can still detect at least as well under a
    *mild* deviation. This is the approximate-orbit regime of Section 14.1: the
    relational advantage degrades gradually, not abruptly.

    The claim that keeps this honest is the higher-mode test above, where the
    deviation leaves the fundamental subspace entirely and Model C collapses.
    """
    subset = results[
        (results["scenario"] == "independent_fundamental")
        & (results["m"] == 6)
        & (results["total_length"] == max(LENGTHS) * 2)
    ]
    power = subset.groupby("detector")["power_calibrated"].mean()
    assert power["shared_orbit"] > power["full"]
    assert power["fundamental"] > power["full"]


# --------------------------------------------------------------------------
# Table 8: relative-shift recovery
# --------------------------------------------------------------------------


@pytest.mark.parametrize("m", [2, 4, 6])
def test_relative_shift_recovery_at_the_longest_length(results, m):
    """Manuscript Table 8 at total length 6,400: mean accuracy 1.000, 1.000,
    0.997, 0.9945, 0.989 for m = 2..6. The slight decline with m is expected --
    the detector maximises over more candidate shifts.

    This grid tops out at total length 3,200 rather than 6,400, so the weakest
    effect at m = 6 does not reach the manuscript's accuracy; the assertion is
    made at the largest effect, plus monotonicity in effect size."""
    subset = results[
        (results["scenario"] == "exact_orbit")
        & (results["m"] == m)
        & (results["detector"] == "shared_orbit")
        & (results["total_length"] == max(LENGTHS) * 2)
    ].sort_values("effect")
    assert subset["shift_accuracy"].iloc[-1] > 0.95
    assert subset["shift_accuracy"].mean() > 0.8
    # Accuracy must not degrade as the signal strengthens.
    accuracy = subset["shift_accuracy"].to_numpy()
    assert np.all(np.diff(accuracy) > -0.02)


# --------------------------------------------------------------------------
# null calibration
# --------------------------------------------------------------------------


def test_raw_mdl_rule_is_conservative_for_the_regular_detectors(results):
    """Models A and B pay a penalty that grows with n, so with no change present
    the raw zero-threshold rule fires well below the 5% calibration target."""
    regular = results[results["detector"].isin(["full", "fundamental"])]
    assert regular["null_rate_zero_threshold"].max() <= 0.06
    assert regular["null_rate_zero_threshold"].mean() < 0.02


def test_shared_orbit_raw_rule_is_not_conservative(results):
    """The other side of Model C's constant penalty, and the reason the study
    calibrates every detector to a common 5% null rather than comparing raw
    zero-threshold rules.

    A penalty of ``log(m - 1)`` does not grow with n, so it provides no
    increasing protection -- and at m = 2 it is exactly zero, leaving the raw
    rule with no protection at all. Section 6.3's singular qualification applies
    here: near orbit collapse the relative label is unidentifiable and the null
    score is not well separated from zero."""
    shared = results[results["detector"] == "shared_orbit"]
    assert shared["null_rate_zero_threshold"].max() > 0.1

    at_m2 = shared[(shared["m"] == 2) & (shared["segment_length"] == min(LENGTHS))]
    assert at_m2["null_rate_zero_threshold"].max() > 0.1

    # The calibrated critical values are what make the comparison fair; they are
    # strictly positive exactly where the raw rule over-fires.
    assert (shared.loc[shared["null_rate_zero_threshold"] > 0.1, "critical_value"] > 0).all()


def test_every_configuration_produced_a_finite_score(results):
    assert np.isfinite(results["mean_score"]).all()
    assert np.isfinite(results["critical_value"]).all()
    assert results["power_calibrated"].between(0.0, 1.0).all()


# --------------------------------------------------------------------------
# analysis integrity, on real Monte Carlo output
# --------------------------------------------------------------------------


@pytest.mark.parametrize("detector,scenario", [
    ("full", "higher_mode"),
    ("fundamental", "independent_fundamental"),
    ("shared_orbit", "exact_orbit"),
])
@pytest.mark.parametrize("m", [4, 6])
def test_weighted_and_unweighted_slopes_agree(results, detector, scenario, m):
    """The review's concern that the reported slopes might be artifacts of an
    unweighted fit on few aggregate design points. Monte Carlo-variance
    weighting must not move a slope enough to change its interpretation."""
    subset = results[
        (results["detector"] == detector) & (results["scenario"] == scenario) & (results["m"] == m)
    ]
    ols = score_regression(subset)
    wls = score_regression(subset, weighted=True)
    assert wls["method"] == "wls"
    assert abs(wls["penalty_slope"] - ols["penalty_slope"]) < 0.4
    # Both must still sit on the correct side of the neighbouring predictions.
    predicted = predicted_slope(detector, m)
    assert abs(wls["penalty_slope"] - predicted) < 0.45


def test_design_is_not_badly_conditioned(results):
    """`n * gain` and `log n` are both deterministic in length within an effect,
    so the design's conditioning is worth reporting rather than assuming."""
    summary = score_regression_summary(results)
    assert summary["condition_number"].max() < 5000


@pytest.mark.parametrize("detector,scenario", [
    ("full", "higher_mode"),
    ("fundamental", "independent_fundamental"),
    ("shared_orbit", "exact_orbit"),
])
@pytest.mark.parametrize("m", [4, 6])
def test_penalty_slope_identity_holds_on_real_output(results, detector, scenario, m):
    """``penalty_slope = d/2 - residual_slope`` is an algebraic identity, so it
    must hold to numerical precision on simulated output too. It decomposes each
    empirical slope into the exact penalty we subtract and the raw gain's
    departure from ``n * G``."""
    subset = results[
        (results["detector"] == detector) & (results["scenario"] == scenario) & (results["m"] == m)
    ]
    score = score_regression(subset)["penalty_slope"]
    residual = gain_residual_regression(subset)["residual_slope"]
    assert predicted_slope(detector, m) - residual == pytest.approx(score, abs=1e-8)


@pytest.mark.parametrize("m", [2, 4, 6])
def test_gain_residual_slopes_are_small_for_every_detector(results, m):
    """Theory says the maximised gain is ``n * G + O(1)``, so the residual
    log-length slope should be near zero for all three detectors -- including
    the full detector, whose large penalty slope comes from the exact penalty
    rather than from any drift in its gain."""
    for detector, scenario in (
        ("full", "higher_mode" if m >= 4 else "exact_orbit"),
        ("fundamental", "independent_fundamental"),
        ("shared_orbit", "exact_orbit"),
    ):
        subset = results[
            (results["detector"] == detector)
            & (results["scenario"] == scenario)
            & (results["m"] == m)
        ]
        residual = gain_residual_regression(subset)["residual_slope"]
        assert abs(residual) < 0.5, f"{detector} m={m}: residual slope {residual}"


def test_crossover_intervals_are_produced_and_contain_their_estimates(results):
    boot = crossover_bootstrap(results, n_boot=200)
    internal = boot[boot["status"] == "internal"].dropna(subset=["ci_low", "ci_high"])
    assert len(internal) > 5
    for _, row in internal.iterrows():
        assert row["ci_low"] <= row["crossover_length"] <= row["ci_high"]


def test_shared_orbit_advantage_survives_its_uncertainty_interval(results):
    """Table 6's headline claim, with the uncertainty the manuscript omits: the
    shared/full crossover ratio must stay below one across the bootstrap, not
    merely at the point estimate."""
    boot = crossover_ratio_bootstrap(results, scenario="exact_orbit", n_boot=300)
    boot = boot.set_index("m")
    for m in (4, 6):
        assert boot.loc[m, "shared_orbit/full_boot_n"] > 100
        assert boot.loc[m, "shared_orbit/full_ci_high"] < 1.0


def test_no_optimizer_failures_anywhere_in_the_grid(results):
    assert results["optimizer_failures"].sum() == 0


def test_penalty_slope_is_invariant_to_the_split_fraction():
    """Section 4.2 predicts that rho moves only the bounded term
    ``(d/2) log(rho (1 - rho))``, leaving the ``log n`` coefficient at ``d/2``.
    This is the empirical check the manuscript never ran.

    Run on matched exact-orbit data. The property under test belongs to the
    *penalty*, which does not depend on the scenario, so the scenario should be
    chosen to make the slope estimable rather than to stress the detector: on
    weak-signal misspecified data the estimate is noisy enough at an unbalanced
    split -- where the left segment is only a quarter of the sample -- to
    swamp the effect being measured.
    """
    slopes = {}
    for rho in (0.5, 0.25):
        configs = build_grid(
            groups=(4,), scenarios=("exact_orbit",), effects=EFFECTS,
            segment_lengths=(100, 200, 400, 800, 1600), n_alt=N_ALT, n_null=N_NULL,
            split_fractions=(rho,),
        )
        frame = run_grid(configs, workers=4)
        slopes[rho] = score_regression(frame[frame["detector"] == "full"])["penalty_slope"]
    assert abs(slopes[0.5] - slopes[0.25]) < 0.3
    for slope in slopes.values():
        assert slope == pytest.approx(predicted_slope("full", 4), abs=0.3)

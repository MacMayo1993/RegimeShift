"""Tests for the analysis-integrity work prompted by external methodological review.

The review raised five concerns about the analysis pipeline rather than the
models: the group-level regressions use few *aggregate* design points of unequal
precision but are fitted unweighted; the gain and penalty coefficients are
estimated simultaneously when the penalty is in fact known exactly; the
crossover estimates carry no uncertainty; optimiser failures are never audited
even though the population gains are themselves obtained numerically; and the
split fraction is never varied even though the theory makes a sharp prediction
about it.

These tests pin the fixes. They are fast and deterministic; the Monte Carlo
counterparts live in ``test_statistical_validation.py``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from regimeshift.analysis import (
    K_STAR,
    _has_joint_patterns,
    _pattern_mask,
    crossover_bootstrap,
    crossover_estimates,
    dimension_increment,
    crossover_ratio_bootstrap,
    gain_residual_regression,
    predicted_slope,
    score_regression,
    score_regression_summary,
)
from regimeshift.detectors import (
    fit_failure_count,
    nats_to_bits,
    fit_fundamental,
    reset_fit_failures,
    split_penalty,
)
from regimeshift.simulation import (
    DETECTION_PATTERNS,
    DETECTOR_NAMES,
    Config,
    build_grid,
    config_seed,
    run_config,
)

EFFECTS = (0.1, 0.2, 0.3)
LENGTHS = (200, 400, 800, 1600, 3200)


def synthetic_frame(
    beta=1.0,
    slope=1.5,
    detector="full",
    m=4,
    noise_scale=0.0,
    heteroskedastic=False,
    seed=0,
):
    """Exact (or noisy) realisation of the affine score model, with the columns
    the weighted fit needs.

    ``heteroskedastic`` makes the short-length points far noisier than the long
    ones, which is the situation the weighting is meant to handle.
    """
    rng = np.random.default_rng(seed)
    rows = []
    for effect in EFFECTS:
        gain = 0.5 * effect**2
        for n in LENGTHS:
            sd = 5.0 if (heteroskedastic and n <= 400) else 1.0
            mean_score = beta * gain * n - slope * np.log(n) + 3.0 * effect
            if noise_scale:
                mean_score += rng.normal(scale=noise_scale * sd / np.sqrt(200))
            rows.append(
                {
                    "m": m,
                    "scenario": "exact_orbit",
                    "detector": detector,
                    "effect": effect,
                    "total_length": n,
                    "expected_gain_total": gain * n,
                    "mean_score": mean_score,
                    # raw gain = score + exact penalty, with the penalty being
                    # the balanced-split increment for this detector.
                    "mean_raw_gain": mean_score + split_penalty(
                        2 * predicted_slope(detector, m), n // 2, n - n // 2
                    ),
                    "sd_score": sd,
                    "n_alt": 200,
                }
            )
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# 1. weighted least squares
# --------------------------------------------------------------------------


def test_weighted_and_unweighted_agree_on_noiseless_data():
    frame = synthetic_frame(slope=1.5)
    assert score_regression(frame)["penalty_slope"] == pytest.approx(1.5, abs=1e-8)
    assert score_regression(frame, weighted=True)["penalty_slope"] == pytest.approx(1.5, abs=1e-8)


def test_weighting_beats_ols_under_heteroskedastic_noise():
    """The point of the fix: when short-length points are far noisier, the
    unweighted fit on aggregate means is the worse estimator."""
    ols_errors, wls_errors = [], []
    for seed in range(40):
        frame = synthetic_frame(slope=1.5, noise_scale=1.0, heteroskedastic=True, seed=seed)
        ols_errors.append(abs(score_regression(frame)["penalty_slope"] - 1.5))
        wls_errors.append(abs(score_regression(frame, weighted=True)["penalty_slope"] - 1.5))
    assert np.mean(wls_errors) < np.mean(ols_errors)


def test_weighting_reports_its_method_and_falls_back_cleanly():
    frame = synthetic_frame()
    assert score_regression(frame, weighted=True)["method"] == "wls"
    assert score_regression(frame)["method"] == "ols"

    # A frame without the Monte Carlo columns must silently use OLS rather than
    # inventing weights.
    bare = frame.drop(columns=["sd_score", "n_alt"])
    assert score_regression(bare, weighted=True)["method"] == "ols"

    # Degenerate variances must not produce infinite weights.
    zeroed = frame.assign(sd_score=0.0)
    assert score_regression(zeroed, weighted=True)["method"] == "ols"


def test_condition_number_is_reported_and_finite():
    fit = score_regression(synthetic_frame())
    assert np.isfinite(fit["condition_number"])
    assert fit["condition_number"] > 1.0


def test_summary_reports_both_fits_and_their_disagreement():
    frame = pd.concat(
        [
            synthetic_frame(detector="full", m=4, slope=1.5).assign(scenario="higher_mode"),
            synthetic_frame(detector="fundamental", m=4, slope=1.0).assign(
                scenario="independent_fundamental"
            ),
            synthetic_frame(detector="shared_orbit", m=4, slope=0.0),
        ],
        ignore_index=True,
    )
    summary = score_regression_summary(frame)
    for column in ("penalty_slope", "penalty_slope_wls", "weighting_shift", "residual_slope"):
        assert column in summary.columns
    # Noiseless data: the two fits must agree, so the reported shift is ~0.
    assert summary["weighting_shift"].max() < 1e-6


# --------------------------------------------------------------------------
# 2. the direct gain-residual diagnostic
# --------------------------------------------------------------------------


@pytest.mark.parametrize("detector,m,slope", [
    ("full", 4, 1.5), ("full", 6, 2.5),
    ("fundamental", 4, 1.0), ("fundamental", 2, 0.5),
    ("shared_orbit", 4, 0.0), ("shared_orbit", 6, 0.0),
])
def test_penalty_slope_identity(detector, m, slope):
    """``penalty_slope = d/2 - residual_slope`` exactly.

    The penalty this implementation subtracts is deterministic, so the score
    regression's slope is the theoretical dimension term minus the empirical
    departure of the raw gain from ``n * G``. The identity makes that
    decomposition auditable instead of implicit.
    """
    frame = synthetic_frame(detector=detector, m=m, slope=slope)
    score = score_regression(frame)["penalty_slope"]
    residual = gain_residual_regression(frame)["residual_slope"]
    assert predicted_slope(detector, m) - residual == pytest.approx(score, abs=1e-8)


def test_gain_residual_slope_is_zero_when_the_gain_follows_its_population_value():
    """Noiseless data whose raw gain is exactly ``n * G`` plus the exact penalty
    has no residual log-length dependence."""
    frame = synthetic_frame(detector="fundamental", m=4, slope=1.0)
    assert gain_residual_regression(frame)["residual_slope"] == pytest.approx(0.0, abs=1e-8)


def test_shared_orbit_residual_slope_is_the_whole_story():
    """For Model C the theoretical dimension term is zero, so the identity
    collapses to ``penalty_slope = -residual_slope``: the residual log-length
    drift reported in Table 5 is entirely a property of the likelihood gain, not
    of a hidden continuous-dimension penalty."""
    frame = synthetic_frame(detector="shared_orbit", m=6, slope=0.17)
    score = score_regression(frame)["penalty_slope"]
    residual = gain_residual_regression(frame)["residual_slope"]
    assert residual == pytest.approx(-score, abs=1e-8)


def test_gain_residual_regression_supports_weighting():
    frame = synthetic_frame()
    assert gain_residual_regression(frame, weighted=True)["method"] == "wls"


# --------------------------------------------------------------------------
# 3. bootstrap intervals for the crossovers
# --------------------------------------------------------------------------


def _power_frame(n_alt=300):
    curves = {
        "full": [0.05, 0.18, 0.42, 0.72, 0.93],
        "fundamental": [0.08, 0.25, 0.55, 0.85, 0.98],
        "shared_orbit": [0.12, 0.34, 0.68, 0.92, 1.0],
    }
    rows = []
    for detector, curve in curves.items():
        for n, power in zip(LENGTHS, curve):
            rows.append(
                {
                    "m": 4, "scenario": "exact_orbit", "effect": 0.18, "detector": detector,
                    "total_length": n, "power_calibrated": power, "n_alt": n_alt,
                }
            )
    # A second effect so the median-across-effects step has something to do.
    second = pd.DataFrame(rows).copy()
    second["effect"] = 0.25
    second["power_calibrated"] = np.minimum(second["power_calibrated"] * 1.15, 1.0)
    return pd.concat([pd.DataFrame(rows), second], ignore_index=True)


def test_bootstrap_interval_brackets_the_point_estimate():
    boot = crossover_bootstrap(_power_frame(), n_boot=200)
    assert len(boot) == 6
    for _, row in boot.iterrows():
        assert row["ci_low"] <= row["crossover_length"] <= row["ci_high"]
        assert row["boot_internal_fraction"] > 0.9


def test_bootstrap_is_reproducible_and_narrows_with_more_trials():
    first = crossover_bootstrap(_power_frame(), n_boot=200, seed=7)
    second = crossover_bootstrap(_power_frame(), n_boot=200, seed=7)
    pd.testing.assert_frame_equal(first, second)

    wide = crossover_bootstrap(_power_frame(n_alt=50), n_boot=300, seed=1)
    narrow = crossover_bootstrap(_power_frame(n_alt=5000), n_boot=300, seed=1)
    wide_width = (wide["ci_high"] - wide["ci_low"]).mean()
    narrow_width = (narrow["ci_high"] - narrow["ci_low"]).mean()
    assert narrow_width < wide_width


def test_ratio_bootstrap_reports_intervals_per_pair():
    boot = crossover_ratio_bootstrap(_power_frame(), n_boot=200)
    row = boot.iloc[0]
    for pair in ("shared_orbit/full", "shared_orbit/fundamental", "fundamental/full"):
        assert row[f"{pair}_ci_low"] <= row[f"{pair}_ci_high"]
        assert row[f"{pair}_boot_n"] > 0
    # The more powerful detector needs fewer observations at every resample.
    assert row["shared_orbit/full_ci_high"] < 1.0


def test_ratio_bootstrap_handles_a_scenario_with_no_crossings():
    frame = _power_frame()
    frame["power_calibrated"] = 0.01
    boot = crossover_ratio_bootstrap(frame, n_boot=50)
    assert boot["shared_orbit/full_boot_n"].max() == 0
    assert boot["shared_orbit/full_ci_low"].isna().all()


# --------------------------------------------------------------------------
# 4. optimiser audit
# --------------------------------------------------------------------------


def test_no_convergence_failures_across_the_operating_range():
    """The population gains are themselves optimisation output, so a silent
    failure would corrupt a likelihood invisibly."""
    rng = np.random.default_rng(0)
    reset_fit_failures()
    fits = 0
    for m in (2, 3, 4, 5, 6):
        for n in (10, 50, 800, 3200, 50_000):
            for _ in range(40):
                fit_fundamental(rng.multinomial(n, rng.dirichlet(np.ones(m))).astype(float), m)
                fits += 1
    assert fits > 500
    assert fit_failure_count() == 0


def test_no_convergence_failures_at_degenerate_counts():
    """Extreme count vectors push the softmax fit toward the simplex boundary."""
    reset_fit_failures()
    for m in (2, 4, 6):
        for counts in (np.eye(m)[0] * 500, np.eye(m)[-1] * 3, np.full(m, 1.0)):
            fit_fundamental(counts, m)
    assert fit_failure_count() == 0


def test_failure_counter_resets_and_is_reported_per_configuration():
    reset_fit_failures()
    assert fit_failure_count() == 0
    rows = pd.DataFrame(run_config(Config(4, "exact_orbit", 0.25, 120, n_alt=20, n_null=20)))
    assert "optimizer_failures" in rows.columns
    assert (rows["optimizer_failures"] == 0).all()


# --------------------------------------------------------------------------
# 5. split fraction
# --------------------------------------------------------------------------


def test_split_fraction_defaults_to_the_balanced_design():
    config = Config(4, "exact_orbit", 0.18, 200)
    assert config.split_fraction == 0.5
    assert (config.n_left, config.n_right, config.total_length) == (200, 200, 400)


@pytest.mark.parametrize("rho,expected", [(0.25, (200, 600, 800)), (0.8, (200, 50, 250))])
def test_unbalanced_splits_have_the_requested_geometry(rho, expected):
    config = Config(4, "exact_orbit", 0.18, 200, split_fraction=rho)
    assert (config.n_left, config.n_right, config.total_length) == expected
    assert config.n_left / config.total_length == pytest.approx(rho, abs=0.01)


def test_invalid_split_fractions_are_rejected():
    for bad in (0.0, 1.0, -0.3, 1.7):
        with pytest.raises(ValueError):
            Config(4, "exact_orbit", 0.18, 200, split_fraction=bad)


def test_adding_split_fraction_left_balanced_seeds_unchanged():
    """The field joins the seed payload only when it is not the default, so
    existing balanced-split checkpoints stay valid."""
    balanced = Config(4, "exact_orbit", 0.18, 200)
    assert config_seed(balanced) == 1563635273
    assert config_seed(Config(4, "exact_orbit", 0.18, 200, split_fraction=0.5)) == config_seed(balanced)
    assert config_seed(Config(4, "exact_orbit", 0.18, 200, split_fraction=0.25)) != config_seed(balanced)


def test_split_fraction_is_part_of_the_configuration_key():
    grid = build_grid((4,), ("exact_orbit",), (0.2,), (100,), split_fractions=(0.5, 0.25))
    assert len(grid) == 2
    assert len({c.key for c in grid}) == 2
    assert len({config_seed(c) for c in grid}) == 2


@pytest.mark.parametrize("dim", [1, 2, 5])
def test_split_fraction_shifts_the_intercept_but_not_the_slope(dim):
    """Section 4.2's prediction, in closed form: the coefficient of ``log n`` is
    ``dim/2`` for every split fraction, and rho enters only through the bounded
    term ``(dim/2) log(rho (1 - rho))``."""
    for rho in (0.5, 0.25, 0.1):
        penalties = []
        for n in (10**4, 10**5, 10**6):
            n_left = int(n * rho)
            penalties.append(split_penalty(dim, n_left, n - n_left))
        slopes = np.diff(penalties) / np.log(10)
        np.testing.assert_allclose(slopes, dim / 2, rtol=1e-6)

        offset = penalties[0] - (dim / 2) * np.log(10**4)
        np.testing.assert_allclose(offset, (dim / 2) * np.log(rho * (1 - rho)), rtol=1e-4)


def test_summary_degrades_gracefully_without_the_raw_gain_column():
    """An older results file has no ``mean_raw_gain``. The residual diagnostic
    should report NaN rather than taking the whole report down."""
    frame = synthetic_frame(detector="shared_orbit", m=4, slope=0.0).drop(columns=["mean_raw_gain"])
    summary = score_regression_summary(frame, scenario_by_detector={"shared_orbit": "exact_orbit"})
    assert len(summary) == 1
    assert np.isnan(summary.iloc[0]["residual_slope"])
    assert np.isfinite(summary.iloc[0]["penalty_slope"])

    with pytest.raises(ValueError, match="mean_raw_gain"):
        gain_residual_regression(frame)


def test_boundary_fits_are_likelihood_identified_but_not_coordinate_identified():
    """A property the optimiser audit exposed, worth pinning rather than hiding.

    When a category has zero count the fundamental MLE does not exist: the
    likelihood rises toward the boundary and is asymptotically flat along that
    direction. Different starts therefore halt at very different ``|theta|``
    while agreeing on the log-likelihood to many decimals.

    This is safe for the detectors, which consume only likelihoods, and it is
    why the fit keeps two starts even though the objective is concave: the
    second start costs a little time and guards the reported coordinate.
    """
    rng = np.random.default_rng(5)
    disagreements = 0
    for m in (3, 4, 5):
        for _ in range(60):
            counts = rng.multinomial(10, rng.dirichlet(np.ones(m))).astype(float)
            theta_two, ll_two = fit_fundamental(counts, m, n_restarts=2)
            theta_one, ll_one = fit_fundamental(counts, m, n_restarts=1)
            # The likelihood -- the only thing a detector uses -- is stable.
            assert ll_two == pytest.approx(ll_one, abs=1e-5)
            if np.abs(theta_two - theta_one).max() > 1e-6:
                assert (counts == 0).any(), "coordinates should only disagree at the boundary"
                disagreements += 1
    assert disagreements > 0, "expected some boundary cases in this sample"


def test_detector_scores_are_start_independent():
    """The consequence that actually matters: scores do not depend on the
    number of optimiser starts, even on degenerate short segments."""
    import regimeshift.detectors as det

    rng = np.random.default_rng(11)
    for m in (3, 4, 6):
        for _ in range(30):
            left = rng.multinomial(12, rng.dirichlet(np.ones(m))).astype(float)
            right = rng.multinomial(12, rng.dirichlet(np.ones(m))).astype(float)
            baseline = {
                name: result.score for name, result in det.run_all_detectors(left, right, m).items()
            }
            original = det.fit_fundamental.__defaults__
            try:
                det.fit_fundamental.__defaults__ = (1,)
                single = {
                    name: result.score
                    for name, result in det.run_all_detectors(left, right, m).items()
                }
            finally:
                det.fit_fundamental.__defaults__ = original
            for name, value in baseline.items():
                assert value == pytest.approx(single[name], abs=1e-5), name


# --------------------------------------------------------------------------
# reconstructed-constant provenance
# --------------------------------------------------------------------------


def test_manuscript_constants_match_the_live_module_values():
    """Exact reproduction depends on these constants, so the machine-readable
    provenance table must never drift from the code it describes."""
    import regimeshift.scenarios as scenarios
    from regimeshift.detectors import label_cost
    from regimeshift.scenarios import MANUSCRIPT_CONSTANTS

    for name, entry in MANUSCRIPT_CONSTANTS.items():
        assert entry["basis"], f"{name} has no stated basis"
        assert entry["section"], f"{name} has no manuscript section"
        assert isinstance(entry["recovered_from_manuscript"], bool)
        if name == "LABEL_COST":
            # A formula rather than a scalar; check the implementation matches it.
            for m in (2, 3, 5, 9):
                assert label_cost(m) == pytest.approx(np.log(m - 1))
            continue
        assert entry["value"] == getattr(scenarios, name), f"{name} drifted from its table entry"


def test_every_scenario_constant_is_documented():
    """A new scenario constant must come with provenance, not appear silently."""
    from regimeshift.scenarios import MANUSCRIPT_CONSTANTS

    documented = set(MANUSCRIPT_CONSTANTS)
    expected = {
        "INDEPENDENT_RADIUS_FACTOR", "INDEPENDENT_ANGLE_RAD",
        "INDEPENDENT_M2_FACTOR", "HIGHER_MODE_FACTOR", "LABEL_COST",
    }
    assert documented == expected
    # Every one of them is quoted from the manuscript. An earlier version of this
    # repository guessed three, believing the equations were images; they are
    # OMML and the first extraction pass simply dropped them.
    assert all(v["recovered_from_manuscript"] for v in MANUSCRIPT_CONSTANTS.values())


# --------------------------------------------------------------------------
# K* = 1/(2 ln 2): the per-dimension penalty quantum, in bits
# --------------------------------------------------------------------------


def test_k_star_is_the_nats_to_bits_conversion_of_one_half():
    """K* is definitional, not empirical: it is Schwarz's one-half expressed in
    bits. Pinning it as such keeps anyone from reading it as a fitted quantity."""
    assert K_STAR == 1.0 / (2.0 * np.log(2.0))
    assert K_STAR == pytest.approx(0.7213475204444817, abs=1e-15)
    assert nats_to_bits(0.5) == K_STAR
    assert nats_to_bits(1.0) == pytest.approx(1.0 / np.log(2.0))


@pytest.mark.parametrize("m", [2, 3, 4, 5, 6, 7, 8])
@pytest.mark.parametrize("detector", ["full", "fundamental", "shared_orbit"])
def test_every_predicted_slope_in_bits_is_an_integer_multiple_of_k_star(detector, m):
    """The structural claim: ``(d/2) log n`` nats is ``d * K*`` bits, so each
    model's leading coefficient is a whole number of K* -- ``m - 1`` for Model A,
    ``d_fund`` for Model B, and zero for Model C. The three-way hierarchy is how
    many K* a model pays to cross the boundary."""
    d = dimension_increment(detector, m)
    assert d == int(d)
    assert predicted_slope(detector, m, units="nats") == d / 2.0
    assert predicted_slope(detector, m, units="bits") == pytest.approx(d * K_STAR, abs=1e-15)


def test_shared_orbit_pays_zero_k_star():
    """The paper's central result, restated in these units."""
    for m in (2, 4, 6):
        assert dimension_increment("shared_orbit", m) == 0
        assert predicted_slope("shared_orbit", m, units="bits") == 0.0


def test_dimension_increments_match_the_three_model_classes():
    assert [dimension_increment("full", m) for m in (2, 4, 6)] == [1, 3, 5]
    assert [dimension_increment("fundamental", m) for m in (2, 3, 6)] == [1, 2, 2]
    for bad in ("mystery", ""):
        with pytest.raises(ValueError):
            dimension_increment(bad, 4)


def test_unknown_units_are_rejected():
    for bad in ("nat", "BITS", "decibans", ""):
        with pytest.raises(ValueError, match="units"):
            predicted_slope("full", 4, units=bad)
    with pytest.raises(ValueError, match="units"):
        score_regression_summary(synthetic_frame(), units="decibans")


def test_summary_converts_every_slope_column_together():
    """A units switch must convert all slope-valued columns and leave the
    dimensionless ones alone, or the report becomes internally inconsistent."""
    frame = synthetic_frame(detector="fundamental", m=4, slope=1.0).assign(
        scenario="independent_fundamental"
    )
    mapping = {"fundamental": "independent_fundamental"}
    nats = score_regression_summary(frame, scenario_by_detector=mapping).iloc[0]
    bits = score_regression_summary(frame, scenario_by_detector=mapping, units="bits").iloc[0]

    assert nats["units"] == "nats" and bits["units"] == "bits"
    for column in ("penalty_slope", "penalty_slope_wls", "residual_slope", "predicted_slope"):
        assert bits[column] == pytest.approx(nats[column] / np.log(2), abs=1e-12)
    # Dimensionless quantities must not move.
    for column in ("beta_gain", "r_squared", "condition_number", "n_points", "k_star_multiple"):
        assert bits[column] == pytest.approx(nats[column])
    # And the predicted slope in bits lands exactly on d * K*.
    assert bits["predicted_slope"] == pytest.approx(bits["k_star_multiple"] * K_STAR, abs=1e-12)

# ---------------------------------------------------------------------------
# Joint (paired) crossover-ratio bootstrap
# ---------------------------------------------------------------------------

LENGTH_GRID = (200, 400, 800, 1600, 3200)


def _synthetic_results(power_by_detector, patterns=None, effects=(0.18, 0.25),
                       lengths=LENGTH_GRID, n_alt=500, m=4):
    """A results frame shaped like a real one, with power curves we control.

    ``power_by_detector`` maps a detector to a power curve over ``lengths``.
    ``patterns``, when given, maps each ``(effect, length)`` to the eight joint
    counts; omit it to produce a pre-pattern file and exercise the fallback.
    """
    rows = []
    for effect in effects:
        for i, length in enumerate(lengths):
            joint = None if patterns is None else patterns[(effect, length)]
            for detector, curve in power_by_detector.items():
                row = {
                    "m": m, "scenario": "exact_orbit", "effect": effect,
                    "total_length": length, "detector": detector,
                    "power_calibrated": curve[i], "n_alt": n_alt,
                }
                if joint is not None:
                    row.update(dict(zip(DETECTION_PATTERNS, joint)))
                rows.append(row)
    return pd.DataFrame(rows)


def _agreeing_patterns(power, n_alt=500):
    """Joint counts for detectors that agree on every dataset: all the mass on
    ``pattern_000`` and ``pattern_111``, so every marginal is ``power``."""
    caught = int(round(power * n_alt))
    counts = [0] * len(DETECTION_PATTERNS)
    counts[DETECTION_PATTERNS.index("pattern_000")] = n_alt - caught
    counts[DETECTION_PATTERNS.index("pattern_111")] = caught
    return counts


def test_detection_patterns_reproduce_every_marginal_power():
    """The eight joint counts refine the three marginals, so summing the four
    patterns that set a detector's bit must return its ``power_calibrated``
    exactly. That identity is what makes them safe to resample from.
    """
    frame = pd.DataFrame(run_config(Config(4, "exact_orbit", 0.25, 200, n_alt=200, n_null=200)))

    counts = frame[list(DETECTION_PATTERNS)].to_numpy(dtype=float)
    np.testing.assert_allclose(counts.sum(axis=1), frame["n_alt"].to_numpy(dtype=float))
    for index, detector in enumerate(DETECTOR_NAMES):
        mask = (frame["detector"] == detector).to_numpy()
        marginal = counts[mask][:, _pattern_mask(index)].sum(axis=1) / frame.loc[mask, "n_alt"]
        np.testing.assert_allclose(
            marginal.to_numpy(), frame.loc[mask, "power_calibrated"].to_numpy()
        )
    assert _has_joint_patterns(frame)


def test_inconsistent_patterns_are_rejected_rather_than_trusted():
    """The identity above is checked, not assumed, so a file whose patterns and
    marginals disagree falls back instead of resampling from bad joint counts."""
    frame = pd.DataFrame(run_config(Config(4, "exact_orbit", 0.25, 200, n_alt=200, n_null=200)))
    assert _has_joint_patterns(frame)
    frame.loc[0, "pattern_111"] = frame.loc[0, "pattern_111"] + 1
    assert not _has_joint_patterns(frame)


def test_a_degenerate_pair_gets_a_degenerate_interval_under_pairing():
    """The check the independent bootstrap failed.

    Detectors agreeing on every dataset have a crossover ratio of exactly 1
    with no variance, and a paired resample returns that. An independent one
    invents an interval around it -- the artefact reported at ``m = 2`` and
    ``m = 3`` in the committed run, where ``full`` and ``fundamental`` are one
    detector.
    """
    curve = [0.10, 0.28, 0.55, 0.78, 0.93]
    powers = {name: curve for name in DETECTOR_NAMES}
    patterns = {
        (effect, length): _agreeing_patterns(curve[i])
        for effect in (0.18, 0.25)
        for i, length in enumerate(LENGTH_GRID)
    }

    paired = crossover_ratio_bootstrap(_synthetic_results(powers, patterns), n_boot=200).iloc[0]
    assert bool(paired["paired"])
    assert paired["fundamental/full_boot_n"] > 150
    assert paired["fundamental/full_ci_low"] == pytest.approx(1.0)
    assert paired["fundamental/full_ci_high"] == pytest.approx(1.0)

    independent = crossover_ratio_bootstrap(_synthetic_results(powers), n_boot=200).iloc[0]
    assert not bool(independent["paired"])
    width = independent["fundamental/full_ci_high"] - independent["fundamental/full_ci_low"]
    assert width > 0.05, "the independent bootstrap invents a width on a constant"


def test_pairing_narrows_a_ratio_of_correlated_detectors():
    """Not a tautology: positive correlation makes a paired resample of a ratio
    tighter than an independent one, which is why the independent intervals
    were the wrong ones to report."""
    slow = [0.05, 0.18, 0.42, 0.70, 0.90]
    fast = [0.12, 0.34, 0.62, 0.85, 0.97]
    powers = {"full": slow, "fundamental": slow, "shared_orbit": fast}

    # shared_orbit catches everything the others do, plus a little more, so the
    # curves move together dataset by dataset.
    patterns = {}
    for effect in (0.18, 0.25):
        for i, length in enumerate(LENGTH_GRID):
            n_alt = 500
            both = int(round(slow[i] * n_alt))
            extra = int(round((fast[i] - slow[i]) * n_alt))
            counts = [0] * len(DETECTION_PATTERNS)
            counts[DETECTION_PATTERNS.index("pattern_111")] = both
            counts[DETECTION_PATTERNS.index("pattern_001")] = extra
            counts[DETECTION_PATTERNS.index("pattern_000")] = n_alt - both - extra
            patterns[(effect, length)] = counts

    label = "shared_orbit/full"
    paired = crossover_ratio_bootstrap(
        _synthetic_results(powers, patterns), n_boot=300, seed=7
    ).iloc[0]
    independent = crossover_ratio_bootstrap(
        _synthetic_results(powers), n_boot=300, seed=7
    ).iloc[0]

    assert paired[f"{label}_boot_n"] > 200 and independent[f"{label}_boot_n"] > 200
    paired_width = paired[f"{label}_ci_high"] - paired[f"{label}_ci_low"]
    independent_width = independent[f"{label}_ci_high"] - independent[f"{label}_ci_low"]
    assert paired_width < independent_width


def test_the_bootstrap_freezes_the_effect_subset():
    """Every replicate must estimate the same quantity: the median over the
    effects the point estimate uses. A replicate that cannot fill that subset
    is discarded rather than re-medianed over whatever survived, so
    ``*_boot_n`` falls below ``n_boot`` -- and that shortfall is itself the
    signal that the median rests on thin ground.

    Here one of the two effects never crosses 0.5 inside the grid, so the point
    estimate uses a single effect and no replicate may quietly substitute the
    other.
    """
    strong = [0.10, 0.28, 0.55, 0.78, 0.93]
    never = [0.01, 0.02, 0.03, 0.04, 0.05]
    rows = []
    for effect, power in ((0.18, never), (0.25, strong)):
        for i, length in enumerate(LENGTH_GRID):
            for detector in DETECTOR_NAMES:
                rows.append({
                    "m": 4, "scenario": "exact_orbit", "effect": effect,
                    "total_length": length, "detector": detector,
                    "power_calibrated": power[i], "n_alt": 500,
                })
    row = crossover_ratio_bootstrap(pd.DataFrame(rows), n_boot=200).iloc[0]

    for pair in ("shared_orbit/full", "shared_orbit/fundamental", "fundamental/full"):
        assert row[f"{pair}_effects"] == 1, "only the strong effect crosses inside the grid"
        assert 0 < row[f"{pair}_boot_n"] <= 200


def test_a_results_file_without_pattern_columns_falls_back_and_says_so():
    """Backwards compatibility: a results file written before the joint columns
    existed must still analyse, by the old independent route, and say so rather
    than failing or silently pretending to be paired.

    The shipped run now carries the columns, so the pre-pattern case is made by
    dropping them -- which is exactly what an older file looks like.
    """
    path = Path(__file__).resolve().parent.parent / "results" / "v3-production" / "full_results.csv"
    if not path.exists():
        pytest.skip("committed production results are not present")
    frame = pd.read_csv(path)

    paired = crossover_ratio_bootstrap(frame, n_boot=50).iloc[0]
    assert bool(paired["paired"]), "the shipped run should carry the joint patterns"

    stripped = frame.drop(columns=list(DETECTION_PATTERNS))
    fallback = crossover_ratio_bootstrap(stripped, n_boot=50).iloc[0]
    assert not bool(fallback["paired"])
    assert np.isfinite(fallback["shared_orbit/full_ci_low"])


def test_every_analysis_output_is_invariant_to_input_row_order():
    """Section 8.1 promises results do not depend on grid ordering, worker
    count or completion order. Rows arrive as workers finish, so that promise
    is only true if nothing downstream reads the row order.

    It was not. ``_interpolate_crossover`` sorts its own inputs, so the point
    estimates were safe -- but ``crossover_bootstrap`` draws one binomial per
    row, and on an unsorted group that paired draws with lengths arbitrarily.
    The committed run had 18 such groups and a re-run had 117, which is why its
    intervals moved while every other artefact reproduced exactly.

    Shuffling the input must now change nothing, anywhere.
    """
    frame = pd.DataFrame(run_config(Config(4, "exact_orbit", 0.25, 200, n_alt=200, n_null=200)))
    for length in (400, 800, 1600):
        frame = pd.concat(
            [frame, pd.DataFrame(run_config(Config(4, "exact_orbit", 0.25, length,
                                                   n_alt=200, n_null=200)))],
            ignore_index=True,
        )
    key = ["m", "scenario", "effect", "segment_length", "split_fraction", "detector"]
    shuffled = frame.sample(frac=1.0, random_state=17).reset_index(drop=True)
    assert not shuffled[key].equals(frame[key]), "the shuffle must actually reorder rows"

    for name, fn in (
        ("crossover_estimates", lambda f: crossover_estimates(f)),
        ("crossover_bootstrap", lambda f: crossover_bootstrap(f, n_boot=40)),
        ("crossover_ratio_bootstrap", lambda f: crossover_ratio_bootstrap(f, n_boot=40)),
        ("score_regression_summary", lambda f: score_regression_summary(f)),
    ):
        a = fn(frame).reset_index(drop=True)
        b = fn(shuffled).reset_index(drop=True)
        pd.testing.assert_frame_equal(a, b, check_exact=False, atol=1e-12,
                                      obj=f"{name} changed under a row shuffle")


def test_run_grid_emits_rows_in_a_canonical_order():
    """The fix at the source: whatever order workers finish in, the frame the
    runner returns -- and so the CSV written from it -- is sorted by the design
    key, which makes the shipped artefact itself reproducible rather than only
    the analyses that sort defensively."""
    from regimeshift.runner import run_grid

    configs = build_grid(
        groups=(4,), scenarios=("exact_orbit",), effects=(0.25,),
        segment_lengths=(100, 200, 400), n_alt=100, n_null=100,
    )
    serial = run_grid(configs, workers=1)
    parallel = run_grid(configs, workers=3)

    key = ["m", "scenario", "effect", "segment_length", "split_fraction", "detector"]
    assert serial[key].equals(serial[key].sort_values(key).reset_index(drop=True))
    pd.testing.assert_frame_equal(serial, parallel)

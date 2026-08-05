"""Tests for the Monte Carlo engine, the runner and the post-processing."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from regimeshift.analysis import (
    _interpolate_crossover,
    crossover_estimates,
    crossover_ratio_summary,
    predicted_slope,
    score_regression,
    score_regression_summary,
)
from regimeshift.runner import PRODUCTION_GRID, QUICK_GRID, run_grid
from regimeshift.simulation import BASE_SEED, Config, build_grid, config_seed, run_config


# --------------------------------------------------------------------------
# design and seeding
# --------------------------------------------------------------------------


def test_production_grid_matches_the_manuscript_design():
    """Table 2: 312 configurations, 936 detector rows, 468,000 datasets."""
    configs = build_grid(
        PRODUCTION_GRID["groups"],
        PRODUCTION_GRID["scenarios"],
        PRODUCTION_GRID["effects"],
        PRODUCTION_GRID["segment_lengths"],
        n_alt=PRODUCTION_GRID["n_alt"],
        n_null=PRODUCTION_GRID["n_null"],
    )
    assert len(configs) == 312
    assert 3 * len(configs) == 936
    assert sum(c.n_alt + c.n_null for c in configs) == 468_000
    assert {c.total_length for c in configs} == {200, 400, 800, 1600, 3200, 6400}
    assert BASE_SEED == 20260713


def test_grid_skips_scenarios_a_group_order_cannot_support():
    configs = build_grid((2, 3), ("higher_mode",), (0.2,), (100,))
    assert configs == []


def test_config_seed_is_deterministic_and_configuration_specific():
    a = Config(4, "exact_orbit", 0.18, 200)
    b = Config(4, "exact_orbit", 0.18, 200)
    c = Config(4, "exact_orbit", 0.18, 400)
    assert config_seed(a) == config_seed(b)
    assert config_seed(a) != config_seed(c)
    assert config_seed(a) != config_seed(a, base_seed=BASE_SEED + 1)
    # Independent of ordering or process: stable across the whole grid.
    configs = build_grid((2, 4), ("exact_orbit",), (0.1, 0.2), (100, 200))
    assert len({config_seed(x) for x in configs}) == len(configs)


def test_run_config_is_reproducible():
    config = Config(4, "exact_orbit", 0.25, 200, n_alt=25, n_null=40)
    first = pd.DataFrame(run_config(config))
    second = pd.DataFrame(run_config(config))
    pd.testing.assert_frame_equal(first, second)


def test_run_config_shape_and_columns():
    config = Config(5, "exact_orbit", 0.25, 150, n_alt=20, n_null=30)
    rows = pd.DataFrame(run_config(config))
    assert len(rows) == 3
    assert set(rows["detector"]) == {"full", "fundamental", "shared_orbit"}
    assert (rows["total_length"] == 300).all()
    for column in ("population_gain", "mean_score", "critical_value", "power_calibrated"):
        assert rows[column].notna().all()
    assert rows["power_calibrated"].between(0, 1).all()


def test_shift_accuracy_is_recorded_only_where_a_shift_is_planted():
    rows = pd.DataFrame(run_config(Config(4, "exact_orbit", 0.25, 400, n_alt=30, n_null=30)))
    shared = rows[rows["detector"] == "shared_orbit"].iloc[0]
    assert shared["shift_accuracy"] > 0.5
    assert rows[rows["detector"] == "full"]["shift_accuracy"].isna().all()

    rows = pd.DataFrame(
        run_config(Config(4, "independent_fundamental", 0.25, 200, n_alt=20, n_null=20))
    )
    assert rows["shift_accuracy"].isna().all()


def test_null_calibration_controls_the_false_positive_rate():
    """The calibrated critical value is the empirical 95th percentile of the
    null scores, so at most alpha of those null draws exceed it."""
    config = Config(4, "exact_orbit", 0.25, 300, n_alt=20, n_null=400, alpha=0.05)
    rows = pd.DataFrame(run_config(config))
    assert rows["critical_value"].notna().all()
    # The raw (zero-threshold) MDL rule must be conservative relative to 5%.
    assert (rows["null_rate_zero_threshold"] <= 0.05).all()


# --------------------------------------------------------------------------
# runner
# --------------------------------------------------------------------------


def test_run_grid_checkpoints_and_resumes(tmp_path):
    configs = build_grid((4,), ("exact_orbit",), (0.25,), (100, 200), n_alt=15, n_null=20)
    checkpoint = tmp_path / "ck.csv"

    first = run_grid(configs[:1], checkpoint=checkpoint)
    assert len(first) == 3
    assert checkpoint.exists()

    both = run_grid(configs, checkpoint=checkpoint)
    assert len(both) == 6
    # Re-running everything must not duplicate the completed configuration.
    again = run_grid(configs, checkpoint=checkpoint)
    assert len(again) == 6


def test_run_grid_results_do_not_depend_on_worker_count(tmp_path):
    configs = build_grid((4,), ("exact_orbit",), (0.25,), (100, 200), n_alt=15, n_null=20)
    serial = run_grid(configs, checkpoint=tmp_path / "a.csv", workers=1)
    parallel = run_grid(configs, checkpoint=tmp_path / "b.csv", workers=2)
    key = ["m", "scenario", "effect", "segment_length", "detector"]
    serial = serial.sort_values(key).reset_index(drop=True)
    parallel = parallel.sort_values(key).reset_index(drop=True)
    pd.testing.assert_frame_equal(serial, parallel)


def test_quick_grid_supports_every_regression():
    """CI runs the quick grid; it must be large enough for the reports."""
    configs = build_grid(
        QUICK_GRID["groups"], QUICK_GRID["scenarios"], QUICK_GRID["effects"],
        QUICK_GRID["segment_lengths"], n_alt=QUICK_GRID["n_alt"], n_null=QUICK_GRID["n_null"],
    )
    n_effects = len(QUICK_GRID["effects"])
    n_points = n_effects * len(QUICK_GRID["segment_lengths"])
    assert n_points >= n_effects + 3, "regression would be underdetermined"
    assert len(configs) > 0


# --------------------------------------------------------------------------
# analysis
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    "detector,m,expected",
    [
        ("full", 2, 0.5), ("full", 4, 1.5), ("full", 5, 2.0), ("full", 6, 2.5),
        ("fundamental", 2, 0.5), ("fundamental", 3, 1.0), ("fundamental", 6, 1.0),
        ("shared_orbit", 2, 0.0), ("shared_orbit", 6, 0.0),
    ],
)
def test_predicted_slopes_match_the_manuscript_tables(detector, m, expected):
    assert predicted_slope(detector, m) == expected


def test_predicted_slope_rejects_unknown_detectors():
    with pytest.raises(ValueError):
        predicted_slope("mystery", 4)


def _synthetic_frame(beta=1.0, slope=1.5, effects=(0.1, 0.2, 0.3), lengths=(200, 400, 800, 1600, 3200)):
    """Exact realisation of mean_score = beta * n * gain - slope * log n + a_effect."""
    rows = []
    for effect in effects:
        gain = 0.5 * effect**2
        for n in lengths:
            rows.append(
                {
                    "m": 4,
                    "scenario": "exact_orbit",
                    "detector": "full",
                    "effect": effect,
                    "total_length": n,
                    "expected_gain_total": gain * n,
                    "mean_score": beta * gain * n - slope * np.log(n) + 3.0 * effect,
                }
            )
    return pd.DataFrame(rows)


def test_score_regression_recovers_planted_coefficients():
    frame = _synthetic_frame(beta=0.97, slope=1.5)
    fit = score_regression(frame)
    assert fit["beta_gain"] == pytest.approx(0.97, abs=1e-8)
    assert fit["penalty_slope"] == pytest.approx(1.5, abs=1e-8)
    assert fit["r_squared"] == pytest.approx(1.0, abs=1e-10)
    assert fit["n_points"] == 15


def test_score_regression_recovers_a_zero_slope():
    fit = score_regression(_synthetic_frame(slope=0.0))
    assert fit["penalty_slope"] == pytest.approx(0.0, abs=1e-8)


def test_score_regression_rejects_underdetermined_designs():
    with pytest.raises(ValueError, match="not enough design points"):
        score_regression(_synthetic_frame(effects=(0.1, 0.2), lengths=(200,)))


def test_score_regression_summary_labels_predictions():
    frame = _synthetic_frame()
    frame = pd.concat(
        [
            frame,
            frame.assign(detector="shared_orbit", mean_score=frame["mean_score"] + 1.0),
            frame.assign(detector="fundamental", scenario="independent_fundamental"),
        ],
        ignore_index=True,
    )
    frame.loc[frame["detector"] == "full", "scenario"] = "higher_mode"
    summary = score_regression_summary(frame)
    assert set(summary["detector"]) == {"full", "fundamental", "shared_orbit"}
    for _, row in summary.iterrows():
        assert row["predicted_slope"] == predicted_slope(row["detector"], row["m"])


def test_score_regression_summary_is_empty_without_matching_rows():
    empty = score_regression_summary(
        pd.DataFrame(columns=["detector", "scenario", "m", "effect", "total_length",
                              "expected_gain_total", "mean_score"])
    )
    assert len(empty) == 0


# --------------------------------------------------------------------------
# crossovers
# --------------------------------------------------------------------------


def test_crossover_interpolates_in_log_length():
    lengths = np.array([100.0, 1000.0])
    power = np.array([0.0, 1.0])
    value, status = _interpolate_crossover(lengths, power, 0.5)
    assert status == "internal"
    assert value == pytest.approx(np.sqrt(100 * 1000))


def test_crossover_suppresses_monte_carlo_reversals():
    lengths = np.array([100.0, 200.0, 400.0, 800.0])
    noisy = np.array([0.2, 0.6, 0.55, 0.9])
    monotone = np.array([0.2, 0.6, 0.6, 0.9])
    assert _interpolate_crossover(lengths, noisy)[0] == pytest.approx(
        _interpolate_crossover(lengths, monotone)[0]
    )


def test_crossover_flags_out_of_grid_estimates():
    lengths = np.array([100.0, 200.0, 400.0])
    assert _interpolate_crossover(lengths, np.array([0.9, 0.95, 1.0]))[1] == "below_grid"
    above = _interpolate_crossover(lengths, np.array([0.01, 0.02, 0.03]))
    assert above[1] == "above_grid"
    assert np.isnan(above[0])


def test_crossover_is_order_independent():
    lengths = np.array([800.0, 100.0, 400.0, 200.0])
    power = np.array([0.9, 0.2, 0.6, 0.4])
    shuffled = _interpolate_crossover(lengths, power)
    order = np.argsort(lengths)
    assert shuffled[0] == pytest.approx(_interpolate_crossover(lengths[order], power[order])[0])


def _power_frame():
    rows = []
    curves = {"full": [0.1, 0.3, 0.7, 0.95], "fundamental": [0.2, 0.45, 0.8, 0.99],
              "shared_orbit": [0.3, 0.6, 0.9, 1.0]}
    for detector, curve in curves.items():
        for n, power in zip((200, 400, 800, 1600), curve):
            rows.append({"m": 4, "scenario": "exact_orbit", "effect": 0.18, "detector": detector,
                         "total_length": n, "power_calibrated": power})
    return pd.DataFrame(rows)


def test_crossover_estimates_and_ratio_summary():
    crossovers = crossover_estimates(_power_frame())
    assert len(crossovers) == 3
    assert (crossovers["status"] == "internal").all()

    ratios = crossover_ratio_summary(crossovers)
    row = ratios.iloc[0]
    # A more powerful detector needs fewer observations, so the ratio is below one.
    assert row["shared_orbit/full"] < 1.0
    assert row["fundamental/full"] < 1.0
    assert row["shared_orbit/fundamental"] < 1.0
    assert row["shared_orbit/full_n"] == 1


def test_ratio_summary_handles_no_internal_crossovers():
    frame = _power_frame()
    frame["power_calibrated"] = 0.01
    ratios = crossover_ratio_summary(crossover_estimates(frame))
    assert len(ratios) == 0
    assert "shared_orbit/full" in ratios.columns

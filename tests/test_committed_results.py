"""Integrity checks for the committed production results.

The repository ships one full production run so that the manuscript's numbers
can be read, re-analysed and cited without anyone re-running a multi-hour grid.
That only helps if the shipped files are demonstrably the ones the manifest
describes, so these tests verify the provenance record rather than trusting it.

They skip cleanly when the results are absent, so a fresh checkout or a
partial clone still has a green suite.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from regimeshift.manifest import file_digest
from regimeshift.runner import PRODUCTION_GRID
from regimeshift.simulation import BASE_SEED, build_grid

RESULTS = Path(__file__).resolve().parent.parent / "results" / "v3-production"
MANIFEST = RESULTS / "run_manifest.json"

pytestmark = pytest.mark.skipif(
    not MANIFEST.exists(), reason="committed production results are not present"
)


@pytest.fixture(scope="module")
def manifest() -> dict:
    return json.loads(MANIFEST.read_text())


def test_every_shipped_file_matches_its_recorded_checksum(manifest):
    """The point of shipping a manifest: a hand-edited CSV must be detectable."""
    assert manifest["files"], "manifest lists no files"
    for name, record in manifest["files"].items():
        path = RESULTS / name
        assert path.exists(), f"{name} is in the manifest but missing from disk"
        assert path.stat().st_size == record["bytes"], f"{name} changed size"
        assert file_digest(path) == record["sha256"], f"{name} does not match its checksum"


def test_no_unrecorded_csv_files(manifest):
    on_disk = {p.name for p in RESULTS.glob("*.csv") if p.name != "checkpoint.csv"}
    assert on_disk == set(manifest["files"]), "results directory and manifest disagree"


def test_manifest_records_the_provenance_a_reviewer_needs(manifest):
    assert manifest["git"]["commit"], "no commit recorded"
    assert len(manifest["git"]["commit"]) == 40
    env = manifest["environment"]
    assert env["python"] and env["platform"]
    for package in ("numpy", "scipy", "pandas"):
        assert env["packages"][package], f"no version recorded for {package}"
    assert manifest["base_seed"] == BASE_SEED


def test_manifest_design_matches_the_manuscript_grid(manifest):
    """Table 2: 312 configurations, 936 detector rows, 468,000 datasets."""
    assert manifest["grid"] == "production"
    assert manifest["configurations"] == 312
    assert manifest["detector_rows"] == 936
    assert manifest["simulated_datasets"] == 468_000

    design = manifest["design"]
    for key, expected in PRODUCTION_GRID.items():
        actual = design[key]
        if isinstance(expected, tuple):
            assert list(expected) == list(actual), f"{key} differs from PRODUCTION_GRID"
        else:
            assert expected == actual, f"{key} differs from PRODUCTION_GRID"


def test_shipped_results_have_the_expected_shape():
    results = pd.read_csv(RESULTS / "full_results.csv")
    assert len(results) == 936
    assert set(results["detector"]) == {"full", "fundamental", "shared_orbit"}
    assert set(results["m"]) == {2, 3, 4, 5, 6}
    assert results["optimizer_failures"].sum() == 0
    assert results["mean_score"].notna().all()
    # The design is exactly the one build_grid produces.
    configs = build_grid(
        PRODUCTION_GRID["groups"], PRODUCTION_GRID["scenarios"],
        PRODUCTION_GRID["effects"], PRODUCTION_GRID["segment_lengths"],
        n_alt=PRODUCTION_GRID["n_alt"], n_null=PRODUCTION_GRID["n_null"],
    )
    assert len(configs) * 3 == len(results)


def test_shipped_regression_summary_tracks_the_predictions():
    """A smoke test on the shipped numbers themselves: the full and fundamental
    slopes should sit near their predicted dimensions, and the shared-orbit
    slope far below both."""
    summary = pd.read_csv(RESULTS / "score_regression_summary.csv")
    for _, row in summary.iterrows():
        if row["detector"] in ("full", "fundamental"):
            assert abs(row["penalty_slope"] - row["predicted_slope"]) < 0.35, row.to_dict()
        else:
            assert abs(row["penalty_slope"]) < 0.4, row.to_dict()


# ---------------------------------------------------------------------------
# Claims the manuscript states about this run, pinned against the run itself
# ---------------------------------------------------------------------------


def test_models_a_and_b_are_the_same_detector_at_m_two_and_three():
    """At ``m = 2, 3`` the fundamental component is the whole nontrivial
    tangent space, so ``full`` and ``fundamental`` are one hypothesis.

    In the shipped run their calibrated power is not merely close but *equal*
    at every design point. That is what makes the ``fundamental/full`` crossover
    ratio identically 1 with no variance, and hence what makes the bootstrap
    interval reported for it in Section 9.4 an artefact of resampling the two
    detectors independently rather than a check that passed.
    """
    results = pd.read_csv(RESULTS / "full_results.csv")
    keys = ["m", "scenario", "effect", "total_length"]

    for m in (2, 3):
        wide = results[results["m"] == m].pivot_table(
            index=keys, columns="detector", values="power_calibrated"
        )
        assert (wide["full"] == wide["fundamental"]).all()

        scores = results[results["m"] == m].pivot_table(
            index=keys, columns="detector", values="mean_score"
        )
        assert (scores["full"] - scores["fundamental"]).abs().max() < 1e-9

    # and at m = 4 they genuinely differ, so the check above has content
    wide = results[results["m"] == 4].pivot_table(
        index=keys, columns="detector", values="power_calibrated"
    )
    assert (wide["full"] != wide["fundamental"]).any()


def test_the_independent_ratio_bootstrap_widens_a_constant():
    """Quantifies the artefact the previous test explains.

    ``fundamental/full`` is exactly 1 at ``m = 2, 3``; the committed interval is
    roughly +-9%. Pinned so the manuscript's caveat cannot drift from its size.
    """
    summary = pd.read_csv(RESULTS / "crossover_ratio_summary.csv")
    boot = pd.read_csv(RESULTS / "crossover_ratio_bootstrap.csv")
    exact = summary[summary["scenario"] == "exact_orbit"].set_index("m")
    exact_boot = boot[boot["scenario"] == "exact_orbit"].set_index("m")

    for m in (2, 3):
        assert exact.loc[m, "fundamental/full"] == 1.0
        low = exact_boot.loc[m, "fundamental/full_ci_low"]
        high = exact_boot.loc[m, "fundamental/full_ci_high"]
        assert low < 0.94 and high > 1.06


def test_the_worst_raw_null_rate_is_not_at_m_two():
    """Section 9.6 reported the worst zero-threshold null rates as sitting at
    ``m = 2``. They do not: the worst row is ``m = 3``, and ``m = 2`` has the
    *lowest* mean of any group order. The driver is weak effect on short
    segments -- proximity to orbit collapse -- at every ``m``, and it eases as
    ``m`` grows.
    """
    results = pd.read_csv(RESULTS / "full_results.csv")
    orbit = results[results["detector"] == "shared_orbit"]

    worst = orbit.loc[orbit["null_rate_zero_threshold"].idxmax()]
    assert worst["m"] == 3
    assert worst["effect"] == pytest.approx(0.08)
    assert worst["total_length"] == 200

    by_m = orbit.groupby("m")["null_rate_zero_threshold"].mean()
    assert by_m.idxmin() == 2
    assert by_m.loc[3] > by_m.loc[6]


def test_crossover_ratio_medians_rest_on_few_effects_at_small_m():
    """Section 9.4's medians are taken over the effects where *both* detectors
    cross inside the grid. At ``m = 2, 3`` that is two of four, so the "median
    across effects" is the midpoint of a pair -- and because each column keeps
    its own surviving subset, the columns do not compose.
    """
    summary = pd.read_csv(RESULTS / "crossover_ratio_summary.csv")
    exact = summary[summary["scenario"] == "exact_orbit"].set_index("m")

    assert exact.loc[2, "shared_orbit/full_n"] == 2
    assert exact.loc[3, "shared_orbit/full_n"] == 2
    assert exact.loc[6, "shared_orbit/full_n"] == 4

    composed = (
        exact.loc[5, "shared_orbit/fundamental"] * exact.loc[5, "fundamental/full"]
    )
    assert abs(composed - exact.loc[5, "shared_orbit/full"]) > 1e-3


def test_shared_orbit_residual_slopes_are_significantly_nonzero():
    """Section 9.3 calls these "small group-dependent drift". By the run's own
    standard errors several are many standard errors from the structural
    prediction of zero -- which is the point, since the identity
    ``slope = dd/2 - s`` localises them in the likelihood gain. They are a
    quantified target for the singular analysis, not noise.
    """
    summary = pd.read_csv(RESULTS / "score_regression_summary.csv")
    orbit = summary[summary["detector"] == "shared_orbit"].set_index("m")

    t_wls = (orbit["penalty_slope_wls"] / orbit["penalty_slope_wls_se"]).abs()
    assert t_wls.loc[6] > 5.0
    assert (t_wls > 2.0).sum() >= 3

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

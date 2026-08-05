"""Model selection among the geometries, instead of assuming one is known.

Every efficiency number in the manuscript is an oracle number: the matching
detector is chosen in advance. These tests cover selecting it from the data.

The central correctness property is that selection is done on *absolute* code
lengths, not on detector scores. Detector scores are measured against different
nulls -- Model A against a pooled unrestricted multinomial, Models B/C/D against
a pooled fundamental coordinate -- so they are not mutually comparable. The code
lengths are, and they reproduce the detector scores exactly as differences.
"""

from __future__ import annotations

import collections

import numpy as np
import pytest

from regimeshift.detectors import (
    full_detector,
    fundamental_detector,
    shared_orbit_detector,
)
from regimeshift.fourier import rotation_matrix
from regimeshift.scenarios import build_segments
from regimeshift.selection import (
    ALTERNATIVES,
    CANDIDATES,
    code_lengths,
    generating_family_for_scenario,
    select_model,
)

GROUPS = [2, 3, 4, 5, 6]


def counts(p, n, rng):
    return rng.multinomial(n, p)


def orbit_distance(segments) -> float:
    """Distance from the nearest exact orbit, as a multiple of the effect."""
    m = segments.m
    best = min(
        np.linalg.norm(segments.theta_right - rotation_matrix(m, r) @ segments.theta_left)
        for r in range(1, m)
    )
    return best / segments.effect


# --------------------------------------------------------------------------
# code lengths are consistent with the detectors
# --------------------------------------------------------------------------


@pytest.mark.parametrize("m", GROUPS)
def test_detector_scores_are_differences_of_code_lengths(m):
    """The bridge between the two views, exact to numerical precision. Each
    detector's score is its null's code length minus its own."""
    segments = build_segments(m, "exact_orbit", 0.25)
    rng = np.random.default_rng(m)
    for _ in range(6):
        cL = counts(segments.p_left, 700, rng)
        cR = counts(segments.p_right, 700, rng)
        L = code_lengths(cL, cR, m)
        assert full_detector(cL, cR, m).score == pytest.approx(
            L["null_full"] - L["full"], abs=1e-9
        )
        assert fundamental_detector(cL, cR, m).score == pytest.approx(
            L["null_fundamental"] - L["fundamental"], abs=1e-9
        )
        assert shared_orbit_detector(cL, cR, m).score == pytest.approx(
            L["null_fundamental"] - L["shared_orbit"], abs=1e-9
        )


@pytest.mark.parametrize("m", GROUPS)
def test_code_lengths_are_finite_and_complete(m):
    rng = np.random.default_rng(m + 1)
    cL = counts(np.full(m, 1 / m), 400, rng)
    cR = counts(np.full(m, 1 / m), 400, rng)
    L = code_lengths(cL, cR, m)
    assert set(L) == set(CANDIDATES)
    assert all(np.isfinite(v) for v in L.values())


def test_candidate_subsets_are_respected():
    rng = np.random.default_rng(0)
    cL = counts(np.full(4, 0.25), 300, rng)
    cR = counts(np.full(4, 0.25), 300, rng)
    subset = ("null_fundamental", "shared_orbit")
    L = code_lengths(cL, cR, 4, candidates=subset)
    assert set(L) == set(subset)
    assert select_model(cL, cR, 4, candidates=subset).selected in subset


def test_empty_segments_are_rejected():
    with pytest.raises(ValueError):
        code_lengths(np.zeros(4), np.array([1.0, 1, 1, 1]), 4)


# --------------------------------------------------------------------------
# what the selector picks
# --------------------------------------------------------------------------


@pytest.mark.parametrize("m", [4, 5, 6])
def test_exact_orbit_data_selects_the_shared_orbit_model(m):
    segments = build_segments(m, "exact_orbit", 0.25)
    rng = np.random.default_rng(m * 3)
    picks = collections.Counter(
        select_model(counts(segments.p_left, 1600, rng), counts(segments.p_right, 1600, rng), m).selected
        for _ in range(40)
    )
    assert picks["shared_orbit"] >= 30, picks
    assert generating_family_for_scenario("exact_orbit") == "shared_orbit"


@pytest.mark.parametrize("m", [4, 6])
def test_higher_mode_data_selects_the_full_model_once_it_is_visible(m):
    """The change lies outside the fundamental subspace, so only Model A can
    describe it -- but it takes data to see that."""
    segments = build_segments(m, "higher_mode", 0.25)
    rng = np.random.default_rng(m * 5)
    picks = collections.Counter(
        select_model(counts(segments.p_left, 3200, rng), counts(segments.p_right, 3200, rng), m).selected
        for _ in range(25)
    )
    assert picks["full"] >= 20, picks


@pytest.mark.parametrize("m", GROUPS)
def test_no_change_selects_a_null(m):
    """The selector answers "did anything change" as well as "what kind"."""
    segments = build_segments(m, "exact_orbit", 0.25)
    rng = np.random.default_rng(m + 11)
    declared = [
        select_model(counts(segments.p_left, 1200, rng), counts(segments.p_left, 1200, rng), m).declared_change
        for _ in range(40)
    ]
    assert sum(declared) <= 4, f"{sum(declared)}/40 false changes"


@pytest.mark.parametrize("m", [4, 6])
def test_selection_sharpens_with_sample_size(m):
    segments = build_segments(m, "exact_orbit", 0.25)
    rng = np.random.default_rng(m)
    rates = []
    for n in (150, 1600):
        hits = sum(
            select_model(counts(segments.p_left, n, rng), counts(segments.p_right, n, rng), m).selected
            == "shared_orbit"
            for _ in range(40)
        )
        rates.append(hits / 40)
    assert rates[1] > rates[0]


# --------------------------------------------------------------------------
# the finding: "generated from B" is not "best described by B"
# --------------------------------------------------------------------------


@pytest.mark.parametrize("m,expected_distance", [(3, 1.12), (4, 0.76), (5, 0.53), (6, 0.40)])
def test_independent_fundamental_drifts_toward_an_orbit_as_m_grows(m, expected_distance):
    """A property of the manuscript's scenario, not of any detector.

    Section 8.2 fixes the angular offset at 0.713 rad while the one-step
    rotation is ``2 pi / m``, which *shrinks* as ``m`` grows. The scenario
    therefore slides toward a one-step orbit: at ``m = 3`` the two differ by
    1.4 rad, at ``m = 6`` by only 0.33 rad. Its distance from the nearest exact
    orbit falls from 1.12 effects to 0.40.
    """
    segments = build_segments(m, "independent_fundamental", 0.25)
    assert orbit_distance(segments) == pytest.approx(expected_distance, abs=0.01)


def test_orbit_distance_is_monotone_in_group_order():
    distances = [
        orbit_distance(build_segments(m, "independent_fundamental", 0.25)) for m in (3, 4, 5, 6)
    ]
    assert distances == sorted(distances, reverse=True)


@pytest.mark.parametrize("m,expected", [(3, "fundamental"), (6, "approximate_orbit")])
def test_independent_fundamental_is_best_described_by_different_models_at_different_m(m, expected):
    """The consequence, and the reason ``generating_family_for_scenario`` is
    bookkeeping rather than ground truth: at ``m = 3`` the scenario really is a
    Model B change, but by ``m = 6`` it sits 0.40 effects from an orbit and an
    approximate-orbit code describes it better."""
    segments = build_segments(m, "independent_fundamental", 0.25)
    rng = np.random.default_rng(m * 13)
    picks = collections.Counter(
        select_model(counts(segments.p_left, 3200, rng), counts(segments.p_right, 3200, rng), m).selected
        for _ in range(30)
    )
    assert picks.most_common(1)[0][0] == expected, picks


# --------------------------------------------------------------------------
# the Selection record
# --------------------------------------------------------------------------


def test_selection_reports_a_margin_and_a_shift():
    segments = build_segments(6, "exact_orbit", 0.25)
    rng = np.random.default_rng(2)
    result = select_model(
        counts(segments.p_left, 2000, rng), counts(segments.p_right, 2000, rng), 6
    )
    assert result.selected == "shared_orbit"
    assert result.declared_change is True
    assert result.selected_shift == 1
    assert result.margin > 0
    assert result.code_lengths[result.selected] == min(result.code_lengths.values())


def test_margin_is_small_when_the_geometries_are_hard_to_tell_apart():
    """Short segments should not produce confident geometry claims. The margin
    is the honest summary the oracle comparison cannot report."""
    segments = build_segments(6, "exact_orbit", 0.08)
    rng = np.random.default_rng(4)
    short = np.mean([
        select_model(counts(segments.p_left, 120, rng), counts(segments.p_right, 120, rng), 6).margin
        for _ in range(25)
    ])
    long = np.mean([
        select_model(counts(segments.p_left, 3200, rng), counts(segments.p_right, 3200, rng), 6).margin
        for _ in range(25)
    ])
    assert short < long


def test_declared_change_matches_the_alternative_set():
    rng = np.random.default_rng(9)
    segments = build_segments(5, "exact_orbit", 0.25)
    for _ in range(10):
        result = select_model(
            counts(segments.p_left, 900, rng), counts(segments.p_right, 900, rng), 5
        )
        assert result.declared_change == (result.selected in ALTERNATIVES)


def test_unknown_scenario_has_no_recorded_family():
    with pytest.raises(ValueError, match="generating family"):
        generating_family_for_scenario("not_a_scenario")
    assert generating_family_for_scenario("approximate_orbit", deviation=0.0) == "shared_orbit"
    assert generating_family_for_scenario("approximate_orbit", deviation=0.5) == "approximate_orbit"

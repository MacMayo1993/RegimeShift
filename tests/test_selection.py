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
    subspace change, but by ``m = 6`` it sits 0.40 effects from an orbit and an
    approximate-orbit code describes it better.

    At ``m = 3`` the expected answer arrives as a *tie* with ``full``, because
    there the fundamental component is the whole tangent space and the two are
    the same model."""
    segments = build_segments(m, "independent_fundamental", 0.25)
    rng = np.random.default_rng(m * 13)
    hits = 0
    for _ in range(30):
        result = select_model(
            counts(segments.p_left, 3200, rng), counts(segments.p_right, 3200, rng), m
        )
        hits += expected in result.tied
    assert hits >= 24, f"m={m}: {expected} chosen or tied in {hits}/30"


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


# --------------------------------------------------------------------------
# the fix: hold the orbit distance constant instead of the angle
# --------------------------------------------------------------------------


@pytest.mark.parametrize("m", [2, 3, 4, 5, 6, 7, 8])
def test_fixed_distance_scenario_holds_its_orbit_distance(m):
    """The corrected variant sits exactly ``INDEPENDENT_ORBIT_DISTANCE`` effects
    from the nearest orbit at every group order, where the manuscript's version
    slides from 1.12 down to 0.40."""
    from regimeshift.scenarios import INDEPENDENT_ORBIT_DISTANCE

    segments = build_segments(m, "independent_fundamental_fixed_distance", 0.25)
    assert orbit_distance(segments) == pytest.approx(INDEPENDENT_ORBIT_DISTANCE, abs=1e-9)
    assert segments.planted_shift is None
    for s in range(m):
        assert not np.allclose(segments.p_right, np.roll(segments.p_left, s), atol=1e-6)


def test_fixed_distance_is_flat_where_the_manuscript_version_drifts():
    drifting = [
        orbit_distance(build_segments(m, "independent_fundamental", 0.25)) for m in (3, 4, 5, 6)
    ]
    fixed = [
        orbit_distance(build_segments(m, "independent_fundamental_fixed_distance", 0.25))
        for m in (3, 4, 5, 6)
    ]
    assert max(drifting) - min(drifting) > 0.6
    assert max(fixed) - min(fixed) < 1e-9


@pytest.mark.parametrize("m", [3, 4, 5, 6])
def test_holding_the_distance_requires_leaving_the_left_radius(m):
    """Not a design choice but a geometric fact, and the reason the manuscript's
    fixed radius ratio could not have held the distance either: adjacent orbit
    points are ``2 sin(pi/m)`` apart, so on a circle of the same radius no point
    is further than ``sin(pi/m)`` from all of them -- 0.5 effects at ``m = 6``."""
    from regimeshift.scenarios import INDEPENDENT_ORBIT_DISTANCE

    best_on_circle = np.sin(np.pi / m)
    assert INDEPENDENT_ORBIT_DISTANCE > best_on_circle
    segments = build_segments(m, "independent_fundamental_fixed_distance", 0.25)
    assert np.linalg.norm(segments.theta_right) > np.linalg.norm(segments.theta_left)


@pytest.mark.parametrize("m", [3, 4, 5, 6])
def test_the_fix_removes_the_selection_drift(m):
    """The point of the whole exercise. On the corrected scenario the selector
    recovers ``fundamental`` at every group order, where on the manuscript's it
    switches to ``approximate_orbit`` by ``m = 5``."""
    segments = build_segments(m, "independent_fundamental_fixed_distance", 0.25)
    rng = np.random.default_rng(m * 17)
    hits = 0
    for _ in range(25):
        result = select_model(
            counts(segments.p_left, 2400, rng), counts(segments.p_right, 2400, rng), m
        )
        # At m = 3 the fundamental component spans the whole tangent space, so
        # "fundamental" and "full" are the same model and come back tied.
        hits += "fundamental" in result.tied
    assert hits >= 22, f"m={m}: {hits}/25"


def test_an_unreachable_orbit_distance_is_rejected():
    """Below ``sin(pi/m)`` the midpoint construction has no real radius."""
    import regimeshift.scenarios as scenarios

    original = scenarios.INDEPENDENT_ORBIT_DISTANCE
    try:
        scenarios.INDEPENDENT_ORBIT_DISTANCE = 0.1
        with pytest.raises(ValueError, match="unreachable"):
            scenarios.build_segments(3, "independent_fundamental_fixed_distance", 0.25)
    finally:
        scenarios.INDEPENDENT_ORBIT_DISTANCE = original


@pytest.mark.parametrize("m", [2, 3])
def test_full_and_fundamental_tie_where_they_are_the_same_model(m):
    """At ``m = 2`` and ``m = 3`` the fundamental component spans the whole
    nontrivial tangent space, so the two candidates describe one hypothesis.
    Their code lengths must coincide, and the selector must report the tie
    rather than picking between them on floating-point noise."""
    segments = build_segments(m, "exact_orbit", 0.25)
    rng = np.random.default_rng(m * 23)
    for _ in range(8):
        cL = counts(segments.p_left, 900, rng)
        cR = counts(segments.p_right, 900, rng)
        L = code_lengths(cL, cR, m)
        assert L["full"] == pytest.approx(L["fundamental"], abs=1e-8)
        result = select_model(cL, cR, m)
        if "fundamental" in result.tied or "full" in result.tied:
            assert {"full", "fundamental"} <= set(result.tied)


@pytest.mark.parametrize("m", [4, 5, 6])
def test_full_and_fundamental_are_distinguishable_above_m_three(m):
    segments = build_segments(m, "exact_orbit", 0.25)
    rng = np.random.default_rng(m * 29)
    cL = counts(segments.p_left, 1200, rng)
    cR = counts(segments.p_right, 1200, rng)
    L = code_lengths(cL, cR, m)
    assert abs(L["full"] - L["fundamental"]) > 1.0


def test_ties_break_toward_the_less_structured_candidate():
    """A tie must never become a claim of structure the data cannot support."""
    from regimeshift.selection import CANDIDATES as order

    segments = build_segments(3, "exact_orbit", 0.25)
    rng = np.random.default_rng(1)
    result = select_model(counts(segments.p_left, 800, rng), counts(segments.p_right, 800, rng), 3)
    if len(result.tied) > 1:
        positions = [order.index(name) for name in result.tied]
        assert positions == sorted(positions)
        assert result.selected == result.tied[0]

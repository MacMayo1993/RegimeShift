"""Tests for the data-generating scenarios and the population gains."""

from __future__ import annotations

import numpy as np
import pytest

from regimeshift.fourier import fundamental_dimension, probabilities
from regimeshift.gains import (
    full_gain,
    fundamental_gain,
    population_gains,
    shared_orbit_gain,
    weighted_jensen_shannon,
)
from regimeshift.scenarios import SCENARIOS, build_segments

GROUPS = [2, 3, 4, 5, 6]
EFFECTS = [0.08, 0.12, 0.18, 0.25]


# --------------------------------------------------------------------------
# scenarios
# --------------------------------------------------------------------------


@pytest.mark.parametrize("m", GROUPS)
@pytest.mark.parametrize("effect", EFFECTS)
def test_exact_orbit_segments_are_group_related(m, effect):
    segments = build_segments(m, "exact_orbit", effect)
    np.testing.assert_allclose(segments.p_right, np.roll(segments.p_left, 1), atol=1e-13)
    assert segments.planted_shift == 1


@pytest.mark.parametrize("m", GROUPS)
@pytest.mark.parametrize("effect", EFFECTS)
def test_left_coordinate_has_the_requested_fisher_norm(m, effect):
    segments = build_segments(m, "exact_orbit", effect)
    assert np.linalg.norm(segments.theta_left) == pytest.approx(effect)
    # The group action is an isometry, so the right state has the same norm.
    assert np.linalg.norm(segments.theta_right) == pytest.approx(effect)


@pytest.mark.parametrize("m", [3, 4, 5, 6])
@pytest.mark.parametrize("effect", EFFECTS)
def test_independent_fundamental_is_not_an_exact_orbit(m, effect):
    segments = build_segments(m, "independent_fundamental", effect)
    assert segments.planted_shift is None
    for s in range(m):
        assert not np.allclose(segments.p_right, np.roll(segments.p_left, s), atol=1e-6)


@pytest.mark.parametrize("effect", EFFECTS)
def test_independent_fundamental_at_m2_is_not_an_exact_orbit(effect):
    segments = build_segments(2, "independent_fundamental", effect)
    for s in range(2):
        assert not np.allclose(segments.p_right, np.roll(segments.p_left, s), atol=1e-6)


@pytest.mark.parametrize("m", [4, 5, 6])
@pytest.mark.parametrize("effect", EFFECTS)
def test_higher_mode_leaves_the_fundamental_family(m, effect):
    """The higher-mode scenario must not be representable inside the fundamental
    family -- otherwise it would not be a misspecification test."""
    from regimeshift.detectors import fit_fundamental

    segments = build_segments(m, "higher_mode", effect)
    for p in (segments.p_left, segments.p_right):
        _, ll = fit_fundamental(p, m)
        entropy_ll = float(np.sum(p * np.log(p)))
        # A perfect fit would reach the entropy bound; a genuine departure does not.
        assert ll < entropy_ll - 1e-6


def test_higher_mode_requires_m_at_least_four():
    for m in (2, 3):
        with pytest.raises(ValueError, match="higher_mode"):
            build_segments(m, "higher_mode", 0.2)


@pytest.mark.parametrize("m", GROUPS)
@pytest.mark.parametrize("scenario", SCENARIOS)
def test_segments_are_valid_distributions(m, scenario):
    if scenario == "higher_mode" and m < 4:
        pytest.skip("scenario undefined for this group order")
    segments = build_segments(m, scenario, 0.18)
    for p in (segments.p_left, segments.p_right, segments.p_null):
        assert p.shape == (m,)
        assert np.all(p > 0)
        np.testing.assert_allclose(p.sum(), 1.0, atol=1e-14)


def test_unknown_scenario_and_bad_effect_are_rejected():
    with pytest.raises(ValueError):
        build_segments(4, "not_a_scenario", 0.1)
    with pytest.raises(ValueError):
        build_segments(4, "exact_orbit", 0.0)


# --------------------------------------------------------------------------
# population gains
# --------------------------------------------------------------------------


@pytest.mark.parametrize("m", GROUPS)
def test_gains_vanish_when_there_is_no_change(m):
    p = probabilities(np.full(fundamental_dimension(m), 0.2), m)
    assert full_gain(p, p) == pytest.approx(0.0, abs=1e-12)
    assert fundamental_gain(p, p, m) == pytest.approx(0.0, abs=1e-9)
    gain, _ = shared_orbit_gain(p, p, m)
    assert gain <= 1e-9  # nonidentity alignment cannot help identical segments


@pytest.mark.parametrize("m", GROUPS)
@pytest.mark.parametrize("effect", EFFECTS)
def test_gains_are_positive_on_exact_orbit_data(m, effect):
    segments = build_segments(m, "exact_orbit", effect)
    gains = population_gains(segments.p_left, segments.p_right, m)
    for name, value in gains.items():
        assert value > 0, name


@pytest.mark.parametrize("m", GROUPS)
def test_shared_orbit_gain_identifies_the_planted_shift(m):
    segments = build_segments(m, "exact_orbit", 0.25)
    _, shift = shared_orbit_gain(segments.p_left, segments.p_right, m)
    assert shift == segments.planted_shift


@pytest.mark.parametrize("m", GROUPS)
@pytest.mark.parametrize("effect", EFFECTS)
def test_shared_orbit_gain_equals_the_fundamental_gain_on_exact_orbits(m, effect):
    """On exact-orbit populations the constrained alternative of Model C loses
    nothing relative to Model B: both recover the truth exactly."""
    segments = build_segments(m, "exact_orbit", effect)
    shared, _ = shared_orbit_gain(segments.p_left, segments.p_right, m)
    assert shared == pytest.approx(fundamental_gain(segments.p_left, segments.p_right, m), abs=1e-8)


@pytest.mark.parametrize("m", [4, 5, 6])
def test_shared_orbit_gain_falls_short_on_independent_changes(m):
    """When the segments are not group-related, sharing one state costs signal."""
    segments = build_segments(m, "independent_fundamental", 0.25)
    shared, _ = shared_orbit_gain(segments.p_left, segments.p_right, m)
    assert shared < fundamental_gain(segments.p_left, segments.p_right, m) - 1e-6


@pytest.mark.parametrize("m", [4, 5, 6])
def test_constrained_gains_lose_the_out_of_subspace_signal(m):
    """Section 9.5: under higher-mode misspecification the full gain retains
    signal the constrained families cannot represent."""
    segments = build_segments(m, "higher_mode", 0.25)
    gains = population_gains(segments.p_left, segments.p_right, m)
    assert gains["full"] > gains["fundamental"] > 0
    # Strict at every group order now that the change is confined to the higher
    # mode. Under the earlier reading these coincided at m = 4, because the
    # mode-2 sign flip there is exactly a one-step shift of the fundamental part.
    assert gains["fundamental"] > gains["shared_orbit"] + 1e-9


@pytest.mark.parametrize("m", [2, 3])
def test_full_and_fundamental_gains_coincide_for_small_groups(m):
    rng = np.random.default_rng(1)
    for _ in range(10):
        p_left = rng.dirichlet(np.ones(m))
        p_right = rng.dirichlet(np.ones(m))
        assert full_gain(p_left, p_right) == pytest.approx(
            fundamental_gain(p_left, p_right, m), abs=1e-8
        )


def test_weighted_jsd_matches_the_symmetric_definition():
    p = np.array([0.5, 0.3, 0.2])
    q = np.array([0.2, 0.3, 0.5])
    mix = 0.5 * (p + q)
    expected = 0.5 * np.sum(p * np.log(p / mix)) + 0.5 * np.sum(q * np.log(q / mix))
    assert weighted_jensen_shannon(p, q) == pytest.approx(expected)
    assert weighted_jensen_shannon(p, q) == pytest.approx(weighted_jensen_shannon(q, p))


def test_weighted_jsd_respects_asymmetric_weights():
    p = np.array([0.7, 0.3])
    q = np.array([0.3, 0.7])
    assert weighted_jensen_shannon(p, q, w_left=0.5) > weighted_jensen_shannon(p, q, w_left=0.9)


def test_higher_mode_is_not_an_exact_orbit_at_any_shift():
    """The corrected reading of Section 8.2: the mode-2 sign flip *is* the
    change, so both segments share a fundamental coordinate and differ only
    outside the fundamental subspace.

    An earlier reading also rotated the fundamental coordinate. That left the
    change carrying a full-strength exact-orbit component, so Model C kept about
    half the population signal and stayed competitive -- contradicting Table 7.
    Under this reading the scenario is a genuine misspecification for Model C at
    every group order, and at no shift is it an orbit relation.
    """
    for m in (4, 5, 6):
        segments = build_segments(m, "higher_mode", 0.25)
        for s in range(m):
            assert not np.allclose(segments.p_right, np.roll(segments.p_left, s), atol=1e-6)
        # Both sides share the same fundamental coordinate.
        np.testing.assert_allclose(segments.theta_left, segments.theta_right, atol=1e-14)


@pytest.mark.parametrize("m", [4, 5, 6])
def test_higher_mode_reproduces_the_manuscript_misspecification_pattern(m):
    """Table 7 reports full power 0.962-0.999 against fundamental 0.182-0.399
    and shared-orbit 0.044-0.082 -- the last of which is *below* the 5% nominal
    rate. The population gains must show the same structure: the fundamental
    family keeps only a small share of the full gain, and the shared-orbit gain
    is negative, meaning aligned pooling is worse than not aligning at all.
    """
    segments = build_segments(m, "higher_mode", 0.25)
    gains = population_gains(segments.p_left, segments.p_right, m)
    assert gains["full"] > 0
    assert 0 < gains["fundamental"] < 0.1 * gains["full"]
    assert gains["shared_orbit"] < 0


# ---------------------------------------------------------------------------
# geometry of the fixed-distance scenario (Section 9.8)
# ---------------------------------------------------------------------------


def _nearest_nonidentity_orbit_distance(segments):
    """Distance to the nearest orbit point the shared-orbit fit can reach.

    Model C maximises over *nonidentity* shifts only, so that -- not the
    distance to the full orbit including the identity -- is the signal it
    cannot capture, and the quantity both scenarios are parameterised by.
    """
    from regimeshift.fourier import rotation_matrix

    return min(
        float(np.linalg.norm(segments.theta_right - rotation_matrix(segments.m, s) @ segments.theta_left))
        for s in range(1, segments.m)
    )


@pytest.mark.parametrize("m", [2, 3, 4, 5, 6])
def test_fixed_distance_scenario_hits_its_target_at_every_group_order(m):
    from regimeshift.scenarios import INDEPENDENT_ORBIT_DISTANCE

    segments = build_segments(m, "independent_fundamental_fixed_distance", 1.0)
    assert _nearest_nonidentity_orbit_distance(segments) == pytest.approx(
        INDEPENDENT_ORBIT_DISTANCE, abs=1e-9
    )


def test_original_scenario_slides_toward_an_orbit_as_m_grows():
    """The defect Section 9.8 identifies: a fixed angular offset against a
    shrinking one-step rotation. Pins the table it reports."""
    expected = {3: 1.12, 4: 0.76, 5: 0.53, 6: 0.40}
    for m, distance in expected.items():
        segments = build_segments(m, "independent_fundamental", 1.0)
        assert _nearest_nonidentity_orbit_distance(segments) == pytest.approx(distance, abs=0.005)


@pytest.mark.parametrize("m", [3, 4, 5, 6])
def test_a_same_radius_point_cannot_reach_the_target_distance(m):
    """Why the fixed-distance variant must leave the left coordinate's circle.

    The bound is ``2 sin(pi / 2m)`` -- the chord to the angular midpoint --
    *not* half the adjacent-vertex chord ``2 sin(pi/m)``, which is the distance
    to a point that does not lie on the circle at all. The two agree to 3% at
    ``m = 6``, which is why the looser reading is easy to make; they differ by
    13% at ``m = 3``.
    """
    from regimeshift.scenarios import INDEPENDENT_ORBIT_DISTANCE

    angles = np.linspace(0.0, 2.0 * np.pi, 20001)
    points = np.stack([np.cos(angles), np.sin(angles)], axis=1)
    orbit = np.stack(
        [[np.cos(2 * np.pi * k / m), np.sin(2 * np.pi * k / m)] for k in range(m)]
    )
    reachable = np.min(np.linalg.norm(points[:, None, :] - orbit[None, :, :], axis=2), axis=1).max()

    assert reachable == pytest.approx(2 * np.sin(np.pi / (2 * m)), abs=1e-3)
    assert reachable > np.sin(np.pi / m)          # the bound the text used to give
    assert reachable < INDEPENDENT_ORBIT_DISTANCE  # so a unit radius cannot reach it


@pytest.mark.parametrize("m", [3, 4, 5, 6])
def test_the_midpoint_ray_floor_is_the_constraint_the_code_enforces(m):
    """``sin(pi/m)`` is a real quantity here, just a different one: the minimum
    distance along the angular-midpoint ray, minimised over radius. It is what
    the scenario's discriminant requires the target to clear."""
    from regimeshift.scenarios import build_segments as _build

    radii = np.linspace(0.01, 6.0, 60001)
    midpoint = np.array([np.cos(np.pi / m), np.sin(np.pi / m)])
    distances = np.linalg.norm(radii[:, None] * midpoint[None, :] - np.array([1.0, 0.0]), axis=1)
    assert distances.min() == pytest.approx(np.sin(np.pi / m), abs=1e-4)

    with pytest.raises(ValueError, match="unreachable"):
        import regimeshift.scenarios as sc

        original = sc.INDEPENDENT_ORBIT_DISTANCE
        try:
            sc.INDEPENDENT_ORBIT_DISTANCE = np.sin(np.pi / m) * 0.9
            _build(m, "independent_fundamental_fixed_distance", 1.0)
        finally:
            sc.INDEPENDENT_ORBIT_DISTANCE = original


def test_holding_orbit_distance_fixed_does_not_hold_signal_strength_fixed():
    """The corrected variant removes one confound and introduces a larger one.

    Distance to the nearest orbit is constant by construction, but the size of
    the change is not: the full population gain spans a factor of about 30
    across group orders, against roughly 5 for the scenario it replaces. Any
    cross-``m`` reading of either scenario is partly reading signal strength.
    """
    def full_gain_of(scenario, m):
        segments = build_segments(m, scenario, 0.25)
        return population_gains(segments.p_left, segments.p_right, m)["full"]

    fixed = {m: full_gain_of("independent_fundamental_fixed_distance", m) for m in (2, 3, 4, 5, 6)}
    original = {m: full_gain_of("independent_fundamental", m) for m in (2, 3, 4, 5, 6)}

    assert max(fixed.values()) / min(fixed.values()) > 20.0
    assert max(original.values()) / min(original.values()) < 10.0


def test_higher_mode_residual_gain_is_not_uniformly_three_percent():
    """Section 9.5 quotes "about 3% of the full gain" as though it held across
    group orders. It is about 3% at ``m = 5, 6`` and nearly twice that at
    ``m = 4``, where mode 2 is the one-dimensional sign representation."""
    retained = {}
    for m in (4, 5, 6):
        segments = build_segments(m, "higher_mode", 0.25)
        gains = population_gains(segments.p_left, segments.p_right, m)
        retained[m] = gains["fundamental"] / gains["full"]
        assert gains["shared_orbit"] < 0.0  # aligned pooling is worse than not aligning

    assert retained[4] > 0.05
    assert retained[5] == pytest.approx(0.031, abs=0.004)
    assert retained[6] == pytest.approx(0.030, abs=0.004)

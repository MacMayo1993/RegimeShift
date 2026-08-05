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
    assert gains["fundamental"] >= gains["shared_orbit"] - 1e-12
    if m > 4:
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


def test_higher_mode_at_m4_is_still_an_exact_orbit_in_the_full_space():
    """A property worth pinning down: at m = 4 the mode-2 component is the sign
    representation, which flips under a one-step shift. With the antisymmetric
    placement used by the scenario, the two segments therefore remain exact
    group transforms of one another -- but in the *full* simplex, outside the
    fundamental family. The constrained detectors are still misspecified,
    because neither family can represent the mode-2 component at all."""
    segments = build_segments(4, "higher_mode", 0.25)
    np.testing.assert_allclose(segments.p_right, np.roll(segments.p_left, 1), atol=1e-13)

    # At m = 5 and m = 6 the higher mode does not flip sign, so the orbit
    # relation is genuinely broken as well.
    for m in (5, 6):
        segments = build_segments(m, "higher_mode", 0.25)
        for s in range(m):
            assert not np.allclose(segments.p_right, np.roll(segments.p_left, s), atol=1e-6)

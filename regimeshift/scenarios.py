"""Data-generating scenarios (Section 8.2).

Three scenarios exercise the three levels of geometric constraint:

``exact_orbit``
    The right state is a one-step cyclic transform of the left. Matches Model C
    and is contained in Models A and B.

``independent_fundamental``
    Both regimes lie in the fundamental family but are not related by a cyclic
    shift. Matches Model B, generally violates Model C.

``higher_mode``
    A mode-2 Fourier component is added with opposite signs on the two sides,
    placing signal outside the fundamental component. A misspecification test
    for Models B and C.

The scenario constants below are the reproducible defaults of this
implementation; see ``docs/paper-notes.md`` for which manuscript values were
recoverable from the source document.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .fourier import (
    fundamental_dimension,
    higher_mode_logits,
    probabilities,
    rotation_matrix,
)

__all__ = [
    "SCENARIOS",
    "INDEPENDENT_RADIUS_FACTOR",
    "INDEPENDENT_ANGLE_RAD",
    "INDEPENDENT_M2_FACTOR",
    "HIGHER_MODE_FACTOR",
    "Segments",
    "build_segments",
]

SCENARIOS = ("exact_orbit", "independent_fundamental", "higher_mode")

#: Radius ratio of the right coordinate in the independent-fundamental scenario.
INDEPENDENT_RADIUS_FACTOR = 0.85
#: Angular offset (radians) of the right coordinate, m >= 3.
INDEPENDENT_ANGLE_RAD = 0.713
#: Right/left coordinate ratio in the independent-fundamental scenario at m = 2.
INDEPENDENT_M2_FACTOR = -0.6
#: Amplitude of the higher mode, as a multiple of the effect size.
HIGHER_MODE_FACTOR = 0.8


@dataclass(frozen=True)
class Segments:
    """Population description of one two-segment configuration."""

    m: int
    scenario: str
    effect: float
    p_left: np.ndarray
    p_right: np.ndarray
    theta_left: np.ndarray
    theta_right: np.ndarray
    planted_shift: int | None
    """The true relative group element, or ``None`` when no exact orbit relation
    holds."""

    @property
    def p_null(self) -> np.ndarray:
        """No-change distribution used for null calibration of this scenario."""
        return self.p_left


def _base_theta(m: int, effect: float) -> np.ndarray:
    """Left coordinate with Fisher norm equal to ``effect``."""
    d = fundamental_dimension(m)
    theta = np.zeros(d)
    theta[0] = effect
    return theta


def build_segments(m: int, scenario: str, effect: float) -> Segments:
    """Construct the population distributions for one configuration."""
    if scenario not in SCENARIOS:
        raise ValueError(f"unknown scenario {scenario!r}; expected one of {SCENARIOS}")
    if effect <= 0:
        raise ValueError("effect must be positive")

    theta_left = _base_theta(m, effect)

    if scenario == "exact_orbit":
        theta_right = rotation_matrix(m, 1) @ theta_left
        p_left = probabilities(theta_left, m)
        p_right = probabilities(theta_right, m)
        return Segments(m, scenario, effect, p_left, p_right, theta_left, theta_right, 1)

    if scenario == "independent_fundamental":
        if m == 2:
            theta_right = INDEPENDENT_M2_FACTOR * theta_left
        else:
            c, s = np.cos(INDEPENDENT_ANGLE_RAD), np.sin(INDEPENDENT_ANGLE_RAD)
            rot = np.array([[c, -s], [s, c]])
            theta_right = INDEPENDENT_RADIUS_FACTOR * (rot @ theta_left)
        p_left = probabilities(theta_left, m)
        p_right = probabilities(theta_right, m)
        return Segments(m, scenario, effect, p_left, p_right, theta_left, theta_right, None)

    # higher_mode: exact-orbit base plus an antisymmetric mode-2 component.
    if m < 4:
        raise ValueError(
            "the higher_mode scenario needs m >= 4: at m = 2 mode 2 is trivial and at "
            "m = 3 mode 2 is the conjugate of mode 1, so it lies in the fundamental "
            "component and is not a misspecification"
        )
    theta_right = rotation_matrix(m, 1) @ theta_left
    extra = higher_mode_logits(m, HIGHER_MODE_FACTOR * effect, mode=2)
    p_left = probabilities(theta_left, m, extra_logits=extra)
    p_right = probabilities(theta_right, m, extra_logits=-extra)
    return Segments(m, scenario, effect, p_left, p_right, theta_left, theta_right, None)

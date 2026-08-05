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

``approximate_orbit``
    Section 14.1: ``eta_R = R eta_L + delta``, a one-step orbit plus a
    controllable deviation *inside* the fundamental subspace. ``deviation = 0``
    is the exact orbit; growing it sweeps continuously toward an independent
    subspace change, which is how far a relational detector can be pushed
    before its advantage disappears.

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
    "MANUSCRIPT_CONSTANTS",
    "Segments",
    "build_segments",
]

SCENARIOS = ("exact_orbit", "independent_fundamental", "higher_mode", "approximate_orbit")

#: Radius ratio of the right coordinate in the independent-fundamental scenario.
INDEPENDENT_RADIUS_FACTOR = 0.72
#: Angular offset (radians) of the right coordinate, m >= 3.
INDEPENDENT_ANGLE_RAD = 0.713
#: Right/left coordinate ratio in the independent-fundamental scenario at m = 2.
INDEPENDENT_M2_FACTOR = -0.55
#: Amplitude of the higher mode, as a multiple of the effect size.
HIGHER_MODE_FACTOR = 0.85

#: Machine-readable provenance for each scenario constant taken from the source
#: manuscript: value, manuscript section, and how it was obtained.
#:
#: All of these are now quoted directly from the document. An earlier version of
#: this file guessed three of them, on the mistaken belief that the manuscript
#: rendered its equations as images. It does not -- the equations are Office Math
#: (OMML), and the first extraction pass simply dropped them by reading only
#: ``w:t`` elements. See ``docs/paper-notes.md``.
#:
#: A test asserts each entry matches the live module constant, so the table
#: cannot drift from the code.
MANUSCRIPT_CONSTANTS = {
    "INDEPENDENT_RADIUS_FACTOR": {
        "value": INDEPENDENT_RADIUS_FACTOR,
        "section": "8.2",
        "recovered_from_manuscript": True,
        "basis": (
            "Section 8.2: 'the right coordinate had radius 0.72 times the left radius "
            "and angle 0.713 radians'."
        ),
    },
    "INDEPENDENT_ANGLE_RAD": {
        "value": INDEPENDENT_ANGLE_RAD,
        "section": "8.2",
        "recovered_from_manuscript": True,
        "basis": "Stated in readable body text.",
    },
    "INDEPENDENT_M2_FACTOR": {
        "value": INDEPENDENT_M2_FACTOR,
        "section": "8.2",
        "recovered_from_manuscript": True,
        "basis": (
            "Section 8.2: 'For g = 2, the right coordinate was -0.55 times the left "
            "coordinate.' The sign matters: -1 would BE the exact orbit."
        ),
    },
    "HIGHER_MODE_FACTOR": {
        "value": HIGHER_MODE_FACTOR,
        "section": "8.2",
        "recovered_from_manuscript": True,
        "basis": (
            "Section 8.2: 'a mode-2 Fourier component with amplitude 0.85 times the "
            "effect was added with opposite signs on the two sides of the boundary'."
        ),
    },
    "LABEL_COST": {
        "value": "log(m - 1)",
        "section": "3.3, 7.3",
        "recovered_from_manuscript": True,
        "basis": (
            "Section 3.3: 'Under a uniform two-part label code, its cost is log(g - 1) "
            "nats.' Section 7.3 confirms log 1 = 0 at g = 2. Still one legitimate "
            "two-part code rather than a uniquely determined MDL constant."
        ),
    },
}


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
    deviation: float = 0.0
    """Fisher-norm size of the departure from an exact orbit, as a multiple of
    ``effect``. Nonzero only for the ``approximate_orbit`` scenario."""

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


def build_segments(m: int, scenario: str, effect: float, deviation: float = 0.0) -> Segments:
    """Construct the population distributions for one configuration.

    ``deviation`` applies only to ``approximate_orbit``, where it sets the size
    of the departure from an exact orbit as a multiple of ``effect``.
    """
    if scenario not in SCENARIOS:
        raise ValueError(f"unknown scenario {scenario!r}; expected one of {SCENARIOS}")
    if effect <= 0:
        raise ValueError("effect must be positive")
    if deviation < 0:
        raise ValueError("deviation must be non-negative")
    if deviation and scenario != "approximate_orbit":
        raise ValueError(f"deviation applies only to approximate_orbit, not {scenario!r}")

    theta_left = _base_theta(m, effect)

    if scenario == "approximate_orbit":
        # One-step orbit, displaced perpendicular to the rotated state so the
        # deviation is a pure departure from the orbit rather than a rescaling.
        rotated = rotation_matrix(m, 1) @ theta_left
        if fundamental_dimension(m) == 1:
            direction = np.array([1.0])
        else:
            unit = rotated / np.linalg.norm(rotated)
            direction = np.array([-unit[1], unit[0]])
        theta_right = rotated + deviation * effect * direction
        p_left = probabilities(theta_left, m)
        p_right = probabilities(theta_right, m)
        planted = 1 if deviation == 0 else None
        return Segments(m, scenario, effect, p_left, p_right,
                        theta_left, theta_right, planted, deviation)

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

    # higher_mode: the mode-2 sign flip *is* the change. Both segments share the
    # same fundamental coordinate, which only sets the operating point away from
    # uniform; the two sides differ solely by an antisymmetric mode-2 component.
    #
    # Section 8.2 says a mode-2 component "was added with opposite signs on the
    # two sides of the boundary". An earlier reading of this implementation also
    # rotated the fundamental coordinate, so the change carried a full-strength
    # exact-orbit component *plus* the higher mode -- which left Model C about
    # half the population signal and kept it competitive, contradicting Table 7.
    # With the change confined to the higher mode, the population gains
    # reproduce Table 7's pattern: the fundamental family retains only a few
    # percent of the full gain and the shared-orbit gain goes negative, matching
    # the reported below-nominal shared-orbit power of 0.044-0.082.
    if m < 4:
        raise ValueError(
            "the higher_mode scenario needs m >= 4: at m = 2 mode 2 is trivial and at "
            "m = 3 mode 2 is the conjugate of mode 1, so it lies in the fundamental "
            "component and is not a misspecification"
        )
    theta_right = theta_left
    extra = higher_mode_logits(m, HIGHER_MODE_FACTOR * effect, mode=2)
    p_left = probabilities(theta_left, m, extra_logits=extra)
    p_right = probabilities(theta_right, m, extra_logits=-extra)
    return Segments(m, scenario, effect, p_left, p_right, theta_left, theta_right, None)

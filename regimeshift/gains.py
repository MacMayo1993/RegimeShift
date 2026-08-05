"""Detector-specific population gains (Section 6.1).

The relevant signal strength is the expected log-likelihood advantage *within
the model being fitted*. For the unrestricted multinomial family the pooled
null is the arithmetic mixture and the gain is the weighted Jensen-Shannon
divergence. For the constrained families the pooled null is a KL projection
rather than the mixture, so ordinary JSD must not be reused; each gain is
computed by optimising the corresponding population log-likelihood.

All gains are per observation, in nats. Multiplying by the total length ``n``
gives the leading term of the expected raw gain.
"""

from __future__ import annotations

import numpy as np

from .detectors import fit_fundamental

__all__ = [
    "weighted_jensen_shannon",
    "full_gain",
    "fundamental_gain",
    "shared_orbit_gain",
    "population_gains",
]


def _kl(p: np.ndarray, q: np.ndarray) -> float:
    p = np.asarray(p, dtype=float)
    q = np.asarray(q, dtype=float)
    nz = p > 0
    return float(np.sum(p[nz] * np.log(p[nz] / q[nz])))


def weighted_jensen_shannon(p_left: np.ndarray, p_right: np.ndarray, w_left: float = 0.5) -> float:
    """Weighted Jensen-Shannon divergence, the Model A population gain."""
    w_right = 1.0 - w_left
    p_bar = w_left * np.asarray(p_left, float) + w_right * np.asarray(p_right, float)
    return w_left * _kl(p_left, p_bar) + w_right * _kl(p_right, p_bar)


def full_gain(p_left: np.ndarray, p_right: np.ndarray, w_left: float = 0.5) -> float:
    """Model A population gain per observation."""
    return weighted_jensen_shannon(p_left, p_right, w_left)


def fundamental_gain(p_left: np.ndarray, p_right: np.ndarray, m: int, w_left: float = 0.5) -> float:
    """Model B population gain per observation.

    Alternative: each segment is KL-projected onto the fundamental family
    independently. Null: one coordinate is fitted to the weighted mixture.
    """
    w_right = 1.0 - w_left
    p_left = np.asarray(p_left, float)
    p_right = np.asarray(p_right, float)
    _, ll_left = fit_fundamental(p_left, m)
    _, ll_right = fit_fundamental(p_right, m)
    _, ll_null = fit_fundamental(w_left * p_left + w_right * p_right, m)
    return float(w_left * ll_left + w_right * ll_right - ll_null)


def shared_orbit_gain(p_left: np.ndarray, p_right: np.ndarray, m: int, w_left: float = 0.5):
    """Model C population gain per observation and the maximising shift.

    Returns ``(gain, shift)``. The maximisation runs over nonidentity shifts
    only, matching the detector.
    """
    w_right = 1.0 - w_left
    p_left = np.asarray(p_left, float)
    p_right = np.asarray(p_right, float)
    _, ll_null = fit_fundamental(w_left * p_left + w_right * p_right, m)

    best_ll, best_shift = -np.inf, None
    for s in range(1, m):
        pooled = w_left * p_left + w_right * np.roll(p_right, -s)
        _, ll = fit_fundamental(pooled, m)
        if ll > best_ll:
            best_ll, best_shift = ll, s
    return float(best_ll - ll_null), best_shift


def population_gains(p_left: np.ndarray, p_right: np.ndarray, m: int, w_left: float = 0.5) -> dict[str, float]:
    """All three population gains for one pair of segment distributions."""
    shared, _ = shared_orbit_gain(p_left, p_right, m, w_left)
    return {
        "full": full_gain(p_left, p_right, w_left),
        "fundamental": fundamental_gain(p_left, p_right, m, w_left),
        "shared_orbit": shared,
    }

"""The three known-boundary detectors (Models A, B and C of the manuscript).

Model A -- full independent change: each segment gets an unrestricted
multinomial parameter. Continuous-dimension increment ``m - 1``.

Model B -- independent fundamental-subspace change: each segment gets its own
parameter inside the fundamental invariant subspace. Increment ``d_fund``.

Model C -- shared exact-orbit transition: both segments share one continuous
state and differ only by a relative cyclic shift. Continuous-dimension
increment zero; the alternative pays only a discrete relative-label code.

All three detectors are scored as ``maximised log-likelihood gain minus an
explicit complexity increment``, in nats, at a *known* boundary. No location
cost is applied to any detector (Section 4.3).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize

from .fourier import (
    fourier_design_matrix,
    full_dimension,
    fundamental_dimension,
)

__all__ = [
    "DetectorResult",
    "split_penalty",
    "label_cost",
    "multinomial_loglik",
    "fit_fundamental",
    "fundamental_loglik",
    "full_detector",
    "fundamental_detector",
    "shared_orbit_detector",
    "run_all_detectors",
    "DETECTORS",
]


@dataclass(frozen=True)
class DetectorResult:
    """Outcome of scoring one two-segment dataset with one detector."""

    name: str
    raw_gain: float
    """Maximised log-likelihood advantage of the alternative, in nats."""
    penalty: float
    """Explicit complexity increment, in nats."""
    score: float
    """``raw_gain - penalty``. The raw MDL rule declares a change when > 0."""
    dimension_increment: float
    """Continuous-dimension increment of the alternative."""
    selected_shift: int | None = None
    """Relative group element chosen by the shared-orbit detector."""


def split_penalty(dim: int, n_left: int, n_right: int) -> float:
    """Exact known-split regular complexity increment (Section 4.2).

    ``(dim / 2) * log(n_left * n_right / n)``. The coefficient of ``log n`` is
    ``dim / 2``; the split fraction only affects the bounded term.
    """
    if n_left <= 0 or n_right <= 0:
        raise ValueError("both segments must be non-empty")
    n = n_left + n_right
    return 0.5 * dim * (np.log(n_left) + np.log(n_right) - np.log(n))


def label_cost(m: int) -> float:
    """Two-part code length for a *nonidentity* relative group element, in nats.

    The alternative of Model C ranges over the ``m - 1`` nonidentity elements of
    C_m, so a uniform label code costs ``log(m - 1)``. For ``m = 2`` there is a
    single nonidentity shift and the cost is zero. This is constant in ``n`` --
    the defining property of Model C.
    """
    if m < 2:
        raise ValueError("group order must be >= 2")
    return float(np.log(m - 1))


def multinomial_loglik(counts: np.ndarray) -> float:
    """Maximised multinomial log-likelihood (kernel, excluding the multinomial
    coefficient, which cancels between null and alternative)."""
    counts = np.asarray(counts, dtype=float)
    n = counts.sum()
    if n <= 0:
        return 0.0
    nz = counts > 0
    return float(np.sum(counts[nz] * (np.log(counts[nz]) - np.log(n))))


def _neg_loglik_and_grad(theta: np.ndarray, counts: np.ndarray, B: np.ndarray):
    logits = B @ theta
    logits = logits - logits.max()
    w = np.exp(logits)
    z = w.sum()
    p = w / z
    n = counts.sum()
    ll = float(counts @ np.log(p))
    grad = B.T @ (counts - n * p)
    return -ll, -grad


def fit_fundamental(counts: np.ndarray, m: int, n_restarts: int = 2) -> tuple[np.ndarray, float]:
    """MLE of the fundamental-family coordinate for a count vector.

    Returns ``(theta_hat, max_loglik)``. Optimisation is L-BFGS-B with analytic
    gradients in Cartesian Fourier coordinates; the objective is concave in the
    logits, so restarts only guard against pathological line searches.
    """
    counts = np.asarray(counts, dtype=float)
    if counts.shape != (m,):
        raise ValueError(f"counts must have shape ({m},), got {counts.shape}")
    if counts.sum() <= 0:
        return np.zeros(fundamental_dimension(m)), 0.0

    B = fourier_design_matrix(m)
    d = B.shape[1]
    # Warm start from the first-order (Fisher-orthonormal) projection of the
    # empirical frequencies onto the fundamental component.
    freq = counts / counts.sum()
    start = (B - B.mean(axis=0, keepdims=True)).T @ freq
    starts = [start, np.zeros(d)][:max(1, n_restarts)]

    best_theta, best_ll = None, -np.inf
    for x0 in starts:
        res = minimize(
            _neg_loglik_and_grad,
            x0,
            args=(counts, B),
            jac=True,
            method="L-BFGS-B",
            options={"maxiter": 500, "ftol": 1e-14, "gtol": 1e-10},
        )
        if -res.fun > best_ll:
            best_ll, best_theta = -res.fun, res.x
    return best_theta, float(best_ll)


def fundamental_loglik(theta: np.ndarray, counts: np.ndarray, m: int) -> float:
    """Log-likelihood of ``counts`` under the fundamental family at ``theta``."""
    B = fourier_design_matrix(m)
    return -_neg_loglik_and_grad(np.asarray(theta, dtype=float), np.asarray(counts, float), B)[0]


def full_detector(counts_left: np.ndarray, counts_right: np.ndarray, m: int) -> DetectorResult:
    """Model A: unrestricted multinomial, BIC-scored known-split increment."""
    cL = np.asarray(counts_left, dtype=float)
    cR = np.asarray(counts_right, dtype=float)
    gain = multinomial_loglik(cL) + multinomial_loglik(cR) - multinomial_loglik(cL + cR)
    dim = full_dimension(m)
    pen = split_penalty(dim, int(cL.sum()), int(cR.sum()))
    return DetectorResult("full", gain, pen, gain - pen, dim)


def fundamental_detector(counts_left: np.ndarray, counts_right: np.ndarray, m: int) -> DetectorResult:
    """Model B: independently fitted coordinates inside the fundamental subspace."""
    cL = np.asarray(counts_left, dtype=float)
    cR = np.asarray(counts_right, dtype=float)
    _, ll_null = fit_fundamental(cL + cR, m)
    _, ll_left = fit_fundamental(cL, m)
    _, ll_right = fit_fundamental(cR, m)
    gain = ll_left + ll_right - ll_null
    dim = fundamental_dimension(m)
    pen = split_penalty(dim, int(cL.sum()), int(cR.sum()))
    return DetectorResult("fundamental", gain, pen, gain - pen, dim)


def shared_orbit_detector(counts_left: np.ndarray, counts_right: np.ndarray, m: int) -> DetectorResult:
    """Model C: one shared continuous state plus a relative cyclic shift.

    For each nonidentity shift ``s`` the right counts are aligned by ``g^{-s}``,
    pooled with the left counts, and a single shared coordinate is fitted. The
    alternative takes the shift with the largest shared-state likelihood. The
    penalty is the constant label cost -- no location cost and no
    continuous-dimension increment.
    """
    cL = np.asarray(counts_left, dtype=float)
    cR = np.asarray(counts_right, dtype=float)
    _, ll_null = fit_fundamental(cL + cR, m)

    best_ll, best_shift = -np.inf, None
    for s in range(1, m):
        pooled = cL + np.roll(cR, -s)
        _, ll = fit_fundamental(pooled, m)
        if ll > best_ll:
            best_ll, best_shift = ll, s

    gain = best_ll - ll_null
    pen = label_cost(m)
    return DetectorResult("shared_orbit", gain, pen, gain - pen, 0.0, best_shift)


DETECTORS = {
    "full": full_detector,
    "fundamental": fundamental_detector,
    "shared_orbit": shared_orbit_detector,
}


def run_all_detectors(counts_left: np.ndarray, counts_right: np.ndarray, m: int) -> dict[str, DetectorResult]:
    """Score one dataset with all three detectors."""
    return {name: fn(counts_left, counts_right, m) for name, fn in DETECTORS.items()}

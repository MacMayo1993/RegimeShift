"""The known-boundary detectors (Models A, B and C of the manuscript, plus D).

Model A -- full independent change: each segment gets an unrestricted
multinomial parameter. Continuous-dimension increment ``m - 1``.

Model B -- independent fundamental-subspace change: each segment gets its own
parameter inside the fundamental invariant subspace. Increment ``d_fund``.

Model C -- shared exact-orbit transition: both segments share one continuous
state and differ only by a relative cyclic shift. Continuous-dimension
increment zero; the alternative pays only a discrete relative-label code.

Model D -- approximate orbit (Section 14.1): a shared state plus a relative
shift plus a *shrunk* deviation, interpolating between C and B. Not part of the
manuscript's production comparison; see ``approximate_orbit_detector``.

Every detector is scored as ``maximised log-likelihood gain minus an explicit
complexity increment``, in nats, at a *known* boundary. No location cost is
applied to any detector (Section 4.3).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize

from .fourier import (
    fourier_design_matrix,
    rotation_matrix as _rotation,
    full_dimension,
    fundamental_dimension,
)

__all__ = [
    "DetectorResult",
    "GRADIENT_TOLERANCE",
    "K_STAR",
    "nats_to_bits",
    "fit_failure_count",
    "reset_fit_failures",
    "split_penalty",
    "label_cost",
    "multinomial_loglik",
    "fit_fundamental",
    "fundamental_loglik",
    "full_detector",
    "fundamental_detector",
    "shared_orbit_detector",
    "approximate_orbit_detector",
    "fit_approximate_orbit",
    "deviation_penalty",
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


K_STAR = 1.0 / (2.0 * np.log(2.0))
"""The per-dimension leading penalty rate, in bits per e-fold: ``1/(2 ln 2)``.

Scoring throughout this package is in nats, where a regular known-split
increment has leading term ``(d/2) log n``. Converting to bits gives

    (d/2) log2(n)  =  d / (2 ln 2) * ln n  =  d * K_STAR * ln n

so on regular strata every leading coefficient in the framework is an integer
multiple of ``K_STAR``: ``(m - 1)`` of them for Model A, ``d_fund`` for Model B,
and *zero* for Model C. The three-way hierarchy is how many ``K_STAR`` a model
pays to cross the boundary.

The "integer multiple" part is conditional, and the condition is not decorative.
It holds because regular BIC counts a whole number of parameters and charges
half a ``log n`` for each. Under singular learning theory the leading
coefficient is a real log canonical threshold, which need not be a half-integer
at all -- and orbit collapse (``eta = 0``, or any state with a nontrivial
stabiliser) is exactly such a singularity. The counting picture describes the
regular part of the problem.

Note also that this constant is definitional, not a discovery -- it is just the
nats-to-bits conversion of Schwarz's one-half. Anything counting half a
parameter per e-fold in base 2 produces it. Appendix D of the manuscript notes
that the same number appears in an East-model inverse-gap asymptotic and is
careful to call the resemblance suggestive rather than explanatory; the
resemblance carries weight only if the dynamical occurrence is *not* similarly
a units artefact.
"""

GRADIENT_TOLERANCE = 1e-6
"""Relative first-order tolerance defining a converged fundamental-family fit."""

_FIT_FAILURES = 0
"""Count of fundamental-family fits whose best optimiser run did not converge.

The population gains used as regressors are themselves obtained by numerical
optimisation, so a silent convergence failure would corrupt a likelihood with no
visible symptom. Every fit checks ``OptimizeResult.success`` and increments this
counter; :func:`~regimeshift.simulation.run_config` reports the per-configuration
total, and the test suite asserts it stays at zero across the grid.

The counter is process-local, which is what the parallel runner needs: each
worker accumulates the failures of the configurations it ran.
"""


def reset_fit_failures() -> None:
    """Zero the convergence-failure counter for this process."""
    global _FIT_FAILURES
    _FIT_FAILURES = 0


def fit_failure_count() -> int:
    """Number of non-converged fits since the last reset in this process."""
    return _FIT_FAILURES


def nats_to_bits(nats: float) -> float:
    """Convert a codelength or codelength coefficient from nats to bits."""
    return nats / np.log(2.0)


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

    n = counts.sum()
    best_theta, best_ll, best_ok = None, -np.inf, False
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
            # Judge convergence on the first-order condition, not on
            # OptimizeResult.success: under our deliberately tight ftol,
            # L-BFGS-B reports ABNORMAL whenever its line search cannot improve
            # at machine precision, which routinely happens *at* the optimum
            # (observed gradient norms ~5e-8). The gradient is
            # B^T (counts - n p), so it scales with n and the test is relative.
            best_ok = bool(res.success) or np.abs(res.jac).max() <= GRADIENT_TOLERANCE * max(1.0, n)

    if not best_ok:
        global _FIT_FAILURES
        _FIT_FAILURES += 1
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

    **The zero increment is a regular-stratum statement.** It requires the shared
    state to lie on a *regular orbit*, i.e. to have trivial stabiliser in the
    cyclic group, so that the ``m`` points ``{R^s eta}`` are distinct and the
    relative shift is identifiable. At ``eta = 0`` the orbit collapses to a
    point, every shift acts trivially, and ordinary BIC dimension counting is
    not the right marginal-likelihood theory there at all -- the singular
    framework of Watanabe and of Drton and Plummer applies instead. The detector
    still runs at such points; what changes is the interpretation of its
    penalty. See Proposition 1 of the v4 manuscript for the assumptions stated
    in full, and its Section 9.6 for what collapse does to the ``m = 2`` null in
    practice.
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


# --------------------------------------------------------------------------
# Model D: approximate orbit (Section 14.1)
# --------------------------------------------------------------------------


def deviation_penalty(
    dim: int, n_left: int, n_right: int, deviation_scale: float
) -> float:
    """Codelength for a shrunk deviation vector, in nats.

    Section 14.1 interpolates between Models B and C with
    ``eta_R = R^r eta_L + delta`` under "a shrinkage prior or code on delta".
    Take that code to be a Gaussian prior ``delta ~ N(0, tau^2 I_d)`` and apply
    the Laplace approximation. The shared state ``eta`` is *not* known: it is
    estimated jointly with ``delta``, and the two are correlated, because moving
    ``eta`` and compensating with ``delta`` leaves the right segment's fit
    unchanged. In Fisher-orthonormal coordinates the joint information in
    ``(eta, delta)`` is, per direction,

        [[L1 + L2,  L2],
         [     L2,  L2]]

    so the information that actually constrains ``delta`` is the Schur
    complement -- the *profile* information left after ``eta`` is optimised
    away:

        J_eff = L2 - L2^2 / (L1 + L2) = L1 * L2 / (L1 + L2),

    the harmonic-style combination of the two segment lengths. The penalty is
    therefore

        (d / 2) * log(1 + tau^2 * L1 * L2 / (L1 + L2)).

    Using ``L2`` alone here would be right only if the shared state were known,
    or if the left segment were infinitely informative; on a balanced split it
    overstates the effective information by a factor of two (``L/2`` against the
    correct ``L/4``), which inflates the penalty by ``(d/2) log 2``. That is a
    bounded error, so it leaves the leading coefficient alone -- but Model D's
    whole contribution is the bounded term, so it is exactly the term that must
    be right.

    This is the isotropic reduction of the general correction
    ``0.5 * logdet(I + tau^2 J_eff)``; it is exact under the local
    identity-information approximation at the Fisher-orthonormal reference
    point, and approximate away from it.

    The two limits are exactly the models being interpolated:

    * ``tau -> 0``  : cost 0, delta pinned at zero -- the exact orbit, Model C.
    * ``tau -> inf``: cost ``(d/2) log(tau^2 L1 L2 / (L1 + L2))``, i.e. a leading
      ``(d/2) log n`` -- the independent-subspace rate of Model B. (The limit is
      improper; the divergent ``d log tau`` is the price of an unbounded prior.)

    Note what this does *not* do. For any **fixed** ``tau > 0`` the leading
    coefficient is ``d/2`` -- Model B's, not something in between. The
    interpolation lives in the bounded term, i.e. in the finite-sample regime,
    which is exactly where "how much deviation can be tolerated" is a question.
    A genuine interpolation of the *leading* coefficient needs ``tau`` shrinking
    with ``n``.
    """
    if deviation_scale < 0:
        raise ValueError(f"deviation_scale must be >= 0, got {deviation_scale}")
    if n_left <= 0 or n_right <= 0:
        raise ValueError("both segments must be non-empty")
    if deviation_scale == 0:
        return 0.0
    effective = (n_left * n_right) / (n_left + n_right)
    return 0.5 * dim * float(np.log1p(effective * deviation_scale**2))


def _neg_joint_and_grad(params, counts_left, counts_right, B, rotation, tau2):
    """Objective for one fixed shift: shared state plus a shrunk deviation."""
    d = B.shape[1]
    theta, delta = params[:d], params[d:]
    right_coord = rotation @ theta + delta

    nll_left, grad_left = _neg_loglik_and_grad(theta, counts_left, B)
    nll_right, grad_right = _neg_loglik_and_grad(right_coord, counts_right, B)

    value = nll_left + nll_right
    grad_theta = grad_left + rotation.T @ grad_right
    grad_delta = grad_right.copy()
    if tau2 > 0:
        value += 0.5 * float(delta @ delta) / tau2
        grad_delta += delta / tau2
    return value, np.concatenate([grad_theta, grad_delta])


def fit_approximate_orbit(counts_left, counts_right, m: int, shift: int, deviation_scale: float):
    """Fit ``eta_R = R^shift eta_L + delta`` with a Gaussian code on ``delta``.

    Returns ``(theta, delta, penalised_loglik)``. With ``deviation_scale == 0``
    the deviation is pinned at zero and this reduces exactly to the shared-orbit
    fit of Model C.
    """
    cL = np.asarray(counts_left, dtype=float)
    cR = np.asarray(counts_right, dtype=float)
    d = fundamental_dimension(m)
    rotation = _rotation(m, shift)

    if deviation_scale == 0:
        # delta == 0: aligned pooling, identical to Model C's fit.
        theta, ll = fit_fundamental(cL + np.roll(cR, -shift), m)
        return theta, np.zeros(d), ll

    B = fourier_design_matrix(m)
    tau2 = deviation_scale**2
    theta0, _ = fit_fundamental(cL + np.roll(cR, -shift), m)
    start = np.concatenate([theta0, np.zeros(d)])

    res = minimize(
        _neg_joint_and_grad,
        start,
        args=(cL, cR, B, rotation, tau2),
        jac=True,
        method="L-BFGS-B",
        options={"maxiter": 1000, "ftol": 1e-14, "gtol": 1e-10},
    )
    n = cL.sum() + cR.sum()
    if not (res.success or np.abs(res.jac).max() <= GRADIENT_TOLERANCE * max(1.0, n)):
        global _FIT_FAILURES
        _FIT_FAILURES += 1
    return res.x[:d], res.x[d:], float(-res.fun)


def approximate_orbit_detector(
    counts_left: np.ndarray, counts_right: np.ndarray, m: int, deviation_scale: float = 0.1
) -> DetectorResult:
    """Model D: one shared state, a relative shift, and a *shrunk* deviation.

    The alternative is ``eta_R = R^r eta_L + delta`` with a Gaussian code on
    ``delta`` of width ``deviation_scale``. This nests the two constrained
    models: ``deviation_scale = 0`` is exactly Model C, and large
    ``deviation_scale`` approaches Model B's rate.

    The shift ranges over nonidentity elements, matching Model C, so the
    ``deviation_scale -> 0`` limit is that detector exactly rather than merely
    asymptotically.
    """
    cL = np.asarray(counts_left, dtype=float)
    cR = np.asarray(counts_right, dtype=float)
    _, ll_null = fit_fundamental(cL + cR, m)

    best = (-np.inf, None, None)
    for s in range(1, m):
        _, delta, ll = fit_approximate_orbit(cL, cR, m, s, deviation_scale)
        if ll > best[0]:
            best = (ll, s, delta)
    ll_alt, best_shift, _ = best

    dim = fundamental_dimension(m)
    pen = label_cost(m) + deviation_penalty(
        dim, int(cL.sum()), int(cR.sum()), deviation_scale
    )
    gain = ll_alt - ll_null
    # Asymptotically a fixed positive scale pays the full dimension; only at
    # exactly zero does the continuous increment vanish.
    increment = 0.0 if deviation_scale == 0 else float(dim)
    return DetectorResult("approximate_orbit", gain, pen, gain - pen, increment, best_shift)

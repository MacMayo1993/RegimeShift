"""The shared-orbit statistic at and near orbit collapse (Section 6.4).

Proposition 1 gives Model C a zero continuous-dimension increment on a *regular*
orbit. It says nothing about what happens as the shared state approaches the
origin, where the ``g`` orbit points merge and the relative shift stops being
identifiable. That is the regime the raw MDL rule actually fails in, and this
module is the limit law that describes it.

The picture is not a singular one. For each fixed nonidentity ``r`` the
alternative ``eta -> (p(eta), p(R^r eta))`` is a perfectly regular
``d_fund``-dimensional family at the origin, with positive-definite quadratic KL
expansion; Model C's parameter space is a finite *union* of ``g - 1`` such
branches, all passing through the same point. A finite union of regular branches
keeps the regular leading exponent -- the branch multiplicity moves into the
bounded term -- so no new ``log L`` coefficient appears and Proposition 1's
``Delta d_C = 0`` survives collapse. What does happen is that the profiled
log-likelihood gain converges to a nondegenerate ``O_p(1)`` law: a maximum, over
the nonidentity shifts, of differences of Gaussian projection energies.

Writing ``U_1, U_2 ~ N(0, I_d)`` for the two segments' normalised scores,
``rho = L_1 / L``, and taking the shared state at the local value
``eta_L = h / sqrt(L)`` with ``w_r = R^r h - h``:

    gain_r = || sqrt(rho) U_1 + sqrt(1-rho) R^-r U_2
                                    - (1-rho) R^-r w_r ||^2 / 2
             - || sqrt(rho) U_1 + sqrt(1-rho) U_2 ||^2 / 2
             + sqrt(1-rho) U_2 . w_r  -  (1-rho) ||w_r||^2 / 2

    W_{g,rho}(h) = max over r != e of gain_r.

``h = 0`` is exact collapse, where every ``w_r`` vanishes and the law reduces to
a difference of two projection energies. Large ``||h||`` is a regular orbit
state, where the ``-||w_r||^2`` term drives the statistic to minus infinity and
the detector stops false-alarming.

**Scaling conventions, stated once.** ``h`` is the local coordinate of the shared
state under ``eta_L = h / sqrt(L)`` with ``L = L_1 + L_2`` the *total* length,
and ``U_1, U_2`` are the two segments' scores normalised to ``N(0, I_d)``. That
is the convention every function here uses.

On a balanced split the same law has a much shorter statement, in the
*per-segment* convention ``eta_n = h_n / sqrt(n)`` with ``n = L_1 = L_2``. Put
``Y_i = h_n + Z_i`` with ``Z_1, Z_2`` iid ``N(0, I_d)``; then

    W_g(h_n) = 1/4 max over r != e of [ ||Y_1 + R^-r Y_2||^2 - ||Y_1 + Y_2||^2 ].

Signal plus noise in each segment, compared aligned against unaligned. It is the
same distribution -- the expansion of the general form differs from it by
``h . R^-r h - h . R^r h``, which is zero for a rotation -- and the two agree to
machine precision under ``h = sqrt(2) * h_n``, which
``tests/test_collapse.py`` checks. The naive extension of the ``Y`` form to
unbalanced splits, weighting ``Y_i`` by ``sqrt(rho_i)``, is *not* the general law:
it happens to be right at ``h = 0`` and is wrong everywhere else.

Two consequences the raw rule cannot survive. The limit does not depend on ``L``,
so the zero-threshold false-positive rate does not vanish as data accumulates;
and it is not stochastically dominated by any constant label cost, so no choice
of ``log(g-1)``, ``log g`` or any other structural constant controls it. The
correct object is a *critical value* read off this law -- which is a frequentist
threshold at a chosen level, not a codelength constant, and the two must not be
conflated.
"""

from __future__ import annotations

import numpy as np

from .fourier import fundamental_dimension, rotation_matrix

__all__ = [
    "collapse_law",
    "critical_value",
    "false_alarm_rate",
    "least_favourable_sweep",
]


def _rotations(m: int):
    return [(rotation_matrix(m, r), rotation_matrix(m, -r)) for r in range(1, m)]


def collapse_law(
    m: int,
    h: np.ndarray | None = None,
    rho: float = 0.5,
    reps: int = 200_000,
    rng: np.random.Generator | None = None,
) -> np.ndarray:
    """Draws from ``W_{g,rho}(h)``, the limit of the shared-orbit raw gain.

    ``h`` is the local coordinate of the shared state, ``eta_L = h / sqrt(L)``;
    ``None`` means exact collapse. Returns ``reps`` draws, in nats, on the same
    scale as :func:`~regimeshift.detectors.shared_orbit_detector`'s ``raw_gain``.
    """
    if not 0.0 < rho < 1.0:
        raise ValueError(f"rho must lie in (0, 1), got {rho}")
    d = fundamental_dimension(m)
    h = np.zeros(d) if h is None else np.atleast_1d(np.asarray(h, dtype=float))
    if h.shape != (d,):
        raise ValueError(f"h must have shape ({d},) for m={m}, got {h.shape}")
    rng = np.random.default_rng() if rng is None else rng

    a, b = np.sqrt(rho), np.sqrt(1.0 - rho)
    U1 = rng.normal(size=(reps, d))
    U2 = rng.normal(size=(reps, d))
    base = np.sum((a * U1 + b * U2) ** 2, axis=1) / 2.0

    best = np.full(reps, -np.inf)
    for R, Rinv in _rotations(m):
        w = R @ h - h
        shifted = a * U1 + b * (U2 @ Rinv.T) - (1.0 - rho) * (Rinv @ w)[None, :]
        gain = (
            np.sum(shifted ** 2, axis=1) / 2.0
            - base
            + b * (U2 @ w)
            - (1.0 - rho) * (w @ w) / 2.0
        )
        best = np.maximum(best, gain)
    return best


def critical_value(
    m: int,
    alpha: float = 0.05,
    rho: float = 0.5,
    h: np.ndarray | None = None,
    reps: int = 200_000,
    rng: np.random.Generator | None = None,
) -> float:
    """The ``1 - alpha`` quantile of the collapse law, in nats.

    This is a **critical value**, not a codelength constant: it is defined by a
    chosen error rate and moves with ``alpha``. It answers "what threshold makes
    the raw shared-orbit rule an asymptotically valid test near collapse", not
    "what should the relative shift cost to encode".
    """
    draws = collapse_law(m, h=h, rho=rho, reps=reps, rng=rng)
    return float(np.quantile(draws, 1.0 - alpha))


def false_alarm_rate(
    m: int,
    threshold: float | None = None,
    rho: float = 0.5,
    h: np.ndarray | None = None,
    reps: int = 200_000,
    rng: np.random.Generator | None = None,
) -> float:
    """Asymptotic rate at which the raw rule fires against a fixed threshold.

    ``threshold=None`` uses the label cost ``log(g-1)`` the two-part code
    charges, which is what the detector applies.
    """
    if threshold is None:
        threshold = float(np.log(m - 1)) if m > 2 else 0.0
    return float(np.mean(collapse_law(m, h=h, rho=rho, reps=reps, rng=rng) > threshold))


def least_favourable_sweep(
    m: int,
    alpha: float = 0.05,
    rho: float = 0.5,
    magnitudes: np.ndarray | None = None,
    n_angles: int = 9,
    reps: int = 200_000,
    seed: int = 20260713,
) -> dict:
    """Is exact collapse the least favourable local null?

    Sweeps ``h`` over magnitudes and, in two dimensions, one orbit sector, and
    compares each ``1 - alpha`` quantile against the one at ``h = 0``. The sweep
    uses **common random numbers** across every ``h``: without them the maximum
    over a few hundred separately-estimated quantiles is biased upward by its own
    Monte Carlo noise, which is large enough to manufacture an excess where there
    is none.

    Common random numbers reduce the noise but do not remove selection bias
    entirely: the reported maximum is still a maximum over a grid of correlated
    estimates. The sweep therefore re-evaluates its own argmax against the origin
    on an **independent** draw, and returns that confirmation separately. A
    result is only worth quoting when both the swept excess and the independent
    re-check come back at zero.

    What this establishes is a numerical finding over the evaluated grid, not
    uniform size control. Whether exact collapse is least favourable for every
    ``h``, every ``rho`` and every ``alpha`` is a conjecture; this function
    supplies evidence for it at the settings actually tested.
    """
    d = fundamental_dimension(m)
    rng = np.random.default_rng(seed + m)
    a, b = np.sqrt(rho), np.sqrt(1.0 - rho)
    U1 = rng.normal(size=(reps, d))
    U2 = rng.normal(size=(reps, d))
    base = np.sum((a * U1 + b * U2) ** 2, axis=1) / 2.0
    rots = [(R, Rinv, U2 @ Rinv.T) for R, Rinv in _rotations(m)]

    def quantile(h):
        best = np.full(reps, -np.inf)
        for R, Rinv, U2R in rots:
            w = R @ h - h
            shifted = a * U1 + b * U2R - (1.0 - rho) * (Rinv @ w)[None, :]
            gain = (np.sum(shifted ** 2, axis=1) / 2.0 - base
                    + b * (U2 @ w) - (1.0 - rho) * (w @ w) / 2.0)
            best = np.maximum(best, gain)
        return float(np.quantile(best, 1.0 - alpha))

    if magnitudes is None:
        magnitudes = np.concatenate([np.linspace(0.0, 1.5, 31), np.linspace(1.6, 6.0, 23)])
    angles = [0.0] if d == 1 else np.linspace(0.0, 2.0 * np.pi / m, n_angles)

    at_zero = quantile(np.zeros(d))
    best_q, best_h = at_zero, np.zeros(d)
    exceed = 0
    for mag in magnitudes:
        candidates = ([np.array([mag]), np.array([-mag])] if d == 1
                      else [mag * np.array([np.cos(t), np.sin(t)]) for t in angles])
        for h in candidates:
            q = quantile(h)
            if q > best_q:
                best_q, best_h = q, h
            if q - at_zero > 0.01:
                exceed += 1

    # Independent re-check at the selected argmax, on fresh draws, so the
    # headline number is not itself the product of grid selection.
    confirm = np.random.default_rng(seed + 7919 + m)
    indep_zero = float(np.quantile(collapse_law(m, None, rho, reps, confirm), 1.0 - alpha))
    confirm = np.random.default_rng(seed + 104729 + m)
    indep_best = float(np.quantile(collapse_law(m, best_h, rho, reps, confirm), 1.0 - alpha))

    return {
        "m": m,
        "alpha": alpha,
        "rho": rho,
        "quantile_at_collapse": at_zero,
        "sweep_maximum": best_q,
        "excess": best_q - at_zero,
        "argmax_norm": float(np.linalg.norm(best_h)),
        "points_exceeding_by_0.01": exceed,
        "independent_excess_at_argmax": indep_best - indep_zero,
        "reps": reps,
    }

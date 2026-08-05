"""Fisher-orthonormal cyclic Fourier geometry for direct categorical models.

Implements Section 5 of the manuscript. The direct model identifies the group
order ``m`` of the cyclic group C_m with the alphabet size, and the group acts
by cyclically permuting category coordinates.

Conventions
-----------
* Tangent vectors live in ``T = {v in R^m : sum_j v_j = 0}``.
* At the uniform distribution ``u = (1/m, ..., 1/m)`` the Fisher inner product
  is ``<v, w>_F = m * sum_j v_j w_j``.
* The fundamental (first) Fourier mode has real dimension ``d = 1`` for ``m = 2``
  (the sign representation) and ``d = 2`` for ``m >= 3``.
* The group element ``g^s`` maps a probability vector ``p`` to ``np.roll(p, s)``,
  i.e. ``(g^s p)_j = p_{(j - s) mod m}``. In fundamental coordinates this is a
  planar rotation by ``2 pi s / m`` (a sign flip when ``m = 2``).
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "fundamental_dimension",
    "full_dimension",
    "fisher_inner_product",
    "fourier_design_matrix",
    "fundamental_tangent_basis",
    "rotation_matrix",
    "probabilities",
    "roll_probabilities",
    "higher_mode_logits",
]


def full_dimension(m: int) -> int:
    """Continuous dimension of the unrestricted ``m``-category simplex."""
    _check_m(m)
    return m - 1


def fundamental_dimension(m: int) -> int:
    """Real dimension of the fundamental invariant component."""
    _check_m(m)
    return 1 if m == 2 else 2


def _check_m(m: int) -> None:
    if int(m) != m or m < 2:
        raise ValueError(f"group order must be an integer >= 2, got {m!r}")


def fisher_inner_product(u: np.ndarray, v: np.ndarray, m: int | None = None) -> float:
    """Fisher inner product of two tangent vectors at the uniform distribution."""
    u = np.asarray(u, dtype=float)
    v = np.asarray(v, dtype=float)
    if m is None:
        m = u.shape[-1]
    return float(m * np.dot(u, v))


def fourier_design_matrix(m: int) -> np.ndarray:
    """Logit design matrix ``B`` of shape ``(m, d)`` for the fundamental family.

    The family is ``p(theta) = softmax(B @ theta)``. The scaling is chosen so
    that ``d p / d theta_a`` at ``theta = 0`` is Fisher-orthonormal, which makes
    ``|theta|_2`` the Fisher norm of the local perturbation.
    """
    _check_m(m)
    j = np.arange(m)
    if m == 2:
        return np.array([[1.0], [-1.0]])
    phi = 2.0 * np.pi * j / m
    scale = np.sqrt(2.0)
    return np.column_stack([scale * np.cos(phi), scale * np.sin(phi)])


def fundamental_tangent_basis(m: int) -> np.ndarray:
    """Fisher-orthonormal tangent basis of the fundamental component.

    Returns an array of shape ``(d, m)`` whose rows are the derivatives of
    :func:`probabilities` at ``theta = 0``.
    """
    B = fourier_design_matrix(m)
    # d p_j / d theta_a = (1/m) * (B[j, a] - mean_j B[j, a]); columns are centred.
    return (B - B.mean(axis=0, keepdims=True)).T / m


def rotation_matrix(m: int, steps: int = 1) -> np.ndarray:
    """Action of ``g^steps`` on fundamental coordinates."""
    _check_m(m)
    if m == 2:
        return np.array([[(-1.0) ** (steps % 2)]])
    phi = 2.0 * np.pi * steps / m
    return np.array([[np.cos(phi), -np.sin(phi)], [np.sin(phi), np.cos(phi)]])


def probabilities(theta: np.ndarray, m: int, extra_logits: np.ndarray | None = None) -> np.ndarray:
    """Softmax probabilities of the fundamental family at ``theta``.

    ``extra_logits`` optionally adds a component outside the fundamental
    subspace (used by the misspecification scenario).
    """
    theta = np.atleast_1d(np.asarray(theta, dtype=float))
    d = fundamental_dimension(m)
    if theta.shape != (d,):
        raise ValueError(f"theta must have shape ({d},) for m={m}, got {theta.shape}")
    logits = fourier_design_matrix(m) @ theta
    if extra_logits is not None:
        logits = logits + np.asarray(extra_logits, dtype=float)
    logits = logits - logits.max()
    w = np.exp(logits)
    return w / w.sum()


def roll_probabilities(p: np.ndarray, steps: int) -> np.ndarray:
    """Apply ``g^steps`` to a probability (or count) vector."""
    return np.roll(np.asarray(p), steps)


def higher_mode_logits(m: int, amplitude: float, mode: int = 2) -> np.ndarray:
    """Logit perturbation along a higher Fourier mode (outside the fundamental).

    For ``m = 4`` and ``mode = 2`` this is the one-dimensional sign
    representation; for ``m = 6`` it is a two-dimensional higher mode.
    """
    _check_m(m)
    if mode <= 0 or mode >= m:
        raise ValueError(f"mode must satisfy 0 < mode < m, got {mode}")
    phi = 2.0 * np.pi * mode * np.arange(m) / m
    v = np.cos(phi)
    return amplitude * (v - v.mean())

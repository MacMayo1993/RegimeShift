"""Block (codon-phase) families: the group acts on phases, not on the alphabet.

Section 11. The direct model used everywhere else identifies the group order
with the alphabet size, which confounds them: raising ``m`` simultaneously
changes the number of candidate shifts, the simplex dimension, the category
sparsity and the geometric separation of neighbouring orbit states. Both
external reviews named this. A block family separates them.

A codon-phase model has ``g = 3`` reading-frame phases and ``a = 4`` nucleotide
symbols. Each phase carries its own distribution over the alphabet, so a regime
is a tuple of ``g`` distributions and the observation for one segment is a
``(g, a)`` count array. The group permutes the phase blocks and leaves the
alphabet alone.

Note the manuscript's notation in Section 11 is the reverse of the direct
model's: there ``g`` is the group order and ``m`` the alphabet size. This module
uses ``g`` and ``a`` to keep them unmistakable.

Dimensions, matching Section 11 exactly for the codon case ``g = 3, a = 4``::

    d_full        = g (a - 1)                    = 3 x 3 = 9
    d_phase_fund  = dim(V_fund) x (a - 1)        = 2 x 3 = 6

where ``V_fund`` is the fundamental isotypic component of the phase
representation: real dimension 1 for ``g = 2`` (the sign representation) and 2
for ``g >= 3``, exactly as in the direct model. The phase representation
decomposes as ``R^g = V_triv + V_fund + (higher modes)``.

The geometry is the direct model's tensored with an ``(a - 1)``-dimensional
alphabet-contrast space: the group acts on the phase-mode index only, identically
in every contrast direction, so the rotation is ``R (x) I``. That is why the
dimensions multiply.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize

from .detectors import (
    DetectorResult,
    GRADIENT_TOLERANCE,
    label_cost,
    split_penalty,
)
from .fourier import fourier_design_matrix, fundamental_dimension, rotation_matrix

__all__ = [
    "BlockGeometry",
    "CODON",
    "validate_block_pair",
    "block_multinomial_loglik",
    "block_probabilities",
    "fit_block_fundamental",
    "block_full_detector",
    "block_fundamental_detector",
    "block_shared_orbit_detector",
    "run_block_detectors",
]


@dataclass(frozen=True)
class BlockGeometry:
    """A cyclic group of order ``g`` permuting ``g`` blocks over an alphabet of
    size ``a``."""

    g: int
    """Group order, i.e. the number of phase blocks."""
    a: int
    """Alphabet size within each block."""

    def __post_init__(self) -> None:
        if int(self.g) != self.g or self.g < 2:
            raise ValueError(f"group order must be an integer >= 2, got {self.g!r}")
        if int(self.a) != self.a or self.a < 2:
            raise ValueError(f"alphabet size must be an integer >= 2, got {self.a!r}")

    @property
    def contrast_dimension(self) -> int:
        """Free dimensions within one block: ``a - 1``."""
        return self.a - 1

    @property
    def mode_dimension(self) -> int:
        """Real dimension of the fundamental phase component: 1 or 2."""
        return fundamental_dimension(self.g)

    @property
    def full_dimension(self) -> int:
        """``g (a - 1)`` -- every block free."""
        return self.g * self.contrast_dimension

    @property
    def fundamental_dimension(self) -> int:
        """``dim(V_fund) (a - 1)`` -- the fundamental phase isotypic component."""
        return self.mode_dimension * self.contrast_dimension

    @property
    def shape(self) -> tuple[int, int]:
        return (self.g, self.a)


#: The manuscript's codon-phase case: three reading frames, four nucleotides.
CODON = BlockGeometry(g=3, a=4)


def _contrast_basis(a: int) -> np.ndarray:
    """Orthonormal basis of the mean-zero contrast space, shape ``(a - 1, a)``.

    Any orthonormal basis of ``{v : sum v = 0}`` works; the Helmert basis is used
    because it is real, deterministic and independent of the group structure, so
    nothing about the alphabet can leak into the phase geometry.
    """
    rows = []
    for k in range(1, a):
        row = np.zeros(a)
        row[:k] = 1.0
        row[k] = -float(k)
        rows.append(row / np.linalg.norm(row))
    return np.array(rows)


def block_probabilities(theta: np.ndarray, geometry: BlockGeometry) -> np.ndarray:
    """Per-block probabilities from fundamental coordinates.

    ``theta`` has shape ``(mode_dimension, a - 1)``. Returns a ``(g, a)`` array
    whose rows are probability vectors.
    """
    theta = np.asarray(theta, dtype=float)
    expected = (geometry.mode_dimension, geometry.contrast_dimension)
    if theta.shape != expected:
        raise ValueError(f"theta must have shape {expected}, got {theta.shape}")
    # phase modes (g, mode_dim) x coordinates (mode_dim, a-1) -> (g, a-1),
    # then lift the contrasts into logits over the alphabet.
    logits = fourier_design_matrix(geometry.g) @ theta @ _contrast_basis(geometry.a)
    logits = logits - logits.max(axis=1, keepdims=True)
    weights = np.exp(logits)
    return weights / weights.sum(axis=1, keepdims=True)


def block_multinomial_loglik(counts: np.ndarray) -> float:
    """Maximised unrestricted log-likelihood: each block fitted independently."""
    counts = np.asarray(counts, dtype=float)
    total = 0.0
    for row in counts:
        n = row.sum()
        if n <= 0:
            continue
        nz = row > 0
        total += float(np.sum(row[nz] * (np.log(row[nz]) - np.log(n))))
    return total


def _neg_block_loglik_and_grad(flat, counts, design, contrasts, shape):
    theta = flat.reshape(shape)
    logits = design @ theta @ contrasts
    logits = logits - logits.max(axis=1, keepdims=True)
    weights = np.exp(logits)
    probabilities = weights / weights.sum(axis=1, keepdims=True)

    loglik = float(np.sum(counts * np.log(probabilities)))
    block_totals = counts.sum(axis=1, keepdims=True)
    residual = counts - block_totals * probabilities
    grad = design.T @ residual @ contrasts.T
    return -loglik, -grad.ravel()


def fit_block_fundamental(counts: np.ndarray, geometry: BlockGeometry):
    """MLE of the fundamental coordinates for one segment's block counts."""
    counts = np.asarray(counts, dtype=float)
    if counts.shape != geometry.shape:
        raise ValueError(f"counts must have shape {geometry.shape}, got {counts.shape}")
    shape = (geometry.mode_dimension, geometry.contrast_dimension)
    if counts.sum() <= 0:
        return np.zeros(shape), 0.0

    design = fourier_design_matrix(geometry.g)
    contrasts = _contrast_basis(geometry.a)
    result = minimize(
        _neg_block_loglik_and_grad,
        np.zeros(int(np.prod(shape))),
        args=(counts, design, contrasts, shape),
        jac=True,
        method="L-BFGS-B",
        options={"maxiter": 1000, "ftol": 1e-14, "gtol": 1e-10},
    )
    n = counts.sum()
    if not (result.success or np.abs(result.jac).max() <= GRADIENT_TOLERANCE * max(1.0, n)):
        import regimeshift.detectors as detectors

        detectors._FIT_FAILURES += 1
    return result.x.reshape(shape), float(-result.fun)


def validate_block_pair(counts_left, counts_right, geometry: BlockGeometry):
    """Validate and coerce a block segment pair, as ``validate_pair`` does for the
    direct model.

    The shape check is the one that matters. ``block_full_detector`` computes its
    likelihood block by block from the array it is handed and its penalty from
    ``geometry``, so a ``(3, 4)`` array scored against ``BlockGeometry(6, 4)``
    used to return a score built from two different geometries -- a wrong number
    rather than an error, while the other two block detectors raised. Every block
    detector now routes its inputs through here so the contract is uniform.
    """
    arrays = []
    for name, counts in (("counts_left", counts_left), ("counts_right", counts_right)):
        c = np.asarray(counts, dtype=float)
        if c.shape != geometry.shape:
            raise ValueError(f"{name} must have shape {geometry.shape}, got {c.shape}")
        if not np.all(np.isfinite(c)):
            raise ValueError(f"{name} must be finite")
        if np.any(c < 0):
            raise ValueError(f"{name} must be nonnegative")
        if c.sum() <= 0:
            raise ValueError(f"{name} must be nonempty (positive total count)")
        arrays.append(c)
    return arrays[0], arrays[1]


def block_full_detector(counts_left, counts_right, geometry: BlockGeometry) -> DetectorResult:
    """Model A on blocks: every phase free on each side. Increment ``g(a-1)``."""
    cL, cR = validate_block_pair(counts_left, counts_right, geometry)
    gain = (
        block_multinomial_loglik(cL)
        + block_multinomial_loglik(cR)
        - block_multinomial_loglik(cL + cR)
    )
    dim = geometry.full_dimension
    penalty = split_penalty(dim, int(cL.sum()), int(cR.sum()))
    return DetectorResult("block_full", gain, penalty, gain - penalty, dim)


def block_fundamental_detector(counts_left, counts_right, geometry: BlockGeometry) -> DetectorResult:
    """Model B on blocks: each side free inside the fundamental phase component."""
    cL, cR = validate_block_pair(counts_left, counts_right, geometry)
    _, ll_null = fit_block_fundamental(cL + cR, geometry)
    _, ll_left = fit_block_fundamental(cL, geometry)
    _, ll_right = fit_block_fundamental(cR, geometry)
    gain = ll_left + ll_right - ll_null
    dim = geometry.fundamental_dimension
    penalty = split_penalty(dim, int(cL.sum()), int(cR.sum()))
    return DetectorResult("block_fundamental", gain, penalty, gain - penalty, dim)


def block_shared_orbit_detector(counts_left, counts_right, geometry: BlockGeometry) -> DetectorResult:
    """Model C on blocks: one shared state, the right side phase-shifted.

    The shift permutes phase blocks and leaves the alphabet untouched, which is
    the whole point of the block family.
    """
    cL, cR = validate_block_pair(counts_left, counts_right, geometry)
    _, ll_null = fit_block_fundamental(cL + cR, geometry)

    best_ll, best_shift = -np.inf, None
    for shift in range(1, geometry.g):
        aligned = cL + np.roll(cR, -shift, axis=0)
        _, ll = fit_block_fundamental(aligned, geometry)
        if ll > best_ll:
            best_ll, best_shift = ll, shift

    gain = best_ll - ll_null
    penalty = label_cost(geometry.g)
    return DetectorResult("block_shared_orbit", gain, penalty, gain - penalty, 0.0, best_shift)


def run_block_detectors(counts_left, counts_right, geometry: BlockGeometry) -> dict[str, DetectorResult]:
    """Score one block dataset with all three detectors."""
    return {
        "block_full": block_full_detector(counts_left, counts_right, geometry),
        "block_fundamental": block_fundamental_detector(counts_left, counts_right, geometry),
        "block_shared_orbit": block_shared_orbit_detector(counts_left, counts_right, geometry),
    }


def block_rotation(theta: np.ndarray, geometry: BlockGeometry, steps: int = 1) -> np.ndarray:
    """Act on fundamental coordinates by ``g^steps``: ``R (x) I``."""
    return rotation_matrix(geometry.g, steps) @ np.asarray(theta, dtype=float)

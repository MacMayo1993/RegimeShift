"""Model selection among the geometries, rather than assuming one (oracle) is known.

Every efficiency number in the manuscript is an *oracle* number: the detector
matching the data-generating geometry is chosen in advance. Both external
reviews flagged this. An analyst does not know whether a change is unrestricted,
confined to the invariant subspace, or an exact orbit -- deciding that is part of
the problem.

Doing it properly needs absolute code lengths, not the detector scores. A
detector score is ``gain - penalty`` measured against *its own* null, and those
nulls differ: Model A pools an unrestricted multinomial, while Models B, C and D
pool a fundamental-family coordinate. Comparing the scores directly would
compare quantities with different origins. What is comparable is the total
description length of the same data under each hypothesis, which is what
:func:`code_lengths` returns.

The relationship to the detectors is exact, and asserted in the tests::

    score_A = L(null_full)        - L(full)
    score_B = L(null_fundamental) - L(fundamental)
    score_C = L(null_fundamental) - L(shared_orbit)

so nothing here contradicts the three-model comparison; it re-expresses it on a
common scale and adds the null hypotheses as candidates in their own right.

Parameter costs are the same BIC-style ``(d/2) log n`` used throughout: a block
of ``d`` continuous parameters fitted from ``n`` observations costs
``(d/2) log n`` nats. As elsewhere in this package these are regular
approximations, not exact universal codes.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .detectors import (
    deviation_penalty,
    fit_approximate_orbit,
    fit_fundamental,
    label_cost,
    multinomial_loglik,
)
from .fourier import full_dimension, fundamental_dimension

__all__ = [
    "CANDIDATES",
    "ALTERNATIVES",
    "Selection",
    "code_lengths",
    "select_model",
    "generating_family_for_scenario",
]

#: Candidate hypotheses, in increasing order of structure. The order is also the
#: tie-break order in :func:`select_model`: less structured first, so a tie never
#: becomes a claim of structure the data cannot support.
CANDIDATES = (
    "null_full",
    "null_fundamental",
    "full",
    "fundamental",
    "shared_orbit",
    "approximate_orbit",
)

#: The candidates that assert a change occurred.
ALTERNATIVES = ("full", "fundamental", "shared_orbit", "approximate_orbit")

#: Which candidate family a scenario is *generated* from.
#:
#: This is not the same as which candidate gives the shortest description. The
#: manuscript's ``independent_fundamental`` scenario is generated with free
#: coordinates inside the fundamental family, yet at ``m = 5`` and ``m = 6`` it
#: sits close enough to a one-step orbit that ``approximate_orbit`` is the better
#: code -- see ``docs/paper-notes.md``. Use this for bookkeeping, not as ground
#: truth for whether a selection was "correct".
_GENERATING_FAMILY = {
    "exact_orbit": "shared_orbit",
    "independent_fundamental": "fundamental",
    "higher_mode": "full",
    "approximate_orbit": "approximate_orbit",
}


def generating_family_for_scenario(scenario: str, deviation: float = 0.0) -> str:
    """The candidate family a scenario is generated from.

    ``approximate_orbit`` with zero deviation is an exact orbit, so it maps to
    ``shared_orbit`` rather than to the interpolating model.

    Read the caveat on ``_GENERATING_FAMILY``: being generated from a family
    does not make that family the shortest description of the result.
    """
    if scenario == "approximate_orbit" and deviation == 0.0:
        return "shared_orbit"
    try:
        return _GENERATING_FAMILY[scenario]
    except KeyError:
        raise ValueError(f"no generating family recorded for scenario {scenario!r}") from None


@dataclass(frozen=True)
class Selection:
    """Outcome of selecting among the candidate geometries."""

    selected: str
    """Candidate with the shortest total description length."""
    code_lengths: dict[str, float]
    """Total description length of the data under each candidate, in nats."""
    declared_change: bool
    """Whether the selected candidate is an alternative rather than a null."""
    selected_shift: int | None
    """Relative group element, when the selected candidate has one."""
    tied: tuple[str, ...] = ()
    """Candidates whose code length is indistinguishable from the selected one.

    Usually just the selection itself. It has more than one entry when two
    candidates describe the *same* hypothesis: at ``m = 2`` and ``m = 3`` the
    fundamental component spans the whole nontrivial tangent space, so ``full``
    and ``fundamental`` are the same model and their code lengths agree to
    numerical precision. Selecting between them would be reading floating-point
    noise, so the tie is reported instead."""

    @property
    def margin(self) -> float:
        """Nats by which the selected candidate beats the runner-up.

        A small margin means the data does not distinguish the geometries, which
        is information the oracle comparison throws away.
        """
        ordered = sorted(self.code_lengths.values())
        return float(ordered[1] - ordered[0]) if len(ordered) > 1 else float("inf")


def code_lengths(
    counts_left: np.ndarray,
    counts_right: np.ndarray,
    m: int,
    deviation_scale: float = 0.05,
    candidates: tuple[str, ...] = CANDIDATES,
) -> dict[str, float]:
    """Total description length of the two segments under each candidate, in nats.

    Lower is better. The multinomial coefficient is omitted throughout, which is
    a constant of the data and so cancels in every comparison.
    """
    cL = np.asarray(counts_left, dtype=float)
    cR = np.asarray(counts_right, dtype=float)
    nL, nR = cL.sum(), cR.sum()
    if nL <= 0 or nR <= 0:
        raise ValueError("both segments must be non-empty")
    pooled = cL + cR
    n = nL + nR

    d_full = full_dimension(m)
    d_fund = fundamental_dimension(m)
    lengths: dict[str, float] = {}

    if "null_full" in candidates:
        lengths["null_full"] = -multinomial_loglik(pooled) + 0.5 * d_full * np.log(n)
    if "null_fundamental" in candidates:
        _, ll = fit_fundamental(pooled, m)
        lengths["null_fundamental"] = -ll + 0.5 * d_fund * np.log(n)
    if "full" in candidates:
        lengths["full"] = (
            -multinomial_loglik(cL) - multinomial_loglik(cR)
            + 0.5 * d_full * (np.log(nL) + np.log(nR))
        )
    if "fundamental" in candidates:
        _, ll_left = fit_fundamental(cL, m)
        _, ll_right = fit_fundamental(cR, m)
        lengths["fundamental"] = (
            -ll_left - ll_right + 0.5 * d_fund * (np.log(nL) + np.log(nR))
        )
    if "shared_orbit" in candidates:
        best = max(fit_fundamental(cL + np.roll(cR, -s), m)[1] for s in range(1, m))
        lengths["shared_orbit"] = -best + 0.5 * d_fund * np.log(n) + label_cost(m)
    if "approximate_orbit" in candidates:
        best = max(
            fit_approximate_orbit(cL, cR, m, s, deviation_scale)[2] for s in range(1, m)
        )
        lengths["approximate_orbit"] = (
            -best
            + 0.5 * d_fund * np.log(n)
            + label_cost(m)
            + deviation_penalty(d_fund, int(nR), deviation_scale)
        )
    return lengths


def _selected_shift(name: str, counts_left, counts_right, m, deviation_scale) -> int | None:
    if name not in ("shared_orbit", "approximate_orbit"):
        return None
    cL = np.asarray(counts_left, dtype=float)
    cR = np.asarray(counts_right, dtype=float)
    if name == "shared_orbit":
        scores = [fit_fundamental(cL + np.roll(cR, -s), m)[1] for s in range(1, m)]
    else:
        scores = [fit_approximate_orbit(cL, cR, m, s, deviation_scale)[2] for s in range(1, m)]
    return int(np.argmax(scores)) + 1


def select_model(
    counts_left: np.ndarray,
    counts_right: np.ndarray,
    m: int,
    deviation_scale: float = 0.05,
    candidates: tuple[str, ...] = CANDIDATES,
    tie_tolerance: float = 1e-8,
) -> Selection:
    """Pick the candidate geometry with the shortest total description length.

    This answers both questions at once -- whether a change occurred, and what
    kind -- without being told the answer to the second.

    Candidates within ``tie_tolerance`` nats of the best are treated as tied and
    reported in :attr:`Selection.tied`. Ties are broken toward the *less*
    structured candidate, in :data:`CANDIDATES` order, so a tie never becomes a
    claim of structure the data cannot support. This matters at ``m = 2`` and
    ``m = 3``, where ``full`` and ``fundamental`` are literally the same model.
    """
    lengths = code_lengths(counts_left, counts_right, m, deviation_scale, candidates)
    best = min(lengths.values())
    tied = tuple(name for name in candidates if name in lengths
                 and lengths[name] - best <= tie_tolerance)
    selected = tied[0]
    return Selection(
        selected=selected,
        code_lengths=lengths,
        declared_change=selected in ALTERNATIVES,
        selected_shift=_selected_shift(selected, counts_left, counts_right, m, deviation_scale),
        tied=tied,
    )

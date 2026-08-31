"""Monte Carlo engine: null calibration, alternative power, raw-score records.

The experiment separates two analyses (Section 8.3):

* **Raw MDL score analysis** -- the mean uncalibrated score is regressed on
  ``n * gain`` and ``log n`` to recover the theoretical penalty coefficient.
* **Calibrated power analysis** -- an additive critical value is set to the
  empirical 95th percentile of the null scores at each configuration, and power
  is estimated at that common nominal false-positive target.

Every configuration draws from a deterministic configuration-specific seed, so
results are reproducible and independent of execution order or worker count.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from typing import Sequence

import numpy as np

from .detectors import fit_failure_count, reset_fit_failures, run_all_detectors
from .gains import population_gains
from .scenarios import Segments, build_segments

__all__ = [
    "BASE_SEED",
    "DETECTION_PATTERNS",
    "Config",
    "config_seed",
    "sample_counts",
    "run_config",
    "build_grid",
]

#: Base random seed of the production run reported in the manuscript.
BASE_SEED = 20260713

DETECTOR_NAMES = ("full", "fundamental", "shared_orbit")

#: Column names for the joint calibrated-detection pattern counts.
#:
#: Every alternative dataset is scored by *all three* detectors, so each one
#: yields a triple of calibrated detection indicators. ``power_calibrated``
#: keeps only the three marginals, which is enough to report power but not to
#: resample the detectors together: a bootstrap built from the marginals alone
#: treats detectors as independent when they are strongly positively correlated,
#: and that inflates every interval on a *ratio* of two of them.
#:
#: These eight counts are the joint distribution's sufficient statistic. Bit
#: ``i`` of the suffix is detection by ``DETECTOR_NAMES[i]``, so ``pattern_101``
#: counts datasets the full and shared-orbit detectors caught and the
#: fundamental detector missed. They sum to ``n_alt``, and summing the four
#: patterns whose bit ``i`` is set recovers ``power_calibrated`` for that
#: detector exactly -- which :mod:`regimeshift.analysis` asserts before using
#: them.
DETECTION_PATTERNS = tuple(
    "pattern_" + format(code, "03b") for code in range(2 ** len(DETECTOR_NAMES))
)


@dataclass(frozen=True)
class Config:
    """One point of the Monte Carlo design."""

    m: int
    scenario: str
    effect: float
    segment_length: int
    """Length of the *left* segment. With the default balanced split the right
    segment matches it and the total is twice this."""
    n_alt: int = 500
    n_null: int = 1000
    alpha: float = 0.05
    split_fraction: float = 0.5
    """Fraction ``rho`` of the total length falling in the left segment.

    Section 4.2 predicts that ``rho`` shifts only the bounded term of the
    complexity increment, leaving the coefficient of ``log n`` at ``d/2``. The
    default 0.5 reproduces the manuscript's balanced design exactly."""

    def __post_init__(self) -> None:
        if not 0.0 < self.split_fraction < 1.0:
            raise ValueError(f"split_fraction must lie in (0, 1), got {self.split_fraction}")
        if self.n_right <= 0:
            raise ValueError("split_fraction leaves the right segment empty")

    @property
    def n_left(self) -> int:
        return self.segment_length

    @property
    def total_length(self) -> int:
        return int(round(self.segment_length / self.split_fraction))

    @property
    def n_right(self) -> int:
        return self.total_length - self.segment_length

    @property
    def key(self) -> tuple:
        return (self.m, self.scenario, self.effect, self.segment_length, self.split_fraction)


def config_seed(config: Config, base_seed: int = BASE_SEED) -> int:
    """Deterministic per-configuration seed, independent of grid ordering."""
    payload = f"{base_seed}|{config.m}|{config.scenario}|{config.effect!r}|{config.segment_length}"
    # The split fraction joins the payload only when it is not the balanced
    # default, so adding this field left every existing balanced-split seed --
    # and therefore every existing checkpoint -- unchanged.
    if config.split_fraction != 0.5:
        payload += f"|{config.split_fraction!r}"
    # SHA-256 rather than hash(): stable across processes and Python runs.
    digest = hashlib.sha256(payload.encode()).digest()
    return int.from_bytes(digest[:8], "big") % (2**32 - 1)


def sample_counts(p: np.ndarray, n: int, rng: np.random.Generator) -> np.ndarray:
    """Draw a multinomial count vector of ``n`` categorical observations."""
    return rng.multinomial(n, np.asarray(p, dtype=float))


def run_config(config: Config, base_seed: int = BASE_SEED, segments: Segments | None = None) -> list[dict]:
    """Run one configuration and return one result row per detector.

    Null samples are drawn from the scenario's no-change distribution on both
    sides; the empirical ``1 - alpha`` quantile (higher interpolation) of each
    detector's null scores becomes that detector's calibrated critical value.
    Alternative samples are then scored against those fixed thresholds.
    """
    if segments is None:
        segments = build_segments(config.m, config.scenario, config.effect)
    rng = np.random.default_rng(config_seed(config, base_seed))
    m, nL, nR = config.m, config.n_left, config.n_right
    reset_fit_failures()

    null_scores = {name: np.empty(config.n_null) for name in DETECTOR_NAMES}
    p_null = segments.p_null
    for i in range(config.n_null):
        cL = sample_counts(p_null, nL, rng)
        cR = sample_counts(p_null, nR, rng)
        for name, res in run_all_detectors(cL, cR, m).items():
            null_scores[name][i] = res.score

    critical = {
        name: float(np.quantile(vals, 1.0 - config.alpha, method="higher"))
        for name, vals in null_scores.items()
    }

    alt_scores = {name: np.empty(config.n_alt) for name in DETECTOR_NAMES}
    alt_gain = {name: np.empty(config.n_alt) for name in DETECTOR_NAMES}
    shift_hits = np.zeros(config.n_alt)
    for i in range(config.n_alt):
        cL = sample_counts(segments.p_left, nL, rng)
        cR = sample_counts(segments.p_right, nR, rng)
        results = run_all_detectors(cL, cR, m)
        for name, res in results.items():
            alt_scores[name][i] = res.score
            alt_gain[name][i] = res.raw_gain
        if segments.planted_shift is not None:
            shift_hits[i] = float(results["shared_orbit"].selected_shift == segments.planted_shift)

    gains = population_gains(segments.p_left, segments.p_right, m, w_left=nL / (nL + nR))
    failures = fit_failure_count()

    # Joint calibrated-detection patterns, retained so that a ratio bootstrap can
    # resample the detectors together rather than independently. See
    # DETECTION_PATTERNS for the encoding.
    detected = np.stack(
        [alt_scores[name] > critical[name] for name in DETECTOR_NAMES], axis=1
    )
    codes = detected @ (1 << np.arange(len(DETECTOR_NAMES) - 1, -1, -1))
    pattern_counts = np.bincount(codes, minlength=len(DETECTION_PATTERNS))

    rows = []
    for name in DETECTOR_NAMES:
        rows.append(
            {
                "m": m,
                "scenario": config.scenario,
                "effect": config.effect,
                "segment_length": config.segment_length,
                "split_fraction": config.split_fraction,
                "n_left": nL,
                "n_right": nR,
                "total_length": config.total_length,
                "detector": name,
                "population_gain": gains[name],
                "expected_gain_total": gains[name] * config.total_length,
                "mean_raw_gain": float(alt_gain[name].mean()),
                "mean_score": float(alt_scores[name].mean()),
                "sd_score": float(alt_scores[name].std(ddof=1)),
                "critical_value": critical[name],
                "null_rate_zero_threshold": float((null_scores[name] > 0).mean()),
                "power_zero_threshold": float((alt_scores[name] > 0).mean()),
                "power_calibrated": float((alt_scores[name] > critical[name]).mean()),
                "shift_accuracy": (
                    float(shift_hits.mean())
                    if (name == "shared_orbit" and segments.planted_shift is not None)
                    else float("nan")
                ),
                "n_alt": config.n_alt,
                "n_null": config.n_null,
                "alpha": config.alpha,
                "optimizer_failures": failures,
                "seed": config_seed(config, base_seed),
                **{
                    column: int(count)
                    for column, count in zip(DETECTION_PATTERNS, pattern_counts)
                },
            }
        )
    return rows


def build_grid(
    groups: Sequence[int],
    scenarios: Sequence[str],
    effects: Sequence[float],
    segment_lengths: Sequence[int],
    n_alt: int = 500,
    n_null: int = 1000,
    alpha: float = 0.05,
    split_fractions: Sequence[float] = (0.5,),
) -> list[Config]:
    """Cartesian design, skipping combinations a scenario does not support."""
    configs: list[Config] = []
    for m in groups:
        for scenario in scenarios:
            if scenario == "higher_mode" and m < 4:
                continue
            for effect in effects:
                for length in segment_lengths:
                    for rho in split_fractions:
                        configs.append(
                            Config(m, scenario, float(effect), int(length),
                                   n_alt, n_null, alpha, float(rho))
                        )
    return configs

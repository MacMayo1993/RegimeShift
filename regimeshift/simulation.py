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

from .detectors import run_all_detectors
from .gains import population_gains
from .scenarios import Segments, build_segments

__all__ = [
    "BASE_SEED",
    "Config",
    "config_seed",
    "sample_counts",
    "run_config",
    "build_grid",
]

#: Base random seed of the production run reported in the manuscript.
BASE_SEED = 20260713

DETECTOR_NAMES = ("full", "fundamental", "shared_orbit")


@dataclass(frozen=True)
class Config:
    """One point of the Monte Carlo design."""

    m: int
    scenario: str
    effect: float
    segment_length: int
    """Length of *each* side; total length is twice this."""
    n_alt: int = 500
    n_null: int = 1000
    alpha: float = 0.05

    @property
    def total_length(self) -> int:
        return 2 * self.segment_length

    @property
    def key(self) -> tuple:
        return (self.m, self.scenario, self.effect, self.segment_length)


def config_seed(config: Config, base_seed: int = BASE_SEED) -> int:
    """Deterministic per-configuration seed, independent of grid ordering."""
    payload = f"{base_seed}|{config.m}|{config.scenario}|{config.effect!r}|{config.segment_length}"
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
    m, nL = config.m, config.segment_length
    nR = config.segment_length

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

    gains = population_gains(segments.p_left, segments.p_right, m)

    rows = []
    for name in DETECTOR_NAMES:
        rows.append(
            {
                "m": m,
                "scenario": config.scenario,
                "effect": config.effect,
                "segment_length": config.segment_length,
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
                "seed": config_seed(config, base_seed),
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
) -> list[Config]:
    """Cartesian design, skipping combinations a scenario does not support."""
    configs: list[Config] = []
    for m in groups:
        for scenario in scenarios:
            if scenario == "higher_mode" and m < 4:
                continue
            for effect in effects:
                for length in segment_lengths:
                    configs.append(Config(m, scenario, float(effect), int(length), n_alt, n_null, alpha))
    return configs

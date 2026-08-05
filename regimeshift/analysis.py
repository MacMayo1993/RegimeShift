"""Post-processing: raw-score regressions, power crossovers, ratio summaries.

Section 8.3 and Appendix B. The raw-score regression is the primary test of the
theoretical penalty coefficients; the calibrated crossovers answer the separate
practical question of required sample length at a common false-positive target,
and their slopes are *not* expected to equal the raw MDL coefficient.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

__all__ = [
    "PREDICTED_SLOPES",
    "predicted_slope",
    "score_regression",
    "score_regression_summary",
    "crossover_estimates",
    "crossover_ratio_summary",
]


def predicted_slope(detector: str, m: int) -> float:
    """Theoretical coefficient of ``log n`` in the complexity increment."""
    if detector == "full":
        return (m - 1) / 2.0
    if detector == "fundamental":
        return (1 if m == 2 else 2) / 2.0
    if detector == "shared_orbit":
        return 0.0
    raise ValueError(f"unknown detector {detector!r}")


#: Convenience mapping used by the tests and the reports.
PREDICTED_SLOPES = {
    detector: {m: predicted_slope(detector, m) for m in range(2, 9)}
    for detector in ("full", "fundamental", "shared_orbit")
}


def score_regression(frame: pd.DataFrame) -> dict:
    """Regress mean score on ``n * gain``, ``log n`` and effect indicators.

    The design is ``mean_score ~ beta * (n * gain) - c * log n + a_effect``.
    Returns ``beta``, its standard error, the empirical penalty slope ``c``
    (defined as minus the ``log n`` coefficient), its standard error, and R^2.
    """
    frame = frame.dropna(subset=["mean_score", "expected_gain_total", "total_length"])
    n_points = len(frame)
    effects = np.sort(frame["effect"].unique())
    if n_points < len(effects) + 3:
        raise ValueError(f"not enough design points ({n_points}) to fit the regression")

    x_gain = frame["expected_gain_total"].to_numpy(dtype=float)
    x_log = np.log(frame["total_length"].to_numpy(dtype=float))
    indicators = np.column_stack([(frame["effect"].to_numpy() == e).astype(float) for e in effects])
    X = np.column_stack([x_gain, x_log, indicators])
    y = frame["mean_score"].to_numpy(dtype=float)

    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ coef
    dof = n_points - np.linalg.matrix_rank(X)
    sigma2 = float(resid @ resid) / dof if dof > 0 else np.nan
    cov = sigma2 * np.linalg.pinv(X.T @ X)
    se = np.sqrt(np.diag(cov))
    ss_tot = float(((y - y.mean()) ** 2).sum())
    r2 = 1.0 - float(resid @ resid) / ss_tot if ss_tot > 0 else np.nan

    return {
        "beta_gain": float(coef[0]),
        "beta_gain_se": float(se[0]),
        "penalty_slope": float(-coef[1]),
        "penalty_slope_se": float(se[1]),
        "r_squared": r2,
        "n_points": int(n_points),
    }


def score_regression_summary(results: pd.DataFrame, scenario_by_detector: dict[str, str] | None = None) -> pd.DataFrame:
    """Run :func:`score_regression` for each detector and group order.

    ``scenario_by_detector`` selects the *matched* scenario for each detector:
    the full detector against the full-space higher-mode change (defined for
    ``m >= 4``), and each constrained detector against the scenario its
    hypothesis describes.
    """
    scenario_by_detector = scenario_by_detector or {
        "full": "higher_mode",
        "fundamental": "independent_fundamental",
        "shared_orbit": "exact_orbit",
    }
    rows = []
    for detector, scenario in scenario_by_detector.items():
        subset = results[(results["detector"] == detector) & (results["scenario"] == scenario)]
        for m, group in subset.groupby("m"):
            try:
                fit = score_regression(group)
            except (ValueError, np.linalg.LinAlgError) as exc:
                rows.append({"detector": detector, "scenario": scenario, "m": int(m), "error": str(exc)})
                continue
            fit.update(
                {
                    "detector": detector,
                    "scenario": scenario,
                    "m": int(m),
                    "predicted_slope": predicted_slope(detector, int(m)),
                }
            )
            rows.append(fit)
    columns = [
        "detector", "scenario", "m", "beta_gain", "beta_gain_se",
        "penalty_slope", "penalty_slope_se", "predicted_slope", "r_squared", "n_points",
    ]
    if not rows:
        return pd.DataFrame(columns=columns)
    frame = pd.DataFrame(rows)
    ordered = [c for c in columns if c in frame.columns] + [c for c in frame.columns if c not in columns]
    return frame[ordered].sort_values(["detector", "m"]).reset_index(drop=True)


def _interpolate_crossover(lengths: np.ndarray, power: np.ndarray, target: float = 0.5):
    """Interpolate the ``target``-power length, linearly in log total length.

    The power sequence is first replaced by its cumulative maximum to suppress
    Monte Carlo reversals. Returns ``(length, status)`` where status is one of
    ``"internal"``, ``"below_grid"`` or ``"above_grid"``.
    """
    order = np.argsort(lengths)
    lengths = np.asarray(lengths, float)[order]
    power = np.maximum.accumulate(np.asarray(power, float)[order])

    if power[0] >= target:
        return float(lengths[0]), "below_grid"
    if power[-1] < target:
        return float("nan"), "above_grid"
    idx = int(np.argmax(power >= target))
    x0, x1 = np.log(lengths[idx - 1]), np.log(lengths[idx])
    y0, y1 = power[idx - 1], power[idx]
    if y1 == y0:
        return float(np.exp(x1)), "internal"
    return float(np.exp(x0 + (target - y0) * (x1 - x0) / (y1 - y0))), "internal"


def crossover_estimates(results: pd.DataFrame, power_column: str = "power_calibrated", target: float = 0.5) -> pd.DataFrame:
    """Estimate the ``target``-power total length for every detector and effect."""
    rows = []
    keys = ["m", "scenario", "effect", "detector"]
    for (m, scenario, effect, detector), group in results.groupby(keys):
        length, status = _interpolate_crossover(
            group["total_length"].to_numpy(), group[power_column].to_numpy(), target
        )
        rows.append(
            {
                "m": int(m),
                "scenario": scenario,
                "effect": float(effect),
                "detector": detector,
                "power_column": power_column,
                "target": target,
                "crossover_length": length,
                "status": status,
            }
        )
    return pd.DataFrame(rows).sort_values(keys).reset_index(drop=True)


def crossover_ratio_summary(
    crossovers: pd.DataFrame,
    scenario: str = "exact_orbit",
    pairs: tuple[tuple[str, str], ...] = (
        ("shared_orbit", "full"),
        ("shared_orbit", "fundamental"),
        ("fundamental", "full"),
    ),
    internal_only: bool = True,
) -> pd.DataFrame:
    """Median crossover-length ratios per group order. A ratio below one favours
    the numerator (fewer observations required)."""
    subset = crossovers[crossovers["scenario"] == scenario]
    if internal_only:
        subset = subset[subset["status"] == "internal"]
    columns = ["m"] + [f"{num}/{den}{suffix}" for num, den in pairs for suffix in ("", "_n")]
    if subset.empty:
        return pd.DataFrame(columns=columns)
    wide = subset.pivot_table(index=["m", "effect"], columns="detector", values="crossover_length")

    rows = []
    for m, group in wide.groupby(level="m"):
        row = {"m": int(m)}
        for num, den in pairs:
            if num in group.columns and den in group.columns:
                ratio = (group[num] / group[den]).dropna()
                row[f"{num}/{den}"] = float(ratio.median()) if len(ratio) else float("nan")
                row[f"{num}/{den}_n"] = int(len(ratio))
            else:
                row[f"{num}/{den}"] = float("nan")
                row[f"{num}/{den}_n"] = 0
        rows.append(row)
    return pd.DataFrame(rows).sort_values("m").reset_index(drop=True)

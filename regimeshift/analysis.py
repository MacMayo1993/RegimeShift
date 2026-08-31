"""Post-processing: raw-score regressions, power crossovers, ratio summaries.

Section 8.3 and Appendix B. The raw-score regression is the primary test of the
theoretical penalty coefficients; the calibrated crossovers answer the separate
practical question of required sample length at a common false-positive target,
and their slopes are *not* expected to equal the raw MDL coefficient.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .detectors import K_STAR, nats_to_bits
from .simulation import BASE_SEED

__all__ = [
    "PREDICTED_SLOPES",
    "K_STAR",
    "UNITS",
    "dimension_increment",
    "predicted_slope",
    "score_regression",
    "gain_residual_regression",
    "score_regression_summary",
    "crossover_estimates",
    "crossover_bootstrap",
    "crossover_ratio_bootstrap",
    "crossover_ratio_summary",
]


UNITS = ("nats", "bits")

#: Slope-valued report columns, i.e. the ones that carry units of
#: codelength-per-e-fold and must be converted together.
_SLOPE_COLUMNS = (
    "penalty_slope", "penalty_slope_se",
    "penalty_slope_wls", "penalty_slope_wls_se",
    "residual_slope", "residual_slope_se",
    "predicted_slope", "weighting_shift",
)


def dimension_increment(detector: str, m: int) -> int:
    """Continuous-dimension increment ``d`` of a detector's alternative."""
    if detector == "full":
        return m - 1
    if detector == "fundamental":
        return 1 if m == 2 else 2
    if detector == "shared_orbit":
        return 0
    raise ValueError(f"unknown detector {detector!r}")


def predicted_slope(detector: str, m: int, units: str = "nats") -> float:
    """Theoretical coefficient of ``log n`` in the complexity increment.

    In nats this is ``d/2``. In bits it is ``d * K_STAR`` exactly, since
    ``(d/2) log2(n) = d/(2 ln 2) * ln n`` -- so every leading coefficient in the
    framework is an integer multiple of ``K_STAR = 1/(2 ln 2)``, and the
    three-way hierarchy is how many of them a model pays: ``m - 1`` for Model A,
    ``d_fund`` for Model B, and zero for Model C.
    """
    if units not in UNITS:
        raise ValueError(f"units must be one of {UNITS}, got {units!r}")
    nats = dimension_increment(detector, m) / 2.0
    return nats if units == "nats" else nats_to_bits(nats)


#: Convenience mapping used by the tests and the reports.
PREDICTED_SLOPES = {
    detector: {m: predicted_slope(detector, m) for m in range(2, 9)}
    for detector in ("full", "fundamental", "shared_orbit")
}


def _design(frame: pd.DataFrame, response: str):
    """Build the shared design matrix: ``n * gain``, ``log n``, effect indicators."""
    frame = frame.dropna(subset=[response, "expected_gain_total", "total_length"])
    n_points = len(frame)
    effects = np.sort(frame["effect"].unique())
    if n_points < len(effects) + 3:
        raise ValueError(f"not enough design points ({n_points}) to fit the regression")

    x_gain = frame["expected_gain_total"].to_numpy(dtype=float)
    x_log = np.log(frame["total_length"].to_numpy(dtype=float))
    indicators = np.column_stack([(frame["effect"].to_numpy() == e).astype(float) for e in effects])
    X = np.column_stack([x_gain, x_log, indicators])
    return frame, X, frame[response].to_numpy(dtype=float)


def _monte_carlo_weights(frame: pd.DataFrame) -> np.ndarray | None:
    """Inverse Monte Carlo variance of each aggregate mean, when recoverable.

    Each design point is a mean over ``n_alt`` trials, so its sampling variance
    is ``sd_score^2 / n_alt`` -- known, and generally *unequal* across the grid.
    Unweighted least squares on the aggregate means therefore misstates the
    standard errors and can bias the fitted slope. Returns ``None`` when the
    frame lacks the columns (e.g. synthetic fixtures), so callers fall back to
    ordinary least squares.
    """
    if not {"sd_score", "n_alt"}.issubset(frame.columns):
        return None
    variance = frame["sd_score"].to_numpy(dtype=float) ** 2 / frame["n_alt"].to_numpy(dtype=float)
    if not np.all(np.isfinite(variance)) or np.any(variance <= 0):
        return None
    return 1.0 / variance


def _fit(X: np.ndarray, y: np.ndarray, weights: np.ndarray | None) -> dict:
    """Least squares, optionally weighted, with standard errors and R^2."""
    if weights is None:
        Xw, yw = X, y
    else:
        root = np.sqrt(weights)[:, None]
        Xw, yw = X * root, y * root.ravel()

    coef, *_ = np.linalg.lstsq(Xw, yw, rcond=None)
    resid = yw - Xw @ coef
    dof = len(y) - np.linalg.matrix_rank(Xw)
    # Scaling by the residual mean square rather than trusting the weights'
    # absolute scale keeps the standard errors honest if the variance model is
    # only approximately right.
    sigma2 = float(resid @ resid) / dof if dof > 0 else np.nan
    se = np.sqrt(np.diag(sigma2 * np.linalg.pinv(Xw.T @ Xw)))

    raw_resid = y - X @ coef
    ss_tot = float(((y - y.mean()) ** 2).sum())
    r2 = 1.0 - float(raw_resid @ raw_resid) / ss_tot if ss_tot > 0 else np.nan
    return {"coef": coef, "se": se, "r_squared": r2, "condition_number": float(np.linalg.cond(Xw))}


def score_regression(frame: pd.DataFrame, weighted: bool = False) -> dict:
    """Regress mean score on ``n * gain``, ``log n`` and effect indicators.

    The design is ``mean_score ~ beta * (n * gain) - c * log n + a_effect``.
    Returns ``beta``, its standard error, the empirical penalty slope ``c``
    (defined as minus the ``log n`` coefficient), its standard error, R^2, and
    the design condition number.

    Set ``weighted=True`` for Monte Carlo-variance-weighted least squares. Each
    group-level regression has only ``n_effects * n_lengths`` aggregate design
    points -- far fewer than the number of simulated datasets suggests -- and
    those points have unequal precision, so the unweighted fit is not
    automatically the right one. Reporting both is the check that the estimated
    slope is not an artifact of the specification.
    """
    frame, X, y = _design(frame, "mean_score")
    weights = _monte_carlo_weights(frame) if weighted else None
    fit = _fit(X, y, weights)
    return {
        "beta_gain": float(fit["coef"][0]),
        "beta_gain_se": float(fit["se"][0]),
        "penalty_slope": float(-fit["coef"][1]),
        "penalty_slope_se": float(fit["se"][1]),
        "r_squared": fit["r_squared"],
        "condition_number": fit["condition_number"],
        "method": "wls" if weights is not None else "ols",
        "n_points": int(len(y)),
    }


def gain_residual_regression(frame: pd.DataFrame, weighted: bool = False) -> dict:
    """Regress ``mean_raw_gain - n * gain`` on ``log n``: the direct diagnostic.

    The penalty this implementation subtracts is computed exactly and is
    deterministic, so estimating the gain coefficient and the penalty
    coefficient simultaneously is avoidable. Writing the fitted mean raw gain as
    ``n * G + a + s * log n`` and the exact penalty as ``(d/2) log n + c``, the
    score regression's penalty slope satisfies the identity

        penalty_slope = d/2 - s

    where ``s`` -- returned here as ``residual_slope`` -- is the *only*
    empirical quantity in it. Theory says the raw gain is ``n * G + O(1)``, so
    ``s`` should be near zero for all three detectors. This isolates the
    departure of the likelihood gain from its population value, which for the
    shared-orbit detector (``d = 0``) is the entire source of the residual
    log-length slope reported in Table 5.
    """
    if "mean_raw_gain" not in frame.columns:
        raise ValueError(
            "the gain-residual diagnostic needs a 'mean_raw_gain' column; "
            "results written before this diagnostic existed do not carry one"
        )
    frame = frame.copy()
    frame["gain_residual"] = frame["mean_raw_gain"] - frame["expected_gain_total"]
    frame, X, y = _design(frame, "gain_residual")
    weights = _monte_carlo_weights(frame) if weighted else None
    fit = _fit(X, y, weights)
    return {
        "residual_slope": float(fit["coef"][1]),
        "residual_slope_se": float(fit["se"][1]),
        "r_squared": fit["r_squared"],
        "condition_number": fit["condition_number"],
        "method": "wls" if weights is not None else "ols",
        "n_points": int(len(y)),
    }


def score_regression_summary(
    results: pd.DataFrame,
    scenario_by_detector: dict[str, str] | None = None,
    units: str = "nats",
) -> pd.DataFrame:
    """Run :func:`score_regression` for each detector and group order.

    ``scenario_by_detector`` selects the *matched* scenario for each detector:
    the full detector against the full-space higher-mode change (defined for
    ``m >= 4``), and each constrained detector against the scenario its
    hypothesis describes.

    ``units`` converts the slope-valued columns. The fits themselves are done in
    nats; ``"bits"`` divides them by ``ln 2``, which makes each predicted slope
    an exact multiple of ``K_STAR`` and adds a ``k_star_multiple`` column
    reporting that multiple (which is just ``d``).
    """
    if units not in UNITS:
        raise ValueError(f"units must be one of {UNITS}, got {units!r}")
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
                weighted = score_regression(group, weighted=True)
            except (ValueError, np.linalg.LinAlgError) as exc:
                rows.append({"detector": detector, "scenario": scenario, "m": int(m), "error": str(exc)})
                continue
            fit.update(
                {
                    "detector": detector,
                    "scenario": scenario,
                    "m": int(m),
                    "predicted_slope": predicted_slope(detector, int(m)),
                    "penalty_slope_wls": weighted["penalty_slope"],
                    "penalty_slope_wls_se": weighted["penalty_slope_se"],
                    "beta_gain_wls": weighted["beta_gain"],
                    "weighting_shift": abs(weighted["penalty_slope"] - fit["penalty_slope"]),
                }
            )
            # The residual diagnostic needs a column that older results files do
            # not carry; its absence must not take the whole report down.
            try:
                residual = gain_residual_regression(group)
                fit["residual_slope"] = residual["residual_slope"]
                fit["residual_slope_se"] = residual["residual_slope_se"]
            except (ValueError, np.linalg.LinAlgError):
                fit["residual_slope"] = float("nan")
                fit["residual_slope_se"] = float("nan")
            rows.append(fit)
    columns = [
        "detector", "scenario", "m", "beta_gain", "beta_gain_se",
        "penalty_slope", "penalty_slope_se", "penalty_slope_wls", "penalty_slope_wls_se",
        "weighting_shift", "residual_slope", "residual_slope_se",
        "predicted_slope", "k_star_multiple", "units",
        "r_squared", "condition_number", "n_points",
    ]
    if not rows:
        return pd.DataFrame(columns=columns)
    frame = pd.DataFrame(rows)
    if "detector" in frame.columns and "m" in frame.columns:
        frame["k_star_multiple"] = [
            dimension_increment(d, int(m)) for d, m in zip(frame["detector"], frame["m"])
        ]
    frame["units"] = units
    if units == "bits":
        for column in _SLOPE_COLUMNS:
            if column in frame.columns:
                frame[column] = frame[column].astype(float).map(nats_to_bits)
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


def crossover_bootstrap(
    results: pd.DataFrame,
    n_boot: int = 500,
    power_column: str = "power_calibrated",
    target: float = 0.5,
    ci: float = 0.95,
    seed: int = BASE_SEED,
) -> pd.DataFrame:
    """Bootstrap confidence intervals for the 50%-power crossover lengths.

    The point estimate passes through a cumulative-maximum stabilisation, an
    interpolation over a handful of lengths, and (for the ratios) a median
    across effects. None of that carries an uncertainty statement on its own, so
    the whole pipeline is resampled: each design point's power is redrawn as
    ``Binomial(n_alt, power_hat) / n_alt`` and the crossover recomputed.

    Two limitations are deliberate and worth stating rather than hiding. The
    resampling captures binomial power noise only -- not the extra variability
    of the empirical 95th-percentile critical value, which is estimated from
    ``n_null`` draws and has roughly ``alpha * n_null`` observations in the
    relevant tail. And detectors are resampled independently although they score
    the *same* simulated datasets and are positively correlated, which makes the
    ratio intervals conservative (too wide) rather than optimistic.

    **How wide, measurably.** The second limitation is not a mild conservatism,
    and the grid contains its own calibration. At ``m = 2`` and ``m = 3`` the
    fundamental component spans the whole nontrivial tangent space, so ``full``
    and ``fundamental`` are the same detector scoring the same datasets: in the
    committed run their calibrated power agrees *exactly* at every design point.
    The ``fundamental/full`` ratio is therefore identically 1 with zero
    variance -- and :func:`crossover_ratio_bootstrap`, which resamples the same
    way, reports [0.916, 1.098] at ``m = 2`` and [0.926, 1.094] at ``m = 3``,
    because it draws two independent binomials from one shared power curve.
    A second symptom in the same row: at ``m = 2``
    ``shared_orbit/full`` gets [0.695, 0.814] while ``shared_orbit/fundamental``
    gets [0.701, 0.821], two intervals for what is mathematically one number.

    Read that +-9% as the measured inflation on the rows that carry real
    content, and do not read the ``m = 2, 3`` rows as a passed check. Removing
    it needs a *joint* resample of per-dataset detector outcomes, which
    :mod:`regimeshift.simulation` does not currently retain; that is the natural
    next revision of this analysis, and it is a change to the runner rather than
    to this function. ``tests/test_committed_results.py`` pins both the
    degeneracy and the width of the interval it produces, so the artefact
    cannot be mistaken for a result.
    """
    rng = np.random.default_rng(seed)
    lo_q, hi_q = (1.0 - ci) / 2.0, 1.0 - (1.0 - ci) / 2.0
    rows = []
    for (m, scenario, effect, detector), group in results.groupby(
        ["m", "scenario", "effect", "detector"]
    ):
        lengths = group["total_length"].to_numpy(dtype=float)
        power = group[power_column].to_numpy(dtype=float)
        trials = group["n_alt"].to_numpy(dtype=int)
        point, status = _interpolate_crossover(lengths, power, target)

        draws = []
        for _ in range(n_boot):
            resampled = rng.binomial(trials, np.clip(power, 0.0, 1.0)) / trials
            value, boot_status = _interpolate_crossover(lengths, resampled, target)
            if boot_status == "internal" and np.isfinite(value):
                draws.append(value)

        rows.append(
            {
                "m": int(m),
                "scenario": scenario,
                "effect": float(effect),
                "detector": detector,
                "crossover_length": point,
                "status": status,
                "ci_low": float(np.quantile(draws, lo_q)) if draws else float("nan"),
                "ci_high": float(np.quantile(draws, hi_q)) if draws else float("nan"),
                "boot_internal_fraction": len(draws) / n_boot if n_boot else float("nan"),
                "n_boot": n_boot,
            }
        )
    return pd.DataFrame(rows).sort_values(["m", "scenario", "effect", "detector"]).reset_index(drop=True)


def crossover_ratio_bootstrap(
    results: pd.DataFrame,
    scenario: str = "exact_orbit",
    pairs: tuple[tuple[str, str], ...] = (
        ("shared_orbit", "full"),
        ("shared_orbit", "fundamental"),
        ("fundamental", "full"),
    ),
    n_boot: int = 500,
    power_column: str = "power_calibrated",
    target: float = 0.5,
    ci: float = 0.95,
    seed: int = BASE_SEED,
) -> pd.DataFrame:
    """Bootstrap intervals for the median crossover-length *ratios* per group.

    The full pipeline -- resample power, interpolate crossovers, take the median
    ratio across effects -- is repeated, so the interval covers interpolation
    and median-across-effects variability as well as binomial noise. See
    :func:`crossover_bootstrap` for what it still does not cover, including the
    independent-resampling artefact that its ``m = 2, 3`` rows measure directly.

    Two further properties of the *point* estimates belong with any reading of
    them, and :func:`crossover_ratio_summary` reports the counts that expose
    both.

    A ratio is formed only at effects where *both* detectors cross inside the
    grid, so the median can rest on very few points -- two of the four effect
    levels at ``m = 2`` and ``m = 3``, where a "median across effects" is just
    the midpoint of a pair. And because each column keeps its own surviving
    subset, the columns are not mutually consistent: at ``m = 5`` the reported
    ``shared_orbit/full`` is 0.630 while ``(shared_orbit/fundamental) *
    (fundamental/full)`` is 0.624.

    The same filter runs inside the loop below, so a replicate whose crossover
    falls out of grid drops that effect from *its* median. The bootstrap
    distribution therefore mixes medians taken over different effect subsets
    rather than resampling one fixed estimator. Fixing that means freezing the
    subset to the one the point estimate uses and discarding replicates that
    cannot fill it -- which changes the committed intervals, so it is left as a
    stated defect rather than applied silently here.
    """
    rng = np.random.default_rng(seed)
    lo_q, hi_q = (1.0 - ci) / 2.0, 1.0 - (1.0 - ci) / 2.0
    subset = results[results["scenario"] == scenario]
    rows = []

    for m, group in subset.groupby("m"):
        effects = np.sort(group["effect"].unique())
        curves = {
            (float(effect), detector): sub.sort_values("total_length")
            for (effect, detector), sub in group.groupby(["effect", "detector"])
        }
        boot = {pair: [] for pair in pairs}
        for _ in range(n_boot):
            crossings: dict[tuple[float, str], float] = {}
            for key, sub in curves.items():
                lengths = sub["total_length"].to_numpy(dtype=float)
                power = np.clip(sub[power_column].to_numpy(dtype=float), 0.0, 1.0)
                trials = sub["n_alt"].to_numpy(dtype=int)
                value, status = _interpolate_crossover(
                    lengths, rng.binomial(trials, power) / trials, target
                )
                if status == "internal" and np.isfinite(value):
                    crossings[key] = value
            for num, den in pairs:
                ratios = [
                    crossings[(float(e), num)] / crossings[(float(e), den)]
                    for e in effects
                    if (float(e), num) in crossings and (float(e), den) in crossings
                ]
                if ratios:
                    boot[(num, den)].append(float(np.median(ratios)))

        row = {"m": int(m), "scenario": scenario, "n_boot": n_boot}
        for pair in pairs:
            draws = boot[pair]
            label = f"{pair[0]}/{pair[1]}"
            row[f"{label}_ci_low"] = float(np.quantile(draws, lo_q)) if draws else float("nan")
            row[f"{label}_ci_high"] = float(np.quantile(draws, hi_q)) if draws else float("nan")
            row[f"{label}_boot_n"] = len(draws)
        rows.append(row)
    return pd.DataFrame(rows).sort_values("m").reset_index(drop=True)


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

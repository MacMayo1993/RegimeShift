"""Command line entry point: ``python -m regimeshift``.

Examples
--------
Quick smoke run (seconds)::

    python -m regimeshift run --grid quick --out results/quick

Full production grid (hours; parallelise it)::

    python -m regimeshift run --grid production --workers 16 --out results/v3

Re-analyse an existing results file::

    python -m regimeshift analyse --results results/v3/full_results.csv --out results/v3
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from .analysis import crossover_estimates, crossover_ratio_summary, score_regression_summary
from .runner import PRODUCTION_GRID, QUICK_GRID, run_grid
from .simulation import BASE_SEED, build_grid

GRIDS = {"production": PRODUCTION_GRID, "quick": QUICK_GRID}


def _write_reports(results: pd.DataFrame, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    results.to_csv(out_dir / "full_results.csv", index=False)

    regressions = score_regression_summary(results)
    regressions.to_csv(out_dir / "score_regression_summary.csv", index=False)

    crossovers = crossover_estimates(results)
    crossovers.to_csv(out_dir / "crossover_estimates.csv", index=False)

    scenarios = [s for s in results["scenario"].unique()]
    ratios = pd.concat(
        [crossover_ratio_summary(crossovers, scenario=s).assign(scenario=s) for s in scenarios],
        ignore_index=True,
    )
    ratios.to_csv(out_dir / "crossover_ratio_summary.csv", index=False)

    print(f"\nwrote 4 files to {out_dir}\n")
    print("Raw-score regressions (penalty slope vs prediction):")
    cols = [c for c in ["detector", "m", "beta_gain", "penalty_slope", "penalty_slope_se", "predicted_slope", "r_squared"] if c in regressions]
    print(regressions[cols].to_string(index=False))
    print("\nMedian calibrated crossover ratios:")
    print(ratios.to_string(index=False))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="regimeshift", description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    run = sub.add_parser("run", help="run a Monte Carlo grid and write reports")
    run.add_argument("--grid", choices=sorted(GRIDS), default="quick")
    run.add_argument("--out", default="results/quick", type=Path)
    run.add_argument("--workers", type=int, default=1)
    run.add_argument("--base-seed", type=int, default=BASE_SEED)
    run.add_argument("--checkpoint", type=Path, default=None, help="defaults to <out>/checkpoint.csv")
    run.add_argument("--no-progress", action="store_true")

    analyse = sub.add_parser("analyse", help="re-run the reports on an existing results file")
    analyse.add_argument("--results", required=True, type=Path)
    analyse.add_argument("--out", required=True, type=Path)

    args = parser.parse_args(argv)

    if args.command == "analyse":
        _write_reports(pd.read_csv(args.results), args.out)
        return 0

    spec = GRIDS[args.grid]
    configs = build_grid(
        spec["groups"], spec["scenarios"], spec["effects"], spec["segment_lengths"],
        n_alt=spec["n_alt"], n_null=spec["n_null"],
    )
    datasets = sum(c.n_alt + c.n_null for c in configs)
    print(f"grid={args.grid}: {len(configs)} configurations, "
          f"{3 * len(configs)} detector rows, {datasets:,} simulated datasets")

    checkpoint = args.checkpoint or (args.out / "checkpoint.csv")
    results = run_grid(
        configs, checkpoint=checkpoint, workers=args.workers,
        base_seed=args.base_seed, progress=not args.no_progress,
    )
    _write_reports(results, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

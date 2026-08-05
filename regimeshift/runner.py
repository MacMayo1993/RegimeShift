"""Deterministic parallel production runner with per-configuration checkpointing.

Configurations are independent and carry their own seeds, so results do not
depend on worker count or completion order. Completed configurations are
appended to the checkpoint file and skipped on resume.
"""

from __future__ import annotations

import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Iterable, Sequence

import pandas as pd

from .simulation import BASE_SEED, Config, run_config

__all__ = ["PRODUCTION_GRID", "QUICK_GRID", "run_grid", "load_checkpoint"]

#: The design reported in the manuscript (Table 2): 312 configurations,
#: 936 detector rows, 468,000 simulated two-segment datasets.
PRODUCTION_GRID = {
    "groups": (2, 3, 4, 5, 6),
    "scenarios": ("exact_orbit", "independent_fundamental", "higher_mode"),
    "effects": (0.08, 0.12, 0.18, 0.25),
    "segment_lengths": (100, 200, 400, 800, 1600, 3200),
    "n_alt": 500,
    "n_null": 1000,
}

#: A small grid for CI and smoke runs: same structure, far fewer draws.
QUICK_GRID = {
    "groups": (2, 4),
    "scenarios": ("exact_orbit", "independent_fundamental", "higher_mode"),
    "effects": (0.08, 0.12, 0.18),
    "segment_lengths": (100, 200, 400, 800),
    "n_alt": 100,
    "n_null": 200,
}

_COLUMN_KEYS = ["m", "scenario", "effect", "segment_length", "split_fraction", "detector"]


def load_checkpoint(path: str | os.PathLike | None) -> pd.DataFrame:
    """Load previously completed rows, or an empty frame."""
    if path is None or not Path(path).exists():
        return pd.DataFrame()
    return pd.read_csv(path)


def _pending(configs: Sequence[Config], done: pd.DataFrame) -> list[Config]:
    if done.empty:
        return list(configs)
    # Checkpoints written before split_fraction existed are all balanced.
    fractions = (
        done["split_fraction"].astype(float)
        if "split_fraction" in done.columns
        else pd.Series(0.5, index=done.index)
    )
    seen = set(
        zip(done["m"], done["scenario"], done["effect"].astype(float),
            done["segment_length"], fractions)
    )
    return [c for c in configs if c.key not in seen]


def run_grid(
    configs: Iterable[Config],
    checkpoint: str | os.PathLike | None = None,
    workers: int = 1,
    base_seed: int = BASE_SEED,
    progress: bool = False,
) -> pd.DataFrame:
    """Run every configuration and return the combined detector-level frame."""
    configs = list(configs)
    done = load_checkpoint(checkpoint)
    pending = _pending(configs, done)
    frames = [done] if not done.empty else []

    def _emit(rows):
        frame = pd.DataFrame(rows)
        frames.append(frame)
        if checkpoint is not None:
            path = Path(checkpoint)
            path.parent.mkdir(parents=True, exist_ok=True)
            header = not path.exists()
            frame.to_csv(path, mode="a", header=header, index=False)

    if workers <= 1:
        for i, config in enumerate(pending, 1):
            _emit(run_config(config, base_seed))
            if progress:
                print(f"[{i}/{len(pending)}] {config.key}", flush=True)
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = {pool.submit(run_config, c, base_seed): c for c in pending}
            for i, future in enumerate(as_completed(futures), 1):
                _emit(future.result())
                if progress:
                    print(f"[{i}/{len(pending)}] {futures[future].key}", flush=True)

    if not frames:
        return pd.DataFrame()
    combined = pd.concat(frames, ignore_index=True)
    return combined.drop_duplicates(subset=_COLUMN_KEYS, keep="last").reset_index(drop=True)

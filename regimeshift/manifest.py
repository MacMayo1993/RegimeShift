"""Provenance capture for a results directory.

Both external reviews asked for the same thing: a permanent record tying the
numbers to an exact commit, an environment, and checksums. Every run written by
:mod:`regimeshift.cli` drops a ``run_manifest.json`` beside its CSVs so the
results carry that record with them rather than depending on the surrounding
repository state.
"""

from __future__ import annotations

import hashlib
import json
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

__all__ = ["file_digest", "git_commit", "environment", "write_manifest"]


def file_digest(path: Path) -> str:
    """SHA-256 of a file, as a hex string."""
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git_commit() -> dict:
    """The repository's current commit and cleanliness, if git is available."""
    def run(*args):
        return subprocess.run(
            args, cwd=Path(__file__).resolve().parent, capture_output=True, text=True, timeout=10
        )

    try:
        head = run("git", "rev-parse", "HEAD")
        if head.returncode != 0:
            return {"commit": None, "note": "not a git checkout"}
        status = run("git", "status", "--porcelain")
        return {
            "commit": head.stdout.strip(),
            "dirty": bool(status.stdout.strip()),
        }
    except (OSError, subprocess.SubprocessError):
        return {"commit": None, "note": "git unavailable"}


def environment() -> dict:
    """Python, platform and dependency versions."""
    versions = {}
    for name in ("numpy", "scipy", "pandas"):
        try:
            versions[name] = __import__(name).__version__
        except Exception:  # pragma: no cover - a missing dependency cannot reach here
            versions[name] = None
    import regimeshift

    return {
        "regimeshift": regimeshift.__version__,
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "machine": platform.machine(),
        "packages": versions,
    }


def write_manifest(out_dir: Path, *, grid: str, spec: dict, base_seed: int,
                   n_configs: int, n_datasets: int, workers: int,
                   elapsed_seconds: float | None = None, units: str = "nats",
                   n_boot: int = 0) -> Path:
    """Write ``run_manifest.json`` describing how a results directory was made."""
    out_dir = Path(out_dir)
    files = sorted(
        p for p in out_dir.glob("*.csv") if p.name != "checkpoint.csv"
    )
    manifest = {
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "grid": grid,
        "design": {k: list(v) if isinstance(v, tuple) else v for k, v in spec.items()},
        "base_seed": base_seed,
        "configurations": n_configs,
        "detector_rows": 3 * n_configs,
        "simulated_datasets": n_datasets,
        "workers": workers,
        "elapsed_seconds": None if elapsed_seconds is None else round(elapsed_seconds, 1),
        "report_units": units,
        "bootstrap_replicates": n_boot,
        "git": git_commit(),
        "environment": environment(),
        "files": {p.name: {"sha256": file_digest(p), "bytes": p.stat().st_size} for p in files},
    }
    path = out_dir / "run_manifest.json"
    path.write_text(json.dumps(manifest, indent=2) + "\n")
    return path

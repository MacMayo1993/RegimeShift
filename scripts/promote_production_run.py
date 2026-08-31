"""Verify a fresh production run against the committed one, then promote it.

The committed run in ``results/v3-production/`` is what the manuscript quotes,
so replacing it is not something to do casually. This script refuses to unless
the new run reproduces every column the old one already had, *bit for bit*.

The point of a re-run is the joint detection-pattern columns
(:data:`~regimeshift.simulation.DETECTION_PATTERNS`), which let the
crossover-ratio bootstrap resample detectors together instead of independently.
Those are new; nothing else may move. If a shared column differs -- because the
environment drifted, or a dependency changed a numerical path -- the right
outcome is a report and a stop, not a silent substitution of every number in
Section 9.

Usage::

    python scripts/promote_production_run.py <new-run-dir> [--promote]

Without ``--promote`` it only reports. The environment the committed run used is
recorded in its ``run_manifest.json``; match it before running.
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

import pandas as pd

from regimeshift.simulation import DETECTION_PATTERNS

COMMITTED = Path(__file__).resolve().parent.parent / "results" / "v3-production"
KEY = ["m", "scenario", "effect", "segment_length", "split_fraction", "detector"]


def compare(old: pd.DataFrame, new: pd.DataFrame) -> tuple[bool, list[str]]:
    """Do the columns the committed run already had reproduce exactly?"""
    notes: list[str] = []

    added = [c for c in new.columns if c not in old.columns]
    removed = [c for c in old.columns if c not in new.columns]
    if removed:
        notes.append(f"columns LOST: {removed}")
    expected_new = set(DETECTION_PATTERNS)
    unexpected = sorted(set(added) - expected_new)
    if unexpected:
        notes.append(f"unexpected new columns: {unexpected}")
    notes.append(f"columns added: {sorted(added)}")

    if len(old) != len(new):
        notes.append(f"row count differs: {len(old)} -> {len(new)}")
        return False, notes

    old_sorted = old.sort_values(KEY).reset_index(drop=True)
    new_sorted = new.sort_values(KEY).reset_index(drop=True)

    identical = True
    for column in old.columns:
        if column not in new.columns:
            identical = False
            continue
        a, b = old_sorted[column], new_sorted[column]
        if a.dtype.kind in "fc" or b.dtype.kind in "fc":
            same = ((a == b) | (a.isna() & b.isna())).all()
        else:
            same = a.equals(b)
        if not same:
            identical = False
            diff = (a != b) & ~(a.isna() & b.isna())
            worst = ""
            if a.dtype.kind == "f":
                worst = f", max |delta| = {(a - b).abs().max():.3e}"
            notes.append(f"  DIFFERS: {column} ({int(diff.sum())} of {len(a)} rows{worst})")
    notes.append("all shared columns identical" if identical else "SHARED COLUMNS DIFFER")
    return identical, notes


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("new_run", type=Path, help="directory holding the fresh run")
    parser.add_argument("--promote", action="store_true",
                        help="copy the new run over results/v3-production/ if it verifies")
    args = parser.parse_args()

    old = pd.read_csv(COMMITTED / "full_results.csv")
    new = pd.read_csv(args.new_run / "full_results.csv")

    identical, notes = compare(old, new)
    for note in notes:
        print(note)

    if not identical:
        print("\nREFUSING to promote: the re-run does not reproduce the committed columns.")
        print("Report the difference rather than replacing the manuscript's numbers.")
        return 1

    print("\nVerified: the re-run reproduces every committed column exactly.")
    if not args.promote:
        print("Re-run with --promote to install it.")
        return 0

    for path in sorted(args.new_run.glob("*")):
        if path.name == "checkpoint.csv":
            continue
        shutil.copy2(path, COMMITTED / path.name)
        print(f"  installed {path.name}")
    print("\nPromoted. Re-check the manuscript tables that quote bootstrap intervals.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

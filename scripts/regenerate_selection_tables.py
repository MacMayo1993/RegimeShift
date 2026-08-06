"""Regenerate the Section 10 and Section 11 tables of the v4 manuscript.

Both tables concern *model selection* rather than detector scores, so neither is
produced by the main production grid. This script is their provenance: it is
deterministic given the seeds below, and the manuscript quotes its output.

    python scripts/regenerate_selection_tables.py

Section 10 reports how often the shortest code is the generating family.
Section 11 reports description-length savings against one common reference,
    L(null_full) - L(M),
which is comparable across model classes in a way detector scores are not --
see the manuscript's Section 10 for why, and for the O(1) convention these
lengths are stated under.
"""

from __future__ import annotations

import numpy as np

from regimeshift.fourier import fundamental_dimension, probabilities, rotation_matrix
from regimeshift.scenarios import build_segments
from regimeshift.selection import code_lengths, select_model

SEED = 20260713
GROUP = 6
EFFECT = 0.25
TAU = 0.05

MODELS = ("full", "fundamental", "shared_orbit", "approximate_orbit")
LABEL = {
    "full": "A",
    "fundamental": "B",
    "shared_orbit": "C",
    "approximate_orbit": "D",
}


def section_10(trials: int = 200) -> None:
    """Recovery of the generating family by shortest total code length."""
    lengths = (200, 800, 3200)
    print(f"\n## Section 10 -- selection without an oracle "
          f"(g={GROUP}, effect={EFFECT}, {trials} trials)\n")
    print("| generated from | " + " | ".join(f"{n}/side" for n in lengths) + " |")
    print("|---|" + "---:|" * len(lengths))

    cases = [
        ("exact orbit -> recovers shared orbit", "exact_orbit", "shared_orbit"),
        ("higher mode -> recovers full", "higher_mode", "full"),
        ("no change -> false-change rate", None, None),
    ]
    for caption, scenario, target in cases:
        rates = []
        for n_side in lengths:
            rng = np.random.default_rng(SEED)
            if scenario is None:
                p_left = p_right = probabilities(
                    np.array([EFFECT, 0.0])[: fundamental_dimension(GROUP)], GROUP
                )
            else:
                segments = build_segments(GROUP, scenario, EFFECT)
                p_left, p_right = segments.p_left, segments.p_right

            hits = 0
            for _ in range(trials):
                left = rng.multinomial(n_side, p_left).astype(float)
                right = rng.multinomial(n_side, p_right).astype(float)
                chosen = select_model(left, right, GROUP, deviation_scale=TAU)
                hits += (
                    chosen.declared_change
                    if target is None
                    else chosen.selected == target
                )
            rates.append(100.0 * hits / trials)
        print(f"| {caption} | " + " | ".join(f"{r:.0f}%" for r in rates) + " |")


def section_11(trials: int = 250, n_side: int = 1200) -> None:
    """Deviation sweep, on a common reference, with paired uncertainty."""
    d = fundamental_dimension(GROUP)
    rotation = rotation_matrix(GROUP, 1)
    theta_left = np.array([EFFECT, 0.0])[:d]
    sweep = (0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 3.0)

    print(f"\n## Section 11 -- deviation sweep "
          f"(g={GROUP}, effect={EFFECT}, {n_side}/side, tau={TAU}, {trials} trials)")
    print("   reported as L(null_full) - L(M), in nats\n")
    print("| deviation | A full | B fundamental | C exact orbit | D approximate | best |")
    print("|---:|---:|---:|---:|---:|---|")

    paired = {}
    for deviation in sweep:
        rng = np.random.default_rng(SEED + int(deviation * 1000))
        saving = {name: [] for name in MODELS}
        advantage = []
        for _ in range(trials):
            delta = rng.normal(0.0, deviation * EFFECT, size=d)
            theta_right = rotation @ theta_left + delta
            left = rng.multinomial(n_side, probabilities(theta_left, GROUP)).astype(float)
            right = rng.multinomial(n_side, probabilities(theta_right, GROUP)).astype(float)

            lengths = code_lengths(left, right, GROUP, deviation_scale=TAU)
            for name in MODELS:
                saving[name].append(lengths["null_full"] - lengths[name])
            rival = min(lengths[n] for n in MODELS if n != "approximate_orbit")
            advantage.append(rival - lengths["approximate_orbit"])

        means = {name: float(np.mean(saving[name])) for name in MODELS}
        best = max(means, key=means.__getitem__)
        row = " | ".join(f"{means[name]:.2f}" for name in MODELS)
        print(f"| {deviation:.2f} | {row} | {LABEL[best]} |")
        paired[deviation] = np.asarray(advantage)

    print("\n| deviation | mean paired advantage | s.e. | datasets where D is shortest |")
    print("|---:|---:|---:|---:|")
    for deviation in sweep[:5]:
        a = paired[deviation]
        se = a.std(ddof=1) / np.sqrt(a.size)
        print(f"| {deviation:.2f} | {a.mean():+.2f} | {se:.2f} | {100 * (a > 0).mean():.0f}% |")


if __name__ == "__main__":
    section_10()
    section_11()

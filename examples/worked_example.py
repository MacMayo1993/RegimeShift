"""A worked example: one interpretable picture of what the three detectors do.

Runs the three detectors over growing sample sizes on two datasets that differ
only in whether the cyclic structure is real:

* **matched** -- an exact-orbit change, the hypothesis Model C encodes;
* **misspecified** -- a higher-mode change that leaves the fundamental subspace.

and writes ``docs/figures/worked-example.svg``.

All three detectors are calibrated to a common 5% false-positive rate first, so
the panels compare detection ability rather than threshold generosity. The
point of the figure is the contrast between them: when the cyclic structure is
real the constrained detectors need markedly less data, and when it is wrong
that advantage narrows or disappears.

Run it with::

    python examples/worked_example.py

It writes SVG directly rather than pulling in a plotting library, so the
package keeps its NumPy/SciPy/pandas-only dependency set and the output is a
diffable text file.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from regimeshift import build_segments, run_all_detectors

DETECTORS = ("full", "fundamental", "shared_orbit")
COLOURS = {"full": "#ef4444", "fundamental": "#3b82f6", "shared_orbit": "#10b981"}
LABELS = {"full": "A. full", "fundamental": "B. fundamental", "shared_orbit": "C. shared orbit"}

M = 6
EFFECT = 0.25
LENGTHS = (25, 50, 100, 200, 400, 800, 1600)
TRIALS = 300
NULL_TRIALS = 300
SEED = 20260713


def detection_rates(scenario: str, seed: int = SEED) -> dict[str, list[float]]:
    """Calibrated detection rate for each detector, by segment length.

    Every detector is first calibrated to the *same* 5% false-positive rate: at
    each length, null samples with no change are scored, and that detector's
    empirical 95th percentile becomes its critical value. Detection rates on the
    alternative are then measured against those thresholds.

    The calibration is not a formality. The raw ``score > 0`` rule is not
    comparable across these three detectors: Model C's penalty is constant in
    ``n`` (and zero at ``m = 2``), so it fires far more readily under the null
    and looks strongest even where its hypothesis is false. Building this figure
    on the raw rule puts the constrained detectors ahead in *both* panels --
    which is a fact about the thresholds, not about the geometry.
    """
    segments = build_segments(M, scenario, EFFECT)
    rng = np.random.default_rng(seed)
    rates = {name: [] for name in DETECTORS}

    for length in LENGTHS:
        null_scores = {name: np.empty(NULL_TRIALS) for name in DETECTORS}
        for i in range(NULL_TRIALS):
            left = rng.multinomial(length, segments.p_null)
            right = rng.multinomial(length, segments.p_null)
            for name, result in run_all_detectors(left, right, M).items():
                null_scores[name][i] = result.score
        critical = {
            name: float(np.quantile(scores, 0.95, method="higher"))
            for name, scores in null_scores.items()
        }

        hits = {name: 0 for name in DETECTORS}
        for _ in range(TRIALS):
            left = rng.multinomial(length, segments.p_left)
            right = rng.multinomial(length, segments.p_right)
            for name, result in run_all_detectors(left, right, M).items():
                hits[name] += result.score > critical[name]
        for name in DETECTORS:
            rates[name].append(hits[name] / TRIALS)
    return rates


def _panel(x0: int, y0: int, width: int, height: int, title: str, subtitle: str,
           rates: dict[str, list[float]]) -> str:
    """One axes-with-curves panel of the figure."""
    logs = np.log(np.array(LENGTHS, dtype=float) * 2)
    lo, hi = logs.min(), logs.max()

    def sx(log_n: float) -> float:
        return x0 + (log_n - lo) / (hi - lo) * width

    def sy(rate: float) -> float:
        return y0 + height - rate * height

    parts = [
        f'<text class="hd" x="{x0 + width / 2:.0f}" y="{y0 - 34}" text-anchor="middle">{title}</text>',
        f'<text class="sub" x="{x0 + width / 2:.0f}" y="{y0 - 16}" text-anchor="middle">{subtitle}</text>',
        f'<line class="axis" x1="{x0}" y1="{y0 + height}" x2="{x0 + width}" y2="{y0 + height}"/>',
        f'<line class="axis" x1="{x0}" y1="{y0}" x2="{x0}" y2="{y0 + height}"/>',
        f'<line class="grid" x1="{x0}" y1="{sy(0.5):.1f}" x2="{x0 + width}" y2="{sy(0.5):.1f}"/>',
        f'<text class="tick" x="{x0 - 8}" y="{sy(0.5) + 4:.1f}" text-anchor="end">50%</text>',
        f'<text class="tick" x="{x0 - 8}" y="{sy(1.0) + 4:.1f}" text-anchor="end">100%</text>',
        f'<text class="tick" x="{x0 - 8}" y="{sy(0.0) + 4:.1f}" text-anchor="end">0%</text>',
    ]
    for length in (LENGTHS[0], LENGTHS[len(LENGTHS) // 2], LENGTHS[-1]):
        x = sx(np.log(length * 2))
        parts.append(f'<text class="tick" x="{x:.1f}" y="{y0 + height + 20}" text-anchor="middle">{length * 2}</text>')
    parts.append(
        f'<text class="tick" x="{x0 + width / 2:.0f}" y="{y0 + height + 40}" text-anchor="middle">total length n (log scale)</text>'
    )

    for name in DETECTORS:
        points = " ".join(
            f"{sx(np.log(length * 2)):.1f},{sy(rate):.1f}"
            for length, rate in zip(LENGTHS, rates[name])
        )
        parts.append(f'<polyline class="curve" style="stroke:{COLOURS[name]}" points="{points}"/>')
        for length, rate in zip(LENGTHS, rates[name]):
            parts.append(
                f'<circle cx="{sx(np.log(length * 2)):.1f}" cy="{sy(rate):.1f}" r="3.2" fill="{COLOURS[name]}"/>'
            )
    return "\n  ".join(parts)


def build_svg(matched: dict[str, list[float]], misspecified: dict[str, list[float]]) -> str:
    legend = []
    for i, name in enumerate(DETECTORS):
        x = 150 + i * 200
        legend.append(f'<line class="curve" style="stroke:{COLOURS[name]}" x1="{x}" y1="452" x2="{x + 26}" y2="452"/>')
        legend.append(f'<text class="lgd" x="{x + 34}" y="456">{LABELS[name]}</text>')

    return f"""<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 820 480" width="820" height="480" role="img" aria-labelledby="t d">
  <title id="t">Detection rate versus sample size, with and without correct cyclic structure</title>
  <desc id="d">All detectors are calibrated to a common five percent false-positive rate. On an exact-orbit change the constrained detectors reach a given detection rate at shorter samples than the full detector. On a higher-mode change outside the fundamental subspace the ordering reverses and the full detector leads at every length.</desc>
  <style>
    .hd    {{ font: 600 15px system-ui, sans-serif; fill: #8b949e; }}
    .sub   {{ font: 12px system-ui, sans-serif; fill: #8b949e; }}
    .tick  {{ font: 11px system-ui, sans-serif; fill: #8b949e; }}
    .lgd   {{ font: 13px system-ui, sans-serif; fill: #8b949e; }}
    .axis  {{ stroke: #8b949e; stroke-width: 1.5; }}
    .grid  {{ stroke: #8b949e; stroke-width: 1; stroke-dasharray: 4 4; opacity: 0.5; }}
    .curve {{ fill: none; stroke-width: 2.5; }}
  </style>
  {_panel(70, 60, 300, 320, "Structure is real", "exact orbit, m = 6", matched)}
  {_panel(480, 60, 300, 320, "Structure is wrong", "higher mode, outside the subspace", misspecified)}
  {chr(10) + "  ".join(legend)}
</svg>
"""


def main() -> None:
    print(f"m = {M}, effect = {EFFECT}, {TRIALS} trials per length, seed {SEED}")
    print("scenario: exact_orbit ...", flush=True)
    matched = detection_rates("exact_orbit")
    print("scenario: higher_mode ...", flush=True)
    misspecified = detection_rates("higher_mode")

    out = Path(__file__).resolve().parent.parent / "docs" / "figures" / "worked-example.svg"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(build_svg(matched, misspecified))
    print(f"\nwrote {out}\n")

    header = f"{'n':>7}  " + "  ".join(f"{LABELS[n]:>16}" for n in DETECTORS)
    for name, rates in (("matched (exact orbit)", matched), ("misspecified (higher mode)", misspecified)):
        print(name)
        print(header)
        for i, length in enumerate(LENGTHS):
            row = "  ".join(f"{rates[d][i]:>16.3f}" for d in DETECTORS)
            print(f"{length * 2:>7}  {row}")
        print()


if __name__ == "__main__":
    main()

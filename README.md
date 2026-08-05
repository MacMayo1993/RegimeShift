# RegimeShift

**Geometry of Regime Shifts** — reference implementation and test harness for
*Geometric Complexity in Cyclic Regime Changes: Full, Fundamental-Subspace, and
Shared-Orbit Models under Minimum Description Length* (Mac Mayo).

The manuscript is in [`docs/`](docs/), together with
[`docs/paper-notes.md`](docs/paper-notes.md), which records exactly which
manuscript quantities this code reproduces and which had to be reconstructed.

## The claim being tested

A categorical regime change can be modelled at three levels of structural
constraint. At a **known boundary** with total length `n`, they carry different
minimum-description-length complexity increments:

| Model | Alternative relation between segments | Continuous increment | Leading penalty |
|---|---|---|---|
| **A. Full independent** | arbitrary separate parameters | `m - 1` | `((m-1)/2) log n` |
| **B. Independent fundamental** | separate parameters inside the invariant subspace | `d_fund` | `(d_fund/2) log n` |
| **C. Shared exact orbit** | one shared state plus a relative group element | `0` | `log(m-1)`, constant in `n` |

with `d_fund = 1` for `m = 2` (the sign representation) and `d_fund = 2` for
`m >= 3`. Model B buys efficiency through *dimension reduction*; Model C buys it
through *parameter sharing*, and that is why its leading logarithmic coefficient
is zero rather than merely small.

## Install

```bash
pip install -e ".[test]"
```

Requires Python 3.10+, NumPy, SciPy and pandas.

## Use

```python
import numpy as np
from regimeshift import build_segments, run_all_detectors

m = 6
segments = build_segments(m, "exact_orbit", effect=0.25)
rng = np.random.default_rng(0)
left = rng.multinomial(1600, segments.p_left)
right = rng.multinomial(1600, segments.p_right)

for name, result in run_all_detectors(left, right, m).items():
    print(f"{name:13s} gain={result.raw_gain:8.3f}  penalty={result.penalty:6.3f}  "
          f"score={result.score:8.3f}")
# full          gain=  23.629  penalty=16.712  score=   6.917
# fundamental   gain=  23.481  penalty= 6.685  score=  16.796
# shared_orbit  gain=  22.832  penalty= 1.609  score=  21.223
```

A positive score is the raw MDL rule's declaration of a change. All three
detectors are scored at a known boundary with no location cost — applying one to
a single detector would confound changepoint multiplicity with model dimension
(Section 4.3).

## Reproducing the Monte Carlo study

```bash
# ~1 minute: small grid, used by CI
python -m regimeshift run --grid quick --out results/quick --workers 4

# the manuscript design: 312 configurations, 936 detector rows,
# 468,000 simulated two-segment datasets, base seed 20260713
python -m regimeshift run --grid production --out results/v3 --workers 16
```

Each configuration carries a deterministic, content-derived seed, so results do
not depend on worker count or completion order, and runs resume from the
checkpoint file. Four reports are written:

| File | Contents |
|---|---|
| `full_results.csv` | every detector-level result row |
| `score_regression_summary.csv` | raw-score regression coefficients vs predictions |
| `crossover_estimates.csv` | estimated 50%-power crossover lengths |
| `crossover_ratio_summary.csv` | median practical sample-length ratios |

## Testing

The test suite has two layers.

**Structural tests** (fast, deterministic, no Monte Carlo) verify the properties
the theory asserts exactly — these are the six validations of Section 7.4 plus
the surrounding invariants:

```bash
pytest -m "not slow"           # ~2 minutes
```

- Fisher orthonormality of the Fourier basis, and that the softmax scaling makes
  `|theta|` the Fisher norm
- exact cyclic equivariance `p(R^s theta) = g^s p(theta)`
- the three continuous-dimension increments, and the penalty hierarchy
- equality of Models A and B at `m = 2, 3`, where the fundamental component
  spans the whole nontrivial tangent space
- recovery of planted relative shifts
- Model C's penalty is constant in `n` while the regular penalties grow at
  `d/2` per `log n` — the contrast that gives the claim its content
- the local Jensen–Shannon coefficient `(1 - cos(2 pi / m)) / 4`
- population gains: zero under no change, equal for B and C on exact orbits,
  strictly smaller for C on independent changes, and degraded for both
  constrained families under higher-mode misspecification

**Statistical validation** (slow) runs a real Monte Carlo grid and checks the
empirical penalty slopes against theory, the calibrated sample-length advantage,
and the misspecification reversal:

```bash
pytest -m slow                 # ~5 minutes on 4 cores
pytest                         # everything
```

These are the claims that can only fail probabilistically, so their tolerances
are stated explicitly in `tests/test_statistical_validation.py` and every run is
seeded.

## What the results show

* The full and fundamental coefficients closely track their predicted
  dimensions — this is the clearest evidence that an independently fitted
  fundamental family obeys a different complexity law from the unrestricted
  multinomial.
* The shared-orbit detector has a **near-zero** leading logarithmic coefficient,
  not an exactly constant finite-sample score. Residual drift is expected from
  maximisation over nonidentity shifts, finite-sample MLE bias, and the singular
  orbit-collapse point at the uniform state (Sections 6.3 and 10.3). The tests
  assert the qualified claim, not exact constancy.
* Under a common 5% null calibration the constrained detectors need fewer
  observations on matched data — and *lose* under higher-mode misspecification.
  The advantage is conditional on structural correctness, and the test suite
  asserts both directions.
* Model C's penalty does not grow with `n`, which cuts both ways: its raw
  zero-threshold rule is *not* conservative under the null (at `m = 2` the label
  cost is exactly zero), which is why every comparison here is made at a common
  calibrated 5%. See `docs/paper-notes.md`.

## Layout

```
regimeshift/
  fourier.py      Fisher-orthonormal cyclic Fourier geometry (Section 5)
  detectors.py    the three detectors and their penalties (Sections 3, 4, 7)
  gains.py        detector-specific population gains (Section 6)
  scenarios.py    the three data-generating scenarios (Section 8.2)
  simulation.py   Monte Carlo engine, calibration, seeding (Section 8)
  analysis.py     score regressions and power crossovers (Section 8.3, App. B)
  runner.py       deterministic parallel runner with checkpointing
  cli.py          python -m regimeshift
tests/            structural tests and statistical validation
docs/             the manuscript and reproduction notes
```

## Scope

This is an offline, known-boundary model comparison over independent
categorical observations. Unknown-boundary scanning, sequential stopping rules,
Markov or hidden-state extensions, and the block (codon-phase) family in which
the group acts on phases rather than on the alphabet are not implemented here;
see Sections 11 and 13 of the manuscript.

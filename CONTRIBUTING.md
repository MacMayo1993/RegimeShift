# Contributing

## Setup

```bash
pip install -e ".[test]"
```

## Before you push

```bash
pytest -m "not slow"      # fast structural suite, ~1 minute
pytest -m slow            # seeded Monte Carlo validation, ~5-10 minutes
```

Both run in CI. The structural suite must stay deterministic and fast — it is
the one that catches real regressions in the geometry and the penalties.

## Where a change belongs

* **New geometry** (different group, different representation) →
  `regimeshift/fourier.py`, with an orthonormality and an equivariance test.
* **A new detector or penalty** → `regimeshift/detectors.py`. Every detector
  must report its `dimension_increment` and pay an explicit, documented
  complexity increment. If your detector's penalty depends on `n`, add a test
  pinning its coefficient; if it does not, add a test pinning that it is
  constant.
* **A new data-generating scenario** → `regimeshift/scenarios.py`, with named
  module-level constants (never inline magic numbers) and a test asserting what
  the scenario is *supposed* to violate.
* **A statistical claim** → `tests/test_statistical_validation.py`, marked
  `slow`, with the tolerance justified in the docstring and the manuscript's
  reference value quoted.

## Ground rules

* **No location cost on a single detector.** Known-boundary comparisons must
  treat all three detectors identically; adding a `log n` location term to one
  confounds changepoint multiplicity with model dimension (Section 4.3).
* **Seed everything.** Monte Carlo configurations derive their seed from their
  own contents via `config_seed`, so results never depend on grid ordering or
  worker count. There is a test for this — keep it passing.
* **Do not silently widen a claim.** The shared-orbit result is a near-zero
  *leading* coefficient with finite-sample residual drift, not exact constancy.
  Tests should assert the qualified claim.
* Reconstructed constants from the manuscript are listed in
  `docs/paper-notes.md`. If you change one, update that table.

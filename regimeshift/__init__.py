"""RegimeShift -- geometric complexity in cyclic regime changes.

Reference implementation of the three known-boundary detectors described in
*Geometric Complexity in Cyclic Regime Changes: Full, Fundamental-Subspace, and
Shared-Orbit Models under Minimum Description Length* (docs/).

  Model A (full)          increment (m - 1)/2 * log n
  Model B (fundamental)   increment d_fund/2 * log n,  d_fund = 1 (m=2) else 2
  Model C (shared orbit)  no continuous increment; constant label cost log(m-1)
"""

from .analysis import (
    crossover_estimates,
    crossover_ratio_summary,
    predicted_slope,
    score_regression,
    score_regression_summary,
)
from .detectors import (
    DetectorResult,
    full_detector,
    fundamental_detector,
    label_cost,
    run_all_detectors,
    shared_orbit_detector,
    split_penalty,
)
from .fourier import (
    fisher_inner_product,
    fourier_design_matrix,
    full_dimension,
    fundamental_dimension,
    fundamental_tangent_basis,
    probabilities,
    rotation_matrix,
)
from .gains import population_gains
from .scenarios import MANUSCRIPT_CONSTANTS, Segments, build_segments
from .selection import CANDIDATES, Selection, code_lengths, select_model
from .simulation import BASE_SEED, Config, build_grid, run_config

__version__ = "3.1.0"

__all__ = [
    "__version__",
    "BASE_SEED",
    "Config",
    "DetectorResult",
    "CANDIDATES",
    "MANUSCRIPT_CONSTANTS",
    "Selection",
    "code_lengths",
    "select_model",
    "Segments",
    "build_grid",
    "build_segments",
    "crossover_estimates",
    "crossover_ratio_summary",
    "fisher_inner_product",
    "fourier_design_matrix",
    "full_detector",
    "full_dimension",
    "fundamental_detector",
    "fundamental_dimension",
    "fundamental_tangent_basis",
    "label_cost",
    "population_gains",
    "predicted_slope",
    "probabilities",
    "rotation_matrix",
    "run_all_detectors",
    "run_config",
    "score_regression",
    "score_regression_summary",
    "shared_orbit_detector",
    "split_penalty",
]

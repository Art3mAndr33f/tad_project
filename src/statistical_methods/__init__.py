from __future__ import annotations

"""Statistical methods for biologically-informed TAD boundary consensus.

Modules
-------
weighted_consensus   — per-algorithm weights from ChIP-seq center_ratio (Method 3)
permutation_chipseq  — permutation test for boundary enrichment (Method 4)
changepoint_tad      — PELT change-point on Insulation Score (Method 1, PoC)
"""

from src.statistical_methods.weighted_consensus import (
    load_weights,
    compute_weighted_consensus,
)

__all__ = [
    "load_weights",
    "compute_weighted_consensus",
]
from src.statistical_methods.probabilistic_boundary import (
    compute_probabilistic_consensus,
    run_em,
)

__all__ += [
    "compute_probabilistic_consensus",
    "run_em",
]

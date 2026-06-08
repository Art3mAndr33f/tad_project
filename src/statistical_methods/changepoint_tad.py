from __future__ import annotations

"""PELT change-point detection on Insulation Score for TAD boundary calling.

Method 1 (Priority 3 — proof-of-concept for thesis Chapter 5).

Pipeline
--------
1. Load Hi-C contact matrix via src/data_prep.py (4-level fallback).
2. Compute Insulation Score (IS) profile using a sliding diamond window.
3. Run PELT (ruptures library) over a grid of penalty values beta.
4. A boundary bin is retained if it is detected at >= min_beta_support
   values of beta (stability filter).
5. Convert surviving boundary bins to genomic coordinates (BED3 format).

IS formula (standard diamond window, Crane et al. 2015):
    IS[i] = mean(M[i-w : i, i : i+w])
where w = window_bins.  Normalised to z-score for PELT stability.

Complies with run_<algorithm> contract (§5 rules.md):
    - Returns pd.DataFrame(columns=["chrom","start","end"]) on any error
    - Never returns None
    - All parameters from cfg dict (no hardcoding)
    - logging instead of print
    - seed=42 (not used here — PELT is deterministic)

Dependencies
------------
    ruptures >= 1.1.7   (pip install ruptures)
    numpy, pandas, scipy (already in environment)
"""

import logging
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

# Module-level import so unittest.mock.patch can target
# "src.statistical_methods.changepoint_tad.load_hic_matrix" correctly.
# Aliased from get_matrix (the real function in src/data_prep.py).
# Graceful fallback: if src.data_prep is unavailable (e.g. in isolated tests),
# load_hic_matrix is set to None and run_pelt_tad returns an empty DataFrame.
try:
    from src.data_prep import get_matrix as load_hic_matrix  # noqa: E402
except Exception:  # noqa: BLE001
    load_hic_matrix = None  # type: ignore[assignment]

_EMPTY = pd.DataFrame(columns=["chrom", "start", "end"])


# ---------------------------------------------------------------------------
# Insulation Score
# ---------------------------------------------------------------------------

def compute_insulation_score(
    matrix: np.ndarray,
    window_bins: int,
) -> np.ndarray:
    """Compute the Insulation Score profile from a Hi-C contact matrix.

    Uses the standard sliding diamond window (Crane et al. 2015).
    IS[i] = mean of the off-diagonal block M[i-w:i, i:i+w].
    Bins within window_bins of the matrix edge receive NaN.

    The output is z-score normalised (ignoring NaN) so that PELT penalty
    values are comparable across chromosomes and resolutions.

    Parameters
    ----------
    matrix : (N, N) symmetric numpy array of contact counts.
    window_bins : int — half-width of the diamond window.

    Returns
    -------
    is_profile : (N,) float64 array  (NaN at edges)
    """
    n = matrix.shape[0]
    is_profile = np.full(n, np.nan, dtype=np.float64)

    for i in range(window_bins, n - window_bins):
        block = matrix[i - window_bins : i, i : i + window_bins]
        is_profile[i] = block.mean()

    # z-score normalisation on valid (non-NaN) bins
    valid = ~np.isnan(is_profile)
    if valid.sum() > 1:
        mu  = is_profile[valid].mean()
        std = is_profile[valid].std()
        if std > 0:
            is_profile[valid] = (is_profile[valid] - mu) / std

    return is_profile


# ---------------------------------------------------------------------------
# PELT boundary detection
# ---------------------------------------------------------------------------

def _run_pelt_single_beta(
    signal: np.ndarray,
    beta: float,
) -> np.ndarray:
    """Run PELT with a single penalty beta; return change-point bin indices.

    Parameters
    ----------
    signal : 1-D float array (NaN already replaced with 0.0 before call).
    beta   : PELT penalty (larger → fewer breakpoints).

    Returns
    -------
    breakpoints : sorted int array of change-point bin indices
                  (ruptures convention: last index == len(signal) excluded)
    """
    try:
        import ruptures as rpt  # noqa: PLC0415  (deferred import)
    except ImportError as exc:
        raise ImportError(
            "ruptures is required for PELT change-point detection. "
            "Install with: pip install ruptures>=1.1.7"
        ) from exc

    algo = rpt.Pelt(model="rbf", min_size=2, jump=1)
    algo.fit(signal.reshape(-1, 1).astype(np.float64))
    # ruptures returns breakpoints including len(signal) as sentinel
    bkps = algo.predict(pen=beta)
    # Remove sentinel (last element == len(signal))
    return np.array(sorted(b for b in bkps if b < len(signal)), dtype=np.int64)


def _stability_filter(
    all_breakpoints: list[np.ndarray],
    n_bins: int,
    min_beta_support: int,
    tolerance_bins: int = 1,
) -> np.ndarray:
    """Keep bins that appear as change-points in >= min_beta_support beta values.

    Two breakpoints at bins i and j are considered the same boundary if
    |i - j| <= tolerance_bins.

    Parameters
    ----------
    all_breakpoints : list of breakpoint arrays, one per beta value.
    n_bins          : total number of bins (for bounds checking).
    min_beta_support: minimum number of beta values a boundary must appear in.
    tolerance_bins  : neighbourhood for matching across beta runs.

    Returns
    -------
    stable_bins : sorted int array of stable boundary bin indices.
    """
    if not all_breakpoints:
        return np.array([], dtype=np.int64)

    # Aggregate all breakpoints into a single array with beta-run labels
    # Count support for each unique (after tolerance merging) boundary
    all_bins: list[int] = []
    for bkps in all_breakpoints:
        all_bins.extend(bkps.tolist())

    if not all_bins:
        return np.array([], dtype=np.int64)

    all_bins_arr = np.array(sorted(set(all_bins)), dtype=np.int64)

    # For each unique candidate bin, count how many beta runs contributed
    # a breakpoint within tolerance_bins
    support = np.zeros(len(all_bins_arr), dtype=np.int64)

    for bkps in all_breakpoints:
        if len(bkps) == 0:
            continue
        # Broadcasting: all_bins_arr (C,) vs bkps (K,)
        diff = np.abs(all_bins_arr[:, None] - bkps[None, :])   # (C, K)
        hit  = (diff <= tolerance_bins).any(axis=1)             # (C,)
        support += hit.astype(np.int64)

    stable_mask = support >= min_beta_support
    return all_bins_arr[stable_mask]


# ---------------------------------------------------------------------------
# Public API — run_<algorithm> contract
# ---------------------------------------------------------------------------

def run_pelt_tad(
    chrom: str,
    resolution: int,
    data_path: str,
    cfg: Optional[dict],
    **kwargs,
) -> pd.DataFrame:
    """Detect TAD boundaries using PELT change-point detection on IS.

    Complies with the run_<algorithm> contract (rules.md §5):
        - Returns pd.DataFrame(columns=["chrom","start","end"]) on any error
        - Never returns None
        - Parameters from cfg (no hardcoding)

    Parameters
    ----------
    chrom      : chromosome with "chr" prefix (e.g. "chr17")
    resolution : bin size in bp (e.g. 25_000)
    data_path  : path to Hi-C RAWobserved data directory
    cfg        : config dict; reads cfg["changepoint"] sub-dict
    **kwargs   : passed through (ignored)

    Returns
    -------
    pd.DataFrame with columns ["chrom", "start", "end"]
    """
    if cfg is None:
        cfg = {}

    cp_cfg = cfg.get("changepoint", {})
    window_bins     = int(cp_cfg.get("window_bins",     5))
    beta_min        = float(cp_cfg.get("beta_min",      0.5))
    beta_max        = float(cp_cfg.get("beta_max",      5.0))
    n_beta          = int(cp_cfg.get("n_beta",          10))
    min_beta_support= int(cp_cfg.get("min_beta_support", 3))
    tolerance_bins  = int(cfg.get("consensus", {}).get("tolerance_bins", 1))

    logger.info(
        "[pelt] %s@%dbp | window=%d bins | beta=[%.1f,%.1f] n=%d | "
        "min_support=%d",
        chrom, resolution, window_bins,
        beta_min, beta_max, n_beta, min_beta_support,
    )

    # ── 1. Load contact matrix ────────────────────────────────────────────────
    try:
        matrix = load_hic_matrix(cfg, chrom, resolution)
    except Exception as exc:  # noqa: BLE001
        logger.warning("[pelt] %s: failed to load matrix — %s", chrom, exc)
        return _EMPTY.copy()

    if matrix is None or matrix.size == 0:
        logger.warning("[pelt] %s: empty matrix", chrom)
        return _EMPTY.copy()

    n_bins = matrix.shape[0]
    logger.debug("[pelt] %s: matrix shape %s", chrom, matrix.shape)

    # ── 2. Insulation Score ───────────────────────────────────────────────────
    try:
        is_profile = compute_insulation_score(matrix, window_bins)
    except Exception as exc:  # noqa: BLE001
        logger.warning("[pelt] %s: IS computation failed — %s", chrom, exc)
        return _EMPTY.copy()

    # Replace NaN with 0.0 for PELT (edge bins)
    signal = np.where(np.isnan(is_profile), 0.0, is_profile)

    if np.all(signal == 0.0):
        logger.warning("[pelt] %s: IS signal is all-zero — no boundaries", chrom)
        return _EMPTY.copy()

    # ── 3. PELT over beta grid ────────────────────────────────────────────────
    beta_grid = np.linspace(beta_min, beta_max, n_beta)
    all_breakpoints: list[np.ndarray] = []

    for beta in beta_grid:
        try:
            bkps = _run_pelt_single_beta(signal, beta)
            all_breakpoints.append(bkps)
            logger.debug(
                "[pelt] %s beta=%.2f → %d breakpoints", chrom, beta, len(bkps)
            )
        except Exception as exc:  # noqa: BLE001
            logger.warning(
                "[pelt] %s beta=%.2f failed — %s", chrom, beta, exc
            )

    if not all_breakpoints:
        logger.warning("[pelt] %s: no PELT runs succeeded", chrom)
        return _EMPTY.copy()

    # ── 4. Stability filter ───────────────────────────────────────────────────
    stable_bins = _stability_filter(
        all_breakpoints, n_bins, min_beta_support, tolerance_bins
    )

    logger.info(
        "[pelt] %s@%dbp: %d stable boundaries (min_beta_support=%d/%d)",
        chrom, resolution, len(stable_bins), min_beta_support, n_beta,
    )

    if len(stable_bins) == 0:
        return _EMPTY.copy()

    # ── 5. Bin indices → genomic coordinates (BED3) ───────────────────────────
    # Each change-point bin i is expressed as [i*res, (i+1)*res)
    # We report the full bin as the boundary interval (consistent with
    # other algorithms in the pipeline).
    starts = stable_bins * resolution
    ends   = starts + resolution

    result = pd.DataFrame({
        "chrom": chrom,
        "start": starts,
        "end":   ends,
    })
    return result.reset_index(drop=True)
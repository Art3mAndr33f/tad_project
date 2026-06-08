from __future__ import annotations

"""Permutation test for ChIP-seq enrichment at TAD boundaries.

Method 4 (Priority 2).

For each boundary b_j:
    obs_density = number of ChIP-seq peaks overlapping [b_j - window_bp, b_j + window_bp]
                  normalised by window size (peaks / bp)

Null distribution: N random positions on the same chromosome
    (excluding ±telomere_margin_bp from both ends).

p-value  = fraction of random draws with density >= obs_density  (one-sided, right)
p_adj    = Benjamini-Hochberg FDR correction
significant = p_adj < fdr_threshold

Vectorised counting reuses the numpy broadcasting pattern from
src/validation.py::_count_ctcf_overlaps  (~100x faster than Python loop).

Output
------
CSV : results/stats/permutation_{track}_{chrom}_{res}bp.csv
BED : results/consensus/permutation_filtered_{track}_{chrom}_{res}bp.bed
      (BED3, only significant boundaries)
"""

import logging
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _count_peaks_vectorised(
    positions: np.ndarray,   # centre positions in bp, shape (N,)
    peak_mids: np.ndarray,   # ChIP-seq peak midpoints in bp, shape (M,)
    window_bp: int,
) -> np.ndarray:
    """Count peaks within ±window_bp of each position.

    Vectorised via numpy broadcasting; identical pattern to
    src/validation.py::_count_ctcf_overlaps.

    Parameters
    ----------
    positions : (N,) array of query positions
    peak_mids : (M,) array of ChIP-seq peak midpoints
    window_bp : half-width of counting window

    Returns
    -------
    counts : (N,) int array
    """
    if len(peak_mids) == 0 or len(positions) == 0:
        return np.zeros(len(positions), dtype=np.int64)

    # Split into chunks to avoid OOM on large chromosomes
    # chunk of 2 000 positions × 44 000 peaks ≈ 352 MB float64 → safe
    chunk = 2_000
    counts = np.empty(len(positions), dtype=np.int64)

    for lo in range(0, len(positions), chunk):
        hi  = min(lo + chunk, len(positions))
        pos_chunk = positions[lo:hi]                            # (C,)
        diff = np.abs(pos_chunk[:, None] - peak_mids[None, :]) # (C, M)
        counts[lo:hi] = (diff <= window_bp).sum(axis=1)

    return counts


def _bh_correction(pvalues: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR correction.

    Returns adjusted p-values (same length as input).
    Equivalent to statsmodels.stats.multitest.multipletests(..., method='fdr_bh')
    but without the statsmodels dependency.
    """
    n = len(pvalues)
    if n == 0:
        return np.array([], dtype=np.float64)

    order   = np.argsort(pvalues)
    ranks   = np.empty(n, dtype=np.float64)
    ranks[order] = np.arange(1, n + 1, dtype=np.float64)

    # p_adj_i = p_i * n / rank_i  (then take cumulative min from largest)
    p_adj = pvalues * (n / ranks)

    # Enforce monotonicity: cumulative min from the largest rank down
    p_adj_ordered = p_adj[order]
    for i in range(n - 2, -1, -1):
        p_adj_ordered[i] = min(p_adj_ordered[i], p_adj_ordered[i + 1])

    result = np.empty(n, dtype=np.float64)
    result[order] = p_adj_ordered
    return np.clip(result, 0.0, 1.0)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def compute_boundary_enrichment(
    boundaries_df: pd.DataFrame,
    chipseq_df: pd.DataFrame,
    chrom: str,
    chrom_size: int,
    window_bp: int,
    n_permutations: int,
    rng: np.random.Generator,
    telomere_margin_bp: int,
    fdr_threshold: float = 0.05,
) -> pd.DataFrame:
    """Permutation test for ChIP-seq enrichment at TAD boundaries.

    Parameters
    ----------
    boundaries_df : pd.DataFrame
        Columns: chrom, start, end.  Boundary centre = (start + end) // 2.
    chipseq_df : pd.DataFrame
        Columns: chrom, start, end.  Peak midpoint = (start + end) // 2.
    chrom : str
        Chromosome (e.g. "chr1").  Only rows matching this chrom are used.
    chrom_size : int
        Total chromosome length in bp (from hg19 chrom.sizes or config).
    window_bp : int
        Half-width of ChIP-seq counting window around each boundary centre.
    n_permutations : int
        Number of random positions drawn for the null distribution.
    rng : np.random.Generator
        Seeded generator — always np.random.default_rng(42) in production.
    telomere_margin_bp : int
        Random positions are drawn from
        [telomere_margin_bp, chrom_size - telomere_margin_bp].
    fdr_threshold : float
        Benjamini-Hochberg FDR cutoff for the 'significant' flag.

    Returns
    -------
    pd.DataFrame
        Original columns (chrom, start, end) plus:
            boundary_centre  — midpoint in bp
            obs_count        — raw ChIP-seq peak count in ±window_bp
            obs_density      — obs_count / (2 * window_bp)  [peaks / bp]
            pvalue           — one-sided empirical p-value
            pvalue_adj       — BH-corrected p-value
            significant      — bool, pvalue_adj < fdr_threshold
        Empty DataFrame with the same columns if no boundaries on chrom.
    """
    _EXTRA_COLS = [
        "boundary_centre", "obs_count", "obs_density",
        "pvalue", "pvalue_adj", "significant",
    ]
    _ALL_COLS = ["chrom", "start", "end"] + _EXTRA_COLS
    _EMPTY = pd.DataFrame(columns=_ALL_COLS)

    # ── 1. Filter to target chromosome ───────────────────────────────────────
    bnd = boundaries_df[boundaries_df["chrom"] == chrom].copy()
    cs  = chipseq_df[chipseq_df["chrom"] == chrom].copy()

    if bnd.empty:
        logger.info(
            "[permutation] %s: no boundaries — returning empty DataFrame", chrom
        )
        return _EMPTY

    n_bnd = len(bnd)
    logger.info(
        "[permutation] %s: %d boundaries, %d ChIP-seq peaks, window=±%d bp",
        chrom, n_bnd, len(cs), window_bp,
    )

    # ── 2. Boundary centres & peak midpoints ─────────────────────────────────
    centres  = ((bnd["start"].values + bnd["end"].values) // 2).astype(np.int64)
    if cs.empty:
        peak_mids = np.array([], dtype=np.int64)
    else:
        peak_mids = ((cs["start"].values + cs["end"].values) // 2).astype(np.int64)

    # ── 3. Observed counts ───────────────────────────────────────────────────
    obs_counts = _count_peaks_vectorised(centres, peak_mids, window_bp)
    obs_density = obs_counts.astype(np.float64) / (2.0 * window_bp)

    # ── 4. Null distribution via permutation ─────────────────────────────────
    rand_lo = telomere_margin_bp
    rand_hi = chrom_size - telomere_margin_bp

    if rand_hi <= rand_lo:
        logger.warning(
            "[permutation] %s: chrom_size=%d is too small for telomere_margin=%d "
            "— using full chromosome for random draws",
            chrom, chrom_size, telomere_margin_bp,
        )
        rand_lo, rand_hi = 0, chrom_size

    # null_counts shape: (n_permutations,)
    # Draw n_permutations random positions; count peaks for each.
    rand_positions = rng.integers(rand_lo, rand_hi, size=n_permutations)
    null_counts    = _count_peaks_vectorised(rand_positions, peak_mids, window_bp)
    null_density   = null_counts.astype(np.float64) / (2.0 * window_bp)

    logger.debug(
        "[permutation] %s: null distribution — mean=%.4f  std=%.4f  "
        "max=%d  (N=%d)",
        chrom,
        null_density.mean(), null_density.std(),
        null_counts.max() if len(null_counts) > 0 else 0,
        n_permutations,
    )

    # ── 5. Empirical p-values  (one-sided: P(null >= observed)) ──────────────
    # Broadcasting: obs_density (N,) vs null_density (P,)
    # p_i = #{null_density >= obs_density_i} / n_permutations
    pvalues = (
        (null_density[None, :] >= obs_density[:, None])
        .sum(axis=1)
        .astype(np.float64)
        / n_permutations
    )

    # ── 6. BH FDR correction ─────────────────────────────────────────────────
    pvalues_adj = _bh_correction(pvalues)
    significant = pvalues_adj < fdr_threshold

    n_sig = significant.sum()
    logger.info(
        "[permutation] %s: %d/%d boundaries significant "
        "(FDR<%.2f, window=±%d bp)",
        chrom, n_sig, n_bnd, fdr_threshold, window_bp,
    )

    # ── 7. Assemble result ────────────────────────────────────────────────────
    result = bnd[["chrom", "start", "end"]].copy().reset_index(drop=True)
    result["boundary_centre"] = centres
    result["obs_count"]       = obs_counts
    result["obs_density"]     = np.round(obs_density, 8)
    result["pvalue"]          = np.round(pvalues, 8)
    result["pvalue_adj"]      = np.round(pvalues_adj, 8)
    result["significant"]     = significant

    return result


# ---------------------------------------------------------------------------
# Convenience I/O
# ---------------------------------------------------------------------------

def save_permutation_results(
    df: pd.DataFrame,
    csv_path: Path,
    bed_path: Path,
) -> None:
    """Save full stats to CSV and significant boundaries to BED3.

    Parameters
    ----------
    df       : output of compute_boundary_enrichment()
    csv_path : results/stats/permutation_{track}_{chrom}_{res}bp.csv
    bed_path : results/consensus/permutation_filtered_{track}_{chrom}_{res}bp.bed
    """
    csv_path = Path(csv_path)
    bed_path = Path(bed_path)
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    bed_path.parent.mkdir(parents=True, exist_ok=True)

    # Full CSV with header
    df.to_csv(csv_path, index=False)
    logger.info("Saved permutation stats → %s (%d rows)", csv_path, len(df))

    # Significant boundaries only → BED3, no header
    sig = df[df["significant"]][["chrom", "start", "end"]]
    sig.to_csv(bed_path, sep="\t", header=False, index=False)
    logger.info(
        "Saved filtered BED → %s (%d significant boundaries)",
        bed_path, len(sig),
    )
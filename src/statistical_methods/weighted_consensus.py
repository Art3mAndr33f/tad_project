from __future__ import annotations

"""Weighted TAD boundary consensus weighted by ChIP-seq center_ratio.

Method 3 (Priority 1).

For each candidate boundary position j (in bin units):
    B_j = Σ(w_i · I[algo_i has boundary at j ± tolerance_bins]) / Σ(w_i)
Boundary is retained if B_j >= theta.

Weights w_i = center_ratio(algo_i, track, chr1, 100kb)
from results/stats/chipseq_validation_summary.csv.
Algorithms with excluded verdicts (NO_DATA, SHIFTED) → w_i = 0.0.

Output: BED4  chrom  start  end  weighted_support
        results/consensus/weighted_consensus_{chrom}_{res}bp.bed
"""

import logging
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# BED reader  (auto-detects header/no-header — same pattern as _load_consensus_bed)
# ---------------------------------------------------------------------------

def _read_bed3(path: Path) -> pd.DataFrame:
    """Read a 3- (or more) column BED file; auto-detect header by first char."""
    empty = pd.DataFrame(columns=["chrom", "start", "end"])
    if not path.exists():
        logger.warning("BED file not found: %s", path)
        return empty

    with path.open() as fh:
        first = fh.readline().strip()

    # Header present when first field is not a chromosome token
    first_field = first.split("\t")[0] if "\t" in first else first.split()[0]
    has_header = not (first_field.startswith("chr") or first_field.lstrip("-").isdigit())

    try:
        df = pd.read_csv(
            path,
            sep="\t",
            header=0 if has_header else None,
            usecols=[0, 1, 2],
        )
        df.columns = ["chrom", "start", "end"]
        df["start"] = pd.to_numeric(df["start"], errors="coerce").astype("Int64")
        df["end"]   = pd.to_numeric(df["end"],   errors="coerce").astype("Int64")
        df = df.dropna(subset=["start", "end"]).copy()
        df["start"] = df["start"].astype(int)
        df["end"]   = df["end"].astype(int)
        return df.reset_index(drop=True)
    except Exception as exc:  # noqa: BLE001
        logger.warning("Failed to read BED %s: %s", path, exc)
        return empty


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Weight transforms
# ---------------------------------------------------------------------------

def _apply_transform(
    raw_weights: dict[str, float],
    transform: str,
    params: dict,
) -> dict[str, float]:
    """Apply a nonlinear transform to non-zero raw weights.

    The transform is applied AFTER exclude_verdicts filtering, so excluded
    algorithms (weight=0.0) are never transformed — they stay at 0.0.

    Supported transforms
    --------------------
    linear   — w = max(0, cr)                          (identity, default)
    power    — w = max(0, cr - baseline) ** gamma
               Handles cr < baseline (e.g. topdom cr=0.798, baseline=1.0)
               by clamping to 0 before exponentiation.
               Negative cr is also safely handled via max(0, ...).
    softmax  — w = exp(alpha * (cr - baseline))
               Applied only to non-excluded algos; excluded stay 0.0.
               Note: does NOT normalise to sum=1 here — normalisation
               happens in compute_weighted_consensus.
    zscore   — w = max(0, (cr - mu) / sigma)
               mu/sigma computed from non-excluded (non-zero) entries only.

    Parameters
    ----------
    raw_weights : {algo: center_ratio} — excluded algos already set to 0.0
    transform   : one of "linear" / "power" / "softmax" / "zscore"
    params      : transform hyperparameters (baseline, gamma, alpha)
    """
    baseline = float(params.get("baseline", 1.0))
    gamma    = float(params.get("gamma",    2.0))
    alpha    = float(params.get("alpha",    5.0))

    # Clamp ALL non-positive weights to 0.0 before transform.
    # This covers three cases:
    #   (a) excluded algos already set to 0.0 by load_weights
    #   (b) pathological cr < 0 (inverted profile) → safe clamp
    #   (c) direct calls to _apply_transform with arbitrary raw_weights
    result = {a: max(0.0, w) for a, w in raw_weights.items()}

    # Non-excluded algorithms: those with positive weight after clamp
    active = {a: w for a, w in result.items() if w > 0.0}

    if not active:
        return result   # all excluded/clamped → nothing to transform

    if transform == "linear":
        # w = max(0, cr)  — already done in load_weights via clamp
        pass

    elif transform == "power":
        # w = max(0, cr - baseline) ** gamma
        # cr < baseline → (cr - baseline) < 0 → max(0, ...) = 0
        # cr < 0        → max(0, negative - baseline) = 0 — also safe
        for algo, cr in active.items():
            excess = cr - baseline
            result[algo] = max(0.0, excess) ** gamma
            logger.debug(
                "_apply_transform power: algo=%s cr=%.4f excess=%.4f w=%.6f",
                algo, cr, excess, result[algo],
            )

    elif transform == "softmax":
        # w = exp(alpha * (cr - baseline))
        # Excluded algos stay 0.0 (not part of softmax).
        # Uses (cr - baseline) as argument so baseline becomes the
        # "neutral" level: cr=baseline → exp(0)=1.
        for algo, cr in active.items():
            result[algo] = float(np.exp(alpha * (cr - baseline)))
            logger.debug(
                "_apply_transform softmax: algo=%s cr=%.4f w=%.6f",
                algo, cr, result[algo],
            )

    elif transform == "zscore":
        # w = max(0, (cr - mu) / sigma)
        # mu/sigma from active (non-excluded) weights only.
        vals = np.array(list(active.values()), dtype=np.float64)
        mu   = float(vals.mean())
        std  = float(vals.std())
        if std < 1e-12:
            # All active weights are identical → uniform after zscore → linear
            logger.warning(
                "_apply_transform zscore: std≈0 for active weights "
                "(all identical center_ratios) — falling back to linear"
            )
        else:
            for algo, cr in active.items():
                result[algo] = max(0.0, (cr - mu) / std)
                logger.debug(
                    "_apply_transform zscore: algo=%s cr=%.4f z=%.4f",
                    algo, cr, result[algo],
                )

    else:
        raise ValueError(
            f"Unknown weight_transform: {transform!r}. "
            "Choose from: linear / power / softmax / zscore"
        )

    return result


def load_weights(
    summary_csv: str,
    track: str = "rad21",
    exclude_verdicts: Optional[list[str]] = None,
    weight_transform: str = "linear",
    transform_params: Optional[dict] = None,
) -> dict[str, float]:
    """Load per-algorithm weights from ChIP-seq validation summary CSV.

    Parameters
    ----------
    summary_csv:
        Path to results/stats/chipseq_validation_summary.csv.
        Required columns: algorithm, center_ratio, verdict, track.
    track:
        ChIP-seq track to use for weighting (e.g. "rad21").
    exclude_verdicts:
        Algorithms with these verdicts receive weight 0.0.
        Defaults to ["NO_DATA", "SHIFTED"].
    weight_transform:
        Nonlinear transform applied to center_ratio after exclusion.
        "linear"  — identity (default, backward-compatible)
        "power"   — max(0, cr - baseline)^gamma
        "softmax" — exp(alpha * (cr - baseline))
        "zscore"  — max(0, (cr - mu) / sigma)
        See _apply_transform() for details.
    transform_params:
        Dict of hyperparameters for the chosen transform.
        Keys: baseline (default 1.0), gamma (default 2.0), alpha (default 5.0).

    Returns
    -------
    dict[str, float]
        {algorithm_name: transformed_weight}.
        Excluded algos → 0.0.  Negative raw weights → 0.0 (clamped).
    """
    if exclude_verdicts is None:
        exclude_verdicts = ["NO_DATA", "SHIFTED"]
    if transform_params is None:
        transform_params = {}

    summary_path = Path(summary_csv)
    if not summary_path.exists():
        raise FileNotFoundError(f"Summary CSV not found: {summary_path}")

    df = pd.read_csv(summary_path)

    required = {"algorithm", "center_ratio", "verdict", "track"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Summary CSV missing columns: {missing}")

    df_track = df[df["track"].str.lower() == track.lower()].copy()

    if df_track.empty:
        logger.warning(
            "No rows found for track='%s' in %s — returning empty weights",
            track, summary_csv,
        )
        return {}

    # ── Step 1: raw weights (center_ratio, clamped; excluded → 0.0) ──────────
    raw: dict[str, float] = {}
    for _, row in df_track.iterrows():
        algo    = str(row["algorithm"])
        verdict = str(row["verdict"]).strip()
        cr      = float(row["center_ratio"])

        if verdict in exclude_verdicts:
            raw[algo] = 0.0
        else:
            raw[algo] = max(cr, 0.0)   # clamp: cr<0 is pathological but safe

    # ── Step 2: nonlinear transform ───────────────────────────────────────────
    weights = _apply_transform(raw, weight_transform, transform_params)

    # ── Step 3: log final weights ─────────────────────────────────────────────
    for algo, w in weights.items():
        verdict = df_track.loc[
            df_track["algorithm"] == algo, "verdict"
        ].values
        v_str = str(verdict[0]) if len(verdict) > 0 else "?"
        cr_val = raw.get(algo, 0.0)
        logger.info(
            "load_weights [%s]: algo='%s' verdict=%s "
            "cr=%.4f → w=%.6f",
            weight_transform, algo, v_str, cr_val, w,
        )

    return weights


# ---------------------------------------------------------------------------

def _boundaries_from_tad_df(
    df: pd.DataFrame,
    chrom: str,
    resolution: int,
) -> np.ndarray:
    """Return unique boundary bin indices (start-bin ∪ end-bin) for chrom."""
    sub = df[df["chrom"] == chrom]
    if sub.empty:
        return np.array([], dtype=np.int64)

    starts = (sub["start"].values // resolution).astype(np.int64)
    ends   = (sub["end"].values   // resolution).astype(np.int64)
    return np.unique(np.concatenate([starts, ends]))


def compute_weighted_consensus(
    tad_files: dict[str, Path],
    weights: dict[str, float],
    chrom: str,
    resolution: int,
    theta: float,
    tolerance_bins: int,
) -> pd.DataFrame:
    """Compute weighted boundary consensus.

    For every candidate boundary bin j:
        B_j = Σ_i [w_i · 1(algo_i has a boundary at j ± tolerance_bins)]
              ─────────────────────────────────────────────────────────────
                                    Σ_i w_i

    Boundary is retained if B_j >= theta.

    Parameters
    ----------
    tad_files:
        {algo_name: path_to_bed}.  Files not found → treated as empty.
    weights:
        {algo_name: weight}.  Missing algos default to 0.0.
    chrom:
        Chromosome with "chr" prefix, e.g. "chr1".
    resolution:
        Bin size in bp (e.g. 100_000).
    theta:
        Inclusion threshold for weighted score in [0, 1].
    tolerance_bins:
        Neighbourhood radius in bins (same as consensus.tolerance_bins).

    Returns
    -------
    pd.DataFrame
        Columns: chrom, start, end, weighted_support.
        Empty DataFrame with the same columns if no boundaries pass threshold.
    """
    _EMPTY = pd.DataFrame(columns=["chrom", "start", "end", "weighted_support"])

    # ── 1. Load boundary bins per algorithm ──────────────────────────────────
    algo_bins: dict[str, np.ndarray] = {}
    for algo, path in tad_files.items():
        df  = _read_bed3(Path(path))
        bins = _boundaries_from_tad_df(df, chrom, resolution)
        algo_bins[algo] = bins
        logger.debug(
            "compute_weighted_consensus: algo=%s chrom=%s → %d boundary bins",
            algo, chrom, len(bins),
        )

    # ── 2. Pool all candidate boundary positions ──────────────────────────────
    all_arrays = [b for b in algo_bins.values() if len(b) > 0]
    if not all_arrays:
        logger.info(
            "[weighted_consensus] %s@%dbp: no boundaries from any algorithm",
            chrom, resolution,
        )
        return _EMPTY

    candidate_bins: np.ndarray = np.unique(np.concatenate(all_arrays))

    # ── 3. Effective total weight (only algos present in tad_files) ───────────
    total_weight = sum(max(weights.get(algo, 0.0), 0.0) for algo in tad_files)
    if total_weight == 0.0:
        logger.warning(
            "[weighted_consensus] %s@%dbp: total_weight=0.0 "
            "(all algorithms excluded or missing weights) — returning empty",
            chrom, resolution,
        )
        return _EMPTY

    # ── 4. Vectorised weighted scoring ───────────────────────────────────────
    # candidate_bins : (N,)
    # For each algorithm, compute a boolean hit-vector via broadcasting.
    scores = np.zeros(len(candidate_bins), dtype=np.float64)

    for algo, bins in algo_bins.items():
        w = max(weights.get(algo, 0.0), 0.0)
        if w == 0.0 or len(bins) == 0:
            continue
        # diff[i, k] = |candidate_bins[i] - bins[k]|
        diff = np.abs(candidate_bins[:, None] - bins[None, :])   # (N, M)
        hit  = (diff <= tolerance_bins).any(axis=1)              # (N,)
        scores += w * hit.astype(np.float64)

    scores /= total_weight   # normalise → [0, 1]

    # ── 5. Apply threshold ────────────────────────────────────────────────────
    mask           = scores >= theta
    selected_bins  = candidate_bins[mask]
    selected_scores = scores[mask]

    logger.info(
        "[weighted_consensus] %s@%dbp: %d/%d boundaries retained "
        "(theta=%.2f, tolerance=%d bins)",
        chrom, resolution,
        len(selected_bins), len(candidate_bins),
        theta, tolerance_bins,
    )

    if len(selected_bins) == 0:
        return _EMPTY

    # ── 6. Bin indices → genomic coordinates ─────────────────────────────────
    result = pd.DataFrame({
        "chrom":            chrom,
        "start":            selected_bins * resolution,
        "end":              selected_bins * resolution + resolution,
        "weighted_support": np.round(selected_scores, 6),
    })
    return result.reset_index(drop=True)


# ---------------------------------------------------------------------------
# Convenience I/O
# ---------------------------------------------------------------------------

def save_weighted_consensus(
    df: pd.DataFrame,
    out_path: Path,
) -> None:
    """Save weighted consensus to BED4 (no header, tab-separated)."""
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, sep="\t", header=False, index=False)
    logger.info("Saved weighted consensus → %s (%d rows)", out_path, len(df))
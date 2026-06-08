from __future__ import annotations

"""Probabilistic TAD boundary model (Method 5).

Crowdsourcing-style latent variable model for boundary consensus.

Motivation
----------
M3 (weighted consensus) hits a ceiling because all algorithms detect
boundaries in overlapping positions — reweighting changes scores but not
the accepted set.  The fundamental issue is that M3 treats sensitivity
(center_ratio) but ignores specificity (false positive rate).

An algorithm with 400 boundaries on chr1 (armatus) contributes noise even
with a decent center_ratio.  A model that accounts for both sens AND fpr
separates signal from noise more cleanly.

Model
-----
For each candidate boundary bin b:
    z_b ∈ {0,1}      — latent: is b a true boundary?
    x_{ib} ∈ {0,1}   — algorithm i detected a boundary at b (±tol bins)

Generative model (Naive Bayes / Dawid-Skene style):
    P(x_{ib}=1 | z_b=1) = sens_i     # sensitivity
    P(x_{ib}=1 | z_b=0) = fpr_i      # false positive rate
    P(z_b=1)            = π_b         # prior (from local Hi-C signal)

Posterior (E-step):
    P(z_b=1 | x_{·b}) ∝
        π_b   · ∏_i sens_i^x_{ib} · (1-sens_i)^(1-x_{ib})
    ────────────────────────────────────────────────────────
    above + (1-π_b) · ∏_i fpr_i^x_{ib}  · (1-fpr_i)^(1-x_{ib})

M-step: re-estimate sens_i and fpr_i from posteriors.

EM runs for max_iter iterations or until convergence (Δ log-likelihood < tol).

Parameter initialisation
------------------------
sens_i = sigmoid(k_sens * max(0, center_ratio_i - 1.0))
         High center_ratio  → high sensitivity.
         cr ≤ 1.0 (random)  → sens ≈ 0.5 (uninformative).

fpr_i  = n_boundaries_i / n_bins_chrom
         Many boundaries relative to chromosome length → high FPR.
         Clamped to [fpr_min, fpr_max] to avoid degenerate likelihoods.

π_b    = sigmoid(k_prior * (local_signal_b / global_mean_signal - 1))
         Empty Hi-C regions → local_signal ≈ 0 → π_b ≈ 0.5 (uninformative).
         Active regions     → local_signal >> mean → π_b → 1.
         If no matrix available: π_b = global_prior (scalar from config).

Output
------
BED4: chrom  start  end  posterior_prob
Only bins with posterior ≥ theta are reported.
results/consensus/probabilistic_boundary_{chrom}_{res}bp.bed

CSV:  results/stats/probabilistic_boundary_{chrom}_{res}bp.csv
Full table: all candidate bins with posterior and whether significant.
"""

import logging
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

_EMPTY = pd.DataFrame(columns=["chrom", "start", "end"])

# ---------------------------------------------------------------------------
# Sigmoid helper
# ---------------------------------------------------------------------------

def _sigmoid(x: np.ndarray | float) -> np.ndarray | float:
    return 1.0 / (1.0 + np.exp(-np.clip(x, -500, 500)))


# ---------------------------------------------------------------------------
# Parameter initialisation
# ---------------------------------------------------------------------------

def init_sensitivity(
    center_ratios: np.ndarray,
    k_sens: float = 8.0,
) -> np.ndarray:
    """Initialise per-algorithm sensitivity from center_ratio.

    sens_i = sigmoid(k_sens * max(0, cr_i - 1.0))

    Properties:
        cr = 1.0  → sigmoid(0)         = 0.500  (uninformative)
        cr = 1.1  → sigmoid(0.8)       = 0.690
        cr = 1.3  → sigmoid(2.4)       = 0.917
        cr = 1.45 → sigmoid(3.6)       = 0.973
        cr ≤ 1.0  → sigmoid(≤0)        ≤ 0.500  (clamped to 0 excess → 0.5)

    Parameters
    ----------
    center_ratios : (n_algos,) array of center_ratio values.
    k_sens        : steepness of sigmoid (larger → more discriminative).

    Returns
    -------
    sens : (n_algos,) array in (0, 1).
    """
    excess = np.maximum(0.0, center_ratios - 1.0)
    return _sigmoid(k_sens * excess).astype(np.float64)


def init_fpr(
    n_boundaries: np.ndarray,
    n_bins: int,
    fpr_min: float = 0.01,
    fpr_max: float = 0.60,
) -> np.ndarray:
    """Initialise per-algorithm FPR from number of boundaries.

    fpr_i = n_boundaries_i / n_bins   (clamped to [fpr_min, fpr_max])

    Rationale: if an algorithm calls B boundaries on a chromosome of N bins,
    then in the absence of any signal it would call each bin with probability
    B/N.  This is a conservative upper bound on FPR.

    Parameters
    ----------
    n_boundaries : (n_algos,) array of boundary counts.
    n_bins       : number of bins in chromosome.
    fpr_min, fpr_max : clamp bounds.

    Returns
    -------
    fpr : (n_algos,) array in [fpr_min, fpr_max].
    """
    raw = n_boundaries.astype(np.float64) / max(n_bins, 1)
    return np.clip(raw, fpr_min, fpr_max).astype(np.float64)


def compute_prior(
    matrix: Optional[np.ndarray],
    n_bins: int,
    k_prior: float = 3.0,
    global_prior: float = 0.1,
) -> np.ndarray:
    """Compute per-bin prior probability of being a true boundary.

    Uses the local Insulation Score signal as a proxy:
        signal_b = mean of diagonal strip M[b-1:b+2, b-1:b+2]
                   (sum of contacts near the diagonal, 3-bin window)
        π_b = sigmoid(k_prior * (signal_b/mean_signal - 1))

    Empty regions (signal ≈ 0) → π_b ≈ sigmoid(-k_prior) ≈ small.
    Dense regions               → π_b → sigmoid(+k_prior).

    If matrix is None → uniform prior = global_prior.

    Parameters
    ----------
    matrix       : (n_bins, n_bins) Hi-C contact matrix or None.
    n_bins       : chromosome bin count.
    k_prior      : sigmoid steepness.
    global_prior : fallback prior when matrix unavailable.

    Returns
    -------
    prior : (n_bins,) array in (0, 1).
    """
    if matrix is None or matrix.size == 0:
        return np.full(n_bins, global_prior, dtype=np.float64)

    N = matrix.shape[0]
    if N != n_bins:
        logger.warning(
            "compute_prior: matrix size %d != n_bins %d — using global prior",
            N, n_bins,
        )
        return np.full(n_bins, global_prior, dtype=np.float64)

    # Local diagonal signal: sum of 3×3 neighbourhood around diagonal
    signal = np.zeros(N, dtype=np.float64)
    for b in range(N):
        lo = max(0, b - 1)
        hi = min(N, b + 2)
        signal[b] = matrix[lo:hi, lo:hi].sum()

    mean_signal = signal.mean()
    if mean_signal < 1e-12:
        return np.full(n_bins, global_prior, dtype=np.float64)

    # Normalised signal: 1.0 = average density
    norm = signal / mean_signal
    prior = _sigmoid(k_prior * (norm - 1.0))

    return prior.astype(np.float64)


# ---------------------------------------------------------------------------
# Observation matrix
# ---------------------------------------------------------------------------

def build_observation_matrix(
    algo_boundaries: dict[str, np.ndarray],
    candidate_bins: np.ndarray,
    tolerance_bins: int,
) -> np.ndarray:
    """Build binary observation matrix X of shape (n_candidates, n_algos).

    X[b, i] = 1 iff algorithm i has a boundary within tolerance_bins of
    candidate_bins[b].

    Parameters
    ----------
    algo_boundaries : {algo_name: sorted int array of boundary bin indices}
    candidate_bins  : (n_candidates,) sorted int array of candidate positions
    tolerance_bins  : matching radius

    Returns
    -------
    X : (n_candidates, n_algos) binary int8 array
    """
    n_cand  = len(candidate_bins)
    n_algos = len(algo_boundaries)
    X       = np.zeros((n_cand, n_algos), dtype=np.int8)

    for j, (_, bins) in enumerate(algo_boundaries.items()):
        if len(bins) == 0:
            continue
        # Vectorised: diff[i, k] = |candidate_bins[i] - bins[k]|
        diff = np.abs(candidate_bins[:, None] - bins[None, :])  # (C, K)
        X[:, j] = (diff <= tolerance_bins).any(axis=1).astype(np.int8)

    return X


# ---------------------------------------------------------------------------
# EM algorithm
# ---------------------------------------------------------------------------

def _log_likelihood(
    posteriors: np.ndarray,
    pi: np.ndarray,
    sens: np.ndarray,
    fpr: np.ndarray,
    X: np.ndarray,
) -> float:
    """Compute expected complete-data log-likelihood (for convergence check)."""
    # Avoid log(0)
    eps = 1e-12

    q = posteriors                          # (B,)
    # log P(x_i | z=1): X*log(s) + (1-X)*log(1-s)  summed over algos → (B,)
    log_p1 = (
        X    * np.log(sens[None, :]   + eps) +
        (1-X)* np.log(1-sens[None,:] + eps)
    ).sum(axis=1)
    # log P(x_i | z=0)
    log_p0 = (
        X    * np.log(fpr[None, :]   + eps) +
        (1-X)* np.log(1-fpr[None,:] + eps)
    ).sum(axis=1)

    ll = (
        q   * (np.log(pi + eps) + log_p1) +
        (1-q)* (np.log(1-pi + eps) + log_p0)
    ).sum()
    return float(ll)


def run_em(
    X: np.ndarray,
    pi: np.ndarray,
    sens_init: np.ndarray,
    fpr_init: np.ndarray,
    max_iter: int = 50,
    tol: float = 1e-4,
    fix_prior: bool = True,
    sens_min: float = 0.50,
    sens_max: float = 0.99,
    fpr_min:  float = 0.01,
    fpr_max:  float = 0.60,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[float]]:
    """Run EM to compute posterior boundary probabilities.

    Parameters
    ----------
    X          : (n_candidates, n_algos) binary observation matrix
    pi         : (n_candidates,) prior P(z_b=1)
    sens_init  : (n_algos,) initial sensitivity
    fpr_init   : (n_algos,) initial FPR
    max_iter   : maximum EM iterations
    tol        : log-likelihood convergence threshold
    fix_prior  : if True, do not update pi (recommended: prior is external)
    sens_min/max, fpr_min/max : clamp bounds for M-step updates

    Returns
    -------
    posteriors : (n_candidates,) P(z_b=1 | x_{·b})
    sens       : (n_algos,) final sensitivity estimates
    fpr        : (n_algos,) final FPR estimates
    ll_history : list of log-likelihoods per iteration
    """
    eps  = 1e-12
    sens = sens_init.copy()
    fpr  = fpr_init.copy()
    ll_history: list[float] = []

    for iteration in range(max_iter):

        # ── E-step: compute posteriors ───────────────────────────────────────
        # log P(x_{·b} | z=1) for each candidate bin b
        log_p1 = (
            X    * np.log(sens[None, :] + eps) +
            (1-X)* np.log(1 - sens[None,:] + eps)
        ).sum(axis=1)   # (B,)

        # log P(x_{·b} | z=0)
        log_p0 = (
            X    * np.log(fpr[None, :] + eps) +
            (1-X)* np.log(1 - fpr[None,:] + eps)
        ).sum(axis=1)   # (B,)

        # Numerically stable log-sum-exp
        log_num = np.log(pi + eps)   + log_p1
        log_den = np.log(1-pi + eps) + log_p0

        # P(z=1 | x) = 1 / (1 + exp(log_den - log_num))
        log_odds = log_num - log_den
        posteriors = _sigmoid(log_odds)   # (B,)

        # ── Convergence check ────────────────────────────────────────────────
        ll = _log_likelihood(posteriors, pi, sens, fpr, X)
        ll_history.append(ll)
        if iteration > 0 and abs(ll - ll_history[-2]) < tol:
            logger.debug(
                "EM converged at iteration %d (Δll=%.2e)",
                iteration, abs(ll - ll_history[-2]),
            )
            break

        # ── M-step: update sens and fpr ──────────────────────────────────────
        q  = posteriors           # (B,) P(z=1)
        q0 = 1.0 - posteriors     # (B,) P(z=0)

        # sens_i = E[x_i | z=1] = Σ_b q_b * x_{ib} / Σ_b q_b
        sum_q  = q.sum()
        sum_q0 = q0.sum()

        if sum_q > eps:
            sens = np.clip(
                (q[:, None] * X).sum(axis=0) / (sum_q + eps),
                sens_min, sens_max,
            )
        if sum_q0 > eps:
            fpr = np.clip(
                (q0[:, None] * X).sum(axis=0) / (sum_q0 + eps),
                fpr_min, fpr_max,
            )

        logger.debug(
            "EM iter %d | ll=%.4f | sens=[%.3f..%.3f] | fpr=[%.3f..%.3f]",
            iteration, ll,
            float(sens.min()), float(sens.max()),
            float(fpr.min()),  float(fpr.max()),
        )

    return posteriors, sens, fpr, ll_history


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def compute_probabilistic_consensus(
    tad_files: dict[str, Path],
    weights_meta: dict[str, dict],
    chrom: str,
    resolution: int,
    cfg: dict,
    matrix: Optional[np.ndarray] = None,
) -> pd.DataFrame:
    """Compute probabilistic boundary consensus.

    Parameters
    ----------
    tad_files    : {algo: path_to_bed}
    weights_meta : {algo: {"center_ratio": float, "n_boundaries": int}}
                   From chipseq_validation_summary.csv.
    chrom        : chromosome string (e.g. "chr1")
    resolution   : bin size in bp
    cfg          : config dict (reads cfg["probabilistic"] sub-dict)
    matrix       : optional Hi-C contact matrix for prior computation

    Returns
    -------
    pd.DataFrame
        Columns: chrom, start, end, posterior_prob
        Only bins with posterior >= theta.
        Empty DataFrame with same columns if no candidates.
    """
    _COLS  = ["chrom", "start", "end", "posterior_prob"]
    _EMPTY = pd.DataFrame(columns=_COLS)

    prob_cfg      = cfg.get("probabilistic", {})
    theta         = float(prob_cfg.get("theta",          0.5))
    tolerance_bins= int(cfg.get("consensus", {}).get("tolerance_bins", 1))
    k_sens        = float(prob_cfg.get("k_sens",         8.0))
    k_prior       = float(prob_cfg.get("k_prior",        3.0))
    global_prior  = float(prob_cfg.get("global_prior",   0.10))
    max_iter      = int(prob_cfg.get("max_iter",         50))
    em_tol        = float(prob_cfg.get("em_tol",         1e-4))
    fpr_min       = float(prob_cfg.get("fpr_min",        0.01))
    fpr_max       = float(prob_cfg.get("fpr_max",        0.60))
    sens_min      = float(prob_cfg.get("sens_min",       0.50))
    sens_max      = float(prob_cfg.get("sens_max",       0.99))
    exclude_from_sources = list(
        cfg.get("weighted_consensus", {}).get(
            "exclude_from_sources", ["weak_consensus", "strong_consensus"]
        )
    )

    # ── 0. Filter exclude_from_sources ───────────────────────────────────────
    tad_files    = {k: v for k, v in tad_files.items()
                   if k not in exclude_from_sources}
    weights_meta = {k: v for k, v in weights_meta.items()
                   if k not in exclude_from_sources}

    if not tad_files:
        logger.warning("[probabilistic] %s: no TAD files after filtering", chrom)
        return _EMPTY

    # ── 1. Load boundary bins per algorithm ──────────────────────────────────
    algo_names: list[str] = []
    algo_bins:  list[np.ndarray] = []
    center_ratios: list[float] = []
    n_boundaries_list: list[int] = []

    for algo, path in tad_files.items():
        path = Path(path)
        if not path.exists():
            logger.warning("[probabilistic] %s: missing BED %s", chrom, path)
            bins = np.array([], dtype=np.int64)
        else:
            try:
                df = pd.read_csv(
                    str(path), sep="\t", header=None,
                    usecols=[0, 1, 2], names=["chrom", "start", "end"],
                )
                df = df[df["chrom"] == chrom]
                starts = (df["start"].values // resolution).astype(np.int64)
                ends   = (df["end"].values   // resolution).astype(np.int64)
                bins   = np.unique(np.concatenate([starts, ends]))
            except Exception as exc:   # noqa: BLE001
                logger.warning("[probabilistic] %s: read error %s — %s",
                               chrom, path, exc)
                bins = np.array([], dtype=np.int64)

        meta = weights_meta.get(algo, {})
        cr   = float(meta.get("center_ratio", 1.0))
        nb   = int(meta.get("n_boundaries", len(bins) // 2))

        algo_names.append(algo)
        algo_bins.append(bins)
        center_ratios.append(cr)
        n_boundaries_list.append(nb)

        logger.debug(
            "[probabilistic] %s algo=%s cr=%.3f n_bnd=%d",
            chrom, algo, cr, nb,
        )

    # ── 2. Candidate bins ────────────────────────────────────────────────────
    all_arrays = [b for b in algo_bins if len(b) > 0]
    if not all_arrays:
        logger.info("[probabilistic] %s: no boundaries from any algorithm", chrom)
        return _EMPTY

    candidate_bins = np.unique(np.concatenate(all_arrays))
    n_cand = len(candidate_bins)

    # Chromosome size estimate for FPR
    if matrix is not None and matrix.size > 0:
        n_bins_chrom = matrix.shape[0]
    else:
        n_bins_chrom = max(
            int(candidate_bins.max()) + 1,
            max(n_boundaries_list) * 3,
        )

    logger.info(
        "[probabilistic] %s@%dbp: %d candidates, %d algos, n_bins_chrom=%d",
        chrom, resolution, n_cand, len(algo_names), n_bins_chrom,
    )

    # ── 3. Parameter initialisation ──────────────────────────────────────────
    cr_arr  = np.array(center_ratios,     dtype=np.float64)
    nb_arr  = np.array(n_boundaries_list, dtype=np.float64)

    sens_0 = init_sensitivity(cr_arr, k_sens=k_sens)
    fpr_0  = init_fpr(nb_arr, n_bins=n_bins_chrom,
                      fpr_min=fpr_min, fpr_max=fpr_max)

    # Prior π_b
    prior = compute_prior(matrix, n_cand, k_prior=k_prior,
                          global_prior=global_prior)
    # Note: compute_prior uses matrix shape; here we pass n_cand as n_bins
    # for the fallback path — when matrix is available, size mismatch is
    # handled inside compute_prior.
    if matrix is not None and matrix.size > 0:
        # Re-compute prior at candidate bin positions only
        full_prior = compute_prior(matrix, matrix.shape[0],
                                   k_prior=k_prior,
                                   global_prior=global_prior)
        # Clamp candidate bins to valid range
        safe_bins = np.clip(candidate_bins, 0, len(full_prior) - 1)
        prior = full_prior[safe_bins]

    logger.debug(
        "[probabilistic] %s: sens_init=[%.3f..%.3f] fpr_init=[%.3f..%.3f]",
        chrom,
        float(sens_0.min()), float(sens_0.max()),
        float(fpr_0.min()),  float(fpr_0.max()),
    )

    # ── 4. Observation matrix ────────────────────────────────────────────────
    algo_bins_dict = dict(zip(algo_names, algo_bins))
    X = build_observation_matrix(algo_bins_dict, candidate_bins, tolerance_bins)

    # ── 5. EM ────────────────────────────────────────────────────────────────
    try:
        posteriors, sens_final, fpr_final, ll_hist = run_em(
            X=X, pi=prior,
            sens_init=sens_0, fpr_init=fpr_0,
            max_iter=max_iter, tol=em_tol,
            sens_min=sens_min, sens_max=sens_max,
            fpr_min=fpr_min,   fpr_max=fpr_max,
        )
    except Exception as exc:   # noqa: BLE001
        logger.error("[probabilistic] %s: EM failed — %s", chrom, exc)
        return _EMPTY

    # ── 6. Log final parameter estimates ─────────────────────────────────────
    for i, algo in enumerate(algo_names):
        logger.info(
            "[probabilistic] %s | algo=%-20s | sens_0=%.3f→%.3f | fpr=%.3f→%.3f",
            chrom, algo,
            float(sens_0[i]),     float(sens_final[i]),
            float(fpr_0[i]),      float(fpr_final[i]),
        )
    logger.info(
        "[probabilistic] %s@%dbp: EM %d iters | ll=[%.2f..%.2f] | "
        "posterior mean=%.3f",
        chrom, resolution,
        len(ll_hist),
        min(ll_hist) if ll_hist else float("nan"),
        max(ll_hist) if ll_hist else float("nan"),
        float(posteriors.mean()),
    )

    # ── 7. Apply threshold ────────────────────────────────────────────────────
    mask           = posteriors >= theta
    selected_bins  = candidate_bins[mask]
    selected_post  = posteriors[mask]

    n_selected = int(mask.sum())
    logger.info(
        "[probabilistic] %s@%dbp: %d/%d boundaries retained (theta=%.2f)",
        chrom, resolution, n_selected, n_cand, theta,
    )

    if n_selected == 0:
        return _EMPTY

    result = pd.DataFrame({
        "chrom":          chrom,
        "start":          selected_bins * resolution,
        "end":            selected_bins * resolution + resolution,
        "posterior_prob": np.round(selected_post, 6),
    })
    return result.reset_index(drop=True)


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def save_probabilistic_consensus(
    df: pd.DataFrame,
    bed_path: Path,
    csv_path: Optional[Path] = None,
) -> None:
    """Save BED4 output (no header) and optional full CSV."""
    bed_path = Path(bed_path)
    bed_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(bed_path, sep="\t", header=False, index=False)
    logger.info("Saved probabilistic BED → %s (%d rows)", bed_path, len(df))

    if csv_path is not None:
        csv_path = Path(csv_path)
        csv_path.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(csv_path, index=False)
        logger.info("Saved probabilistic CSV → %s", csv_path)
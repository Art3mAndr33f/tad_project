"""
validation.py
=============
Биологическая валидация TAD-границ через CTCF ChIP-seq.

1. Enrichment score: observed / mean(random), 1000 пермутаций
2. Z-test p-value
3. Профиль обогащения ±500 kb (deepTools-style через numpy)
4. Анализ консенсусных границ (раздельно для support=2,3,4)
"""

from __future__ import annotations

import logging
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

RNG = np.random.default_rng(42)


# ──────────────────────────────────────────────────────────────────────────────
# Загрузка CTCF peaks
# ──────────────────────────────────────────────────────────────────────────────

def load_ctcf_peaks(bed_path: str) -> pd.DataFrame:
    """Загрузить CTCF BED-файл."""
    df = pd.read_csv(
        bed_path, sep="\t", header=None,
        names=["chrom", "start", "end"],
        usecols=[0, 1, 2],
        dtype={"start": int, "end": int},
    )
    if not df["chrom"].iloc[0].startswith("chr"):
        df["chrom"] = "chr" + df["chrom"].astype(str)
    return df


# ──────────────────────────────────────────────────────────────────────────────
# Окна вокруг границ
# ──────────────────────────────────────────────────────────────────────────────

def _boundary_windows(
    domains_df: pd.DataFrame,
    resolution: int,
    window_bins: int = 1,
) -> pd.DataFrame:
    """
    Создать интервалы (±window_bins * resolution) вокруг каждой границы TAD.
    """
    from src.consensus import extract_boundaries

    bnd = extract_boundaries(domains_df, resolution)
    half = window_bins * resolution
    windows = []
    chrom = domains_df["chrom"].iloc[0] if not domains_df.empty else "chr1"

    for pos in bnd:
        windows.append({
            "chrom": chrom,
            "start": max(0, int(pos) - half),
            "end":   int(pos) + half,
        })
    return pd.DataFrame(windows)


def _count_ctcf_overlaps(
    windows_df: pd.DataFrame,
    ctcf_df: pd.DataFrame,
) -> int:
    """
    Подсчитать число CTCF-пиков, перекрывающихся с хотя бы одним окном.
    Использует numpy без pybedtools.
    """
    if windows_df.empty or ctcf_df.empty:
        return 0

    chrom = windows_df["chrom"].iloc[0]
    ctcf_chrom = ctcf_df[ctcf_df["chrom"] == chrom]
    if ctcf_chrom.empty:
        return 0

    w_starts = windows_df["start"].values
    w_ends   = windows_df["end"].values
    c_starts = ctcf_chrom["start"].values
    c_ends   = ctcf_chrom["end"].values

    count = 0
    for cs, ce in zip(c_starts, c_ends):
        # Перекрытие: c_start < w_end AND c_end > w_start
        if np.any((cs < w_ends) & (ce > w_starts)):
            count += 1
    return count


# ──────────────────────────────────────────────────────────────────────────────
# Enrichment через пермутации
# ──────────────────────────────────────────────────────────────────────────────

def _random_windows(
    n_windows: int,
    window_size: int,
    chrom: str,
    chrom_size: int,
    seed: int = 42,
) -> pd.DataFrame:
    """Генерировать случайные окна той же длины на той же хромосоме."""
    rng = np.random.default_rng(seed)
    starts = rng.integers(0, chrom_size - window_size, size=n_windows)
    return pd.DataFrame({
        "chrom": chrom,
        "start": starts,
        "end":   starts + window_size,
    })


def compute_ctcf_enrichment(
    domains_df: pd.DataFrame,
    ctcf_df: pd.DataFrame,
    chrom: str,
    resolution: int,
    n_permutations: int = 1000,
    window_bins: int = 1,
    seed: int = 42,
) -> dict:
    """
    Вычислить обогащение CTCF на границах TAD.

    Returns
    -------
    dict: observed_count, mean_random, std_random, enrichment_score, p_value, z_score
    """
    from src.statistics import HG19_CHROM_SIZES
    chrom_size = HG19_CHROM_SIZES.get(chrom, 250_000_000)
    window_size = 2 * window_bins * resolution

    if domains_df.empty:
        return {
            "observed_count":  0,
            "mean_random":     np.nan,
            "std_random":      np.nan,
            "enrichment_score": np.nan,
            "p_value":         np.nan,
            "z_score":         np.nan,
        }

    windows_df = _boundary_windows(domains_df, resolution, window_bins)
    n_windows  = len(windows_df)

    if n_windows == 0:
        return {
            "observed_count":  0,
            "mean_random":     np.nan,
            "std_random":      np.nan,
            "enrichment_score": np.nan,
            "p_value":         np.nan,
            "z_score":         np.nan,
        }

    observed = _count_ctcf_overlaps(windows_df, ctcf_df)

    # Пермутации
    random_counts = []
    for perm_i in range(n_permutations):
        rand_windows = _random_windows(n_windows, window_size, chrom, chrom_size,
                                       seed=seed + perm_i)
        cnt = _count_ctcf_overlaps(rand_windows, ctcf_df)
        random_counts.append(cnt)

    random_arr = np.array(random_counts, dtype=float)
    mean_rand  = float(np.mean(random_arr))
    std_rand   = float(np.std(random_arr))

    enrichment = float(observed / (mean_rand + 1e-9))
    z_score    = float((observed - mean_rand) / (std_rand + 1e-9))
    # Односторонний z-test (H1: обогащение)
    from scipy.stats import norm
    p_value = float(1 - norm.cdf(z_score))

    return {
        "observed_count":   int(observed),
        "mean_random":      mean_rand,
        "std_random":       std_rand,
        "enrichment_score": enrichment,
        "p_value":          p_value,
        "z_score":          z_score,
    }


# ──────────────────────────────────────────────────────────────────────────────
# Профиль обогащения CTCF (±500 kb)
# ──────────────────────────────────────────────────────────────────────────────

def compute_ctcf_profile(
    domains_df: pd.DataFrame,
    ctcf_df: pd.DataFrame,
    chrom: str,
    resolution: int,
    profile_range_bp: int = 500_000,
    profile_bin_bp:   int = 10_000,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Построить профиль CTCF-обогащения вокруг TAD-границ.

    Parameters
    ----------
    profile_range_bp : диапазон ±N bp от границы
    profile_bin_bp   : размер бина профиля

    Returns
    -------
    (bin_centers, mean_ctcf_density) — оба 1D array
    """
    from src.statistics import HG19_CHROM_SIZES
    from src.consensus import extract_boundaries

    chrom_size = HG19_CHROM_SIZES.get(chrom, 250_000_000)
    n_bins     = 2 * profile_range_bp // profile_bin_bp
    bin_edges  = np.linspace(-profile_range_bp, profile_range_bp, n_bins + 1)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    ctcf_chrom = ctcf_df[ctcf_df["chrom"] == chrom].copy()
    ctcf_midpoints = ((ctcf_chrom["start"] + ctcf_chrom["end"]) // 2).values

    if len(ctcf_midpoints) == 0:
        return bin_centers, np.zeros(n_bins)

    boundaries = extract_boundaries(domains_df, resolution)
    if len(boundaries) == 0:
        return bin_centers, np.zeros(n_bins)

    # Для каждой границы: смещения CTCF-пиков
    all_offsets = []
    for bnd in boundaries:
        if bnd < profile_range_bp or bnd > chrom_size - profile_range_bp:
            continue
        offsets = ctcf_midpoints - bnd
        mask    = (offsets >= -profile_range_bp) & (offsets < profile_range_bp)
        all_offsets.extend(offsets[mask].tolist())

    if not all_offsets:
        return bin_centers, np.zeros(n_bins)

    offsets_arr = np.array(all_offsets)
    counts, _   = np.histogram(offsets_arr, bins=bin_edges)

    # Нормализация на число границ и размер бина
    density = counts / (len(boundaries) * profile_bin_bp / 1000.0)
    return bin_centers, density


# ──────────────────────────────────────────────────────────────────────────────
# Валидация для консенсусных границ
# ──────────────────────────────────────────────────────────────────────────────

def validate_consensus_boundaries(
    consensus_df: pd.DataFrame,
    ctcf_df: pd.DataFrame,
    chrom: str,
    resolution: int,
    n_permutations: int = 1000,
    seed: int = 42,
) -> Dict[int, dict]:
    """
    Отдельный CTCF-enrichment для консенсусных границ с support=2,3,4.

    Returns
    -------
    {support_level: enrichment_dict}
    """
    if consensus_df.empty:
        return {}

    results = {}
    for support_level in [2, 3, 4]:
        subset = consensus_df[consensus_df["support"] == support_level].copy()
        if subset.empty:
            continue

        # Построить pseudo-domains из позиций
        records = []
        for _, row in subset.iterrows():
            pos = int(row["position"])
            records.append({"chrom": chrom, "start": pos, "end": pos + resolution})
        pseudo_df = pd.DataFrame(records)

        enrich = compute_ctcf_enrichment(
            pseudo_df, ctcf_df, chrom, resolution, n_permutations, seed=seed
        )
        results[support_level] = enrich

    return results


# ──────────────────────────────────────────────────────────────────────────────
# Batch-валидация
# ──────────────────────────────────────────────────────────────────────────────

def run_all_validation(
    all_results: Dict,
    consensus_all: Dict,
    cfg: dict,
) -> pd.DataFrame:
    """
    Запустить CTCF-валидацию для всех алгоритмов/хромосом/разрешений.
    Возвращает DataFrame с метриками обогащения.
    """
    ctcf_df        = load_ctcf_peaks(cfg["paths"]["ctcf_bed"])
    n_perm         = cfg["validation"]["n_permutations"]
    window_bins    = cfg["validation"]["ctcf_window_bp"]
    resolutions    = cfg["resolutions"]
    chromosomes    = cfg["chromosomes"]["all"]

    rows = []
    for res in resolutions:
        for chrom in chromosomes:
            # Валидация алгоритмов
            for algo, chrom_dict in all_results.items():
                df = chrom_dict.get(chrom, {}).get(res)
                if df is None or df.empty:
                    continue
                enrich = compute_ctcf_enrichment(
                    df, ctcf_df, chrom, res, n_perm, window_bins, seed=42
                )
                rows.append({
                    "type":       "algorithm",
                    "algorithm":  algo,
                    "support":    None,
                    "chrom":      chrom,
                    "resolution": res,
                    **enrich,
                })

            # Валидация консенсусных границ
            cons_df = consensus_all.get(chrom, {}).get(res)
            if cons_df is not None and not cons_df.empty:
                cons_results = validate_consensus_boundaries(
                    cons_df, ctcf_df, chrom, res, n_perm, seed=42
                )
                for sup, enrich in cons_results.items():
                    rows.append({
                        "type":       "consensus",
                        "algorithm":  None,
                        "support":    sup,
                        "chrom":      chrom,
                        "resolution": res,
                        **enrich,
                    })

    return pd.DataFrame(rows)
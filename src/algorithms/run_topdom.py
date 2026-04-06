"""
run_topdom.py
=============
Python-порт алгоритма TopDom.

Оригинальная статья: Shin et al. 2016, Nucleic Acids Research.

Алгоритм:
  1. Для каждого бина вычисляется bin-score — среднее взаимодействие
     с соседями в окне w (upstream vs downstream triangle).
  2. Локальные минимумы bin-score — кандидаты в границы TAD.
  3. Статистическая фильтрация через mean±sd (non-gap window).
  4. Из границ формируются домены.

Выбор окна: запускается для w ∈ [3, 5, 10], выбирается w с наибольшим
средним внутридоменным / межддоменным соотношением.
"""

from __future__ import annotations

import logging
from typing import Optional

import numpy as np
import pandas as pd
from scipy.signal import argrelmin

from src.data_prep import get_matrix

logger = logging.getLogger(__name__)


# ──────────────────────────────────────────────────────────────────────────────
# Ядро TopDom
# ──────────────────────────────────────────────────────────────────────────────

def _compute_bin_scores(matrix: np.ndarray, w: int) -> np.ndarray:
    """
    Вычислить bin-score для каждого бина.
    Bin-score[i] = mean(upstream_triangle) - mean(downstream_triangle)
    upstream:   matrix[i-w:i, i-w:i]
    downstream: matrix[i+1:i+w+1, i+1:i+w+1]
    """
    n = matrix.shape[0]
    scores = np.full(n, np.nan, dtype=np.float64)

    for i in range(w, n - w):
        up   = matrix[max(0, i-w):i,   max(0, i-w):i]
        down = matrix[i+1:i+w+1,       i+1:i+w+1]

        mean_up   = np.nanmean(up)   if up.size   > 0 else np.nan
        mean_down = np.nanmean(down) if down.size > 0 else np.nan

        if np.isnan(mean_up) or np.isnan(mean_down):
            scores[i] = np.nan
        else:
            scores[i] = mean_up - mean_down

    return scores


def _find_boundaries(scores: np.ndarray, resolution: int) -> np.ndarray:
    """
    Найти границы TAD как локальные минимумы bin-score
    (фильтрация по mean-sd).
    """
    valid = ~np.isnan(scores)
    if valid.sum() < 3:
        return np.array([], dtype=int)

    mean_s = np.nanmean(scores[valid])
    std_s  = np.nanstd(scores[valid])
    threshold = mean_s - std_s

    # Локальные минимумы
    local_min_idx = argrelmin(np.nan_to_num(scores, nan=0.0), order=1)[0]
    # Фильтр: значение ниже порога
    boundaries = local_min_idx[scores[local_min_idx] < threshold]
    return boundaries


def _boundaries_to_domains(
    boundaries: np.ndarray,
    n_bins: int,
    resolution: int,
    chrom: str,
) -> pd.DataFrame:
    """Преобразовать индексы границ в домены (start, end)."""
    bnd = sorted(set([0] + list(boundaries) + [n_bins - 1]))
    records = []
    for i in range(len(bnd) - 1):
        start = bnd[i] * resolution
        end   = bnd[i + 1] * resolution
        if end > start:
            records.append((chrom, start, end))
    if not records:
        return pd.DataFrame(columns=["chrom", "start", "end"])
    return pd.DataFrame(records, columns=["chrom", "start", "end"])


def _intra_inter_ratio(matrix: np.ndarray, domains_df: pd.DataFrame,
                       resolution: int) -> float:
    """
    Критерий качества окна: среднее внутридоменное / межддоменное взаимодействие.
    """
    if domains_df.empty:
        return 0.0

    intra_vals, inter_vals = [], []
    n = matrix.shape[0]

    for _, row in domains_df.iterrows():
        i0 = int(row["start"] // resolution)
        i1 = int(row["end"]   // resolution)
        i0 = max(0, min(i0, n - 1))
        i1 = max(0, min(i1, n - 1))
        if i1 <= i0:
            continue
        block = matrix[i0:i1, i0:i1]
        intra_vals.append(np.nanmean(block))

    # Межддоменные: берём часть off-diagonal
    rows = domains_df.sample(min(10, len(domains_df)), random_state=42)
    for _, r1 in rows.iterrows():
        for _, r2 in domains_df.iterrows():
            if r1["start"] == r2["start"]:
                continue
            s1, e1 = int(r1["start"] // resolution), int(r1["end"] // resolution)
            s2, e2 = int(r2["start"] // resolution), int(r2["end"] // resolution)
            s1, e1 = max(0, min(s1, n-1)), max(0, min(e1, n-1))
            s2, e2 = max(0, min(s2, n-1)), max(0, min(e2, n-1))
            if e1 > s1 and e2 > s2:
                block = matrix[s1:e1, s2:e2]
                inter_vals.append(np.nanmean(block))

    mean_intra = np.nanmean(intra_vals) if intra_vals else 0.0
    mean_inter = np.nanmean(inter_vals) if inter_vals else 1e-9
    return float(mean_intra / (mean_inter + 1e-9))


def _run_topdom_single_window(
    matrix: np.ndarray,
    chrom: str,
    resolution: int,
    w: int,
) -> pd.DataFrame:
    """Запуск TopDom для одного размера окна."""
    scores     = _compute_bin_scores(matrix, w)
    boundaries = _find_boundaries(scores, resolution)
    domains    = _boundaries_to_domains(boundaries, matrix.shape[0], resolution, chrom)
    return domains


# ──────────────────────────────────────────────────────────────────────────────
# Публичный интерфейс
# ──────────────────────────────────────────────────────────────────────────────

def run_topdom(
    chrom: str,
    resolution: int,
    data_path: str,
    cfg: Optional[dict] = None,
    window_sizes: Optional[list[int]] = None,
    matrix: Optional[np.ndarray] = None,
) -> pd.DataFrame:
    """
    Запустить TopDom с автоматическим выбором окна.

    Parameters
    ----------
    chrom        : хромосома
    resolution   : разрешение в bp
    data_path    : путь к директории processed
    cfg          : конфиг-словарь
    window_sizes : список w для перебора
    matrix       : можно передать напрямую (обходит загрузку)

    Returns
    -------
    pd.DataFrame с колонками chrom, start, end
    """
    if cfg is not None:
        window_sizes = cfg["algorithms"]["topdom"]["window_sizes"]

    window_sizes = window_sizes or [3, 5, 10]

    # Загрузить матрицу
    if matrix is None:
        if cfg is not None:
            matrix = get_matrix(cfg, chrom, resolution)
        else:
            import os
            npy = os.path.join(data_path, f"{chrom}_{resolution}bp.npy")
            if not os.path.exists(npy):
                raise FileNotFoundError(f"Матрица не найдена: {npy}")
            matrix = np.load(npy)

    # Заменить NaN/Inf
    matrix = np.nan_to_num(matrix.astype(np.float64))

    logger.info("[TopDom] %s @ %d bp  | windows=%s", chrom, resolution, window_sizes)

    best_w      = window_sizes[0]
    best_ratio  = -1.0
    best_domains: pd.DataFrame = pd.DataFrame(columns=["chrom", "start", "end"])

    for w in window_sizes:
        domains = _run_topdom_single_window(matrix, chrom, resolution, w)
        ratio   = _intra_inter_ratio(matrix, domains, resolution)
        logger.info("  w=%d → %d TADs | ratio=%.3f", w, len(domains), ratio)
        if ratio > best_ratio:
            best_ratio  = ratio
            best_w      = w
            best_domains = domains

    logger.info("[TopDom] Итого %d TAD (w=%d)", len(best_domains), best_w)
    return best_domains.reset_index(drop=True)
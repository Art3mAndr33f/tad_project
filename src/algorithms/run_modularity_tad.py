"""
src/algorithms/run_modularity_tad.py

TAD detection via Hi-C contact graph modularity maximization.

Метод: 1D-адаптация Newman modularity Q для Hi-C контактных матриц.
  Q_tad(i,j) = mean_B(i,j) = mean[ A[a,b] - k[a]*k[b]/(2m) ]
  для всех пар (a,b) внутри TAD [i,j) с |a-b| <= max_dist_bins.

Нормировка на число пар обеспечивает сравнимость score для TAD
любого размера. Auto-penalty находит диапазон 1–2.5 TAD/Mb.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

_DEFAULT_MAX_DIST_MB = 5.0
_DEFAULT_MIN_TAD_KB  = 100.0
_DEFAULT_MAX_TAD_KB  = 3000.0   # уменьшили с 5000 → меньше гигантских TAD
_SEED = 42

# Сетка penalty от строгого к мягкому
_PENALTY_GRID = [2.0, 1.0, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.001]


# ═══════════════════════════════════════════════════════════════════════════════
def run_modularity_tad(
    chrom: str,
    resolution: int,
    data_path: str,
    cfg: Optional[dict],
    **kwargs,
) -> pd.DataFrame:
    empty = pd.DataFrame(columns=["chrom", "start", "end"])

    # ── 1. Загрузить матрицу ─────────────────────────────────────────────────
    try:
        matrix = _load_matrix(data_path, chrom, resolution)
    except FileNotFoundError as exc:
        logger.error("[modtad] Матрица не найдена %s @ %d: %s", chrom, resolution, exc)
        return empty

    n = matrix.shape[0]
    logger.info("[modtad] %s @ %d — матрица %d×%d загружена", chrom, resolution, n, n)

    if n < 10:
        logger.warning("[modtad] %s @ %d — слишком мало бинов (%d)", chrom, resolution, n)
        return empty

    # ── 2. Параметры из конфига ──────────────────────────────────────────────
    algo_cfg     = (cfg or {}).get("algorithms", {}).get("modularity_tad", {})
    max_dist_mb  = float(algo_cfg.get("max_dist_mb", _DEFAULT_MAX_DIST_MB))
    min_tad_kb   = float(algo_cfg.get("min_tad_kb",  _DEFAULT_MIN_TAD_KB))
    max_tad_kb   = float(algo_cfg.get("max_tad_kb",  _DEFAULT_MAX_TAD_KB))
    penalty_cfg  = algo_cfg.get("penalty", None)   # None = auto

    max_dist_bins = max(5, int(max_dist_mb * 1e6 / resolution))
    min_tad_bins  = max(2, int(min_tad_kb  * 1e3 / resolution))
    max_tad_bins  = max(min_tad_bins + 1, int(max_tad_kb * 1e3 / resolution))

    # Целевой диапазон: 1.0–2.5 TAD/Mb (эмпирика GM12878)
    chrom_mb   = n * resolution / 1e6
    target_min = max(5,  int(chrom_mb * 1.0))
    target_max = max(10, int(chrom_mb * 2.5))

    logger.debug(
        "[modtad] %s @ %d: n=%d, %.1f Mb, target=%d–%d TADs, "
        "max_dist=%d bins, min_tad=%d bins, max_tad=%d bins",
        chrom, resolution, n, chrom_mb, target_min, target_max,
        max_dist_bins, min_tad_bins, max_tad_bins,
    )

    # ── 3. Модулярная матрица B (полосовой формат) ───────────────────────────
    B_band = _compute_B_band(matrix, max_dist_bins)

    # ── 4. TAD-scores: нормированные prefix-суммы ────────────────────────────
    score_compact = _compute_score_compact(B_band, n, max_dist_bins, max_tad_bins)

    # ── 4b. Нормировка score → [0, 1] ────────────────────────────────────────
    # Делим на максимальный score чтобы penalty_grid был независим от
    # абсолютных значений матрицы (RAW counts сильно варьируют между хромосомами)
    sc_max = float(np.max(score_compact))
    if sc_max > 0:
        score_compact = (score_compact / sc_max).astype(np.float32)
        logger.debug("[modtad] Score нормирован на %.2f → диапазон [0, 1]", sc_max)

    # ── Диагностика score-диапазона ──────────────────────────────────────────
    nonzero_scores = score_compact[score_compact != 0]
    if len(nonzero_scores):
        logger.debug(
            "[modtad] Score диапазон: min=%.4f median=%.4f max=%.4f",
            float(np.min(nonzero_scores)),
            float(np.median(nonzero_scores)),
            float(np.max(nonzero_scores)),
        )

    # ── 5. Auto-penalty DP ───────────────────────────────────────────────────
    penalties = [float(penalty_cfg)] if penalty_cfg is not None else _PENALTY_GRID

    best_df  = empty
    best_n   = 0
    best_pen = penalties[0]

    for pen in penalties:
        boundaries = _dp_segment(score_compact, n, pen, min_tad_bins, max_tad_bins)
        n_tads     = len(boundaries)
        logger.debug("[modtad] penalty=%.4f → %d TADs", pen, n_tads)

        if target_min <= n_tads <= target_max:
            best_pen = pen
            best_df  = _to_df(boundaries, chrom, resolution)
            logger.info(
                "[modtad] %s @ %d — penalty=%.4f → %d TADs ✓ (target %d–%d)",
                chrom, resolution, pen, n_tads, target_min, target_max,
            )
            break

        # Запомнить лучший пока не нашли в диапазоне
        if target_min <= n_tads * 2 and n_tads > best_n:
            best_n   = n_tads
            best_pen = pen
            best_df  = _to_df(boundaries, chrom, resolution)

    if len(best_df) < 3:
        logger.warning(
            "[modtad] %s @ %d — не попали в target (%d–%d), "
            "возвращаем %d TADs @ penalty=%.4f",
            chrom, resolution, target_min, target_max,
            len(best_df), best_pen,
        )

    logger.info("[modtad] %s @ %d — итого %d TADs", chrom, resolution, len(best_df))
    return best_df


# ═══════════════════════════════════════════════════════════════════════════════
# Вспомогательные функции
# ═══════════════════════════════════════════════════════════════════════════════

def _load_matrix(data_path: str, chrom: str, resolution: int) -> np.ndarray:
    npy = Path(data_path) / f"{chrom}_{resolution}bp.npy"
    if not npy.exists():
        raise FileNotFoundError(npy)
    m = np.load(str(npy)).astype(np.float64)
    m = np.nan_to_num(m, nan=0.0, posinf=0.0, neginf=0.0)
    return np.maximum(m, m.T)


def _compute_B_band(matrix: np.ndarray, max_dist_bins: int) -> np.ndarray:
    """
    B_band[i, d] = A[i, i+d] - k[i]*k[i+d] / (2m)
    где k[i] = локальная степень (сумма контактов до max_dist_bins).

    shape: (n, max_dist_bins+1)
    """
    n = matrix.shape[0]

    # Локальные степени
    k = np.zeros(n, dtype=np.float64)
    for i in range(n):
        lo = max(0, i - max_dist_bins)
        hi = min(n, i + max_dist_bins + 1)
        k[i] = matrix[i, lo:hi].sum() - matrix[i, i]

    two_m = max(k.sum(), 1.0)

    # Полосовой лапласиан
    B_band = np.zeros((n, max_dist_bins + 1), dtype=np.float64)
    for d in range(1, max_dist_bins + 1):
        diag    = np.diagonal(matrix, offset=d)   # shape (n-d,) — всегда правильно
        B_band[:n - d, d] = diag - k[:n - d] * k[d:n] / two_m

    return B_band


def _compute_score_compact(
    B_band: np.ndarray,
    n: int,
    max_dist_bins: int,
    max_tad_bins: int,
) -> np.ndarray:
    """
    score_compact[i, l] = mean B[a,b] для всех пар (a,b) с i<=a<b<i+l, b-a<=max_dist.
    Нормировка на реальное число пар делает score сравнимым для TAD любого размера.
    shape: (n, max_tad_bins+1), dtype float32
    """
    prefix = np.zeros((n + 1, max_dist_bins + 1), dtype=np.float64)
    for d in range(1, max_dist_bins + 1):
        prefix[1:, d] = np.cumsum(B_band[:, d])

    score_compact = np.zeros((n, max_tad_bins + 1), dtype=np.float32)

    for l in range(1, max_tad_bins + 1):
        # n_pairs = число пар (a, a+d) с 0<=d<=min(l-1, max_dist): = Σ max(0, l-d)
        n_pairs = sum(l - d for d in range(1, min(l, max_dist_bins) + 1))
        if n_pairs == 0:
            continue
        norm  = float(n_pairs)
        max_i = n - l
        if max_i <= 0:
            continue

        i_arr = np.arange(max_i, dtype=np.int32)
        total = np.zeros(max_i, dtype=np.float64)

        for d in range(1, min(l, max_dist_bins) + 1):
            # Пар с этой диагональю внутри [i, i+l): a ∈ [i, i+l-d)
            # Вклад: prefix[i+l-d] - prefix[i]  (но не дальше n)
            end_idx = np.minimum(i_arr + (l - d), n).astype(np.int32)
            total  += prefix[end_idx, d] - prefix[i_arr, d]

        score_compact[:max_i, l] = (total / norm).astype(np.float32)

    return score_compact


def _dp_segment(
    score_compact: np.ndarray,
    n: int,
    penalty: float,
    min_tad_bins: int,
    max_tad_bins: int,
) -> list[tuple[int, int]]:
    """
    DP: найти разбиение [0, n) максимизирующее sum(score) - penalty * n_tads.
    """
    dp     = np.full(n + 1, -1e18, dtype=np.float64)
    parent = np.full(n + 1, -1,    dtype=np.int32)
    dp[0]  = 0.0

    sc_cols = score_compact.shape[1]

    for j in range(min_tad_bins, n + 1):
        i_min = max(0, j - max_tad_bins)
        i_max = j - min_tad_bins + 1

        i_arr = np.arange(i_min, i_max, dtype=np.int32)
        l_arr = j - i_arr

        valid = (l_arr >= min_tad_bins) & (l_arr <= max_tad_bins) & (l_arr < sc_cols)
        i_arr = i_arr[valid]
        l_arr = l_arr[valid]
        if len(i_arr) == 0:
            continue

        dp_prev    = dp[i_arr]
        reachable  = dp_prev > -1e17
        if not reachable.any():
            continue

        q_vals     = score_compact[i_arr, l_arr].astype(np.float64)
        candidates = np.where(reachable, dp_prev + q_vals - penalty, -1e18)

        best = int(np.argmax(candidates))
        if candidates[best] > dp[j]:
            dp[j]     = candidates[best]
            parent[j] = int(i_arr[best])

    # Traceback
    tads: list[tuple[int, int]] = []
    pos = n
    while pos > 0:
        prev = int(parent[pos])
        if prev < 0:
            logger.warning("[modtad] DP traceback прерван на pos=%d", pos)
            return [(0, n)]
        tads.append((prev, pos))
        pos = prev
    tads.reverse()
    return tads


def _to_df(
    boundaries: list[tuple[int, int]],
    chrom: str,
    resolution: int,
) -> pd.DataFrame:
    records = [
        {"chrom": chrom, "start": s * resolution, "end": e * resolution}
        for s, e in boundaries
    ]
    return pd.DataFrame(records, columns=["chrom", "start", "end"])

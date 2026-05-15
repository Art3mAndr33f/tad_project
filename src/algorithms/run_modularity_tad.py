# src/algorithms/run_modularity_tad.py
"""ModularityTAD: DP-сегментация Hi-C на основе intra/inter-ratio score.

score(i, j) = mean_OE_intra(i,j) − mean_OE_flanks(i,j)

O/E нормализация (наблюдаемое / ожидаемое по диагонали) убирает
distance-decay без bias в сторону малых или крупных TAD.

v2.2: O/E нормализация вместо log1p.
"""

from __future__ import annotations

import logging
from typing import Optional

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


# ──────────────────────────────────────────────────────────────────────────────
# O/E нормализация
# ──────────────────────────────────────────────────────────────────────────────

def _oe_normalize(A: np.ndarray) -> np.ndarray:
    """Observed / Expected нормализация по диагоналям.

    Для каждого расстояния d: делим все контакты на среднее по диагонали d.
    После нормировки среднее значение на каждой диагонали = 1.
    Это полностью убирает distance-decay, не создавая bias ни к малым,
    ни к крупным TAD.
    """
    n   = A.shape[0]
    out = A.copy().astype(np.float64)

    for d in range(n):
        idx  = np.arange(n - d)
        diag = A[idx, idx + d]

        # Среднее только ненулевых контактов (нулевые = пропущенные данные)
        nz   = diag[diag > 0]
        mean_d = nz.mean() if nz.size > 0 else 0.0

        if mean_d > 0:
            out[idx,     idx + d] = diag / mean_d
            if d > 0:
                out[idx + d, idx    ] = A[idx + d, idx] / mean_d

    return out


# ──────────────────────────────────────────────────────────────────────────────
# 2D Prefix-суммы
# ──────────────────────────────────────────────────────────────────────────────

def _build_prefix2d(A: np.ndarray) -> np.ndarray:
    """P[i, j] = sum(A[0:i, 0:j]).  Размер (n+1)×(n+1)."""
    n = A.shape[0]
    P = np.zeros((n + 1, n + 1), dtype=np.float64)
    P[1:, 1:] = np.cumsum(np.cumsum(A, axis=0), axis=1)
    return P


def _rect_sum(P: np.ndarray, r0: int, c0: int, r1: int, c1: int) -> float:
    """Сумма A[r0:r1, c0:c1] через prefix-суммы."""
    if r0 >= r1 or c0 >= c1:
        return 0.0
    return float(P[r1, c1] - P[r0, c1] - P[r1, c0] + P[r0, c0])


# ──────────────────────────────────────────────────────────────────────────────
# Score: OE-intra − OE-flanks
# ──────────────────────────────────────────────────────────────────────────────

def _score_intra_minus_flanks(
    P: np.ndarray,
    i: int,
    j: int,
    n: int,
) -> float:
    """score(i,j) = mean_OE_intra(i,j) − mean_OE_flanks(i,j).

    TAD = [i, j)  (0-based, полуоткрытый).
    flanks: [max(0,i−size), i) × [i,j)  и  [i,j) × [j, min(n,j+size)).
    """
    size = j - i
    if size <= 0:
        return -np.inf

    intra_sum  = _rect_sum(P, i, i, j, j)
    intra_mean = intra_sum / (size * size)

    lf_start = max(0, i - size)
    rf_end   = min(n, j + size)

    left_sum  = _rect_sum(P, lf_start, i, i, j)
    right_sum = _rect_sum(P, i, j, j, rf_end)

    left_n  = (i - lf_start) * size
    right_n = size * (rf_end - j)
    flank_n = left_n + right_n

    flank_mean = (left_sum + right_sum) / flank_n if flank_n > 0 else 0.0

    return intra_mean - flank_mean


# ──────────────────────────────────────────────────────────────────────────────
# DP-сегментация
# ──────────────────────────────────────────────────────────────────────────────

def _dp_segment(
    score_matrix: np.ndarray,
    penalty: float,
    min_tad_bins: int,
    max_tad_bins: int,
    n: int,
) -> list[tuple[int, int]]:
    """DP: обязательное полное покрытие [0, n) TAD-сегментами.

    Ключевые отличия от предыдущей версии:
    - Нет pass-through: dp[j] недостижимо если нет валидного TAD [i,j).
    - Нет merge-шага: каждый сегмент в backtrace — настоящий TAD.
    - Если n не достижимо — принудительно расширяем последний TAD до n.

    dp[j] = max суммарный score для полного покрытия [0, j).
    back[j] = начало последнего TAD, заканчивающегося в j.
    """
    NEG_INF = -1e18
    dp   = np.full(n + 1, NEG_INF, dtype=np.float64)
    back = np.full(n + 1, -1,      dtype=np.int32)
    dp[0] = 0.0

    for j in range(min_tad_bins, n + 1):
        i_lo = max(0, j - max_tad_bins)
        i_hi = j - min_tad_bins          # включительно

        for i in range(i_lo, i_hi + 1):
            if dp[i] < NEG_INF / 2:
                continue
            s = float(score_matrix[i, j - 1])
            if s <= NEG_INF / 2:         # невалидная ячейка
                continue
            val = dp[i] + s - penalty
            if val > dp[j]:
                dp[j]   = val
                back[j] = i

    # Если n недостижимо — расширить последний достижимый TAD до n
    if dp[n] < NEG_INF / 2 or back[n] < 0:
        # Найти ближайший достижимый j*, расширить [back[j*], n)
        for j_last in range(n - 1, min_tad_bins - 1, -1):
            if dp[j_last] > NEG_INF / 2 and back[j_last] >= 0:
                back[n] = back[j_last]
                dp[n]   = dp[j_last]    # score не пересчитываем
                # Обрезать backtrace: j_last больше не финальная точка
                # Перестроить: вместо (back[j_last], j_last) сделать (back[j_last], n)
                # Достаточно просто выставить back[n] и оборвать цепочку на j_last
                # путём установки back[j_last] в "занятое" значение
                dp[j_last] = NEG_INF    # исключить j_last как конечную точку
                break
        else:
            # Полный fallback: один TAD на весь хром
            back[n] = 0
            dp[n]   = float(score_matrix[0, n - 1]) - penalty

    # Backtrace — только реальные TAD-сегменты
    segments: list[tuple[int, int]] = []
    j = n
    while j > 0:
        i = int(back[j])
        if i < 0:
            # Докрыть остаток одним сегментом
            segments.append((0, j))
            break
        segments.append((i, j))
        j = i
        if j == 0:
            break

    segments.reverse()
    return segments

# ──────────────────────────────────────────────────────────────────────────────
# Auto-penalty sweep
# ──────────────────────────────────────────────────────────────────────────────

def _auto_penalty_sweep(
    score_matrix: np.ndarray,
    min_tad_bins: int,
    max_tad_bins: int,
    n: int,
    chrom_mb: float,
) -> tuple[float, list[tuple[int, int]]]:
    """Подобрать penalty: target 0.8–2.5 TAD/Mb.

    После O/E нормализации score(хороший TAD) ≈ 0.3–1.5,
    поэтому sweep идёт в диапазоне [0.01, 2.0].
    """
    target_lo  = max(5,   int(chrom_mb * 0.8))
    target_hi  = min(200, int(chrom_mb * 2.5))
    target_mid = (target_lo + target_hi) / 2.0

    # Медиана положительных score как anchor
    valid  = score_matrix[score_matrix > 0]
    median = float(np.median(valid)) if valid.size > 0 else 1.0
    logger.info("[modtad] O/E score stats: median=%.4f, max=%.4f",
                median, float(score_matrix.max()))

    # Абсолютные penalties (O/E score не зависит от данных, только от структуры)
    penalties = [2.0, 1.5, 1.0, 0.8, 0.6, 0.5, 0.4, 0.3, 0.2, 0.15, 0.1,
                 0.08, 0.05, 0.03, 0.01]

    best_pen  = penalties[len(penalties) // 2]
    best_segs: list[tuple[int, int]] = []
    best_dist = np.inf

    for pen in penalties:
        segs = _dp_segment(score_matrix, pen, min_tad_bins, max_tad_bins, n)
        n_t  = len(segs)
        dist = abs(n_t - target_mid)
        logger.info("[modtad] penalty=%.3f → %d TADs (target %d–%d)",
                    pen, n_t, target_lo, target_hi)

        if target_lo <= n_t <= target_hi:
            logger.info("[modtad] ✓ penalty=%.3f hits target: %d TADs", pen, n_t)
            return pen, segs

        if dist < best_dist:
            best_dist = dist
            best_pen  = pen
            best_segs = segs

    logger.warning("[modtad] sweep done, best=%d TADs (target %d–%d)",
                   len(best_segs), target_lo, target_hi)
    return best_pen, best_segs


# ──────────────────────────────────────────────────────────────────────────────
# Главная функция
# ──────────────────────────────────────────────────────────────────────────────

def run_modularity_tad(
    chrom: str,
    resolution: int,
    data_path: str,
    cfg: Optional[dict] = None,
    **kwargs,
) -> pd.DataFrame:
    """Запустить ModularityTAD для одной хромосомы."""
    _empty = pd.DataFrame(columns=["chrom", "start", "end"])

    algo_cfg = (cfg or {}).get("algorithms", {}).get("modularity_tad", {})

    max_dist_mb  = float(algo_cfg.get("max_dist_mb",  5.0))
    penalty_cfg  = algo_cfg.get("penalty", None)
    min_tad_kb   = float(algo_cfg.get("min_tad_kb",  100.0))
    max_tad_kb   = float(algo_cfg.get("max_tad_kb", 3000.0))

    min_tad_bins  = max(2, int(min_tad_kb  * 1_000 / resolution))
    max_tad_bins  = max(min_tad_bins + 1,
                        int(max_tad_kb * 1_000 / resolution))
    max_dist_bins = int(max_dist_mb * 1_000_000 / resolution)

    # ── Загрузка матрицы ────────────────────────────────────────────────────
    try:
        from src.data_prep import get_matrix  # type: ignore
        matrix = get_matrix(cfg, chrom, resolution)
    except Exception as exc:
        logger.error("[modtad] get_matrix failed: %s", exc, exc_info=True)
        return _empty

    if matrix is None or matrix.size == 0:
        logger.error("[modtad] Empty matrix for %s @ %d", chrom, resolution)
        return _empty

    A = np.asarray(matrix, dtype=np.float64)
    n = A.shape[0]

    # Маскировать дальние контакты
    for d in range(max_dist_bins, n):
        idx = np.arange(n - d)
        A[idx, idx + d] = 0.0
        A[idx + d, idx] = 0.0

    # Симметризация
    A = (A + A.T) / 2.0

    chrom_mb = n * resolution / 1_000_000
    logger.info("[modtad] %s @ %d: %dx%d (%.1f Mb), min=%d max=%d bins",
                chrom, resolution, n, n, chrom_mb, min_tad_bins, max_tad_bins)

    # ── O/E нормализация ────────────────────────────────────────────────────
    # Убирает distance-decay без bias к размеру TAD.
    # После нормировки среднее на каждой диагонали = 1.
    A_oe = _oe_normalize(A)

    # ── Prefix-суммы и score-матрица ────────────────────────────────────────
    P = _build_prefix2d(A_oe)

    score_matrix = np.full((n, n), -np.inf, dtype=np.float32)
    for i in range(n):
        j_min = i + min_tad_bins
        j_max = min(n, i + max_tad_bins)
        for j in range(j_min, j_max + 1):
            score_matrix[i, j - 1] = _score_intra_minus_flanks(P, i, j, n)

    # ── DP + penalty ────────────────────────────────────────────────────────
    if penalty_cfg is not None:
        pen  = float(penalty_cfg)
        segs = _dp_segment(score_matrix, pen, min_tad_bins, max_tad_bins, n)
        logger.info("[modtad] fixed penalty=%.4f → %d TADs", pen, len(segs))
    else:
        pen, segs = _auto_penalty_sweep(
            score_matrix, min_tad_bins, max_tad_bins, n, chrom_mb
        )

    if not segs:
        logger.warning("[modtad] DP returned 0 segments for %s", chrom)
        return _empty

    # ── Конвертация в bp ────────────────────────────────────────────────────
    tads = []
    for s_bin, e_bin in segs:
        start_bp = s_bin * resolution
        end_bp   = e_bin * resolution
        if (end_bp - start_bp) / 1_000 >= min_tad_kb:
            tads.append({"chrom": chrom, "start": start_bp, "end": end_bp})

    if not tads:
        logger.warning("[modtad] 0 TADs after bp-filter on %s", chrom)
        return _empty

    df_out = pd.DataFrame(tads, columns=["chrom", "start", "end"])
    logger.info("[modtad] FINAL: %d TADs on %s @ %d, median=%.0fkb",
                len(df_out), chrom, resolution,
                (df_out.end - df_out.start).median() / 1000)
    return df_out
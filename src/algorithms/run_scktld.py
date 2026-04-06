"""
run_scktld.py
=============
Адаптация scKTLD для bulk Hi-C данных.

scKTLD (Single-Cell Kernel TAD Learning with Domains) адаптирован
для работы с dense numpy-матрицами одной хромосомы (RAW, без нормализации).

Ограничения памяти (из конфига):
  10 kb  → chr21, chr22
  25 kb  → chr17–chr22
  50 kb  → chr1–chr22
  100 kb → chr1–chr22

Алгоритм:
  1. KNN-граф на основе контактной матрицы (ядро RBF)
  2. Спектральное вложение (dimension=128)
  3. Сегментация через динамическое программирование с штрафом
  4. penalty: автоподбор через grid search по внутренней оценке
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigsh
from sklearn.preprocessing import normalize

from src.data_prep import get_matrix

logger = logging.getLogger(__name__)

RNG = np.random.default_rng(42)


# ──────────────────────────────────────────────────────────────────────────────
# Ядерные функции
# ──────────────────────────────────────────────────────────────────────────────

def _rbf_kernel_matrix(matrix: np.ndarray, sigma: Optional[float] = None) -> np.ndarray:
    """
    Построить RBF-матрицу сходства из контактной матрицы.
    K[i,j] = exp(- ||row_i - row_j||^2 / (2*sigma^2))
    """
    # log1p для стабилизации
    M = np.log1p(matrix.astype(np.float64))
    n = M.shape[0]

    # Вычисляем расстояния через squared norms
    sq_norms = np.sum(M ** 2, axis=1, keepdims=True)
    dists_sq = sq_norms + sq_norms.T - 2.0 * (M @ M.T)
    np.clip(dists_sq, 0, None, out=dists_sq)

    if sigma is None:
        sigma = np.sqrt(np.median(dists_sq[dists_sq > 0])) if np.any(dists_sq > 0) else 1.0

    K = np.exp(-dists_sq / (2.0 * sigma ** 2))
    np.fill_diagonal(K, 0.0)
    return K


def _spectral_embedding(K: np.ndarray, dimension: int) -> np.ndarray:
    """
    Спектральное вложение: нормализованный Лапласиан + eigsh.
    Возвращает матрицу (n, dimension).
    """
    n = K.shape[0]
    dimension = min(dimension, n - 2)

    d = K.sum(axis=1)
    d_inv_sqrt = np.where(d > 0, 1.0 / np.sqrt(d), 0.0)
    # Нормализованный Лапласиан: L_sym = D^{-1/2} K D^{-1/2}
    L_sym = d_inv_sqrt[:, None] * K * d_inv_sqrt[None, :]

    # Разреженный формат для эффективности
    L_sparse = csr_matrix(L_sym)

    try:
        eigenvalues, eigenvectors = eigsh(L_sparse, k=dimension + 1, which="LM")
    except Exception as exc:
        logger.warning("eigsh failed: %s, falling back to eigh", exc)
        eigenvalues, eigenvectors = np.linalg.eigh(L_sym)
        eigenvalues = eigenvalues[-dimension-1:]
        eigenvectors = eigenvectors[:, -dimension-1:]

    # Сортируем по убыванию
    idx = np.argsort(eigenvalues)[::-1]
    eigenvectors = eigenvectors[:, idx[1:dimension+1]]  # пропускаем тривиальный

    embedding = normalize(eigenvectors, norm="l2")
    return embedding.astype(np.float32)


# ──────────────────────────────────────────────────────────────────────────────
# Сегментация через DP
# ──────────────────────────────────────────────────────────────────────────────

def _segment_cost(embedding: np.ndarray, i: int, j: int) -> float:
    """
    Стоимость сегмента [i, j] — среднее квадратичное отклонение внутри.
    """
    if j <= i:
        return 0.0
    seg = embedding[i:j+1]
    centroid = seg.mean(axis=0)
    return float(np.sum((seg - centroid) ** 2))


def _dynamic_programming_segmentation(
    embedding: np.ndarray,
    penalty: float,
    min_size: int = 3,
) -> list[int]:
    """
    Оптимальная сегментация через DP (PELT-like, упрощённая).
    Возвращает список индексов границ (включая 0 и n).
    """
    n = embedding.shape[0]
    # dp[i] = минимальная стоимость сегментации [0, i)
    dp   = np.full(n + 1, np.inf)
    prev = np.full(n + 1, -1, dtype=int)
    dp[0] = 0.0

    # Кэш стоимостей для ускорения
    cost_cache: dict[tuple[int, int], float] = {}

    def cost(a: int, b: int) -> float:
        if (a, b) not in cost_cache:
            cost_cache[(a, b)] = _segment_cost(embedding, a, b - 1)
        return cost_cache[(a, b)]

    for j in range(1, n + 1):
        for i in range(max(0, j - n), j - min_size + 1):
            if dp[i] + cost(i, j) + penalty < dp[j]:
                dp[j] = dp[i] + cost(i, j) + penalty
                prev[j] = i

    # Трассировка
    boundaries = []
    pos = n
    while pos > 0:
        boundaries.append(pos)
        pos = prev[pos]
    boundaries.reverse()
    return boundaries


def _auto_penalty(
    embedding: np.ndarray,
    penalty_grid: Optional[list[float]] = None,
    min_size: int = 3,
) -> float:
    """
    Автоподбор penalty: выбирается значение, при котором число TAD
    стабилизируется (минимальное изменение при увеличении penalty).
    """
    if penalty_grid is None:
        n = embedding.shape[0]
        # Логарифмическая сетка
        penalty_grid = list(np.logspace(
            np.log10(max(0.01, n * 0.001)),
            np.log10(n * 2.0),
            num=15,
        ))

    counts: list[tuple[float, int]] = []
    for p in penalty_grid:
        bnd = _dynamic_programming_segmentation(embedding, p, min_size)
        n_domains = len(bnd) - 1
        counts.append((p, n_domains))
        logger.debug("  penalty=%.3f → %d TADs", p, n_domains)

    # Elbow: наибольшее изменение наклона
    n_counts = np.array([c[1] for c in counts], dtype=float)
    diffs = np.abs(np.diff(n_counts))
    # Выбираем точку после крутого падения
    elbow_idx = int(np.argmax(diffs)) + 1
    best_p = counts[elbow_idx][0]
    logger.info("auto_penalty: выбрано %.4f (elbow @ idx=%d, TADs=%d)",
                best_p, elbow_idx, counts[elbow_idx][1])
    return best_p


# ──────────────────────────────────────────────────────────────────────────────
# Публичный интерфейс
# ──────────────────────────────────────────────────────────────────────────────

def run_scktld(
    chrom: str,
    resolution: int,
    data_path: str,
    cfg: Optional[dict] = None,
    dimension: int = 128,
    penalty: Optional[float] = None,
    matrix: Optional[np.ndarray] = None,
) -> pd.DataFrame:
    """
    Запустить scKTLD (адаптация для bulk Hi-C).

    Parameters
    ----------
    chrom      : хромосома
    resolution : разрешение в bp
    data_path  : путь к директории processed
    cfg        : конфиг-словарь
    dimension  : размерность спектрального вложения
    penalty    : штраф DP (None → автоподбор)
    matrix     : можно передать напрямую (RAW, без нормализации)

    Returns
    -------
    pd.DataFrame с колонками chrom, start, end
    """
    # Проверка ограничений памяти
    if cfg is not None:
        limits = cfg["chromosomes"]["scktld_limits"].get(resolution, [])
        if limits and chrom not in limits:
            logger.warning(
                "[scKTLD] %s @ %d bp пропущена (ограничение памяти, разрешены: %s)",
                chrom, resolution, limits,
            )
            return pd.DataFrame(columns=["chrom", "start", "end"])
        dimension = cfg["algorithms"]["scktld"]["dimension"]

    # Загрузить матрицу (RAW, balance=False)
    if matrix is None:
        if cfg is not None:
            matrix = get_matrix(cfg, chrom, resolution)
        else:
            import os
            npy = os.path.join(data_path, f"{chrom}_{resolution}bp.npy")
            if not os.path.exists(npy):
                raise FileNotFoundError(f"Матрица не найдена: {npy}")
            matrix = np.load(npy)

    matrix = matrix.astype(np.float64)
    n = matrix.shape[0]
    logger.info("[scKTLD] %s @ %d bp  | n_bins=%d dim=%d", chrom, resolution, n, dimension)

    if n < 10:
        logger.warning("[scKTLD] Слишком мало бинов (%d) для %s", n, chrom)
        return pd.DataFrame(columns=["chrom", "start", "end"])

    # 1. RBF-ядро
    K = _rbf_kernel_matrix(matrix)

    # 2. Спектральное вложение
    embedding = _spectral_embedding(K, dimension)

    # 3. Penalty
    if penalty is None:
        penalty = _auto_penalty(embedding)

    # 4. Сегментация
    boundaries = _dynamic_programming_segmentation(
        embedding, penalty, min_size=max(3, resolution // 50000 + 1)
    )

    # 5. Формирование доменов
    records = []
    for k in range(len(boundaries) - 1):
        start = boundaries[k - 1] * resolution if k > 0 else 0
        # Исправленная трассировка: используем сами значения из boundaries
        # как правые края сегментов
        pass

    # Корректная трассировка из DP
    records = []
    prev_end = 0
    for bnd in boundaries[1:]:  # первый элемент = n (полная длина)
        end = bnd * resolution
        if end > prev_end:
            records.append((chrom, prev_end, end))
        prev_end = end

    if not records:
        return pd.DataFrame(columns=["chrom", "start", "end"])

    df = pd.DataFrame(records, columns=["chrom", "start", "end"])
    df = df[df["end"] > df["start"]].reset_index(drop=True)
    logger.info("[scKTLD] Итого %d TAD", len(df))
    return df
"""
run_scktld.py
=============
Адаптация scKTLD для bulk Hi-C данных.

Алгоритм:
  1. KNN-sparse RBF-граф на контактной матрице (без нормализации, RAW)
  2. Нормализованный Лапласиан → eigsh → спектральное вложение
  3. DP-сегментация с penalty (cumsum-ускорение)
  4. _auto_penalty: logspace-сетка + фильтр min_tads/max_tads + elbow

Ограничения памяти (из конфига):
  10 kb  → chr21, chr22
  25 kb  → chr17–chr22
  50 kb  → chr1–chr22
  100 kb → chr1–chr22
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, diags
from scipy.sparse.linalg import eigsh
from sklearn.preprocessing import normalize

from src.data_prep import get_matrix

# ──────────────────────────────────────────────────────────────────────────────
# GPU-детектирование (P2: GPU-адаптация для сервера)
# ──────────────────────────────────────────────────────────────────────────────

def _get_device() -> str:
    """Определить устройство: 'cuda' если доступно, иначе 'cpu'."""
    try:
        import torch
        if torch.cuda.is_available():
            # print т.к. вызывается до настройки logging (module-level)
            print(f"[scKTLD] GPU: {torch.cuda.get_device_name(0)}")
            return "cuda"
    except Exception:
        pass
    print("[scKTLD] CPU режим")
    return "cpu"

DEVICE = _get_device()

# Порог: n < 5000 бинов → плотная матрица ~100 MB, безопасно для GPU
# chr17-chr22 @ 25kb: ~1900–3250 бинов → всегда GPU-путь
_GPU_DENSE_THRESHOLD = 5000

logger = logging.getLogger(__name__)

# ──────────────────────────────────────────────────────────────────────────────
# Константы (переопределяются через cfg)
# ──────────────────────────────────────────────────────────────────────────────
_DEFAULT_DIMENSION = 32
_DEFAULT_KNN_K     = 20
_AUTO_MIN_TADS     = 20
_AUTO_MAX_TADS     = 200


# ──────────────────────────────────────────────────────────────────────────────
# 1. KNN-sparse RBF-ядро
# ──────────────────────────────────────────────────────────────────────────────

def _rbf_kernel_matrix(
    matrix: np.ndarray,
    sigma: Optional[float] = None,
    knn_k: int = _DEFAULT_KNN_K,
) -> csr_matrix:
    """
    KNN-sparse RBF: для каждого бина оставляем knn_k ближайших соседей
    внутри скользящего окна (Hi-C локален — дальние контакты шумные).

    Возвращает симметричную csr_matrix (не dense).
    """
    M = np.log1p(matrix.astype(np.float32))
    n = M.shape[0]
    knn_k  = min(knn_k, n - 1)
    WINDOW = min(500, n - 1)   # ≤500 бинов = ≤12.5 Mb @ 25 kb

    # Оцениваем sigma по выборке из 200 бинов
    sample_idx = np.linspace(0, n - 1, min(200, n), dtype=int)
    sample = M[sample_idx]
    sq  = np.sum(sample ** 2, axis=1, keepdims=True)
    d_s = sq + sq.T - 2.0 * (sample @ sample.T)
    np.clip(d_s, 0, None, out=d_s)
    if sigma is None:
        nz = d_s[d_s > 0]
        sigma = float(np.sqrt(np.median(nz))) if len(nz) > 0 else 1.0
    inv_2s2 = 1.0 / (2.0 * sigma ** 2 + 1e-12)

    rows, cols, vals = [], [], []
    for i in range(n):
        lo = max(0, i - WINDOW)
        hi = min(n, i + WINDOW + 1)
        chunk = M[lo:hi]
        d2 = np.sum((chunk - M[i]) ** 2, axis=1)
        k_local = min(knn_k, len(d2) - 1)
        top_k   = np.argpartition(d2, k_local)[:k_local + 1]
        for t in top_k:
            j = lo + t
            if j == i:
                continue
            w = float(np.exp(-d2[t] * inv_2s2))
            rows.append(i);  cols.append(j);  vals.append(w)
            rows.append(j);  cols.append(i);  vals.append(w)  # симметрия

    return csr_matrix((vals, (rows, cols)), shape=(n, n))


# ──────────────────────────────────────────────────────────────────────────────
# 2. Спектральное вложение
# ──────────────────────────────────────────────────────────────────────────────

def _spectral_embedding(K: csr_matrix, dimension: int) -> np.ndarray:
    """
    Нормализованный Лапласиан → eigendecomposition → вложение.

    Ветки:
      GPU (n < 5000): torch.linalg.eigh на CUDA (~100 MB плотная матрица)
      GPU (n >= 5000): fallback на CPU eigsh (слишком много памяти)
      CPU: scipy eigsh (текущая реализация)
    """
    if not hasattr(K, "toarray"):
        K = csr_matrix(K)

    n = K.shape[0]
    dimension = min(dimension, n - 2)

    # Нормализованный Лапласиан (общий для всех веток)
    d = np.asarray(K.sum(axis=1)).ravel()
    d_inv_sqrt = np.where(d > 0, 1.0 / np.sqrt(d), 0.0)
    D_inv_sqrt = diags(d_inv_sqrt)
    L_sym = D_inv_sqrt @ K @ D_inv_sqrt

    # ── GPU-ветка ─────────────────────────────────────────────────────────────
    if DEVICE == "cuda" and n < _GPU_DENSE_THRESHOLD:
        try:
            import torch
            # Плотная матрица: n=3248 → ~42 MB float32 на GPU
            L_dense = torch.tensor(
                L_sym.toarray(), dtype=torch.float32, device="cuda"
            )
            # eigh возвращает eigenvalues в порядке возрастания
            eigenvalues, eigenvectors = torch.linalg.eigh(L_dense)

            # Берём dimension наибольших (последние в возрастающем порядке)
            # Пропускаем последний (наибольший = 1 для нормализованного L)
            total = eigenvectors.shape[1]
            idx_start = max(0, total - dimension - 1)
            idx_end   = total - 1                        # исключаем наибольший
            vecs = eigenvectors[:, idx_start:idx_end]    # (n, dimension)

            # Перевернуть: хотим убывающий порядок eigenvalues
            vecs = torch.flip(vecs, dims=[1])

            result = vecs.cpu().numpy().astype(np.float32)
            del L_dense, eigenvalues, eigenvectors, vecs
            torch.cuda.empty_cache()

            logger.info(
                "[scKTLD] GPU eigh: n=%d, dim=%d", n, dimension
            )
            return normalize(result, norm="l2")

        except Exception as exc:
            logger.warning(
                "[scKTLD] GPU eigh failed (%s), fallback → CPU eigsh", exc
            )
            # fallthrough к CPU-ветке

    # ── CPU-ветка (sparse eigsh) ──────────────────────────────────────────────
    try:
        eigenvalues, eigenvectors = eigsh(
            L_sym, k=dimension + 1, which="LM",
            maxiter=1000, tol=1e-4,
        )
    except Exception as exc:
        k_fb = min(16, n - 2)
        logger.warning("[scKTLD] eigsh failed (%s), fallback k=%d", exc, k_fb)
        eigenvalues, eigenvectors = eigsh(
            L_sym, k=k_fb + 1, which="LM", tol=1e-3,
        )

    idx          = np.argsort(eigenvalues)[::-1]
    eigenvectors = eigenvectors[:, idx[1: dimension + 1]]
    logger.info("[scKTLD] CPU eigsh: n=%d, dim=%d", n, dimension)
    return normalize(eigenvectors, norm="l2").astype(np.float32)

# ──────────────────────────────────────────────────────────────────────────────
# 3. DP-сегментация (cumsum, O(n²))
# ──────────────────────────────────────────────────────────────────────────────

def _dp_segmentation_fast(
    embedding: np.ndarray,
    penalty: float,
    min_size: int = 3,
) -> list[int]:
    """
    Оптимальная сегментация через DP с cumsum-ускорением.

    Возвращает список правых границ сегментов (в бинах, не включая 0),
    последний элемент всегда == n:
        [b1, b2, ..., n]
    Сегменты: [0, b1), [b1, b2), ..., [b_{k-1}, n)
    """
    n = embedding.shape[0]
    E   = embedding.astype(np.float64)
    # Cumsum для cost(a,b) без явного Python-цикла внутри cost()
    cs  = np.zeros((n + 1, E.shape[1]))
    cs2 = np.zeros(n + 1)
    cs[1:]  = np.cumsum(E, axis=0)
    cs2[1:] = np.cumsum(np.sum(E ** 2, axis=1))

    def cost(a: int, b: int) -> float:
        """Внутрисегментная дисперсия × length для [a, b)."""
        length = b - a
        if length <= 0:
            return 0.0
        seg_sum = cs[b] - cs[a]
        seg_sq  = cs2[b] - cs2[a]
        return float(seg_sq - np.dot(seg_sum, seg_sum) / length)

    dp   = np.full(n + 1, np.inf)
    prev = np.full(n + 1, -1, dtype=int)
    dp[0] = 0.0

    for j in range(min_size, n + 1):
        i_max = j - min_size + 1
        for i in range(0, i_max):
            if dp[i] == np.inf:
                continue
            c = dp[i] + cost(i, j) + penalty
            if c < dp[j]:
                dp[j] = c
                prev[j] = i

    # Трассировка
    boundaries: list[int] = []
    pos = n
    while pos > 0:
        boundaries.append(pos)
        pos = prev[pos]
    boundaries.reverse()
    return boundaries  # [b1, b2, ..., n]


# ──────────────────────────────────────────────────────────────────────────────
# 4. Автоподбор penalty
# ──────────────────────────────────────────────────────────────────────────────

def _auto_penalty(
    embedding: np.ndarray,
    penalty_grid: Optional[list[float]] = None,
    min_size: int = 3,
    min_tads: int = _AUTO_MIN_TADS,
    max_tads: int = _AUTO_MAX_TADS,
) -> float:
    """
    Автоподбор penalty по логарифмической сетке.

    Стратегия:
    1. Просчитать число TAD для каждого penalty в сетке.
    2. Оставить только penalty с min_tads ≤ count ≤ max_tads (valid zone).
    3. Среди валидных найти elbow (максимальное |Δcount|).
    4. Если valid zone пуста — выбрать penalty, ближайшую к середине
       диапазона [min_tads, max_tads].

    Parameters
    ----------
    embedding    : спектральное вложение (n, d)
    penalty_grid : список значений penalty (None → авто)
    min_size     : минимальный размер TAD в бинах
    min_tads     : нижняя граница допустимого числа TAD
    max_tads     : верхняя граница допустимого числа TAD

    Returns
    -------
    float : выбранный penalty
    """
    n = embedding.shape[0]

    if penalty_grid is None:
        # Сетка: от "много TAD" до "мало TAD"
        # Нижняя граница — penalty даёт ~max_tads TAD
        # Верхняя граница — penalty даёт ~min_tads TAD
        # Эмпирически: cost одного бина ≈ trace(cov(embedding))
        typical_cost = float(np.trace(np.cov(embedding.T))) if n > 2 else 1.0
        p_lo = max(typical_cost * 0.05, n * 0.001)
        p_hi = typical_cost * 20.0
        penalty_grid = list(np.logspace(
            np.log10(p_lo),
            np.log10(max(p_hi, p_lo * 10)),
            num=20,
        ))

    # ── Шаг 1: просчёт ────────────────────────────────────────────────────
    counts: list[tuple[float, int]] = []
    for p in penalty_grid:
        bnd     = _dp_segmentation_fast(embedding, p, min_size)
        n_tads  = max(0, len(bnd) - 1)
        counts.append((p, n_tads))
        logger.debug("[scKTLD] auto_penalty  p=%.4f → %d TADs", p, n_tads)

    penalties  = np.array([c[0] for c in counts])
    n_tads_arr = np.array([c[1] for c in counts], dtype=float)

    # ── Шаг 2: valid zone ─────────────────────────────────────────────────
    valid_mask = (n_tads_arr >= min_tads) & (n_tads_arr <= max_tads)
    valid_idx  = np.where(valid_mask)[0]

    if len(valid_idx) == 0:
        # Нет ни одного penalty в допустимом диапазоне →
        # выбираем ближайший к target_tads (середина диапазона)
        target = (min_tads + max_tads) / 2.0
        closest = int(np.argmin(np.abs(n_tads_arr - target)))
        best_p  = float(penalties[closest])
        logger.warning(
            "[scKTLD] auto_penalty: valid zone пуста (min=%d max=%d), "
            "выбран ближайший к target=%.0f: p=%.4f → %d TADs",
            min_tads, max_tads, target, best_p, int(n_tads_arr[closest]),
        )
        return best_p

    # ── Шаг 3: elbow среди валидных ───────────────────────────────────────
    valid_penalties = penalties[valid_idx]
    valid_counts    = n_tads_arr[valid_idx]

    if len(valid_idx) == 1:
        best_p = float(valid_penalties[0])
        logger.info(
            "[scKTLD] auto_penalty: единственный валидный p=%.4f → %d TADs",
            best_p, int(valid_counts[0]),
        )
        return best_p

    # Elbow = индекс максимального абсолютного перепада
    diffs     = np.abs(np.diff(valid_counts))
    elbow_pos = int(np.argmax(diffs)) + 1   # +1: elbow после перепада

    # Берём точку ПОСЛЕ elbow (где count стабилизировался)
    elbow_pos = min(elbow_pos, len(valid_idx) - 1)
    best_p    = float(valid_penalties[elbow_pos])

    logger.info(
        "[scKTLD] auto_penalty: выбрано p=%.4f (elbow_pos=%d, TADs=%d) "
        "[valid zone: %d–%d TADs, %d точек]",
        best_p, elbow_pos, int(valid_counts[elbow_pos]),
        int(valid_counts.min()), int(valid_counts.max()), len(valid_idx),
    )
    return best_p


# ──────────────────────────────────────────────────────────────────────────────
# Публичный интерфейс
# ──────────────────────────────────────────────────────────────────────────────

def run_scktld(
    chrom: str,
    resolution: int,
    data_path: str,
    cfg: Optional[dict] = None,
    dimension: int = _DEFAULT_DIMENSION,
    penalty: Optional[float] = None,
    matrix: Optional[np.ndarray] = None,
) -> pd.DataFrame:
    """
    Запустить scKTLD (адаптация для bulk Hi-C).

    Parameters
    ----------
    chrom      : хромосома ('chr17', ...)
    resolution : разрешение в bp (25000, 50000, 100000)
    data_path  : путь к директории data/processed/
    cfg        : конфиг-словарь (None → использовать значения по умолчанию)
    dimension  : размерность спектрального вложения (из cfg → 32)
    penalty    : штраф DP (None → _auto_penalty)
    matrix     : RAW-матрица (None → загружается через get_matrix)

    Returns
    -------
    pd.DataFrame с колонками ['chrom', 'start', 'end']
    При ошибке / пропуске — пустой DataFrame.
    """
    _empty = pd.DataFrame(columns=["chrom", "start", "end"])

    # ── Параметры из конфига ───────────────────────────────────────────────
    knn_k     = _DEFAULT_KNN_K
    min_tads  = _AUTO_MIN_TADS
    max_tads  = _AUTO_MAX_TADS

    if cfg is not None:
        # Проверка ограничений памяти
        limits = cfg["chromosomes"]["scktld_limits"].get(resolution, [])
        if limits and chrom not in limits:
            logger.warning(
                "[scKTLD] %s @ %d bp пропущена (ограничение памяти). "
                "Разрешены: %s",
                chrom, resolution, limits,
            )
            return _empty

        algo_cfg  = cfg["algorithms"]["scktld"]
        dimension = algo_cfg.get("dimension", dimension)
        knn_k     = algo_cfg.get("knn_k",     knn_k)

    try:
        # ── Загрузка матрицы (RAW, balance=False) ─────────────────────────
        if matrix is None:
            if cfg is not None:
                matrix = get_matrix(cfg, chrom, resolution)
            else:
                import os
                npy = os.path.join(data_path, f"{chrom}_{resolution}bp.npy")
                if not os.path.exists(npy):
                    raise FileNotFoundError(f"Матрица не найдена: {npy}")
                matrix = np.load(npy)

        matrix = np.nan_to_num(matrix.astype(np.float64), nan=0.0, posinf=0.0)
        n      = matrix.shape[0]

        logger.info(
            "[scKTLD] %s @ %d bp | n_bins=%d dim=%d knn_k=%d",
            chrom, resolution, n, dimension, knn_k,
        )

        if n < 10:
            logger.warning("[scKTLD] Слишком мало бинов (%d) для %s", n, chrom)
            return _empty

        # ── 1. RBF-ядро ───────────────────────────────────────────────────
        K = _rbf_kernel_matrix(matrix, knn_k=knn_k)

        # ── 2. Спектральное вложение ───────────────────────────────────────
        embedding = _spectral_embedding(K, dimension)

        # ── 3. Penalty ────────────────────────────────────────────────────
        min_size = max(3, 100_000 // resolution)   # минимум ~100 kb в бинах
        if penalty is None:
            penalty = _auto_penalty(
                embedding,
                min_size=min_size,
                min_tads=min_tads,
                max_tads=max_tads,
            )

        # ── 4. Сегментация ────────────────────────────────────────────────
        boundaries = _dp_segmentation_fast(embedding, penalty, min_size=min_size)

        # ── 5. Формирование доменов ───────────────────────────────────────
        # boundaries = [b1, b2, ..., n]  — правые края сегментов (в бинах)
        # Сегмент k: [prev_bin, bnd_bin) в бинах → [prev_bin*res, bnd_bin*res) в bp
        records: list[tuple[str, int, int]] = []
        prev_bin = 0
        for bnd_bin in boundaries:
            start_bp = prev_bin  * resolution
            end_bp   = bnd_bin   * resolution
            if end_bp > start_bp:
                records.append((chrom, start_bp, end_bp))
            prev_bin = bnd_bin

        if not records:
            logger.warning("[scKTLD] Нет доменов после сегментации (%s)", chrom)
            return _empty

        df = (
            pd.DataFrame(records, columns=["chrom", "start", "end"])
            .query("end > start")
            .reset_index(drop=True)
        )
        logger.info("[scKTLD] Итого %d TAD  (%s @ %d bp)", len(df), chrom, resolution)
        return df

    except Exception as exc:
        logger.error(
            "[scKTLD] Ошибка %s @ %d: %s", chrom, resolution, exc, exc_info=True,
        )
        return _empty
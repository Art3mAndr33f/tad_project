"""
src/algorithms/run_dihmm.py

TAD detection via Directionality Index + Hidden Markov Model.
Reference: Dixon et al. Nature 2012.

Algorithm:
  1. Compute DI(i) = sign(D-U) * (chi2_D + chi2_U)
     where D = downstream contacts sum, U = upstream contacts sum,
     E = (D+U)/2, chi2_X = (X-E)^2 / E
  2. Fit GaussianHMM(n_components=3) on DI series
  3. Assign states: downstream-biased (+DI), unbiased (~0), upstream-biased (-DI)
  4. TAD = region from state transition "upstream-biased → downstream-biased"
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d

logger = logging.getLogger(__name__)

# ── DI window (в бинах) ──────────────────────────────────────────────────────
_DEFAULT_WINDOW_BINS = 10
_HMM_N_ITER = 100
_HMM_N_STATES = 3
_MIN_TAD_BINS = 3   # минимальный TAD = 3 бина (не фильтруем по размеру в bp здесь,
                    # run_dihmm применяет min_tad_kb из конфига ниже)
_SEED = 42


# ═══════════════════════════════════════════════════════════════════════════════
# Публичный интерфейс
# ═══════════════════════════════════════════════════════════════════════════════

def run_dihmm(
    chrom: str,
    resolution: int,
    data_path: str,
    cfg: Optional[dict],
    **kwargs,
) -> pd.DataFrame:
    """
    Детекция TAD через Directionality Index + GaussianHMM.

    Parameters
    ----------
    chrom      : 'chr1', 'chrX', …
    resolution : в bp (25000, 50000, 100000)
    data_path  : путь к data/processed/
    cfg        : конфиг-словарь (может быть None)

    Returns
    -------
    pd.DataFrame(columns=['chrom','start','end'])
    """
    empty = pd.DataFrame(columns=["chrom", "start", "end"])

    # ── 1. Загрузить матрицу ─────────────────────────────────────────────────
    try:
        matrix = _load_matrix(data_path, chrom, resolution)
    except FileNotFoundError as exc:
        logger.error("[dihmm] Матрица не найдена %s @ %d: %s", chrom, resolution, exc)
        return empty

    n = matrix.shape[0]
    logger.info("[dihmm] %s @ %d — матрица %d×%d загружена", chrom, resolution, n, n)

    if n < 2 * _DEFAULT_WINDOW_BINS + 1:
        logger.warning("[dihmm] %s @ %d — слишком мало бинов (%d), пропуск", chrom, resolution, n)
        return empty

    # ── 2. Параметры из конфига ──────────────────────────────────────────────
    algo_cfg = (cfg or {}).get("algorithms", {}).get("dihmm", {})
    window_bins: int  = int(algo_cfg.get("window_bins", _DEFAULT_WINDOW_BINS))
    n_states: int     = int(algo_cfg.get("n_states", _HMM_N_STATES))
    n_iter: int       = int(algo_cfg.get("n_iter", _HMM_N_ITER))
    min_tad_kb: float = float(algo_cfg.get("min_tad_kb", 100.0))
    seed: int         = int((cfg or {}).get("project", {}).get("seed", _SEED))

    min_tad_bins = max(_MIN_TAD_BINS, int(min_tad_kb * 1_000 / resolution))

    # ── 3. Directionality Index ──────────────────────────────────────────────
    di = _compute_di(matrix, window_bins=window_bins)
    logger.debug("[dihmm] DI вычислен, диапазон [%.3f, %.3f]", di.min(), di.max())

    # ── 4. HMM ──────────────────────────────────────────────────────────────
    try:
        states = _fit_hmm(di, n_states=n_states, n_iter=n_iter, seed=seed)
    except Exception as exc:
        logger.error("[dihmm] HMM упал %s @ %d: %s", chrom, resolution, exc, exc_info=True)
        return empty

    # ── 5. Извлечь TAD из состояний ─────────────────────────────────────────
    boundaries = _states_to_tads(states, di, n_bins=n, min_tad_bins=min_tad_bins)

    # ── 6. Конвертировать бины → координаты bp ──────────────────────────────
    records = []
    for (b_start, b_end) in boundaries:
        start_bp = int(b_start) * resolution
        end_bp   = int(b_end)   * resolution
        records.append({"chrom": chrom, "start": start_bp, "end": end_bp})

    df = pd.DataFrame(records, columns=["chrom", "start", "end"])
    logger.info(
        "[dihmm] %s @ %d — найдено %d TAD (window=%d bins, min_tad=%.0f kb)",
        chrom, resolution, len(df), window_bins, min_tad_kb,
    )
    return df


# ═══════════════════════════════════════════════════════════════════════════════
# Внутренние функции
# ═══════════════════════════════════════════════════════════════════════════════

def _load_matrix(data_path: str, chrom: str, resolution: int) -> np.ndarray:
    """Загрузить кэшированную .npy матрицу."""
    npy_path = Path(data_path) / f"{chrom}_{resolution}bp.npy"
    if not npy_path.exists():
        raise FileNotFoundError(npy_path)
    matrix = np.load(str(npy_path)).astype(np.float64)
    matrix = np.nan_to_num(matrix, nan=0.0, posinf=0.0, neginf=0.0)
    # Симметризовать (на случай если матрица верхнетреугольная)
    matrix = np.maximum(matrix, matrix.T)
    return matrix


def _compute_di(matrix: np.ndarray, window_bins: int = 10) -> np.ndarray:
    """
    Directionality Index (Dixon et al. 2012).

    DI(i) = sign(D-U) * [ (D-E)²/E + (U-E)²/E ]
    где E = (D+U)/2,
        D = сумма контактов i→[i+1, i+window] (downstream),
        U = сумма контактов i→[i-window, i-1] (upstream).

    Граничные бины (первые/последние window_bins) заполняются 0.
    """
    n = matrix.shape[0]
    di = np.zeros(n, dtype=np.float64)

    for i in range(window_bins, n - window_bins):
        # Upstream: [i-window, i-1]
        U = np.sum(matrix[i, i - window_bins: i])
        # Downstream: [i+1, i+window]
        D = np.sum(matrix[i, i + 1: i + window_bins + 1])

        total = U + D
        if total == 0.0:
            di[i] = 0.0
            continue

        E = total / 2.0
        # chi² компоненты
        chi2 = ((D - E) ** 2 / E + (U - E) ** 2 / E)
        sign = 1.0 if D > U else (-1.0 if U > D else 0.0)
        di[i] = sign * chi2

    # Лёгкое сглаживание для стабильности HMM
    di = gaussian_filter1d(di, sigma=1.0)
    return di


def _fit_hmm(
    di: np.ndarray,
    n_states: int = 3,
    n_iter: int = 100,
    seed: int = 42,
) -> np.ndarray:
    """
    Обучить GaussianHMM на DI-серии и вернуть последовательность состояний.

    Состояния после ремаппинга:
      0 → upstream-biased   (DI < 0, начало TAD)
      1 → unbiased          (DI ≈ 0)
      2 → downstream-biased (DI > 0, конец TAD)
    """
    from hmmlearn.hmm import GaussianHMM  # type: ignore

    obs = di.reshape(-1, 1)

    # Инициализация средних через квантили DI
    q_low  = np.quantile(di, 0.15)
    q_mid  = np.quantile(di, 0.50)
    q_high = np.quantile(di, 0.85)

    model = GaussianHMM(
        n_components=n_states,
        covariance_type="diag",
        n_iter=n_iter,
        random_state=seed,
        init_params="stc",   # самоинициализировать startprob, transmat, covars
        params="stmc",
    )
    # Ручная инициализация средних для воспроизводимости
    model.means_init = np.array([[q_low], [q_mid], [q_high]])

    model.fit(obs)
    raw_states = model.predict(obs)

    # ── Ремаппинг: упорядочить состояния по среднему DI ──────────────────────
    means = model.means_.flatten()
    order = np.argsort(means)          # [idx_min, idx_mid, idx_max]
    remap = {old: new for new, old in enumerate(order)}
    # После ремаппинга: 0=upstream-biased, 1=unbiased, 2=downstream-biased
    states = np.vectorize(remap.__getitem__)(raw_states)

    logger.debug(
        "[dihmm] HMM средние до ремаппинга: %s → после: upstream<mid<downstream",
        np.sort(means).round(3),
    )
    return states


def _states_to_tads(
    states: np.ndarray,
    di: np.ndarray,
    n_bins: int,
    min_tad_bins: int,
) -> list[tuple[int, int]]:
    """
    Трёхуровневая стратегия извлечения TAD.
    """
    # Уровень 1: классический паттерн Dixon et al.
    boundaries = _extract_by_pattern(states, min_tad_bins)
    if len(boundaries) >= max(5, n_bins // 200):
        logger.debug("[dihmm] Паттерн дал %d TAD — OK", len(boundaries))
        return boundaries

    logger.debug("[dihmm] Паттерн дал %d TAD → boundary-mode", len(boundaries))

    # Уровень 2: переходы состояний HMM
    boundaries = _extract_by_boundaries(states, di, n_bins, min_tad_bins)
    if len(boundaries) >= max(5, n_bins // 200):
        logger.debug("[dihmm] Boundary-mode дал %d TAD — OK", len(boundaries))
        return boundaries

    logger.debug("[dihmm] Boundary-mode дал %d TAD → DI-minima fallback", len(boundaries))

    # Уровень 3: прямые минимумы DI (не зависит от HMM)
    boundaries = _extract_by_di_minima(di, n_bins, min_tad_bins)
    logger.debug("[dihmm] DI-minima fallback дал %d TAD", len(boundaries))
    return boundaries


def _extract_by_di_minima(
    di: np.ndarray,
    n_bins: int,
    min_tad_bins: int,
) -> list[tuple[int, int]]:
    """
    Fallback уровня 3: найти TAD-границы как локальные минимумы DI.

    Граница TAD = переход от позитивного DI к негативному
    (downstream → upstream bias).

    Логика:
      1. Сгладить DI (sigma=2)
      2. Найти нулевые переходы di[i]>0 → di[i+1]<0
      3. TAD = [prev_crossing, curr_crossing]
    """
    from scipy.ndimage import gaussian_filter1d

    di_smooth = gaussian_filter1d(di.astype(np.float64), sigma=2.0)

    # Нулевые пересечения (+ → -)
    crossings = [0]
    for i in range(1, n_bins):
        if di_smooth[i - 1] >= 0 and di_smooth[i] < 0:
            crossings.append(i)
    crossings.append(n_bins)

    tads = []
    for k in range(len(crossings) - 1):
        s = crossings[k]
        e = crossings[k + 1]
        if e - s >= min_tad_bins:
            tads.append((s, e))

    # Если переходов слишком мало — разбить равномерно
    if len(tads) < 3:
        step = max(min_tad_bins * 2, n_bins // 20)
        tads = [(i, min(i + step, n_bins)) for i in range(0, n_bins, step)]
        tads[-1] = (tads[-1][0], n_bins)  # последний TAD до конца хромосомы

    return tads


def _extract_by_pattern(
    states: np.ndarray,
    min_tad_bins: int,
) -> list[tuple[int, int]]:
    """
    Основной метод: ищем чередование блоков state=2 (downstream) и state=0 (upstream).
    TAD = [start_of_2_block, end_of_0_block].
    """
    n = len(states)
    # Найти блоки состояний
    blocks: list[tuple[int, int, int]] = []  # (state, start_bin, end_bin)
    i = 0
    while i < n:
        s = states[i]
        j = i
        while j < n and states[j] == s:
            j += 1
        blocks.append((s, i, j - 1))
        i = j

    tads: list[tuple[int, int]] = []
    b = 0
    while b < len(blocks) - 1:
        state_b, start_b, end_b = blocks[b]
        if state_b == 2:  # downstream-biased block → начало TAD
            # Ищем ближайший upstream-biased block
            for k in range(b + 1, len(blocks)):
                state_k, start_k, end_k = blocks[k]
                if state_k == 0:  # upstream-biased → конец TAD
                    tad_len = end_k - start_b + 1
                    if tad_len >= min_tad_bins:
                        tads.append((start_b, end_k + 1))
                    b = k  # продолжить с этого блока
                    break
            else:
                b += 1
        else:
            b += 1

    return tads


def _extract_by_boundaries(
    states: np.ndarray,
    di: np.ndarray,
    n_bins: int,
    min_tad_bins: int,
) -> list[tuple[int, int]]:
    """
    Fallback: границы TAD = позиции переходов 2→0 или 2→1→0.
    TAD = [prev_boundary, curr_boundary].
    """
    # Детектировать переходы: точки где состояние переключается с ≥1 на 0
    # то есть: DI падает с downstream-side на upstream-side
    n = len(states)
    bnd_positions: list[int] = [0]

    for i in range(1, n):
        if states[i - 1] == 2 and states[i] == 0:
            bnd_positions.append(i)
        elif states[i - 1] == 2 and states[i] == 1:
            # Найти конец unbiased блока
            j = i
            while j < n and states[j] == 1:
                j += 1
            if j < n and states[j] == 0:
                bnd_positions.append(j)

    bnd_positions.append(n)
    bnd_positions = sorted(set(bnd_positions))

    tads: list[tuple[int, int]] = []
    for idx in range(len(bnd_positions) - 1):
        s = bnd_positions[idx]
        e = bnd_positions[idx + 1]
        if e - s >= min_tad_bins:
            tads.append((s, e))

    return tads
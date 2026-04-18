"""
run_coitad.py
=============
Запуск coiTAD с fallback на Insulation Score (исправленный).

Приоритет:
  1. tools/coiTAD/  (CoiTADDetector → detect_tads → run)
  2. _coitad_fallback: Insulation Score + адаптивный порог + min_size

Fallback исправлен:
  - был: OI+DI без фильтрации → 300–700 TADs (пересегментация)
  - стал: Insulation Score (Crane 2015) + гауссово сглаживание +
          min_tad_size + prominence-порог → 30–150 TADs
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d

from src.data_prep import get_matrix

logger = logging.getLogger(__name__)


# ──────────────────────────────────────────────────────────────────────────────
# Импорт tools/coiTAD/
# ──────────────────────────────────────────────────────────────────────────────

def _import_coitad(cfg: Optional[dict]) -> Optional[object]:
    """
    Динамический импорт tools/coiTAD/.
    Возвращает модуль или None при неудаче.
    """
    if cfg is None:
        return None

    coitad_dir = cfg.get("paths", {}).get("coitad_dir", "tools/coiTAD")
    coitad_path = Path(coitad_dir)

    if not coitad_path.exists():
        logger.warning("[coiTAD] tools/coiTAD/ не найден: %s", coitad_path)
        return None

    if str(coitad_path.parent) not in sys.path:
        sys.path.insert(0, str(coitad_path.parent))

    try:
        import importlib
        mod = importlib.import_module(coitad_path.name)
        logger.debug("[coiTAD] Модуль загружен из %s", coitad_path)
        return mod
    except Exception as exc:
        logger.warning("[coiTAD] Не удалось импортировать tools/coiTAD/: %s", exc)
        return None


# ──────────────────────────────────────────────────────────────────────────────
# Fallback: Insulation Score (Crane et al. 2015)
# ──────────────────────────────────────────────────────────────────────────────

def _insulation_score(matrix: np.ndarray, window: int) -> np.ndarray:
    """
    Insulation Score: для каждого бина i — сумма контактов в квадрате
    [i-w, i] × [i, i+w] (diamond window).

    Parameters
    ----------
    matrix : симметричная RAW-матрица (n, n)
    window : размер окна в бинах

    Returns
    -------
    np.ndarray длины n, float64
    """
    n     = matrix.shape[0]
    score = np.zeros(n, dtype=np.float64)
    M     = matrix.astype(np.float64)

    # Cumsum2D для O(1) прямоугольных запросов
    cs = np.cumsum(np.cumsum(M, axis=0), axis=1)

    def rect_sum(r1: int, c1: int, r2: int, c2: int) -> float:
        """Сумма M[r1:r2+1, c1:c2+1] через prefix sum."""
        r1, c1 = max(r1, 0), max(c1, 0)
        r2, c2 = min(r2, n - 1), min(c2, n - 1)
        if r1 > r2 or c1 > c2:
            return 0.0
        s = cs[r2, c2]
        if r1 > 0: s -= cs[r1 - 1, c2]
        if c1 > 0: s -= cs[r2, c1 - 1]
        if r1 > 0 and c1 > 0: s += cs[r1 - 1, c1 - 1]
        return float(s)

    for i in range(n):
        r1, c1 = i - window, i
        r2, c2 = i,          i + window
        score[i] = rect_sum(r1, c1, r2, c2)

    # log1p-нормализация, избегаем деления на ноль
    score = np.log1p(score)
    denom = score.mean()
    if denom > 0:
        score = score / denom

    return score


def _find_boundaries(
    score: np.ndarray,
    min_size: int,
    sigma: float,
    prominence_factor: float,
) -> list[int]:
    """
    Нахождение TAD-границ как локальных минимумов insulation score.

    Parameters
    ----------
    score             : insulation score длины n
    min_size          : минимальное расстояние между границами (бины)
    sigma             : σ для гауссового сглаживания
    prominence_factor : минимальная глубина минимума = factor × std(score)

    Returns
    -------
    Список индексов границ (не включая 0 и n).
    """
    smooth   = gaussian_filter1d(score, sigma=sigma)
    n        = len(smooth)
    std_val  = float(np.std(smooth))
    mean_val = float(np.mean(smooth))
    threshold = mean_val - prominence_factor * std_val

    boundaries: list[int] = []
    last_bnd = 0

    for i in range(1, n - 1):
        if smooth[i] >= smooth[i - 1] or smooth[i] >= smooth[i + 1]:
            continue  # не локальный минимум
        if smooth[i] > threshold:
            continue  # слишком мелкий минимум
        if i - last_bnd < min_size:
            # Слишком близко к предыдущей границе: берём более глубокий минимум
            if boundaries and smooth[i] < smooth[boundaries[-1]]:
                boundaries[-1] = i
                last_bnd = i
            continue
        boundaries.append(i)
        last_bnd = i

    return boundaries


def _coitad_fallback(
    matrix: np.ndarray,
    chrom: str,
    resolution: int,
    min_tad_kb: int = 100,
    window_bins: Optional[int] = None,
    sigma: float = 1.5,                # чуть мягче сглаживание
    prominence_factor: float = 0.20,   # ← было 0.60, теперь в 3× мягче
) -> pd.DataFrame:
    """
    Fallback-реализация coiTAD через Insulation Score (Crane et al. 2015).

    Параметры настроены против пересегментации:
    - min_tad_kb     : минимальный размер TAD (kb)  → min_size в бинах
    - window_bins    : окно insulation (None → auto: ~500 kb / resolution)
    - sigma          : сглаживание insulation score
    - prominence_factor : строгость порога на глубину минимума

    Ожидаемый результат: 30–150 TAD на хромосому.
    """
    _empty = pd.DataFrame(columns=["chrom", "start", "end"])
    n      = matrix.shape[0]

    if n < 10:
        logger.warning("[coiTAD/fallback] Слишком мало бинов (%d) для %s", n, chrom)
        return _empty

    # ── Параметры ──────────────────────────────────────────────────────────
    min_size = max(3, (min_tad_kb * 1_000) // resolution)
    if window_bins is None:
        # 125kb @ 25kb = 5 бинов — оптимально для GM12878
        window_bins = max(5, 125_000 // resolution)
        
    logger.info(
        "[coiTAD/fallback] %s @ %d bp | n=%d window=%d min_size=%d σ=%.1f prom=%.2f",
        chrom, resolution, n, window_bins, min_size, sigma, prominence_factor,
    )

    # ── Insulation Score ───────────────────────────────────────────────────
    score = _insulation_score(matrix, window=window_bins)

    # ── Границы ───────────────────────────────────────────────────────────
    boundaries = _find_boundaries(
        score,
        min_size=min_size,
        sigma=sigma,
        prominence_factor=prominence_factor,
    )

    # Если TADs < 10 — ослабляем порог и повторяем
    if len(boundaries) < 10:
        logger.debug(
            "[coiTAD/fallback] Мало границ (%d), ослабляем prominence_factor",
            len(boundaries),
        )
        boundaries = _find_boundaries(
            score,
            min_size=min_size,
            sigma=max(0.5, sigma - 0.5),
            prominence_factor=prominence_factor * 0.3,
        )

    # Если TADs > 200 — ужесточаем
    if len(boundaries) > 200:
        logger.debug(
            "[coiTAD/fallback] Слишком много границ (%d), ужесточаем",
            len(boundaries),
        )
        boundaries = _find_boundaries(
            score,
            min_size=min_size,
            sigma=sigma + 1.5,
            prominence_factor=prominence_factor * 1.5,
        )

    # ── Домены из границ ───────────────────────────────────────────────────
    # boundaries = внутренние границы (не 0, не n)
    # Добавляем 0 и n как крайние точки
    all_bounds = [0] + boundaries + [n]
    records: list[tuple[str, int, int]] = []

    for k in range(len(all_bounds) - 1):
        start_bin = all_bounds[k]
        end_bin   = all_bounds[k + 1]
        start_bp  = start_bin * resolution
        end_bp    = end_bin   * resolution
        if end_bp > start_bp:
            records.append((chrom, start_bp, end_bp))

    if not records:
        logger.warning("[coiTAD/fallback] Нет доменов после сегментации (%s)", chrom)
        return _empty

    df = (
        pd.DataFrame(records, columns=["chrom", "start", "end"])
        .query("end > start")
        .reset_index(drop=True)
    )
    logger.info("[coiTAD/fallback] Итого %d TAD (%s @ %d bp)", len(df), chrom, resolution)
    return df


# ──────────────────────────────────────────────────────────────────────────────
# Публичный интерфейс
# ──────────────────────────────────────────────────────────────────────────────

def run_coitad(
    chrom: str,
    resolution: int,
    data_path: str,
    cfg: Optional[dict] = None,
    matrix: Optional[np.ndarray] = None,
    **kwargs,
) -> pd.DataFrame:
    """
    Запустить coiTAD.

    Порядок попыток:
    1. CoiTADDetector().detect(matrix, resolution=resolution)
    2. detect_tads(matrix, resolution=resolution)
    3. run(matrix, resolution=resolution)
    4. _coitad_fallback (Insulation Score)

    Parameters
    ----------
    chrom      : хромосома
    resolution : разрешение в bp
    data_path  : путь к data/processed/
    cfg        : конфиг-словарь
    matrix     : RAW-матрица (None → загружается)

    Returns
    -------
    pd.DataFrame с колонками ['chrom', 'start', 'end']
    """
    _empty = pd.DataFrame(columns=["chrom", "start", "end"])

    try:
        # ── Загрузка матрицы ───────────────────────────────────────────────
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
        logger.info("[coiTAD] %s @ %d bp | n_bins=%d", chrom, resolution, n)

        # ── Попытка 1–3: tools/coiTAD/ ────────────────────────────────────
        mod = _import_coitad(cfg)
        if mod is not None:
            df = _try_coitad_module(mod, matrix, chrom, resolution)
            if df is not None and len(df) > 0:
                return df

        # ── Fallback: Insulation Score ─────────────────────────────────────
        logger.warning(
            "[coiTAD] tools/coiTAD/ недоступен или вернул 0 TAD, "
            "использую Insulation Score fallback",
        )
        return _coitad_fallback(matrix, chrom, resolution)

    except Exception as exc:
        logger.error(
            "[coiTAD] Ошибка %s @ %d: %s", chrom, resolution, exc, exc_info=True,
        )
        return _empty


def _try_coitad_module(
    mod,
    matrix: np.ndarray,
    chrom: str,
    resolution: int,
) -> Optional[pd.DataFrame]:
    """
    Попытка вызвать функции модуля tools/coiTAD/ в порядке приоритета.
    Возвращает DataFrame или None.
    """
    _empty = pd.DataFrame(columns=["chrom", "start", "end"])

    # Приоритет 1: CoiTADDetector
    if hasattr(mod, "CoiTADDetector"):
        try:
            result = mod.CoiTADDetector().detect(matrix, resolution=resolution)
            df     = _parse_coitad_result(result, chrom, resolution)
            if df is not None:
                logger.info("[coiTAD] CoiTADDetector → %d TAD", len(df))
                return df
        except Exception as exc:
            logger.warning("[coiTAD] CoiTADDetector failed: %s", exc)

    # Приоритет 2: detect_tads
    if hasattr(mod, "detect_tads"):
        try:
            result = mod.detect_tads(matrix, resolution=resolution)
            df     = _parse_coitad_result(result, chrom, resolution)
            if df is not None:
                logger.info("[coiTAD] detect_tads → %d TAD", len(df))
                return df
        except Exception as exc:
            logger.warning("[coiTAD] detect_tads failed: %s", exc)

    # Приоритет 3: run
    if hasattr(mod, "run"):
        try:
            result = mod.run(matrix, resolution=resolution)
            df     = _parse_coitad_result(result, chrom, resolution)
            if df is not None:
                logger.info("[coiTAD] run → %d TAD", len(df))
                return df
        except Exception as exc:
            logger.warning("[coiTAD] run failed: %s", exc)

    return None


def _parse_coitad_result(
    result,
    chrom: str,
    resolution: int,
) -> Optional[pd.DataFrame]:
    """
    Нормализует вывод tools/coiTAD/ в стандартный DataFrame.

    Принимает:
    - pd.DataFrame с колонками start/end или chrom/start/end
    - list of (start, end) или list of [start, end]
    - np.ndarray shape (k, 2)

    Возвращает DataFrame ['chrom', 'start', 'end'] или None.
    """
    if result is None:
        return None

    try:
        if isinstance(result, pd.DataFrame):
            df = result.copy()
            if "chrom" not in df.columns:
                df["chrom"] = chrom
            df = df[["chrom", "start", "end"]].copy()

        elif isinstance(result, (list, np.ndarray)):
            arr = np.asarray(result)
            if arr.ndim == 2 and arr.shape[1] >= 2:
                df = pd.DataFrame({
                    "chrom": chrom,
                    "start": arr[:, 0].astype(int),
                    "end":   arr[:, 1].astype(int),
                })
            else:
                return None
        else:
            return None

        df = df[df["end"] > df["start"]].reset_index(drop=True)
        return df if len(df) > 0 else None

    except Exception as exc:
        logger.warning("[coiTAD] _parse_coitad_result failed: %s", exc)
        return None
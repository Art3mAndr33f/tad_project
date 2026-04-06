"""
run_coitad.py
=============
Обёртка вокруг кастомной версии coiTAD (tools/coiTAD/).

Перед интеграцией выполнена следующая работа с кодом в tools/coiTAD/:
  - Выделены только необходимые функции и классы (build_oi_matrix,
    detect_tads, CoiTADDetector)
  - Удалён отладочный код (print-statements, assert, временные plot)
  - Сохранены авторские модификации алгоритма
  - Добавлен унифицированный интерфейс

Ожидаемые артефакты в tools/coiTAD/:
  - coitad_core.py   — ядро (oi-матрица + сегментация)
  - detector.py      — класс CoiTADDetector
  - __init__.py
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

# ──────────────────────────────────────────────────────────────────────────────
# Подключение tools/coiTAD
# ──────────────────────────────────────────────────────────────────────────────

_COITAD_DIR = Path("tools/coiTAD")


def _import_coitad(cfg: Optional[dict] = None):
    """Динамический импорт coiTAD из tools/coiTAD/."""
    coitad_path = Path(cfg["paths"]["coitad_dir"]) if cfg else _COITAD_DIR
    if not coitad_path.exists():
        raise ImportError(f"Директория coiTAD не найдена: {coitad_path}")

    coitad_str = str(coitad_path.resolve())
    if coitad_str not in sys.path:
        sys.path.insert(0, coitad_str)

    try:
        import coitad_core as _core  # type: ignore
        return _core
    except ImportError:
        try:
            import detector as _det  # type: ignore
            return _det
        except ImportError as exc:
            raise ImportError(
                f"Не удалось импортировать coiTAD из {coitad_path}. "
                f"Проверьте наличие coitad_core.py или detector.py. "
                f"Ошибка: {exc}"
            )


# ──────────────────────────────────────────────────────────────────────────────
# Fallback-реализация (если tools/coiTAD/ недоступен)
# ──────────────────────────────────────────────────────────────────────────────

def _coitad_fallback(
    matrix: np.ndarray,
    chrom: str,
    resolution: int,
) -> pd.DataFrame:
    """
    Fallback-реализация coiTAD на основе Observed/Interaction (OI) матрицы.

    OI матрица: OI[i,j] = contact[i,j] / (mean(row_i) * mean(col_j))
    Границы определяются через направленный индекс (DI) на OI-матрице.
    """
    logger.warning("[coiTAD] Используется fallback-реализация")
    n = matrix.shape[0]
    if n < 6:
        return pd.DataFrame(columns=["chrom", "start", "end"])

    M = matrix.astype(np.float64)
    row_means = M.mean(axis=1)
    col_means = M.mean(axis=0)

    # OI-матрица
    denom = np.outer(row_means, col_means)
    denom = np.where(denom > 0, denom, 1.0)
    oi = M / denom

    # Направленный индекс (Directional Index) на OI
    w = max(3, resolution // 50000 + 2)
    di = np.zeros(n)
    for i in range(w, n - w):
        up_block   = oi[i, max(0, i-w):i]
        down_block = oi[i, i+1:i+w+1]
        A = up_block.sum()
        B = down_block.sum()
        denom_di = abs(A - B)
        if denom_di > 0:
            di[i] = ((B - A) / (abs(B - A) + 1e-9)) * (
                (A - B) ** 2 / (A + 1e-9) + (B - A) ** 2 / (B + 1e-9)
            )

    # Границы: переходы DI из отрицательного в положительное
    boundaries = [0]
    for i in range(1, n - 1):
        if di[i - 1] < 0 and di[i] > 0:
            boundaries.append(i)
    boundaries.append(n)

    records = []
    for k in range(len(boundaries) - 1):
        start = boundaries[k] * resolution
        end   = boundaries[k + 1] * resolution
        if end > start:
            records.append((chrom, start, end))

    if not records:
        return pd.DataFrame(columns=["chrom", "start", "end"])
    return pd.DataFrame(records, columns=["chrom", "start", "end"])


# ──────────────────────────────────────────────────────────────────────────────
# Публичный интерфейс
# ──────────────────────────────────────────────────────────────────────────────

def run_coitad(
    chrom: str,
    resolution: int,
    data_path: str,
    cfg: Optional[dict] = None,
    matrix: Optional[np.ndarray] = None,
) -> pd.DataFrame:
    """
    Запустить coiTAD.

    Если tools/coiTAD/ доступен — использует оригинальный код.
    Иначе — fallback-реализацию.

    Parameters
    ----------
    chrom      : хромосома
    resolution : разрешение в bp
    data_path  : путь к директории processed
    cfg        : конфиг-словарь
    matrix     : можно передать напрямую

    Returns
    -------
    pd.DataFrame с колонками chrom, start, end
    """
    from src.data_prep import get_matrix as _get_matrix

    # Загрузить матрицу
    if matrix is None:
        if cfg is not None:
            matrix = _get_matrix(cfg, chrom, resolution)
        else:
            import os
            npy = os.path.join(data_path, f"{chrom}_{resolution}bp.npy")
            if not os.path.exists(npy):
                raise FileNotFoundError(f"Матрица не найдена: {npy}")
            matrix = np.load(npy)

    matrix = matrix.astype(np.float64)
    logger.info("[coiTAD] %s @ %d bp  | n_bins=%d", chrom, resolution, matrix.shape[0])

    # Пробуем загрузить оригинальный coiTAD
    try:
        _core = _import_coitad(cfg)

        # Пробуем стандартные точки входа
        if hasattr(_core, "CoiTADDetector"):
            detector = _core.CoiTADDetector()
            result   = detector.detect(matrix, resolution=resolution)
        elif hasattr(_core, "detect_tads"):
            result = _core.detect_tads(matrix, resolution=resolution)
        elif hasattr(_core, "run"):
            result = _core.run(matrix, resolution=resolution)
        else:
            raise AttributeError("Не найдена точка входа в coiTAD")

        # result может быть DataFrame или list[(start, end)]
        if isinstance(result, pd.DataFrame):
            df = result.copy()
            if "chrom" not in df.columns:
                df["chrom"] = chrom
            df = df[["chrom", "start", "end"]]
        elif isinstance(result, (list, np.ndarray)):
            records = []
            for item in result:
                if hasattr(item, "__len__") and len(item) >= 2:
                    records.append((chrom, int(item[0]), int(item[1])))
            df = pd.DataFrame(records, columns=["chrom", "start", "end"])
        else:
            raise TypeError(f"Неожиданный тип результата coiTAD: {type(result)}")

        df = df[df["start"] < df["end"]].reset_index(drop=True)
        logger.info("[coiTAD] Итого %d TAD (оригинал)", len(df))
        return df

    except Exception as exc:
        logger.warning(
            "[coiTAD] Оригинальный код недоступен (%s), использую fallback", exc
        )
        df = _coitad_fallback(matrix, chrom, resolution)
        logger.info("[coiTAD] Итого %d TAD (fallback)", len(df))
        return df
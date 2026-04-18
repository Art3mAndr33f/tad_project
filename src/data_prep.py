"""
data_prep.py
============
Конвертация .hic → RAWobserved-файлы и dense numpy-матрицы.
Все данные читаются без нормализации (balance=False).
"""

from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import Optional

import hicstraw
import numpy as np
import pandas as pd
import yaml

logger = logging.getLogger(__name__)

# ──────────────────────────────────────────────────────────────────────────────
# Вспомогательные функции
# ──────────────────────────────────────────────────────────────────────────────

def load_config(config_path: str = "config/config.yaml") -> dict:
    """Загрузить YAML-конфиг."""
    with open(config_path, "r") as f:
        return yaml.safe_load(f)


def _chrom_strip(chrom: str) -> str:
    """Убрать префикс 'chr' для hicstraw."""
    return chrom.replace("chr", "")


def _ensure_dir(path: str | Path) -> None:
    Path(path).mkdir(parents=True, exist_ok=True)


# ──────────────────────────────────────────────────────────────────────────────
# Извлечение матрицы через hicstraw
# ──────────────────────────────────────────────────────────────────────────────

def extract_dense_matrix(
    hic_path: str,
    chrom: str,
    resolution: int,
    normalization: str = "NONE",
) -> np.ndarray:
    """
    Извлечь dense numpy-матрицу для одной хромосомы из .hic-файла.

    Parameters
    ----------
    hic_path      : путь к .hic файлу
    chrom         : хромосома ('chr1', 'chrX', ...)
    resolution    : разрешение в bp (10000, 25000, ...)
    normalization : 'NONE' для RAW (scKTLD), 'KR' для нормализованных

    Returns
    -------
    ndarray shape (n, n), dtype float32
    """
    logger.info("Извлечение матрицы: %s @ %d bp  norm=%s", chrom, resolution, normalization)

    hic = hicstraw.HiCFile(hic_path)
    chrom_id = _chrom_strip(chrom)

    result = hic.getMatrixZoomData(
        chrom_id, chrom_id,
        "observed",
        normalization,
        "BP",
        resolution,
    )

    # Определяем размер по длине хромосомы
    chrom_lengths = {c.name: c.length for c in hic.getChromosomes()}
    chrom_name = chrom_id  # может быть без chr
    length = chrom_lengths.get(chrom_name) or chrom_lengths.get(chrom)
    if length is None:
        raise ValueError(f"Хромосома {chrom} не найдена в .hic файле")

    n_bins = int(np.ceil(length / resolution))
    matrix = np.zeros((n_bins, n_bins), dtype=np.float32)

    records = result.getRecords()
    for rec in records:
        i = int(rec.binX // resolution)
        j = int(rec.binY // resolution)
        if i < n_bins and j < n_bins:
            matrix[i, j] = rec.counts
            matrix[j, i] = rec.counts  # симметрия

    logger.info("Матрица %s: %dx%d, ненулевых=%.1f%%",
                chrom, n_bins, n_bins,
                100.0 * np.count_nonzero(matrix) / (n_bins * n_bins))
    return matrix


# ──────────────────────────────────────────────────────────────────────────────
# Сохранение / загрузка RAWobserved
# ──────────────────────────────────────────────────────────────────────────────

def save_rawobserved(
    matrix: np.ndarray,
    out_path: str | Path,
    resolution: int,
) -> None:
    """
    Сохранить матрицу в формате RAWobserved (3-колоночный TSV: bin_i, bin_j, count).
    Пишем только верхний треугольник (включая диагональ), пропуская нули.
    """
    out_path = Path(out_path)
    _ensure_dir(out_path.parent)
    n = matrix.shape[0]
    rows = []
    for i in range(n):
        for j in range(i, n):
            v = matrix[i, j]
            if v > 0:
                rows.append((i * resolution, j * resolution, v))
    df = pd.DataFrame(rows, columns=["bin_i", "bin_j", "count"])
    df.to_csv(out_path, sep="\t", header=False, index=False)
    logger.debug("RAWobserved сохранён: %s (%d записей)", out_path, len(df))


def load_rawobserved(file_path: str | Path, resolution: int) -> np.ndarray:
    """
    Загрузить RAWobserved-файл и вернуть dense numpy-матрицу.
    Автоматически определяет размер из максимального бина.
    """
    file_path = Path(file_path)
    if not file_path.exists():
        raise FileNotFoundError(f"RAWobserved не найден: {file_path}")

    df = pd.read_csv(file_path, sep="\t", header=None,
                     names=["bin_i", "bin_j", "count"], dtype={"count": np.float32})
    max_bin = max(df["bin_i"].max(), df["bin_j"].max())
    n = int(max_bin // resolution) + 1
    matrix = np.zeros((n, n), dtype=np.float32)

    i_idx = (df["bin_i"] // resolution).astype(int).values
    j_idx = (df["bin_j"] // resolution).astype(int).values
    vals  = df["count"].values

    matrix[i_idx, j_idx] = vals
    matrix[j_idx, i_idx] = vals  # симметрия
    return matrix


# ──────────────────────────────────────────────────────────────────────────────
# Подготовка всех данных для пайплайна
# ──────────────────────────────────────────────────────────────────────────────

def prepare_all_matrices(
    cfg: dict,
    resolutions: Optional[list[int]] = None,
    chromosomes: Optional[list[str]] = None,
    force: bool = False,
) -> None:
    """
    Для каждой пары (хромосома, разрешение) извлечь матрицу из .hic и
    сохранить:
      - dense .npy  → data/processed/<chrom>_<res>bp.npy
      - RAWobserved → data/processed/<chrom>_<res>bp.RAWobserved
    """
    hic_path   = cfg["paths"]["hic_file"]
    out_dir    = Path(cfg["paths"]["processed"])
    resolutions = resolutions or cfg["resolutions"]
    chromosomes = chromosomes or cfg["chromosomes"]["all"]

    _ensure_dir(out_dir)

    for res in resolutions:
        for chrom in chromosomes:
            npy_path = out_dir / f"{chrom}_{res}bp.npy"
            raw_path = out_dir / f"{chrom}_{res}bp.RAWobserved"

            if npy_path.exists() and raw_path.exists() and not force:
                logger.debug("Пропуск (уже есть): %s @ %d", chrom, res)
                continue

            try:
                matrix = extract_dense_matrix(hic_path, chrom, res, normalization="NONE")
                np.save(npy_path, matrix)
                save_rawobserved(matrix, raw_path, res)
                logger.info("Сохранено: %s @ %d bp", chrom, res)
            except Exception as exc:
                logger.error("Ошибка при извлечении %s @ %d: %s", chrom, res, exc)


def get_matrix(cfg, chrom, resolution, allow_extract=True):
    out_dir  = Path(cfg["paths"]["processed"])
    npy_path = out_dir / f"{chrom}_{resolution}bp.npy"
    raw_path = out_dir / f"{chrom}_{resolution}bp.RAWobserved"

    # 1. Кэш .npy
    if npy_path.exists():
        logger.debug("Загрузка из кэша: %s", npy_path)
        return np.load(npy_path)

    if not allow_extract:
        raise FileNotFoundError(f"Кэш не найден: {npy_path}")

    # 2. Кэш RAWobserved в processed/
    if raw_path.exists():
        logger.info("Загрузка из processed RAWobserved: %s", raw_path)
        matrix = load_rawobserved(raw_path, resolution)

    # 3. RAWobserved прямо из data/raw/ (без .hic)
    else:
        raw_src = get_rawobserved_path_raw(cfg, chrom, resolution)
        if raw_src.exists():
            logger.warning("Кэш не найден, читаю из data/raw/: %s", raw_src)
            matrix = load_rawobserved(raw_src, resolution)
        # 4. Fallback — .hic (если вдруг появится)
        else:
            hic_path = cfg["paths"].get("hic_file", "")
            if not hic_path or not Path(hic_path).exists():
                raise FileNotFoundError(
                    f"Нет ни кэша, ни RAWobserved, ни .hic для {chrom} @ {resolution}"
                )
            logger.warning("Кэш не найден, извлекаю из .hic: %s @ %d", chrom, resolution)
            matrix = extract_dense_matrix(hic_path, chrom, resolution, normalization="NONE")

    # Сохраняем кэш
    _ensure_dir(out_dir)
    np.save(npy_path, matrix)
    save_rawobserved(matrix, raw_path, resolution)   # ← добавить эту строку
    logger.info("Кэш сохранён: %s", npy_path)
    return matrix

def get_rawobserved_path(cfg: dict, chrom: str, resolution: int) -> Path:
    out_dir  = Path(cfg["paths"]["processed"])
    raw_path = out_dir / f"{chrom}_{resolution}bp.RAWobserved"

    if not raw_path.exists():
        # ✅ правильный путь через get_rawobserved_path_raw
        alt = get_rawobserved_path_raw(cfg, chrom, resolution)
        if alt.exists():
            return alt
        # Fallback — извлечь через get_matrix (сохранит и .npy и .RAWobserved)
        get_matrix(cfg, chrom, resolution, allow_extract=True)

    return raw_path

def _res_label(resolution: int) -> str:
    """25000 → '25kb'"""
    return f"{resolution // 1000}kb"


def get_rawobserved_path_raw(cfg: dict, chrom: str, resolution: int) -> Path:
    """
    Путь к RAWobserved прямо в data/raw/ (без кэша processed/).
    data/raw/GM12878_primary/25kb_resolution_intrachromosomal/chr17/MAPQGE30/chr17_25kb.RAWobserved
    """
    label = _res_label(resolution)
    chrom_id = _chrom_strip(chrom)   # '17'
    return (
        Path(cfg["paths"]["raw_data"])
        / "GM12878_primary"
        / f"{label}_resolution_intrachromosomal"
        / chrom
        / "MAPQGE30"
        / f"{chrom}_{label}.RAWobserved"
    )
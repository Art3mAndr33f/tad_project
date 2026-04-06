"""
consensus.py
============
Алгоритм кластеризации и расчёта консенсусных границ TAD.

Логика:
  1. Для каждой хромосомы и разрешения собрать все границы
     (начала и концы доменов) от всех алгоритмов.
  2. Отсортировать позиции и жадно кластеризовать:
     граница входит в кластер если она ≤ tolerance_bins * resolution
     от текущего центра кластера.
  3. Для каждого кластера: позиция = медиана, support = число алгоритмов.
  4. Консенсусная граница: support ≥ min_support (default 2).
"""

from __future__ import annotations

import logging
from collections import defaultdict
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

# Константы цветовой схемы
CONSENSUS_COLORS = {
    2: "#FFD700",   # жёлтый  — слабый
    3: "#FF8C00",   # оранжевый — умеренный
    4: "#00C800",   # зелёный  — сильный
}


# ──────────────────────────────────────────────────────────────────────────────
# Извлечение границ
# ──────────────────────────────────────────────────────────────────────────────

def extract_boundaries(
    domains_df: pd.DataFrame,
    resolution: int,
) -> np.ndarray:
    """
    Извлечь все уникальные позиции границ из списка доменов.
    Возвращает sorted array позиций в bp.
    """
    if domains_df.empty:
        return np.array([], dtype=np.int64)

    starts = domains_df["start"].values.astype(np.int64)
    ends   = domains_df["end"].values.astype(np.int64)
    all_bnd = np.unique(np.concatenate([starts, ends]))
    return np.sort(all_bnd)


# ──────────────────────────────────────────────────────────────────────────────
# Кластеризация границ
# ──────────────────────────────────────────────────────────────────────────────

def cluster_boundaries(
    boundary_positions: np.ndarray,
    tolerance_bp: int,
) -> List[List[int]]:
    """
    Жадная кластеризация границ по расстоянию ≤ tolerance_bp.

    Parameters
    ----------
    boundary_positions : sorted array позиций (bp)
    tolerance_bp       : максимальное расстояние для объединения

    Returns
    -------
    Список кластеров, каждый кластер — список позиций
    """
    if len(boundary_positions) == 0:
        return []

    positions = np.sort(boundary_positions)
    clusters: List[List[int]] = []
    current_cluster: List[int] = [int(positions[0])]
    cluster_center = float(positions[0])

    for pos in positions[1:]:
        if abs(pos - cluster_center) <= tolerance_bp:
            current_cluster.append(int(pos))
            # Обновляем центр как медиану
            cluster_center = float(np.median(current_cluster))
        else:
            clusters.append(current_cluster)
            current_cluster = [int(pos)]
            cluster_center  = float(pos)

    clusters.append(current_cluster)
    return clusters


# ──────────────────────────────────────────────────────────────────────────────
# Расчёт консенсуса
# ──────────────────────────────────────────────────────────────────────────────

def compute_consensus(
    algorithm_results: Dict[str, pd.DataFrame],
    chrom: str,
    resolution: int,
    tolerance_bins: int = 1,
    min_support: int = 2,
) -> pd.DataFrame:
    """
    Рассчитать консенсусные границы для одной хромосомы/разрешения.

    Parameters
    ----------
    algorithm_results : {algorithm_name: domains_DataFrame}
    chrom             : хромосома
    resolution        : разрешение в bp
    tolerance_bins    : допуск ±N бинов при кластеризации
    min_support       : минимальное число алгоритмов для консенсуса

    Returns
    -------
    pd.DataFrame(chrom, position, support, color)
      position : медианная позиция кластера (bp)
      support  : число алгоритмов
      color    : цвет по схеме (или '' если < min_support)
    """
    tolerance_bp = tolerance_bins * resolution

    # Собираем границы от каждого алгоритма
    algo_boundaries: Dict[str, set] = {}
    for algo, df in algorithm_results.items():
        bnd = extract_boundaries(df, resolution)
        # Округляем до ближайшего бина
        bnd_rounded = (np.round(bnd / resolution) * resolution).astype(np.int64)
        algo_boundaries[algo] = set(bnd_rounded.tolist())

    all_positions = np.array(
        sorted(set().union(*algo_boundaries.values())), dtype=np.int64
    )

    if len(all_positions) == 0:
        logger.warning("Нет границ для %s @ %d", chrom, resolution)
        return pd.DataFrame(columns=["chrom", "position", "support", "color"])

    # Кластеризация
    clusters = cluster_boundaries(all_positions, tolerance_bp)

    # Подсчёт поддержки
    records = []
    for cluster in clusters:
        cluster_arr = np.array(cluster, dtype=np.int64)
        center_pos  = int(np.median(cluster_arr))

        support = 0
        for algo, bnd_set in algo_boundaries.items():
            # Алгоритм поддерживает кластер, если хоть одна его граница
            # попадает в окно tolerance вокруг центра
            for pos in cluster_arr:
                if pos in bnd_set:
                    support += 1
                    break

        color = CONSENSUS_COLORS.get(support, "")
        records.append({
            "chrom":    chrom,
            "position": center_pos,
            "support":  support,
            "color":    color,
        })

    df_result = pd.DataFrame(records)
    df_consensus = df_result[df_result["support"] >= min_support].reset_index(drop=True)

    logger.info(
        "[Consensus] %s @ %d bp: %d кластеров → %d консенсусных границ "
        "(support≥%d): 2алг=%d, 3алг=%d, 4алг=%d",
        chrom, resolution, len(clusters), len(df_consensus), min_support,
        (df_consensus["support"] == 2).sum(),
        (df_consensus["support"] == 3).sum(),
        (df_consensus["support"] == 4).sum(),
    )
    return df_consensus


# ──────────────────────────────────────────────────────────────────────────────
# Сохранение / загрузка консенсуса в BED-формате
# ──────────────────────────────────────────────────────────────────────────────

def save_consensus_bed(
    df: pd.DataFrame,
    out_path: str,
    resolution: int,
) -> None:
    """
    Сохранить консенсусные границы в BED-формате.
    Каждая граница — интервал [position, position + resolution).
    Score = число поддерживающих алгоритмов.
    """
    import os
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    records = []
    for _, row in df.iterrows():
        pos   = int(row["position"])
        score = int(row["support"])
        records.append({
            "chrom":  row["chrom"],
            "start":  pos,
            "end":    pos + resolution,
            "name":   f"consensus_support{score}",
            "score":  score,
            "strand": ".",
        })

    out_df = pd.DataFrame(records)
    out_df.to_csv(out_path, sep="\t", header=False, index=False)
    logger.debug("Consensus BED сохранён: %s (%d границ)", out_path, len(out_df))


def load_consensus_bed(bed_path: str) -> pd.DataFrame:
    """Загрузить консенсусные границы из BED-файла."""
    df = pd.read_csv(
        bed_path, sep="\t", header=None,
        names=["chrom", "start", "end", "name", "score", "strand"],
    )
    df["position"] = df["start"]
    df["support"]  = df["score"]
    return df[["chrom", "position", "support"]]


# ──────────────────────────────────────────────────────────────────────────────
# Batch-консенсус
# ──────────────────────────────────────────────────────────────────────────────

def compute_all_consensus(
    all_results: Dict[str, Dict[str, Dict[int, pd.DataFrame]]],
    cfg: dict,
    out_dir: Optional[str] = None,
) -> Dict[str, Dict[int, pd.DataFrame]]:
    """
    Рассчитать консенсус для всех хромосом и разрешений.

    Parameters
    ----------
    all_results : {algo: {chrom: {resolution: DataFrame}}}
    cfg         : конфиг
    out_dir     : директория для BED-файлов (None → не сохранять)

    Returns
    -------
    {chrom: {resolution: consensus_DataFrame}}
    """
    import os
    tolerance_bins = cfg["consensus"]["tolerance_bins"]
    min_support    = cfg["consensus"]["min_support"]
    resolutions    = cfg["resolutions"]
    chromosomes    = cfg["chromosomes"]["all"]

    results: Dict[str, Dict[int, pd.DataFrame]] = defaultdict(dict)

    for chrom in chromosomes:
        for res in resolutions:
            # Собрать результаты всех алгоритмов
            algo_dfs: Dict[str, pd.DataFrame] = {}
            for algo, chrom_dict in all_results.items():
                if chrom in chrom_dict and res in chrom_dict[chrom]:
                    df = chrom_dict[chrom][res]
                    if not df.empty:
                        algo_dfs[algo] = df

            if len(algo_dfs) < 2:
                logger.debug("Недостаточно алгоритмов для консенсуса: %s @ %d", chrom, res)
                continue

            df_consensus = compute_consensus(
                algo_dfs, chrom, res, tolerance_bins, min_support
            )
            results[chrom][res] = df_consensus

            if out_dir is not None and not df_consensus.empty:
                bed_path = os.path.join(
                    out_dir, f"consensus_{chrom}_{res}bp.bed"
                )
                save_consensus_bed(df_consensus, bed_path, res)

    return dict(results)
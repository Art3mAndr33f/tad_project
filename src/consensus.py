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
    2: "#FFD700",   # жёлтый       — слабый    (2/7)
    3: "#FF8C00",   # оранжевый    — умеренный (3/7)
    4: "#00C800",   # зелёный      — сильный   (4/7)
    5: "#008000",   # тёмно-зелёный             (5/7)
    6: "#0000CD",   # синий                     (6/7)
    7: "#8B008B",   # фиолетовый   — максимум  (7/7)
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

        color = CONSENSUS_COLORS.get(min(support, 7), "")
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
        "(support≥%d): 2алг=%d, 3алг=%d, ≥4алг=%d",
        chrom, resolution, len(clusters), len(df_consensus), min_support,
        (df_consensus["support"] == 2).sum(),
        (df_consensus["support"] == 3).sum(),
        (df_consensus["support"] >= 4).sum(),
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
# P2: Консенсус по перекрытию целых TAD-доменов (Jaccard-based)
# ──────────────────────────────────────────────────────────────────────────────

def _jaccard(s_i: int, e_i: int, s_j: int, e_j: int) -> float:
    """Jaccard-индекс двух интервалов."""
    inter = max(0, min(e_i, e_j) - max(s_i, s_j))
    if inter == 0:
        return 0.0
    union = max(e_i, e_j) - min(s_i, s_j)
    return inter / union if union > 0 else 0.0


def compute_tad_consensus(
    algorithm_results: Dict[str, pd.DataFrame],
    chrom: str,
    resolution: int,
    jaccard_threshold: float = 0.5,
    min_support: int = 2,
) -> pd.DataFrame:
    """
    Консенсус по перекрытию целых TAD-доменов.

    Два TAD считаются «согласованными» если их Jaccard ≥ jaccard_threshold.
    Кластеризация Union-Find по всем парам TAD из разных алгоритмов.
    Координаты кластера — медиана start/end входящих TAD.
    Support — число уникальных алгоритмов в кластере.

    Parameters
    ----------
    algorithm_results  : {algo_name: domains_DataFrame}
    chrom              : хромосома (с префиксом chr)
    resolution         : разрешение в bp
    jaccard_threshold  : порог Jaccard для объединения (default 0.5)
    min_support        : минимум алгоритмов для консенсуса (default 2)

    Returns
    -------
    pd.DataFrame(chrom, start, end, support, algorithms)
      support    : число уникальных алгоритмов в кластере
      algorithms : строка с именами алгоритмов через запятую
    """
    # ── Собрать все TAD с меткой алгоритма ───────────────────────────────────
    all_tads: List[Tuple[str, int, int]] = []  # (algo, start, end)
    for algo, df in algorithm_results.items():
        if df is None or df.empty:
            continue
        df_c = df[df["chrom"] == chrom] if "chrom" in df.columns else df
        if df_c.empty:
            continue
        for _, row in df_c.iterrows():
            all_tads.append((algo, int(row["start"]), int(row["end"])))

    n = len(all_tads)
    if n < 2:
        logger.debug(
            "[TAD-Consensus] Недостаточно TAD для %s @ %d (%d шт.)", chrom, resolution, n
        )
        return pd.DataFrame(columns=["chrom", "start", "end", "support", "algorithms"])

    # ── Union-Find ───────────────────────────────────────────────────────────
    parent = list(range(n))

    def _find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def _union(x: int, y: int) -> None:
        parent[_find(x)] = _find(y)

    # Попарное сравнение TAD из разных алгоритмов
    for i in range(n):
        algo_i, s_i, e_i = all_tads[i]
        for j in range(i + 1, n):
            algo_j, s_j, e_j = all_tads[j]
            if algo_i == algo_j:
                continue  # не объединяем TAD одного алгоритма
            if _jaccard(s_i, e_i, s_j, e_j) >= jaccard_threshold:
                _union(i, j)

    # ── Группировка кластеров ────────────────────────────────────────────────
    from collections import defaultdict as _dd
    clusters: Dict[int, List[int]] = _dd(list)
    for i in range(n):
        clusters[_find(i)].append(i)

    records = []
    for root, members in clusters.items():
        algos = {all_tads[m][0] for m in members}
        support = len(algos)
        if support < min_support:
            continue

        starts = np.array([all_tads[m][1] for m in members], dtype=np.int64)
        ends   = np.array([all_tads[m][2] for m in members], dtype=np.int64)
        records.append({
            "chrom":      chrom,
            "start":      int(np.median(starts)),
            "end":        int(np.median(ends)),
            "support":    support,
            "algorithms": ",".join(sorted(algos)),
        })

    if not records:
        return pd.DataFrame(columns=["chrom", "start", "end", "support", "algorithms"])

    df_out = pd.DataFrame(records).sort_values("start").reset_index(drop=True)
    logger.info(
        "[TAD-Consensus] %s @ %d: %d TAD-кластеров (support≥%d, Jaccard≥%.2f)",
        chrom, resolution, len(df_out), min_support, jaccard_threshold,
    )
    return df_out


def save_tad_consensus_bed(
    df: pd.DataFrame,
    out_path: str,
) -> None:
    """Сохранить TAD-консенсус в BED-формате с колонкой support."""
    import os
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    out = df[["chrom", "start", "end"]].copy()
    out["name"]   = "tad_consensus_support" + df["support"].astype(str)
    out["score"]  = df["support"]
    out["strand"] = "."
    out.to_csv(out_path, sep="\t", header=False, index=False)
    logger.debug("TAD-Consensus BED: %s (%d доменов)", out_path, len(out))


# ──────────────────────────────────────────────────────────────────────────────
# Строгий TAD-консенсус (обе границы совпали)
# ──────────────────────────────────────────────────────────────────────────────

def compute_strong_tad_consensus(
    algorithm_results: Dict[str, pd.DataFrame],
    chrom: str,
    resolution: int,
    tolerance_bins: int = 1,
    min_support: int = 2,
) -> pd.DataFrame:
    """
    Строгий TAD-консенсус: TAD включается только если ОБЕ границы
    (start И end) совпадают у ≥ min_support алгоритмов с точностью
    до tolerance_bins бинов.

    Критерий совпадения пары (i, j) из РАЗНЫХ алгоритмов:
        |start_i - start_j| <= tolerance_bins * resolution  AND
        |end_i   - end_j  | <= tolerance_bins * resolution

    Алгоритм:
        1. Собрать все TAD как записи (algo, start, end).
        2. Для каждой пары из разных алгоритмов проверить критерий → Union-Find.
        3. Для каждого кластера:
           - representative start = median start всех элементов кластера,
             округлённый до бина
           - representative end   = median end, округлённый до бина
           - support = число РАЗЛИЧНЫХ алгоритмов в кластере
        4. Вернуть кластеры с support >= min_support, отсортированные по start.

    Parameters
    ----------
    algorithm_results : {algo_name: DataFrame(chrom, start, end)}
    chrom             : "chr1" и т.д. (с префиксом)
    resolution        : разрешение в bp (25000 / 50000 / 100000)
    tolerance_bins    : допуск в бинах (default=1 из config)
    min_support       : минимальное число алгоритмов в кластере

    Returns
    -------
    pd.DataFrame с колонками: chrom, start, end, support, algorithms
    """
    from collections import defaultdict as _dd

    tol_bp = tolerance_bins * resolution

    # ── 1. Собрать все TAD данного chrom ──────────────────────────────────
    all_tads: List[Tuple[str, int, int]] = []  # (algo, start, end)
    for algo, df in algorithm_results.items():
        if df is None or df.empty:
            continue
        df_c = df[df["chrom"] == chrom] if "chrom" in df.columns else df
        if df_c.empty:
            continue
        for _, row in df_c.iterrows():
            all_tads.append((algo, int(row["start"]), int(row["end"])))

    n = len(all_tads)
    if n < 2:
        logger.debug(
            "[StrongConsensus] Недостаточно TAD для %s @ %d (%d шт.)",
            chrom, resolution, n,
        )
        return pd.DataFrame(columns=["chrom", "start", "end", "support", "algorithms"])

    # ── 2. Union-Find ──────────────────────────────────────────────────────
    parent = list(range(n))

    def _find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def _union(x: int, y: int) -> None:
        parent[_find(x)] = _find(y)

    for i in range(n):
        algo_i, s_i, e_i = all_tads[i]
        for j in range(i + 1, n):
            algo_j, s_j, e_j = all_tads[j]
            if algo_i == algo_j:
                continue  # пары внутри одного алгоритма не объединяем
            if abs(s_i - s_j) <= tol_bp and abs(e_i - e_j) <= tol_bp:
                _union(i, j)

    # ── 3. Агрегировать кластеры ──────────────────────────────────────────
    clusters: Dict[int, List[int]] = _dd(list)
    for i in range(n):
        clusters[_find(i)].append(i)

    records = []
    for root, members in clusters.items():
        algos = {all_tads[m][0] for m in members}
        support = len(algos)
        if support < min_support:
            continue

        starts = np.array([all_tads[m][1] for m in members], dtype=np.int64)
        ends   = np.array([all_tads[m][2] for m in members], dtype=np.int64)

        # Округлить медиану до бина
        rep_start = (int(np.median(starts)) // resolution) * resolution
        rep_end   = (int(np.median(ends))   // resolution) * resolution
        if rep_end <= rep_start:
            rep_end = rep_start + resolution

        records.append({
            "chrom":      chrom,
            "start":      rep_start,
            "end":        rep_end,
            "support":    support,
            "algorithms": ",".join(sorted(algos)),
        })

    if not records:
        logger.info(
            "[StrongConsensus] Нет кластеров с support≥%d для %s @ %d bp",
            min_support, chrom, resolution,
        )
        return pd.DataFrame(columns=["chrom", "start", "end", "support", "algorithms"])

    df_out = pd.DataFrame(records).sort_values("start").reset_index(drop=True)
    logger.info(
        "[StrongConsensus] %s @ %d bp → %d строгих TAD (tol=%d бин, support≥%d)",
        chrom, resolution, len(df_out), tolerance_bins, min_support,
    )
    return df_out

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

            # ── TAD-консенсус (Jaccard Union-Find) ──────────────────────────
            jaccard_threshold = cfg["consensus"].get("jaccard_threshold", 0.5)
            tad_df = compute_tad_consensus(
                algo_dfs, chrom, res,
                jaccard_threshold=jaccard_threshold,
                min_support=min_support,
            )
            if not tad_df.empty:
                tad_out_dir = out_dir if out_dir is not None else cfg.get(
                    "paths", {}).get("consensus_out", "results/consensus")
                tad_path = os.path.join(
                    tad_out_dir, f"tad_consensus_{chrom}_{res}bp.bed"
                )
                save_tad_consensus_bed(tad_df, tad_path)
                logger.info(
                    "[TAD-Consensus] сохранён: %s (%d доменов)",
                    tad_path, len(tad_df),
                )

            # ── Строгий TAD-консенсус (обе границы совпали) ─────────────────
            strong_tol = cfg["consensus"].get("strong_boundary_tolerance_bins", 1)
            strong_df = compute_strong_tad_consensus(
                algorithm_results=algo_dfs,
                chrom=chrom,
                resolution=res,
                tolerance_bins=strong_tol,
                min_support=min_support,
            )
            if not strong_df.empty:
                strong_out_dir = out_dir if out_dir is not None else cfg.get(
                    "paths", {}).get("consensus_out", "results/consensus")
                strong_tad_path = os.path.join(
                    strong_out_dir, f"strong_tad_consensus_{chrom}_{res}bp.bed"
                )
                save_tad_consensus_bed(strong_df, strong_tad_path)
                logger.info(
                    "[StrongConsensus] сохранён: %s (%d доменов)",
                    strong_tad_path, len(strong_df),
                )

    return dict(results)
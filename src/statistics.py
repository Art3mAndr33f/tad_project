"""
statistics.py
=============
Все метрики качества TAD-детекции:

1. Базовая статистика доменов (число, размер, покрытие)
2. Попарное сравнение алгоритмов (Jaccard, boundary overlap)
3. Сравнение с эталоном Rao 2014 (precision, recall, F1)
4. Интеграция с консенсусными границами
"""

from __future__ import annotations

import logging
from itertools import combinations
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy import stats as scipy_stats

logger = logging.getLogger(__name__)

_arrowhead_cache: dict[str, Optional[pd.DataFrame]] = {}
_arrowhead_warned: set[str] = set()

def _load_arrowhead_cached(path: str) -> Optional[pd.DataFrame]:
    """Загрузить Arrowhead с кэшированием и однократным предупреждением."""
    if path in _arrowhead_cache:
        return _arrowhead_cache[path]
    try:
        df = pd.read_csv(path, sep="\t", comment="#",
                         names=["chrom1","x1","x2","chrom2","y1","y2"])
        _arrowhead_cache[path] = df
        return df
    except FileNotFoundError:
        if path not in _arrowhead_warned:
            logger.warning("Не удалось загрузить Arrowhead: %s", path)
            _arrowhead_warned.add(path)
        _arrowhead_cache[path] = None
        return None

# ──────────────────────────────────────────────────────────────────────────────
# Размеры хромосом hg19
# ──────────────────────────────────────────────────────────────────────────────

HG19_CHROM_SIZES = {
    "chr1":  249250621, "chr2":  243199373, "chr3":  198022430,
    "chr4":  191154276, "chr5":  180915260, "chr6":  171115067,
    "chr7":  159138663, "chr8":  146364022, "chr9":  141213431,
    "chr10": 135534747, "chr11": 135006516, "chr12": 133851895,
    "chr13": 115169878, "chr14": 107349540, "chr15": 102531392,
    "chr16":  90354753, "chr17":  81195210, "chr18":  78077248,
    "chr19":  59128983, "chr20":  63025520, "chr21":  48129895,
    "chr22":  51304566, "chrX":  155270560,
}


# ──────────────────────────────────────────────────────────────────────────────
# 1. Базовая статистика
# ──────────────────────────────────────────────────────────────────────────────

def compute_basic_stats(
    df: pd.DataFrame,
    chrom: str,
    resolution: int,
    consensus_df: Optional[pd.DataFrame] = None,
    min_consensus_support: int = 2,
) -> dict:
    """
    Вычислить базовую статистику для одного алгоритма/хромосомы/разрешения.

    Returns
    -------
    dict с ключами:
      n_tads, median_size_kb, mean_size_kb, std_size_kb,
      coverage_pct, n_consensus_bnd, pct_consensus_bnd
    """
    chrom_size = HG19_CHROM_SIZES.get(chrom, 0)

    if df.empty:
        return {
            "n_tads":             0,
            "median_size_kb":     np.nan,
            "mean_size_kb":       np.nan,
            "std_size_kb":        np.nan,
            "coverage_pct":       0.0,
            "n_consensus_bnd":    0,
            "pct_consensus_bnd":  0.0,
        }

    sizes_kb = (df["end"] - df["start"]).values / 1000.0
    coverage = (df["end"] - df["start"]).sum()

    # Консенсусные границы
    n_consensus_bnd    = 0
    pct_consensus_bnd  = 0.0
    if consensus_df is not None and not consensus_df.empty:
        from src.consensus import extract_boundaries
        algo_bnd = extract_boundaries(df, resolution)
        consensus_pos = set(
            consensus_df[consensus_df["support"] >= min_consensus_support]["position"]
            .astype(int).tolist()
        )
        tol = resolution
        n_consensus_bnd = sum(
            any(abs(b - cp) <= tol for cp in consensus_pos)
            for b in algo_bnd
        )
        pct_consensus_bnd = (
            100.0 * n_consensus_bnd / len(algo_bnd) if len(algo_bnd) > 0 else 0.0
        )

    return {
        "n_tads":             int(len(df)),
        "median_size_kb":     float(np.median(sizes_kb)),
        "mean_size_kb":       float(np.mean(sizes_kb)),
        "std_size_kb":        float(np.std(sizes_kb)),
        "coverage_pct":       float(100.0 * coverage / chrom_size) if chrom_size else 0.0,
        "n_consensus_bnd":    int(n_consensus_bnd),
        "pct_consensus_bnd":  float(pct_consensus_bnd),
    }


# ──────────────────────────────────────────────────────────────────────────────
# 2. Попарное сравнение
# ──────────────────────────────────────────────────────────────────────────────

def jaccard_domains(df_a: pd.DataFrame, df_b: pd.DataFrame) -> float:
    """
    Jaccard index на уровне доменов.
    Домен считается совпадающим, если start и end идентичны.
    """
    if df_a.empty and df_b.empty:
        return 1.0
    if df_a.empty or df_b.empty:
        return 0.0

    set_a = set(zip(df_a["start"].astype(int), df_a["end"].astype(int)))
    set_b = set(zip(df_b["start"].astype(int), df_b["end"].astype(int)))
    inter = len(set_a & set_b)
    union = len(set_a | set_b)
    return inter / union if union > 0 else 0.0


def boundary_overlap_rate(
    df_a: pd.DataFrame,
    df_b: pd.DataFrame,
    resolution: int,
    tolerance_bins: int = 1,
) -> float:
    """
    % границ алгоритма A, найденных алгоритмом B в пределах ±tolerance_bins.
    """
    from src.consensus import extract_boundaries

    bnd_a = extract_boundaries(df_a, resolution)
    bnd_b = extract_boundaries(df_b, resolution)

    if len(bnd_a) == 0:
        return np.nan
    if len(bnd_b) == 0:
        return 0.0

    tol = tolerance_bins * resolution
    matched = sum(
        any(abs(ba - bb) <= tol for bb in bnd_b)
        for ba in bnd_a
    )
    return 100.0 * matched / len(bnd_a)


def compute_pairwise_matrix(
    algo_results: Dict[str, pd.DataFrame],
    resolution: int,
    tolerance_bins: int = 1,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Вычислить матрицы попарного сравнения (Jaccard и boundary overlap rate).

    Returns
    -------
    (jaccard_df, overlap_df) — оба DataFrame с алгоритмами как строками и столбцами
    """
    algos = list(algo_results.keys())
    n = len(algos)

    jaccard_mat  = np.ones((n, n))
    overlap_mat  = np.full((n, n), np.nan)

    for i, algo_a in enumerate(algos):
        for j, algo_b in enumerate(algos):
            if i == j:
                overlap_mat[i, j] = 100.0
                continue
            jaccard_mat[i, j] = jaccard_domains(
                algo_results[algo_a], algo_results[algo_b]
            )
            overlap_mat[i, j] = boundary_overlap_rate(
                algo_results[algo_a], algo_results[algo_b],
                resolution, tolerance_bins,
            )

    jaccard_df = pd.DataFrame(jaccard_mat, index=algos, columns=algos)
    overlap_df = pd.DataFrame(overlap_mat, index=algos, columns=algos)
    return jaccard_df, overlap_df


# ──────────────────────────────────────────────────────────────────────────────
# 3. Сравнение с эталоном Rao 2014
# ──────────────────────────────────────────────────────────────────────────────

def load_arrowhead(arrowhead_path: str, chrom: Optional[str] = None) -> pd.DataFrame:
    """
    Загрузить Arrowhead domain list из файла.
    Формат: chr1  start  end  ... (TSV, первые 3 колонки)
    """
    df = pd.read_csv(arrowhead_path, sep="\t", header=0)

    # Нормализуем названия колонок
    cols = df.columns.tolist()
    rename_map = {}
    for i, c in enumerate(cols[:3]):
        rename_map[c] = ["chrom", "start", "end"][i]
    df = df.rename(columns=rename_map)

    # Добавить 'chr' если отсутствует
    if not df["chrom"].iloc[0].startswith("chr"):
        df["chrom"] = "chr" + df["chrom"].astype(str)

    if chrom is not None:
        df = df[df["chrom"] == chrom]

    return df[["chrom", "start", "end"]].reset_index(drop=True)


def compare_with_reference(
    pred_df: pd.DataFrame,
    ref_df: pd.DataFrame,
    resolution: int,
    tolerance_bins_list: List[int] = [1, 2],
) -> dict:
    """
    Сравнить результат алгоритма с эталоном Arrowhead.

    Metrics: boundary_recall, boundary_precision, f1, jaccard_domains

    Parameters
    ----------
    pred_df             : домены алгоритма
    ref_df              : домены Arrowhead (фильтрованные по хромосоме)
    resolution          : разрешение в bp
    tolerance_bins_list : список допусков для вычисления

    Returns
    -------
    dict {f"tol{t}": {recall, precision, f1, jaccard}}
    """
    from src.consensus import extract_boundaries

    results = {}
    ref_bnd  = extract_boundaries(ref_df, resolution)
    pred_bnd = extract_boundaries(pred_df, resolution)

    for tol in tolerance_bins_list:
        tol_bp = tol * resolution
        key = f"tol{tol}bin"

        if len(ref_bnd) == 0:
            results[key] = dict(recall=np.nan, precision=np.nan,
                                f1=np.nan, jaccard=np.nan)
            continue

        # Recall: доля reference-границ, найденных алгоритмом
        if len(pred_bnd) > 0:
            recall_count = sum(
                any(abs(rb - pb) <= tol_bp for pb in pred_bnd)
                for rb in ref_bnd
            )
            recall = recall_count / len(ref_bnd)
        else:
            recall = 0.0

        # Precision: доля predicted-границ, совпадающих с reference
        if len(pred_bnd) > 0:
            prec_count = sum(
                any(abs(pb - rb) <= tol_bp for rb in ref_bnd)
                for pb in pred_bnd
            )
            precision = prec_count / len(pred_bnd)
        else:
            precision = 0.0

        # F1
        f1 = (2 * precision * recall / (precision + recall + 1e-9))

        # Jaccard (домены)
        jacc = jaccard_domains(pred_df, ref_df)

        results[key] = {
            "recall":    float(recall),
            "precision": float(precision),
            "f1":        float(f1),
            "jaccard":   float(jacc),
        }

    return results


# ──────────────────────────────────────────────────────────────────────────────
# 4. Сводная статистика
# ──────────────────────────────────────────────────────────────────────────────

def aggregate_statistics(
    all_results: Dict,
    cfg: dict,
    consensus_all: Optional[Dict] = None,
) -> Dict[str, pd.DataFrame]:
    """
    Собрать полную сводную статистику по всем алгоритмам/хромосомам/разрешениям.

    Returns
    -------
    dict с ключами:
      'basic'      : DataFrame базовой статистики
      'pairwise_jaccard'  : DataFrame попарного Jaccard
      'pairwise_overlap'  : DataFrame попарного overlap
      'rao_comparison'    : DataFrame сравнения с Rao
    """
    arrowhead_path = cfg["paths"]["arrowhead_file"]
    resolutions    = cfg["resolutions"]
    chromosomes    = cfg["chromosomes"]["all"]
    algos          = list(all_results.keys())

    basic_rows   = []
    rao_rows     = []
    pairwise_rows_j = []
    pairwise_rows_o = []

    for res in resolutions:
        for chrom in chromosomes:
            # Базовая статистика
            for algo in algos:
                df = all_results.get(algo, {}).get(chrom, {}).get(res)
                if df is None:
                    continue
                cons_df = None
                if consensus_all and chrom in consensus_all and res in consensus_all[chrom]:
                    cons_df = consensus_all[chrom][res]

                stats = compute_basic_stats(df, chrom, res, cons_df)
                basic_rows.append({
                    "algorithm":   algo,
                    "chrom":       chrom,
                    "resolution":  res,
                    **stats,
                })

            # Попарное сравнение
            algo_dfs: Dict[str, pd.DataFrame] = {}
            for algo in algos:
                df = all_results.get(algo, {}).get(chrom, {}).get(res)
                if df is not None and not df.empty:
                    algo_dfs[algo] = df

            if len(algo_dfs) >= 2:
                jac_df, ovl_df = compute_pairwise_matrix(algo_dfs, res)
                for a1, a2 in combinations(list(algo_dfs.keys()), 2):
                    pairwise_rows_j.append({
                        "algo_a": a1, "algo_b": a2,
                        "chrom": chrom, "resolution": res,
                        "jaccard": jac_df.loc[a1, a2],
                    })
                    pairwise_rows_o.append({
                        "algo_a": a1, "algo_b": a2,
                        "chrom": chrom, "resolution": res,
                        "overlap_a_in_b": ovl_df.loc[a1, a2],
                        "overlap_b_in_a": ovl_df.loc[a2, a1],
                    })

            # Сравнение с Rao
            try:
                ref_df = load_arrowhead(arrowhead_path, chrom=chrom)
            except Exception as exc:
                logger.warning("Не удалось загрузить Arrowhead: %s", exc)
                ref_df = pd.DataFrame(columns=["chrom", "start", "end"])

            for algo in algos:
                df = all_results.get(algo, {}).get(chrom, {}).get(res)
                if df is None or df.empty:
                    continue
                rao_metrics = compare_with_reference(
                    df, ref_df, res, [1, 2]
                )
                row = {"algorithm": algo, "chrom": chrom, "resolution": res}
                for tol_key, metrics in rao_metrics.items():
                    for metric, val in metrics.items():
                        row[f"{tol_key}_{metric}"] = val
                rao_rows.append(row)

    return {
        "basic":             pd.DataFrame(basic_rows),
        "pairwise_jaccard":  pd.DataFrame(pairwise_rows_j),
        "pairwise_overlap":  pd.DataFrame(pairwise_rows_o),
        "rao_comparison":    pd.DataFrame(rao_rows),
    }


def save_statistics(stats_dict: Dict[str, pd.DataFrame], out_dir: str) -> None:
    """Сохранить статистику в CSV-файлы."""
    import os
    os.makedirs(out_dir, exist_ok=True)
    for name, df in stats_dict.items():
        path = os.path.join(out_dir, f"{name}.csv")
        df.to_csv(path, index=False)
        logger.info("Статистика сохранена: %s (%d строк)", path, len(df))
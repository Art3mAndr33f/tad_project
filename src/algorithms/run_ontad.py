# src/algorithms/run_ontad.py
"""OnTAD: hierarchical TAD detection (An et al. 2019).

Стратегия v2.3: recursive top-down subdivision.
Начинаем с top-level (depth=1) TAD, рекурсивно заменяем крупные
на детей пока не достигнем target или не исчерпаем иерархию.
"""

from __future__ import annotations

import logging
import os
import subprocess
import tempfile
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


# ──────────────────────────────────────────────────────────────────────────────
# Парсинг .tad файла
# ──────────────────────────────────────────────────────────────────────────────

def _parse_tad_file(tad_path: str) -> pd.DataFrame:
    rows = []
    with open(tad_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 5:
                continue
            try:
                rows.append({
                    "start_bin": int(parts[0]),
                    "end_bin":   int(parts[1]),
                    "depth":     int(parts[2]),
                    "score1":    float(parts[3]),
                    "score2":    float(parts[4]),
                })
            except (ValueError, IndexError):
                pass
    if not rows:
        return pd.DataFrame(
            columns=["start_bin", "end_bin", "depth", "score1", "score2"]
        )
    return pd.DataFrame(rows)


# ──────────────────────────────────────────────────────────────────────────────
# Recursive top-down subdivision
# ──────────────────────────────────────────────────────────────────────────────

def _get_children(
    df: pd.DataFrame,
    parent_s: int,
    parent_e: int,
    parent_depth: int,
    min_tad_bins: int,
) -> pd.DataFrame:
    """Прямые дети: depth == parent_depth+1, координаты строго внутри родителя."""
    return df[
        (df["depth"] == parent_depth + 1) &
        (df["start_bin"] >= parent_s) &
        (df["end_bin"]   <= parent_e) &
        ((df["end_bin"] - df["start_bin"]) >= min_tad_bins)
    ].copy()


def _subdivide(
    current: list[dict],
    df: pd.DataFrame,
    min_tad_bins: int,
    target_lo: int,
    target_hi: int,
    split_bins: int,
) -> list[dict]:
    """Один проход subdivision: заменить TAD крупнее split_bins детьми.

    Возвращает новый список TAD после замены.
    split_bins снижается на каждом уровне рекурсии до тех пор,
    пока не достигнем target или не исчерпаем иерархию.
    """
    if target_lo <= len(current) <= target_hi:
        return current

    new_current: list[dict] = []
    changed = False

    for row in current:
        sz = row["end_bin"] - row["start_bin"]

        if sz <= split_bins:
            new_current.append(row)
            continue

        children = _get_children(
            df,
            row["start_bin"],
            row["end_bin"],
            row["depth"],
            min_tad_bins,
        )

        if len(children) < 2:
            new_current.append(row)
            continue

        new_current.extend(children.to_dict("records"))
        changed = True

    if not changed:
        return current

    logger.debug("_subdivide: split_bins=%d → %d TADs", split_bins, len(new_current))
    return new_current


def _recursive_subdivide(
    df_all: pd.DataFrame,
    min_tad_bins: int,
    target_lo: int,
    target_hi: int,
    n_bins: int,
) -> pd.DataFrame:
    """Top-down recursive subdivision до попадания в [target_lo, target_hi].

    Алгоритм:
    1. Старт: depth=1 TAD >= min_tad_bins.
    2. Рассчитать split_threshold: начинаем с n_bins // target_lo
       (средний ожидаемый размер), снижаем до min_tad_bins * 3.
    3. На каждом уровне: заменять TAD крупнее порога их детьми.
    4. Остановиться при достижении target или исчерпании иерархии.
    """
    if df_all.empty:
        return df_all.iloc[0:0].copy()

    df = df_all[
        (df_all["depth"] > 0) &
        ((df_all["end_bin"] - df_all["start_bin"]) >= min_tad_bins)
    ].copy()

    if df.empty:
        return df_all.iloc[0:0].copy()

    # Стартовый набор: depth=1
    min_d = df["depth"].min()
    current_rows = df[df["depth"] == min_d].to_dict("records")

    if not current_rows:
        return df_all.iloc[0:0].copy()

    logger.debug("recursive_subdivide start: %d TADs at depth=%d",
                 len(current_rows), min_d)

    # Sweep split_threshold от крупного к мелкому
    target_mid = (target_lo + target_hi) / 2.0
    max_depth = int(df["depth"].max())

    # Пороги для subdivision (от крупного к мелкому)
    thresholds = []
    t = n_bins // max(1, target_lo)
    while t >= min_tad_bins * 2:
        thresholds.append(t)
        t = int(t * 0.7)
    if not thresholds:
        thresholds = [min_tad_bins * 4, min_tad_bins * 3, min_tad_bins * 2]

    best_rows = list(current_rows)
    best_dist = abs(len(current_rows) - target_mid)

    for thr in thresholds:
        new_rows = _subdivide(
            list(current_rows), df, min_tad_bins, target_lo, target_hi, thr
        )
        n_new = len(new_rows)
        dist  = abs(n_new - target_mid)

        logger.debug("threshold=%d → %d TADs", thr, n_new)

        if dist < best_dist:
            best_dist = dist
            best_rows = new_rows

        if target_lo <= n_new <= target_hi:
            logger.info("recursive_subdivide: hit target at thr=%d → %d TADs",
                        thr, n_new)
            best_rows = new_rows
            break

        current_rows = new_rows   # продолжить дробить

    result = (
        pd.DataFrame(best_rows)
        .sort_values("start_bin")
        .reset_index(drop=True)
    )
    logger.info("recursive_subdivide FINAL: %d TADs (target %d–%d)",
                len(result), target_lo, target_hi)
    return result


# ──────────────────────────────────────────────────────────────────────────────
# Subprocess
# ──────────────────────────────────────────────────────────────────────────────

def _run_ontad_subprocess(
    matrix_path: str,
    out_prefix: str,
    n_bins: int,
    penalty: float,
    minsz: int,
    maxsz: int,
    use_log2: bool,
    ontad_bin: str,
) -> Optional[str]:
    cmd = [
        ontad_bin, matrix_path,
        "-penalty", str(penalty),
        "-minsz",   str(minsz),
        "-maxsz",   str(max(minsz + 1, min(maxsz, n_bins - 1))),
        "-o",       out_prefix,
    ]
    if use_log2:
        cmd.append("-log2")

    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=600,
        )
        if result.returncode != 0:
            logger.warning("OnTAD exit %d: %s",
                           result.returncode, result.stderr[:300])
            return None

        tad_file = out_prefix + ".tad"
        if not Path(tad_file).exists():
            candidates = list(
                Path(out_prefix).parent.glob(Path(out_prefix).name + "*.tad")
            )
            tad_file = str(candidates[0]) if candidates else None
        return tad_file

    except subprocess.TimeoutExpired:
        logger.error("OnTAD timeout")
        return None
    except Exception as exc:
        logger.error("OnTAD subprocess: %s", exc)
        return None


# ──────────────────────────────────────────────────────────────────────────────
# Python fallback
# ──────────────────────────────────────────────────────────────────────────────

def _python_fallback(
    matrix: np.ndarray,
    chrom: str,
    resolution: int,
    min_tad_bins: int,
) -> pd.DataFrame:
    logger.info("[ontad] IS-fallback for %s @ %d", chrom, resolution)
    from scipy.ndimage import gaussian_filter1d  # type: ignore
    from scipy.signal import argrelmin           # type: ignore

    n = len(matrix)
    window = min(10, n // 4)
    ins = np.array([
        matrix[max(0, i - window):i, i:min(n, i + window)].mean()
        if i >= window and i < n - window else 0.0
        for i in range(n)
    ])
    ins_sm   = gaussian_filter1d(ins, sigma=2)
    ins_norm = ins_sm - np.nanmean(ins_sm)
    mins     = argrelmin(ins_norm, order=max(1, window // 2))[0]

    boundaries = sorted({0} | set(mins.tolist()) | {n})
    tads = [
        {"chrom": chrom,
         "start": boundaries[i] * resolution,
         "end":   boundaries[i + 1] * resolution}
        for i in range(len(boundaries) - 1)
        if boundaries[i + 1] - boundaries[i] >= min_tad_bins
    ]
    df = pd.DataFrame(tads, columns=["chrom", "start", "end"])
    logger.info("[ontad] Fallback: %d TADs", len(df))
    return df


# ──────────────────────────────────────────────────────────────────────────────
# Главная функция
# ──────────────────────────────────────────────────────────────────────────────

def run_ontad(
    chrom: str,
    resolution: int,
    data_path: str,
    cfg: Optional[dict] = None,
    **kwargs,
) -> pd.DataFrame:
    _empty = pd.DataFrame(columns=["chrom", "start", "end"])

    algo_cfg  = (cfg or {}).get("algorithms", {}).get("ontad", {})
    paths_cfg = (cfg or {}).get("paths", {})

    penalty_cfg = algo_cfg.get("penalty", None)
    minsz       = int(algo_cfg.get("minsz",   3))
    maxsz       = int(algo_cfg.get("maxsz", 200))
    use_log2    = bool(algo_cfg.get("log2", True))
    min_tad_kb  = float(algo_cfg.get("min_tad_kb", 100.0))
    ontad_bin   = paths_cfg.get("ontad_bin", "tools/ontad/OnTAD")

    min_tad_bins = max(minsz, int(min_tad_kb * 1_000 / resolution))

    try:
        from src.data_prep import get_matrix  # type: ignore
        matrix = get_matrix(cfg, chrom, resolution)
    except Exception as exc:
        logger.error("[ontad] get_matrix: %s", exc, exc_info=True)
        return _empty

    if matrix is None or matrix.size == 0:
        return _empty

    matrix   = np.asarray(matrix, dtype=np.float64)
    n_bins   = matrix.shape[0]
    chrom_mb = n_bins * resolution / 1_000_000

    target_lo  = max(5,   int(chrom_mb * 0.8))
    target_hi  = min(200, int(chrom_mb * 2.5))
    target_mid = (target_lo + target_hi) / 2.0

    logger.info("[ontad] %s @ %d: %d bins (%.1f Mb), target %d–%d TADs",
                chrom, resolution, n_bins, chrom_mb, target_lo, target_hi)

    ontad_ok = Path(ontad_bin).exists() and os.access(ontad_bin, os.X_OK)
    if not ontad_ok:
        logger.warning("[ontad] Binary not found → fallback")
        return _python_fallback(matrix, chrom, resolution, min_tad_bins)

    # Penalty sweep: начинаем с мелкого penalty (больше иерархии)
    # чтобы recursive_subdivide имел материал для работы
    penalties = [float(penalty_cfg)] if penalty_cfg is not None else \
                [0.1, 0.05, 0.02, 0.01, 0.005]

    with tempfile.TemporaryDirectory() as tmpdir:
        mat_path = os.path.join(tmpdir, f"{chrom}.txt")
        np.savetxt(mat_path, matrix, delimiter="\t", fmt="%.6g")

        best_df: Optional[pd.DataFrame] = None
        best_n   = 0

        for penalty in penalties:
            out_pref = os.path.join(tmpdir, f"{chrom}_p{penalty:.4f}")
            tad_file = _run_ontad_subprocess(
                mat_path, out_pref, n_bins,
                penalty, minsz, maxsz, use_log2, ontad_bin,
            )
            if tad_file is None:
                continue

            df_raw = _parse_tad_file(tad_file)
            if df_raw.empty:
                continue

            depth_counts = df_raw[df_raw.depth > 0]["depth"].value_counts()
            logger.info("[ontad] pen=%.4f: %d rows, depths=%s",
                        penalty, len(df_raw),
                        dict(depth_counts.sort_index()))

            df_sel = _recursive_subdivide(
                df_raw, min_tad_bins, target_lo, target_hi, n_bins
            )
            n_sel = len(df_sel)
            logger.info("[ontad] pen=%.4f → %d TADs after subdivision",
                        penalty, n_sel)

            if target_lo <= n_sel <= target_hi:
                best_df = df_sel
                logger.info("[ontad] ✓ target hit at pen=%.4f: %d TADs",
                            penalty, n_sel)
                break

            if best_df is None or \
               abs(n_sel - target_mid) < abs(best_n - target_mid):
                best_df = df_sel
                best_n  = n_sel

    if best_df is None or best_df.empty:
        logger.warning("[ontad] All penalties failed → fallback")
        return _python_fallback(matrix, chrom, resolution, min_tad_bins)

    # Конвертация бин → bp
    tads = []
    for _, row in best_df.iterrows():
        s_bp = int((row["start_bin"] - 1) * resolution)
        e_bp = min(int(row["end_bin"] * resolution), n_bins * resolution)
        if e_bp - s_bp >= min_tad_bins * resolution:
            tads.append({"chrom": chrom, "start": s_bp, "end": e_bp})

    if not tads:
        return _python_fallback(matrix, chrom, resolution, min_tad_bins)

    df_out = pd.DataFrame(tads, columns=["chrom", "start", "end"])
    logger.info("[ontad] FINAL: %d TADs, median=%.0fkb",
                len(df_out), (df_out.end - df_out.start).median() / 1000)
    return df_out
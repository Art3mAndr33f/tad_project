"""
src/algorithms/run_ontad.py

TAD detection via OnTAD (hierarchical, sliding-average based).
Reference: An et al. Genome Biology 2019.
GitHub: https://github.com/anlin00007/OnTAD

OnTAD принимает N×N матрицу (пробел-разделённую) и возвращает
иерархические TAD. Мы используем только верхний уровень иерархии
(level_1) для консенсусного пайплайна.
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

_DEFAULT_PENALTY  = 0.1
_DEFAULT_MINSZ    = 3       # в бинах
_DEFAULT_MAXSZ    = 200     # в бинах (200 бинов × 25kb = 5 Mb)
_DEFAULT_LOG2     = True
_TIMEOUT_SEC      = 900     # 15 минут на хромосому


# ═══════════════════════════════════════════════════════════════════════════════
# Публичный интерфейс
# ═══════════════════════════════════════════════════════════════════════════════

def run_ontad(
    chrom: str,
    resolution: int,
    data_path: str,
    cfg: Optional[dict],
    **kwargs,
) -> pd.DataFrame:
    empty = pd.DataFrame(columns=["chrom", "start", "end"])

    ontad_bin = _get_ontad_bin(cfg)
    if ontad_bin is None or not Path(ontad_bin).exists():
        logger.error("[ontad] Бинарник OnTAD не найден: %s", ontad_bin)
        return empty

    try:
        matrix = _load_matrix(data_path, chrom, resolution)
    except FileNotFoundError as exc:
        logger.error("[ontad] Матрица не найдена %s @ %d: %s", chrom, resolution, exc)
        return empty

    n = matrix.shape[0]
    logger.info("[ontad] %s @ %d — матрица %d×%d загружена", chrom, resolution, n, n)

    algo_cfg = (cfg or {}).get("algorithms", {}).get("ontad", {})
    minsz    = int(  algo_cfg.get("minsz",  _DEFAULT_MINSZ))
    maxsz    = int(  algo_cfg.get("maxsz",  _DEFAULT_MAXSZ))
    log2flag = bool( algo_cfg.get("log2",   _DEFAULT_LOG2))

    # Целевой диапазон TAD (1.0–2.5 TAD/Mb, как в Armatus)
    chrom_mb   = n * resolution / 1e6
    target_min = max(5,  int(chrom_mb * 1.0))
    target_max = max(10, int(chrom_mb * 2.5))

    logger.debug(
        "[ontad] %s @ %d: %.1f Mb, target=%d–%d TADs",
        chrom, resolution, chrom_mb, target_min, target_max,
    )

    # Penalty sweep: от строгого (мало TAD) к мягкому (много TAD)
    penalty_cfg = algo_cfg.get("penalty", None)
    if penalty_cfg is not None:
        penalties = [float(penalty_cfg)]
    else:
        penalties = [1.0, 0.5, 0.3, 0.2, 0.1, 0.05, 0.02, 0.01]

    best_df  = empty
    best_n   = 0
    best_pen = penalties[0]

    with tempfile.TemporaryDirectory(prefix="ontad_") as tmpdir:
        matrix_path = Path(tmpdir) / f"{chrom}_{resolution}bp.matrix"
        _write_matrix(matrix, matrix_path)

        for pen in penalties:
            out_prefix = Path(tmpdir) / f"output_pen{pen:.3f}"
            cmd = [
                str(ontad_bin), str(matrix_path),
                "-penalty", str(pen),
                "-minsz",   str(minsz),
                "-maxsz",   str(maxsz),
                "-o",       str(out_prefix),
            ]
            if log2flag:
                cmd.append("-log2")

            try:
                proc = subprocess.run(
                    cmd, capture_output=True, text=True, timeout=_TIMEOUT_SEC,
                )
            except subprocess.TimeoutExpired:
                logger.warning("[ontad] Таймаут penalty=%.3f на %s @ %d", pen, chrom, resolution)
                continue
            except Exception as exc:
                logger.warning("[ontad] Subprocess ошибка penalty=%.3f: %s", pen, exc)
                continue

            # Найти output файл
            tad_file = Path(str(out_prefix) + ".tad")
            if not tad_file.exists() or tad_file.stat().st_size == 0:
                logger.debug("[ontad] penalty=%.3f → output не найден", pen)
                continue

            df_pen = _parse_ontad_output(tad_file, chrom, resolution)
            n_tads = len(df_pen)
            logger.debug("[ontad] penalty=%.3f → %d TADs (depth=1)", pen, n_tads)

            if target_min <= n_tads <= target_max:
                best_df  = df_pen
                best_pen = pen
                logger.info(
                    "[ontad] %s @ %d — penalty=%.3f → %d TADs ✓ (target %d–%d)",
                    chrom, resolution, pen, n_tads, target_min, target_max,
                )
                break

            # Запомнить лучший пока не попали в диапазон
            if n_tads > best_n:
                best_n   = n_tads
                best_pen = pen
                best_df  = df_pen

    if best_df.empty or len(best_df) < 3:
        logger.warning(
            "[ontad] %s @ %d — не попали в target (%d–%d), "
            "лучший результат: %d TADs @ penalty=%.3f",
            chrom, resolution, target_min, target_max, len(best_df), best_pen,
        )

    logger.info("[ontad] %s @ %d — найдено %d TAD", chrom, resolution, len(best_df))
    return best_df


# ═══════════════════════════════════════════════════════════════════════════════
# Внутренние функции
# ═══════════════════════════════════════════════════════════════════════════════

def _get_ontad_bin(cfg: Optional[dict]) -> Optional[str]:
    """Найти путь к бинарнику OnTAD."""
    # 1. Из конфига
    bin_path = (cfg or {}).get("paths", {}).get("ontad_bin")
    if bin_path and Path(bin_path).exists():
        return str(bin_path)
    # 2. Дефолтный путь в проекте
    default = Path("tools/ontad/OnTAD")
    if default.exists():
        return str(default)
    # 3. В PATH системы
    import shutil
    found = shutil.which("OnTAD")
    if found:
        return found
    return None


def _load_matrix(data_path: str, chrom: str, resolution: int) -> np.ndarray:
    """Загрузить .npy матрицу."""
    npy_path = Path(data_path) / f"{chrom}_{resolution}bp.npy"
    if not npy_path.exists():
        raise FileNotFoundError(npy_path)
    matrix = np.load(str(npy_path)).astype(np.float64)
    matrix = np.nan_to_num(matrix, nan=0.0, posinf=0.0, neginf=0.0)
    matrix = np.maximum(matrix, matrix.T)
    return matrix


def _write_matrix(matrix: np.ndarray, path: Path) -> None:
    """
    Записать матрицу в формат OnTAD: N строк по N значений, разделитель TAB.
    OnTAD принимает полную NxN матрицу.
    """
    np.savetxt(str(path), matrix, delimiter="\t", fmt="%.6g")


def _parse_ontad_output(
    output_file: Path,
    chrom: str,
    resolution: int,
) -> pd.DataFrame:
    """
    Парсинг .tad формата OnTAD v1.4.

    Формат (TAB-разделённый, 5 колонок, 1-based бины):
        start_bin  end_bin  depth  score1  score2

    depth=0 → root (вся хромосома), пропускаем
    depth=1 → top-level TADs (непересекающиеся) ← используем
    depth>1 → sub-TADs (иерархические) ← пропускаем

    Конвертация: start_bp = (start_bin - 1) * resolution
                 end_bp   = end_bin * resolution
    """
    try:
        lines = [
            ln.strip() for ln in output_file.read_text().splitlines()
            if ln.strip() and not ln.startswith("#")
        ]
    except Exception as exc:
        logger.error("[ontad] Ошибка чтения %s: %s", output_file, exc)
        return pd.DataFrame(columns=["chrom", "start", "end"])

    if not lines:
        return pd.DataFrame(columns=["chrom", "start", "end"])

    records: list[dict] = []

    for line in lines:
        parts = line.split()
        if len(parts) < 3:
            continue
        try:
            start_bin = int(parts[0])
            end_bin   = int(parts[1])
            depth     = int(parts[2])   # ← исправлено: было parts[3]
        except (ValueError, IndexError):
            continue

        # Только верхний уровень иерархии (непересекающиеся TAD)
        if depth != 1:
            continue

        # OnTAD использует 1-based бины
        start_bp = (start_bin - 1) * resolution
        end_bp   = end_bin * resolution

        if end_bp > start_bp:
            records.append({"chrom": chrom, "start": start_bp, "end": end_bp})

    df = pd.DataFrame(records, columns=["chrom", "start", "end"])
    logger.debug("[ontad] Парсинг: %d строк → %d TADs (depth=1)", len(lines), len(df))
    return df

def _ontad_python_fallback(
    matrix: np.ndarray,
    chrom: str,
    resolution: int,
    penalty: float = 0.1,
    minsz: int = 3,
    maxsz: int = 200,
    log2: bool = True,
) -> pd.DataFrame:
    """
    Python-реализация OnTAD через 2D prefix-суммы.
    score(i,j) = mean_intra(i,j) - mean_inter_flanks(i,j)
    """
    n = matrix.shape[0]
    work = np.log2(matrix + 1.0) if log2 else matrix.astype(np.float64)

    # 2D prefix sum: P[i,j] = sum(work[0:i, 0:j])
    P = np.zeros((n + 1, n + 1), dtype=np.float64)
    P[1:, 1:] = np.cumsum(np.cumsum(work, axis=0), axis=1)

    def rsum(r0: int, c0: int, r1: int, c1: int) -> float:
        if r1 <= r0 or c1 <= c0:
            return 0.0
        return float(P[r1, c1] - P[r0, c1] - P[r1, c0] + P[r0, c0])

    def rmean(r0: int, c0: int, r1: int, c1: int) -> float:
        area = (r1 - r0) * (c1 - c0)
        return rsum(r0, c0, r1, c1) / area if area > 0 else 0.0

    # score_table[i, l] для l = minsz..maxsz  (только нужные)
    # Используем векторизацию по i для каждого l
    score_arr = np.zeros((n, maxsz + 1), dtype=np.float32)
    for l in range(minsz, maxsz + 1):
        for i in range(n - l):
            j    = i + l
            intra = rmean(i, i, j, j)
            left  = rmean(i, max(0, i - l), j, i)  if i > 0 else 0.0
            right = rmean(i, j, j, min(n, j + l))   if j < n else 0.0
            inter = 0.5 * (left + right)
            score_arr[i, l] = intra - inter

    # DP
    dp     = np.full(n + 1, -1e18, dtype=np.float64)
    parent = np.full(n + 1, -1,    dtype=np.int32)
    dp[0]  = 0.0

    for j in range(minsz, n + 1):
        i_arr = np.arange(max(0, j - maxsz), j - minsz + 1, dtype=np.int32)
        l_arr = j - i_arr
        valid = (l_arr >= minsz) & (l_arr <= maxsz) & (l_arr < score_arr.shape[1])
        i_arr, l_arr = i_arr[valid], l_arr[valid]
        if not len(i_arr):
            continue
        dp_prev    = dp[i_arr]
        reachable  = dp_prev > -1e17
        if not reachable.any():
            continue
        cands = np.where(reachable, dp_prev + score_arr[i_arr, l_arr].astype(np.float64) - penalty, -1e18)
        best  = int(np.argmax(cands))
        if cands[best] > dp[j]:
            dp[j]     = cands[best]
            parent[j] = int(i_arr[best])

    tads: list[tuple[int, int]] = []
    pos = n
    while pos > 0:
        prev = int(parent[pos])
        if prev < 0:
            return pd.DataFrame(columns=["chrom", "start", "end"])
        tads.append((prev, pos))
        pos = prev
    tads.reverse()

    return pd.DataFrame(
        [{"chrom": chrom, "start": s * resolution, "end": e * resolution} for s, e in tads],
        columns=["chrom", "start", "end"],
    )
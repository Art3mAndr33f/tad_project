"""
run_armatus.py
==============
Subprocess-обёртка вокруг бинарника Armatus.

Алгоритм выбора gamma:
  1. Запустить Armatus для каждого gamma из gamma_values
  2. Для каждой пары (gamma_i, gamma_j) вычислить Jaccard на уровне доменов
  3. gamma с наибольшей суммарной стабильностью (средний Jaccard со всеми
     остальными) считается оптимальной
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

from src.data_prep import get_rawobserved_path

logger = logging.getLogger(__name__)


# ──────────────────────────────────────────────────────────────────────────────
# Внутренние утилиты
# ──────────────────────────────────────────────────────────────────────────────

def _parse_armatus_output(output_file: str, chrom: str) -> pd.DataFrame:
    """
    Распарсить .txt-файл выхода Armatus.
    Формат строки: <chrom>\\t<start>\\t<end>
    """
    records = []
    path = Path(output_file)
    if not path.exists():
        logger.warning("Armatus output не найден: %s", output_file)
        return _empty_df()

    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            try:
                start = int(parts[1])
                end   = int(parts[2])
                records.append((chrom, start, end))
            except ValueError:
                continue

    if not records:
        return _empty_df()
    return pd.DataFrame(records, columns=["chrom", "start", "end"])


def _empty_df() -> pd.DataFrame:
    return pd.DataFrame(columns=["chrom", "start", "end"])


def _jaccard_domains(df_a: pd.DataFrame, df_b: pd.DataFrame) -> float:
    """Jaccard index на уровне доменов (точное совпадение start и end)."""
    if df_a.empty or df_b.empty:
        return 0.0
    set_a = set(zip(df_a["start"], df_a["end"]))
    set_b = set(zip(df_b["start"], df_b["end"]))
    inter = len(set_a & set_b)
    union = len(set_a | set_b)
    return inter / union if union > 0 else 0.0


def _run_single_gamma(
    armatus_bin: str,
    raw_path: str,
    resolution: int,
    gamma: float,
    tmp_dir: str,
    chrom: str,
) -> pd.DataFrame:
    """Запустить Armatus для одного значения gamma."""
    out_prefix = os.path.join(tmp_dir, f"armatus_g{gamma:.3f}")

    cmd = [
        armatus_bin,
        "-i", raw_path,
        "-R",                           # RAWobserved format
        "-r", str(resolution),
        "-g", str(gamma),
        "-o", out_prefix,
        "-z",                           # не выводить в stdout
    ]
    logger.debug("Armatus cmd: %s", " ".join(cmd))

    try:
        proc = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=600,
        )
        if proc.returncode != 0:
            logger.warning(
                "Armatus завершился с кодом %d (gamma=%.1f): %s",
                proc.returncode, gamma, proc.stderr.decode()[:300],
            )
            return _empty_df()
    except subprocess.TimeoutExpired:
        logger.error("Armatus timeout при gamma=%.1f", gamma)
        return _empty_df()
    except FileNotFoundError:
        logger.error("Бинарник Armatus не найден: %s", armatus_bin)
        return _empty_df()

    # Armatus создаёт файл <prefix>.txt или <prefix>_level_0.txt
    candidates = [
        f"{out_prefix}.txt",
        f"{out_prefix}_level_0.txt",
        f"{out_prefix}_domain.txt",
    ]
    for c in candidates:
        if Path(c).exists():
            return _parse_armatus_output(c, chrom)

    logger.warning("Armatus output файл не найден для gamma=%.1f, prefix=%s",
                   gamma, out_prefix)
    return _empty_df()


def _select_best_gamma(
    results: dict[float, pd.DataFrame],
    top_n: int = 3,
) -> float:
    """
    Выбрать наиболее стабильное gamma:
    для каждого gamma вычислить средний Jaccard со всеми остальными gamma.
    """
    gammas = list(results.keys())
    if len(gammas) == 1:
        return gammas[0]

    stability: dict[float, float] = {}
    for g_a in gammas:
        scores = []
        for g_b in gammas:
            if g_a == g_b:
                continue
            scores.append(_jaccard_domains(results[g_a], results[g_b]))
        stability[g_a] = float(np.mean(scores)) if scores else 0.0

    logger.info("Armatus stability: %s",
                {f"{g:.1f}": f"{s:.3f}" for g, s in stability.items()})
    best = max(stability, key=stability.__getitem__)
    logger.info("Выбрано gamma=%.1f (stability=%.3f)", best, stability[best])
    return best


# ──────────────────────────────────────────────────────────────────────────────
# Публичный интерфейс
# ──────────────────────────────────────────────────────────────────────────────

def run_armatus(
    chrom: str,
    resolution: int,
    data_path: str,
    cfg: Optional[dict] = None,
    armatus_bin: str = "tools/armatus/armatus",
    gamma_values: Optional[list[float]] = None,
    raw_path: Optional[str] = None,
) -> pd.DataFrame:
    """
    Запустить Armatus, выбрать оптимальный gamma, вернуть TAD-листинг.

    Parameters
    ----------
    chrom       : хромосома ('chr17')
    resolution  : разрешение в bp
    data_path   : путь к директории processed (для поиска RAWobserved)
    cfg         : конфиг-словарь (опционально)
    armatus_bin : путь к бинарнику
    gamma_values: список значений gamma для перебора
    raw_path    : явный путь к RAWobserved (перекрывает автопоиск)

    Returns
    -------
    pd.DataFrame с колонками chrom, start, end
    """
    if cfg is not None:
        armatus_bin  = cfg["paths"]["armatus_bin"]
        gamma_values = cfg["algorithms"]["armatus"]["gamma_values"]

    gamma_values = gamma_values or [0.1, 0.5, 1.0, 2.0, 5.0]

    # Найти RAWobserved
    if raw_path is None:
        if cfg is not None:
            raw_path = str(get_rawobserved_path(cfg, chrom, resolution))
        else:
            raw_path = str(
                Path(data_path) / f"{chrom}_{resolution}bp.RAWobserved"
            )

    if not Path(raw_path).exists():
        logger.error("RAWobserved не найден: %s", raw_path)
        return _empty_df()

    logger.info("[Armatus] %s @ %d bp  | gammas=%s", chrom, resolution, gamma_values)

    results: dict[float, pd.DataFrame] = {}
    with tempfile.TemporaryDirectory() as tmp_dir:
        for gamma in gamma_values:
            df = _run_single_gamma(
                armatus_bin, raw_path, resolution, gamma, tmp_dir, chrom
            )
            results[gamma] = df
            logger.info("  gamma=%.1f → %d TADs", gamma, len(df))

        best_gamma = _select_best_gamma(results)

    df_best = results[best_gamma].copy()
    df_best = df_best[df_best["start"] < df_best["end"]].reset_index(drop=True)
    logger.info("[Armatus] Итого %d TAD (gamma=%.1f)", len(df_best), best_gamma)
    return df_best
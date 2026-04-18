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
    records = []
    path = Path(output_file)
    if not path.exists():
        logger.warning("Armatus output не найден: %s", output_file)
        return _empty_df()

    with open(path) as f:
        for line in f:
            parts = line.strip().split()
            # ← пропускаем строки где parts[1] или parts[2] не числа
            if len(parts) < 3:
                continue
            try:
                start = int(parts[1])
                end   = int(parts[2])
                records.append((chrom, start, end))
            except ValueError:
                continue   # строки с логом Armatus — пропускаем

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
        "-S",
        "-N",
        "-r", str(resolution),
        "-g", str(gamma),
        "-c", chrom,          # ← добавить, иначе N/A в выводе
        "-o", out_prefix,
        # убрать -j — он меняет логику на "только gamma_max", нам нужна каждая gamma
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
        f"{out_prefix}.consensus.txt",   # ← Armatus пишет сюда
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

# Размеры хромосом hg19 для расчёта ожидаемого числа TAD
_HG19_CHROM_SIZES_MB = {
    "chr1": 249.3, "chr2": 243.2, "chr3": 198.0, "chr4": 191.2,
    "chr5": 180.9, "chr6": 171.1, "chr7": 159.1, "chr8": 146.4,
    "chr9": 141.2, "chr10": 135.5, "chr11": 135.0, "chr12": 133.9,
    "chr13": 115.2, "chr14": 107.3, "chr15": 102.5, "chr16": 90.4,
    "chr17": 81.2, "chr18": 78.1, "chr19": 59.1, "chr20": 63.0,
    "chr21": 48.1, "chr22": 51.3, "chrX": 155.3,
}

# TAD density при 25kb: эмпирически ~1.0–2.0 TAD/Mb (Rao 2014, GM12878)
_TAD_DENSITY_MIN = 0.8   # TAD/Mb → нижняя граница valid zone
_TAD_DENSITY_MAX = 2.5   # TAD/Mb → верхняя граница valid zone


def _select_best_gamma(
    results: dict[float, pd.DataFrame],
    chrom: str = "chr17",
    top_n: int = 3,
) -> float:
    """
    Выбрать оптимальный gamma из перебора.

    Стратегия (аналог _auto_penalty в scKTLD):
    1. Вычислить число TAD и stability (Jaccard) для каждого gamma.
    2. Определить valid zone: min_tads ≤ n_tads ≤ max_tads
       на основе размера хромосомы × эмпирическая плотность TAD.
    3. Среди валидных — выбрать gamma с максимальной stability.
    4. Если valid zone пуста — выбрать gamma с n_tads,
       ближайшим к target_tads (середина диапазона).
    """
    gammas = list(results.keys())
    if len(gammas) == 1:
        return gammas[0]

    # ── Размер valid zone из размера хромосомы ─────────────────────────────
    chrom_mb   = _HG19_CHROM_SIZES_MB.get(chrom, 80.0)
    min_tads   = max(10,  int(chrom_mb * _TAD_DENSITY_MIN))
    max_tads   = max(50,  int(chrom_mb * _TAD_DENSITY_MAX))
    target_tads = (min_tads + max_tads) / 2.0

    logger.debug(
        "[Armatus] valid zone для %s: %d–%d TADs (target=%.0f, chrom=%.1f Mb)",
        chrom, min_tads, max_tads, target_tads, chrom_mb,
    )

    # ── Число TAD на каждый gamma ──────────────────────────────────────────
    n_tads_map: dict[float, int] = {
        g: len(df) for g, df in results.items()
    }

    # ── Stability: средний Jaccard со всеми остальными gamma ───────────────
    stability: dict[float, float] = {}
    for g_a in gammas:
        scores = [
            _jaccard_domains(results[g_a], results[g_b])
            for g_b in gammas if g_b != g_a
        ]
        stability[g_a] = float(np.mean(scores)) if scores else 0.0

    logger.info(
        "Armatus stability: %s",
        {f"{g:.2f}": f"{s:.3f}" for g, s in stability.items()},
    )
    logger.info(
        "Armatus n_tads:    %s | valid zone: %d–%d",
        {f"{g:.2f}": n for g, n in n_tads_map.items()},
        min_tads, max_tads,
    )

    # ── Шаг 1: valid zone ─────────────────────────────────────────────────
    valid_gammas = {
        g: stability[g]
        for g in gammas
        if min_tads <= n_tads_map[g] <= max_tads
    }

    if valid_gammas:
        best = max(valid_gammas, key=valid_gammas.__getitem__)
        logger.info(
            "Выбрано gamma=%.2f (stability=%.3f, TADs=%d) [valid zone]",
            best, stability[best], n_tads_map[best],
        )
        return best

    # ── Шаг 2: valid zone пуста → ближайший к target ──────────────────────
    logger.warning(
        "[Armatus] Valid zone пуста (min=%d max=%d). "
        "Все gammas: %s. Выбираю ближайший к target=%.0f",
        min_tads, max_tads,
        {f"{g:.2f}": n_tads_map[g] for g in gammas},
        target_tads,
    )
    best = min(gammas, key=lambda g: abs(n_tads_map[g] - target_tads))
    logger.info(
        "Выбрано gamma=%.2f (TADs=%d, closest to target=%.0f) [fallback]",
        best, n_tads_map[best], target_tads,
    )
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
            logger.info("  gamma=%.2f → %d TADs", gamma, len(df))

        best_gamma = _select_best_gamma(gamma_results, chrom=chrom, top_n=stability_top_n)

    df_best = results[best_gamma].copy()
    df_best = df_best[df_best["start"] < df_best["end"]].reset_index(drop=True)
    logger.info("[Armatus] Итого %d TAD (gamma=%.1f)", len(df_best), best_gamma)
    return df_best
from __future__ import annotations
import os
"""
visualization_sveta.py
================
Визуализация результатов TAD-детекции в двух стилях:

1) Browser-style:
   - повёрнутая на 45° диагональная Hi-C лента сверху
   - ниже: по строке на каждый алгоритм (TAD-блоки + опц. loop-дуги)
   - консенсусные границы как вертикальные линии

2) Distance-position:
   - прямоугольная карта «позиция × расстояние от диагонали» сверху
   - ниже: по строке на алгоритм (TAD-треугольники-«сталактиты» + опц. loop-точки)
   - консенсусные границы как вертикальные линии

Поддерживаемые форматы Hi-C:
  - .cool / .mcool (cooler)
  - .RAWobserved (Juicer sparse format)
  - Синтетическая матрица (если ничего не указано)

Автономный запуск:
    python visualization_sveta.py \\
        --bed algo1_tads_armatus.bed algo2_tads_topdom.bed \\
        --chrom chr1 --resolution 25000 \\
        --raw chr1_25000bp.RAWobserved \\
        --style both
"""


import logging
import os
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Iterable

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from matplotlib.path import Path as MplPath
from matplotlib.patches import PathPatch, Polygon, Rectangle
from matplotlib.colors import LogNorm
import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix

logger = logging.getLogger(__name__)

# ──────────────────────────────────────────────────────────────────────────────
# Цвета
# ──────────────────────────────────────────────────────────────────────────────

ALGO_COLORS = {
    "armatus":        "#1f77b4",
    "topdom":         "#ff7f0e",
    "scktld":         "#2ca02c",
    "coitad":         "#d62728",
    "dihmm":          "#9467bd",
    "ontad":          "#e377c2",
    "modularity_tad": "#17becf",
    "arrowhead":      "#000000",
}

_EXTRA_PALETTE = [
    "#8c564b", "#bcbd22", "#7f7f7f", "#aec7e8",
    "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
    "#c49c94", "#f7b6d2", "#dbdb8d", "#9edae5",
]

def _algo_color(name: str, index: int = 0) -> str:
    if name in ALGO_COLORS:
        return ALGO_COLORS[name]
    return _EXTRA_PALETTE[index % len(_EXTRA_PALETTE)]

try:
    from src.consensus import CONSENSUS_COLORS
except ImportError:
    CONSENSUS_COLORS = {
        1: "#cccccc",
        2: "#ff9800",
        3: "#e91e63",
        4: "#9c27b0",
        5: "#2196f3",
        6: "#4caf50",
        7: "#f44336",
    }


# ──────────────────────────────────────────────────────────────────────────────
# Загрузка RAWobserved (Juicer sparse format)
# ──────────────────────────────────────────────────────────────────────────────

def load_raw_observed(
    raw_path: str,
    resolution: int,
    chrom_size_bp: Optional[int] = None,
) -> np.ndarray:
    """
    Загрузить Hi-C матрицу из Juicer RAWobserved формата.

    Формат файла (TSV без заголовка):
        pos1_bp  pos2_bp  count
        50000    50000    2.0
        50000    18250000 1.0
        ...

    Parameters
    ----------
    raw_path     : путь к .RAWobserved файлу
    resolution   : разрешение в bp (для конверсии позиций в бины)
    chrom_size_bp: размер хромосомы в bp (если None — определяется из данных)

    Returns
    -------
    Dense numpy matrix (n_bins x n_bins), симметричная
    """
    logger.info("Загружаем RAWobserved: %s", raw_path)

    # Читаем sparse данные
    df = pd.read_csv(
        raw_path,
        sep="\t",
        header=None,
        names=["pos1", "pos2", "count"],
        dtype={"pos1": np.int64, "pos2": np.int64, "count": np.float64},
    )

    if df.empty:
        raise ValueError(f"Пустой файл: {raw_path}")

    # Конвертируем позиции в бины
    df["bin1"] = df["pos1"] // resolution
    df["bin2"] = df["pos2"] // resolution

    # Определяем размер матрицы
    max_bin = max(df["bin1"].max(), df["bin2"].max())

    if chrom_size_bp is not None:
        n_bins = int(np.ceil(chrom_size_bp / resolution))
        n_bins = max(n_bins, max_bin + 1)
    else:
        n_bins = max_bin + 1

    logger.info("  Sparse entries: %d, matrix size: %d x %d",
                len(df), n_bins, n_bins)

    # Строим sparse матрицу
    sparse = coo_matrix(
        (df["count"].values, (df["bin1"].values, df["bin2"].values)),
        shape=(n_bins, n_bins),
        dtype=np.float64,
    )

    # Конвертируем в dense и делаем симметричной
    dense = sparse.toarray()
    matrix = dense + dense.T

    # Диагональ была удвоена — исправляем
    np.fill_diagonal(matrix, np.diag(dense))

    # Статистика
    nonzero = np.count_nonzero(matrix)
    total = matrix.size
    sparsity = 100 * (1 - nonzero / total)
    logger.info("  Non-zero: %d (%.1f%% sparse)", nonzero, sparsity)

    return matrix


def _parse_resolution_from_filename(filepath: str) -> Optional[int]:
    """
    Попытаться извлечь разрешение из имени файла.

    Примеры:
        chr1_25000bp.RAWobserved  →  25000
        chr1_10kb.RAWobserved     →  10000
        matrix_50000.raw          →  50000
    """
    stem = Path(filepath).stem

    # Паттерн: _<NUMBER>bp
    m = re.search(r"_(\d+)bp", stem, re.IGNORECASE)
    if m:
        return int(m.group(1))

    # Паттерн: _<NUMBER>kb
    m = re.search(r"_(\d+)kb", stem, re.IGNORECASE)
    if m:
        return int(m.group(1)) * 1000

    # Паттерн: просто число в конце
    m = re.search(r"_(\d+)$", stem)
    if m:
        val = int(m.group(1))
        # Эвристика: если > 1000, скорее всего это bp
        if val >= 1000:
            return val

    return None


def _parse_chrom_from_filename(filepath: str) -> Optional[str]:
    """
    Попытаться извлечь имя хромосомы из имени файла.

    Примеры:
        chr1_25000bp.RAWobserved  →  chr1
        chromosome2_10kb.raw      →  chr2
    """
    stem = Path(filepath).stem.lower()

    # Паттерн: chr<N> или chromosome<N>
    m = re.match(r"(chr(?:omosome)?(\d+|x|y))", stem, re.IGNORECASE)
    if m:
        chrom = m.group(1)
        # Нормализуем к chr<N>
        if chrom.startswith("chromosome"):
            chrom = "chr" + chrom[10:]
        return chrom

    return None


# ──────────────────────────────────────────────────────────────────────────────
# Извлечение имени алгоритма из имени файла
# ──────────────────────────────────────────────────────────────────────────────

def _algo_name_from_filename(filepath: str) -> str:
    """
    Извлечь имя алгоритма из имени файла.

    Примеры:
        algo1_tads_armatus.bed  →  armatus
        algo2_tads_topdom.bed   →  topdom
        tads_scktld.bed         →  scktld
        anything.bed            →  anything
    """
    stem = Path(filepath).stem

    m = re.search(r"_tads_(\w+)$", stem, re.IGNORECASE)
    if m:
        return m.group(1).lower()

    m = re.match(r"tads_(\w+)$", stem, re.IGNORECASE)
    if m:
        return m.group(1).lower()

    return stem.lower()


# ──────────────────────────────────────────────────────────────────────────────
# Загрузка BED-файлов
# ──────────────────────────────────────────────────────────────────────────────

def load_bed_domains(
    bed_path: str,
    target_chrom: Optional[str] = None,
) -> pd.DataFrame:
    """
    Загрузить TAD-домены из BED-файла (3+ колонки: chrom, start, end).
    """
    df = pd.read_csv(
        bed_path,
        sep="\t",
        header=None,
        comment="#",
        usecols=[0, 1, 2],
        names=["chrom", "start", "end"],
        dtype={"chrom": str},
    )
    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"]   = pd.to_numeric(df["end"],   errors="coerce")
    df = df.dropna(subset=["start", "end"])
    df["start"] = df["start"].astype(int)
    df["end"]   = df["end"].astype(int)

    if not df.empty:
        has_chr = df["chrom"].iloc[0].startswith("chr")
        if target_chrom is not None:
            want_chr = target_chrom.startswith("chr")
            if want_chr and not has_chr:
                df["chrom"] = "chr" + df["chrom"]
            elif not want_chr and has_chr:
                df["chrom"] = df["chrom"].str.replace("^chr", "", regex=True)

    if target_chrom is not None:
        df = df[df["chrom"] == target_chrom].reset_index(drop=True)

    return df


# ──────────────────────────────────────────────────────────────────────────────
# Утилиты
# ──────────────────────────────────────────────────────────────────────────────

def _ensure_dir(path: str) -> None:
    if path:
        Path(path).mkdir(parents=True, exist_ok=True)


def _ensure_parent_dir(file_path: str) -> None:
    parent = os.path.dirname(file_path)
    if parent:
        _ensure_dir(parent)


def _log_matrix(matrix: np.ndarray) -> np.ndarray:
    m = matrix.astype(np.float64)
    m = np.where(m > 0, m, np.nan)
    return np.log1p(m)


def _bp_to_bin(bp: float, resolution: int) -> int:
    return int(np.floor(float(bp) / resolution))


def _format_region_label(chrom: str, start_bin: int, end_bin: int, resolution: int) -> str:
    s = start_bin * resolution / 1e6
    e = end_bin   * resolution / 1e6
    return f"{chrom}:{s:.1f}\u2013{e:.1f} Mb"


def _domain_to_bins(row: pd.Series, resolution: int) -> Tuple[int, int]:
    b0 = int(np.floor(float(row["start"]) / resolution))
    b1 = int(np.ceil(float(row["end"]) / resolution))
    return b0, b1


def _get_vis_params(cfg: Optional[dict], resolution: int, n_bins: int) -> dict:
    vis = cfg.get("visualization", {}) if cfg else {}
    max_dist_bp = int(vis.get("max_distance_bp", 3_000_000))
    panel_bp    = int(vis.get("panel_bp", 10_000_000))
    overlap_bp  = vis.get("panel_overlap_bp", None)
    return {
        "cmap":       vis.get("hic_colormap", "YlOrRd"),
        "max_d_bp":   max_dist_bp,
        "max_d":      max(1, min(n_bins - 1, int(np.ceil(max_dist_bp / resolution)))),
        "panel_bp":   panel_bp,
        "panel_bins": max(2, min(n_bins, int(np.ceil(panel_bp / resolution)))),
        "ov_bp":      overlap_bp,
        "dpi":        int(vis.get("dpi", 300)),
    }


def _auto_overlap(algo_results, resolution, panel_bins, fallback=50):
    lengths = []
    for df in algo_results.values():
        if df is None or df.empty or "start" not in df.columns:
            continue
        lens = (df["end"] - df["start"]).dropna()
        lens = lens[lens > 0] / resolution
        if len(lens):
            lengths.extend(lens.astype(int).tolist())
    if lengths:
        ov = max(20, int(np.percentile(lengths, 95)) + 2)
    else:
        ov = max(20, fallback)
    return min(ov, max(5, panel_bins // 3))


def _make_panels(n_bins, panel_bins, ov) -> List[Tuple[int, int]]:
    panel_bins = min(panel_bins, n_bins)
    ov = min(ov, panel_bins - 1) if panel_bins > 1 else 0
    step = max(1, panel_bins - ov)
    wins = []
    s = 0
    while s < n_bins:
        e = min(n_bins, s + panel_bins)
        if wins and s <= wins[-1][0]:
            break
        wins.append((s, e))
        if e >= n_bins:
            break
        ns = s + step
        if ns + panel_bins > n_bins:
            ns = max(0, n_bins - panel_bins)
        if ns <= s:
            break
        s = ns
    return wins


def _consensus_in_panel(cdf, resolution, w0, w1):
    if cdf is None or cdf.empty:
        return []
    if not {"position", "support"}.issubset(cdf.columns):
        return []
    out = []
    for _, r in cdf.iterrows():
        pos = float(r["position"])
        sup = int(r["support"])
        b = _bp_to_bin(pos, resolution)
        if w0 <= b < w1:
            out.append((pos / 1e6, CONSENSUS_COLORS.get(sup, "#888"), sup))
    return out


# ──────────────────────────────────────────────────────────────────────────────
# Loop helpers
# ──────────────────────────────────────────────────────────────────────────────

def _extract_loop_pairs(loops_df, resolution):
    if loops_df is None or loops_df.empty:
        return []
    cols = set(loops_df.columns)
    pairs = []
    try:
        if {"start1", "end1", "start2", "end2"}.issubset(cols):
            a1 = (pd.to_numeric(loops_df["start1"], errors="coerce") +
                  pd.to_numeric(loops_df["end1"],   errors="coerce")) / 2
            a2 = (pd.to_numeric(loops_df["start2"], errors="coerce") +
                  pd.to_numeric(loops_df["end2"],   errors="coerce")) / 2
        elif {"start1", "start2"}.issubset(cols):
            a1 = pd.to_numeric(loops_df["start1"], errors="coerce")
            a2 = pd.to_numeric(loops_df["start2"], errors="coerce")
        elif {"x1", "x2"}.issubset(cols):
            a1 = pd.to_numeric(loops_df["x1"], errors="coerce")
            a2 = pd.to_numeric(loops_df["x2"], errors="coerce")
        elif {"bin1", "bin2"}.issubset(cols):
            for b1, b2 in zip(
                pd.to_numeric(loops_df["bin1"], errors="coerce"),
                pd.to_numeric(loops_df["bin2"], errors="coerce"),
            ):
                if pd.notna(b1) and pd.notna(b2) and int(b1) != int(b2):
                    pairs.append((min(int(b1), int(b2)), max(int(b1), int(b2))))
            return pairs
        else:
            return []
        for x1, x2 in zip(a1, a2):
            if pd.isna(x1) or pd.isna(x2):
                continue
            b1 = _bp_to_bin(x1, resolution)
            b2 = _bp_to_bin(x2, resolution)
            if b1 != b2:
                pairs.append((min(b1, b2), max(b1, b2)))
    except Exception as exc:
        logger.warning("Loop parse error: %s", exc)
    return pairs


def _draw_loop_arcs(ax, loops_df, resolution, w0, w1, color, alpha=0.6):
    for b1, b2 in _extract_loop_pairs(loops_df, resolution):
        if b1 < w0 or b2 >= w1:
            continue
        x1 = b1 * resolution / 1e6
        x2 = b2 * resolution / 1e6
        span = x2 - x1
        h = min(0.72, 0.18 + 0.25 * np.sqrt(max(span, 1e-6)))
        verts = [(x1, 0.1), ((x1+x2)/2, 0.1+h), (x2, 0.1)]
        codes = [MplPath.MOVETO, MplPath.CURVE3, MplPath.CURVE3]
        ax.add_patch(PathPatch(
            MplPath(verts, codes),
            fc="none", ec=color, lw=0.9, alpha=alpha, zorder=3,
        ))


def _draw_loop_dots(ax, loops_df, resolution, w0, w1, max_d, color, size=18):
    xs, ys = [], []
    for b1, b2 in _extract_loop_pairs(loops_df, resolution):
        if b1 < w0 or b2 >= w1:
            continue
        span = b2 - b1
        if span <= 0 or span > max_d:
            continue
        xs.append((b1+b2)/2 * resolution / 1e6)
        ys.append(span * resolution / 1e6)
    if xs:
        ax.scatter(xs, ys, s=size, c=color, alpha=0.75,
                   edgecolors="white", linewidths=0.35, zorder=4)


# ──────────────────────────────────────────────────────────────────────────────
# TAD drawing primitives
# ──────────────────────────────────────────────────────────────────────────────

def _draw_tad_blocks(ax, df, resolution, w0, w1, color, alpha=0.85):
    """
    Browser-style: TAD = горизонтальный прямоугольник.

    Контрастные границы:
      - белая обводка (зазор между соседними TAD)
      - чёрные вертикальные засечки
      - маленькие треугольные маркеры ▼ сверху на каждой границе
    """
    if df is None or df.empty:
        return

    df_sorted = df.sort_values("start").reset_index(drop=True)
    boundary_positions = set()

    for _, d in df_sorted.iterrows():
        b0, b1 = _domain_to_bins(d, resolution)
        left  = max(b0, w0)
        right = min(b1, w1)
        if right <= left:
            continue

        x0_mb = left  * resolution / 1e6
        x1_mb = right * resolution / 1e6
        w_mb  = x1_mb - x0_mb

        # ── Основной блок
        ax.add_patch(Rectangle(
            (x0_mb, 0.15), w_mb, 0.70,
            facecolor=color,
            edgecolor="white",
            linewidth=2.0,
            alpha=alpha,
            zorder=2,
        ))

        boundary_positions.add(x0_mb)
        boundary_positions.add(x1_mb)

    # ── Засечки и маркеры на границах
    panel_left  = w0 * resolution / 1e6
    panel_right = w1 * resolution / 1e6

    for x_mb in sorted(boundary_positions):
        if panel_left <= x_mb <= panel_right:
            # Вертикальная засечка
            ax.plot(
                [x_mb, x_mb],
                [0.02, 0.98],
                color="black",
                linewidth=0.9,
                alpha=0.75,
                zorder=5,
            )
            # Маркер ▼ сверху
            ax.plot(
                x_mb, 0.98,
                marker="v",
                markersize=3.5,
                color="black",
                alpha=0.8,
                zorder=6,
                clip_on=False,
            )
            # Маркер ▲ снизу
            ax.plot(
                x_mb, 0.02,
                marker="^",
                markersize=3.5,
                color="black",
                alpha=0.8,
                zorder=6,
                clip_on=False,
            )


def _draw_tad_stalactites(ax, df, resolution, w0, w1, max_d, color, lw=1.15, alpha=0.95):
    if df is None or df.empty:
        return
    for _, d in df.iterrows():
        b0, b1 = _domain_to_bins(d, resolution)
        if b0 < w0 or b1 > w1:
            continue
        span = b1 - b0
        if span <= 0 or span > max_d:
            continue
        x0 = b0 * resolution / 1e6
        x1 = b1 * resolution / 1e6
        xm = (b0+b1)/2 * resolution / 1e6
        y  = span * resolution / 1e6
        ax.add_patch(Polygon(
            [[x0, 0], [x1, 0], [xm, y]],
            closed=True, fill=False,
            ec=color, lw=lw, alpha=alpha, zorder=3,
        ))


# ──────────────────────────────────────────────────────────────────────────────
# Hi-C rendering helpers
# ──────────────────────────────────────────────────────────────────────────────

def _plot_rotated_hic(ax, log_sub, w0, w1, resolution, max_d, cmap, vmin, vmax):
    n = log_sub.shape[0]
    edges = np.arange(n + 1, dtype=float)
    I, J = np.meshgrid(edges, edges, indexing="ij")
    X = ((I + J) / 2.0 + w0) * resolution / 1e6
    Y = (J - I) * resolution / 1e6
    ii, jj = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")
    mask = ((jj - ii) < 0) | ((jj - ii) > max_d)
    C = np.ma.array(log_sub, mask=mask)
    im = ax.pcolormesh(X, Y, C, cmap=cmap, vmin=vmin, vmax=vmax,
                       shading="flat", rasterized=True)
    ax.set_xlim(w0 * resolution / 1e6, w1 * resolution / 1e6)
    ax.set_ylim(0, max_d * resolution / 1e6)
    return im


def _build_dist_pos_map(log_sub, max_d):
    n = log_sub.shape[0]
    max_d = min(max_d, n - 1)
    out = np.full((max_d + 1, n), np.nan, dtype=float)
    for d in range(max_d + 1):
        diag = np.diag(log_sub, k=d)
        out[d, :len(diag)] = diag
    return out


# ══════════════════════════════════════════════════════════════════════════════
# ВАРИАНТ 3: Browser-style
# ══════════════════════════════════════════════════════════════════════════════

def _load_consensus_bed(path) -> "pd.DataFrame":
    """Загрузить BED-файл консенсусных доменов (с заголовком support).

    Принимает str или Path. Возвращает DataFrame(chrom,start,end,support)
    или пустой DataFrame при отсутствии файла / некорректном формате.
    """
    from pathlib import Path as _Path
    import pandas as _pd

    p = _Path(path)
    if not p.exists():
        logger.debug("_load_consensus_bed: файл не найден: %s", p)
        return _pd.DataFrame(columns=["chrom", "start", "end", "support"])
    try:
        df = _pd.read_csv(p, sep="\t", comment="#")
    except Exception as exc:
        logger.warning("_load_consensus_bed: ошибка чтения %s: %s", p, exc)
        return _pd.DataFrame(columns=["chrom", "start", "end", "support"])
    for col in ("chrom", "start", "end", "support"):
        if col not in df.columns:
            # BED без заголовка (4-я колонка = name типа tad_consensus_support3)
            if col == "support" and df.shape[1] >= 5:
                try:
                    df["support"] = df.iloc[:, 4].astype(int)
                except Exception:
                    df["support"] = 2
            elif col == "support":
                df["support"] = 2
            else:
                logger.warning("_load_consensus_bed: нет колонки '%s' в %s", col, p)
                return _pd.DataFrame(columns=["chrom", "start", "end", "support"])
    return df[["chrom", "start", "end", "support"]].copy()


def plot_tad_browser_view(
    matrix: np.ndarray,
    algo_results: Dict[str, pd.DataFrame],
    consensus_df: Optional[pd.DataFrame],
    chrom: str,
    resolution: int,
    out_path: str,
    cfg: Optional[dict] = None,
    loop_results: Optional[Dict[str, pd.DataFrame]] = None,
    dpi: int = 300,
    figsize: Optional[Tuple[float, float]] = None,
    weak_consensus_path=None,    # str | None → tad_consensus_*.bed
    strong_consensus_path=None,  # str | None → strong_tad_consensus_*.bed

) -> None:
    _ensure_parent_dir(out_path)
    n_bins = matrix.shape[0]
    if n_bins < 2:
        logger.warning("Matrix too small for browser-view: %s", out_path)
        return

    p = _get_vis_params(cfg, resolution, n_bins)
    panel_bins = p["panel_bins"]
    max_d = p["max_d"]
    cmap = p["cmap"]

    ov = (_auto_overlap(algo_results, resolution, panel_bins, max(20, max_d // 4))
          if p["ov_bp"] is None
          else max(1, int(np.ceil(int(p["ov_bp"]) / resolution))))

    windows = _make_panels(n_bins, panel_bins, ov)
    n_algos = len(algo_results)
    algo_items = list(algo_results.items())

    log_m = _log_matrix(matrix)
    vmin = np.nanpercentile(log_m, 5)
    vmax = np.nanpercentile(log_m, 99)

    if figsize is None:
        fig_w = 18
        panel_h = 2.8 + 0.42 * n_algos
        fig_h = max(4.5, len(windows) * panel_h + 1.5)
        figsize = (fig_w, fig_h)

    fig = plt.figure(figsize=figsize, dpi=dpi)
    outer = fig.add_gridspec(len(windows), 1, hspace=0.35)
    heat_axes = []
    last_im = None

    for pi, (w0, w1) in enumerate(windows):
        sub = log_m[w0:w1, w0:w1]
        subgs = outer[pi].subgridspec(
            1 + n_algos, 1,
            height_ratios=[3.0] + [0.65] * n_algos,
            hspace=0.06,
        )

        ax_hic = fig.add_subplot(subgs[0, 0])
        last_im = _plot_rotated_hic(ax_hic, sub, w0, w1, resolution,
                                     max_d, cmap, vmin, vmax)
        heat_axes.append(ax_hic)
        ax_hic.set_title(_format_region_label(chrom, w0, w1, resolution),
                         loc="left", fontsize=10, fontweight="bold", pad=2)
        ax_hic.set_ylabel("dist, Mb", fontsize=9)
        ax_hic.tick_params(axis="x", labelbottom=False)

        for x_mb, col, _ in _consensus_in_panel(consensus_df, resolution, w0, w1):
            ax_hic.axvline(x=x_mb, color=col, lw=0.9, alpha=0.9, zorder=4)

            # ── weak consensus: вертикальные линии (dashed) ──
            if weak_consensus_path is not None:
                _wdf = _load_consensus_bed(weak_consensus_path)
                if not _wdf.empty:
                    _wdf = _wdf[
                        (_wdf['chrom'] == chrom) &
                        (_wdf['end']   > w0 * 1e6) &
                        (_wdf['start'] < w1 * 1e6)
                    ]
                for _, _row in _wdf.iterrows():
                    _sup = min(int(_row['support']), 7)
                    _col = CONSENSUS_COLORS.get(_sup, '#FFD700')
                    for _bnd_bp in [_row['start'], _row['end']]:
                        ax_hic.axvline(
                            x=_bnd_bp / 1e6, color=_col,
                            linewidth=1.0, alpha=0.70,
                            linestyle='--', zorder=5,
                        )

            # ── strong consensus: заполненные прямоугольники ──
            if strong_consensus_path is not None:
                _sdf = _load_consensus_bed(strong_consensus_path)
                if not _sdf.empty:
                    _sdf = _sdf[
                        (_sdf['chrom'] == chrom) &
                        (_sdf['end']   > w0 * 1e6) &
                        (_sdf['start'] < w1 * 1e6)
                    ]
                for _, _row in _sdf.iterrows():
                    _sup = min(int(_row['support']), 7)
                    _col = CONSENSUS_COLORS.get(_sup, '#FFD700')
                    ax_hic.axvspan(
                        xmin=max(_row['start'] / 1e6, w0),
                        xmax=min(_row['end']   / 1e6, w1),
                        alpha=0.18, color=_col,
                        zorder=2, linewidth=0,
                    )
                    for _bnd_bp in [_row['start'], _row['end']]:
                        if w0 <= _bnd_bp / 1e6 <= w1:
                            ax_hic.axvline(
                                x=_bnd_bp / 1e6, color=_col,
                                linewidth=2.0, alpha=0.90,
                                linestyle='-', zorder=6,
                            )


        for i, (algo, df) in enumerate(algo_items):
            ax = fig.add_subplot(subgs[i + 1, 0], sharex=ax_hic)
            color = _algo_color(algo, i)
            ax.set_ylim(0, 1)
            ax.set_yticks([])
            for sp in ("top", "right", "left"):
                ax.spines[sp].set_visible(False)
            ax.grid(False)

            _draw_tad_blocks(ax, df, resolution, w0, w1, color)
            if loop_results and algo in loop_results:
                _draw_loop_arcs(ax, loop_results[algo], resolution, w0, w1, color)
            for x_mb, c, _ in _consensus_in_panel(consensus_df, resolution, w0, w1):
                ax.axvline(x=x_mb, color=c, lw=0.8, alpha=0.85, zorder=1)

            ax.text(-0.012, 0.5, algo, transform=ax.transAxes,
                    ha="right", va="center", fontsize=8.7,
                    color=color, fontweight="bold")
            is_last = (pi == len(windows) - 1) and (i == n_algos - 1)
            ax.tick_params(axis="x", labelbottom=is_last)
            if is_last:
                ax.set_xlabel(f"Genomic position on {chrom} (Mb)", fontsize=10)

    if last_im is not None:
        cbar = fig.colorbar(last_im, ax=heat_axes, fraction=0.015, pad=0.01)
        cbar.set_label("log(1 + count)", fontsize=9)

    handles = [
        mpatches.Patch(color=_algo_color(a, i), label=a)
        for i, a in enumerate(algo_results.keys())
    ]
    cons_h = [
        Line2D([0], [0], color=c, lw=2, label=f"consensus {n} algo")
        for n, c in sorted(CONSENSUS_COLORS.items())
        if consensus_df is not None and not consensus_df.empty
    ]
    if loop_results:
        handles.append(Line2D([0], [0], color="black", lw=1.2, label="loops (arcs)"))

    fig.legend(handles=handles + cons_h, loc="upper center",
               bbox_to_anchor=(0.5, 0.988),
               ncol=min(5, max(2, len(handles))), frameon=False, fontsize=8)
    fig.suptitle(
        f"Hi-C Browser-style TAD Comparison  \u00b7  {chrom} @ {resolution//1000} kb",
        y=0.998, fontsize=13, fontweight="bold",
    )
    fig.subplots_adjust(left=0.10, right=0.92, top=0.95, bottom=0.04)

    plt.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    logger.info("Browser-view saved: %s", out_path)


# ══════════════════════════════════════════════════════════════════════════════
# ВАРИАНТ 5: Distance-position
# ══════════════════════════════════════════════════════════════════════════════

def plot_tad_distance_view(
    matrix: np.ndarray,
    algo_results: Dict[str, pd.DataFrame],
    consensus_df: Optional[pd.DataFrame],
    chrom: str,
    resolution: int,
    out_path: str,
    cfg: Optional[dict] = None,
    loop_results: Optional[Dict[str, pd.DataFrame]] = None,
    dpi: int = 300,
    figsize: Optional[Tuple[float, float]] = None,
) -> None:
    _ensure_parent_dir(out_path)
    n_bins = matrix.shape[0]
    if n_bins < 2:
        logger.warning("Matrix too small for distance-view: %s", out_path)
        return

    p = _get_vis_params(cfg, resolution, n_bins)
    panel_bins = p["panel_bins"]
    max_d = p["max_d"]
    cmap = p["cmap"]

    ov = (_auto_overlap(algo_results, resolution, panel_bins, max(20, max_d // 4))
          if p["ov_bp"] is None
          else max(1, int(np.ceil(int(p["ov_bp"]) / resolution))))

    windows = _make_panels(n_bins, panel_bins, ov)
    n_algos = len(algo_results)
    algo_items = list(algo_results.items())
    max_d_mb = max_d * resolution / 1e6

    log_m = _log_matrix(matrix)
    vmin = np.nanpercentile(log_m, 5)
    vmax = np.nanpercentile(log_m, 99)

    if figsize is None:
        fig_w = 18
        panel_h = 3.2 + 0.95 * n_algos
        fig_h = max(5, len(windows) * panel_h + 1.5)
        figsize = (fig_w, fig_h)

    fig = plt.figure(figsize=figsize, dpi=dpi)
    outer = fig.add_gridspec(len(windows), 1, hspace=0.35)
    heat_axes = []
    last_im = None

    for pi, (w0, w1) in enumerate(windows):
        sub = log_m[w0:w1, w0:w1]
        dmap = _build_dist_pos_map(sub, max_d)

        subgs = outer[pi].subgridspec(
            1 + n_algos, 1,
            height_ratios=[3.3] + [1.05] * n_algos,
            hspace=0.08,
        )

        x0_mb = w0 * resolution / 1e6
        x1_mb = w1 * resolution / 1e6

        ax_hic = fig.add_subplot(subgs[0, 0])
        last_im = ax_hic.imshow(
            dmap, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax,
            interpolation="nearest", origin="upper",
            extent=[x0_mb, x1_mb, max_d_mb, 0],
        )
        heat_axes.append(ax_hic)
        ax_hic.set_ylabel("dist, Mb", fontsize=9)
        ax_hic.tick_params(axis="x", labelbottom=False)
        ax_hic.set_title(_format_region_label(chrom, w0, w1, resolution),
                         loc="left", fontsize=10, fontweight="bold", pad=2)

        for x_mb, col, _ in _consensus_in_panel(consensus_df, resolution, w0, w1):
            ax_hic.axvline(x=x_mb, color=col, lw=0.9, alpha=0.9, zorder=4)

        for i, (algo, df) in enumerate(algo_items):
            ax = fig.add_subplot(subgs[i + 1, 0], sharex=ax_hic)
            color = _algo_color(algo, i)
            ax.set_xlim(x0_mb, x1_mb)
            ax.set_ylim(0, max_d_mb)
            ax.invert_yaxis()
            ax.set_facecolor("#fbfbfb")
            for sp in ("top", "right"):
                ax.spines[sp].set_visible(False)
            for y in np.linspace(0, max_d_mb, 5):
                ax.axhline(y, color="#ececec", lw=0.6, zorder=0)

            _draw_tad_stalactites(ax, df, resolution, w0, w1, max_d, color)
            if loop_results and algo in loop_results:
                _draw_loop_dots(ax, loop_results[algo], resolution,
                                w0, w1, max_d, color)
            for x_mb, c, _ in _consensus_in_panel(consensus_df, resolution, w0, w1):
                ax.axvline(x=x_mb, color=c, lw=0.8, alpha=0.85, zorder=1)

            ax.set_yticks([])
            ax.text(-0.012, 0.5, algo, transform=ax.transAxes,
                    ha="right", va="center", fontsize=8.7,
                    color=color, fontweight="bold")
            is_last = (pi == len(windows) - 1) and (i == n_algos - 1)
            ax.tick_params(axis="x", labelbottom=is_last)
            if is_last:
                ax.set_xlabel(f"Genomic position on {chrom} (Mb)", fontsize=10)

    if last_im is not None:
        cbar = fig.colorbar(last_im, ax=heat_axes, fraction=0.015, pad=0.01)
        cbar.set_label("log(1 + count)", fontsize=9)

    handles = [
        Line2D([0], [0], color=_algo_color(a, i), lw=2, label=a)
        for i, a in enumerate(algo_results.keys())
    ]
    cons_h = [
        Line2D([0], [0], color=c, lw=2, label=f"consensus {n} algo")
        for n, c in sorted(CONSENSUS_COLORS.items())
        if consensus_df is not None and not consensus_df.empty
    ]
    if loop_results:
        handles.append(Line2D([0], [0], marker="o", color="black",
                              ls="None", ms=5, label="loops"))

    fig.legend(handles=handles + cons_h, loc="upper center",
               bbox_to_anchor=(0.5, 0.988),
               ncol=min(5, max(2, len(handles))), frameon=False, fontsize=8)
    fig.suptitle(
        f"Hi-C Distance-Position TAD Comparison  \u00b7  {chrom} @ {resolution//1000} kb",
        y=0.998, fontsize=13, fontweight="bold",
    )
    fig.subplots_adjust(left=0.10, right=0.92, top=0.95, bottom=0.04)

    plt.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    logger.info("Distance-view saved: %s", out_path)


# ══════════════════════════════════════════════════════════════════════════════
# Обёртка (обратная совместимость)
# ══════════════════════════════════════════════════════════════════════════════

def plot_tad_comparison(
    matrix, algo_results, consensus_df, chrom, resolution, out_path,
    cfg=None, dpi=300, figsize=None, style="browser", loop_results=None,
):
    func = plot_tad_browser_view if style == "browser" else plot_tad_distance_view
    func(matrix=matrix, algo_results=algo_results, consensus_df=consensus_df,
         chrom=chrom, resolution=resolution, out_path=out_path, cfg=cfg,
         loop_results=loop_results, dpi=dpi, figsize=figsize)


# ══════════════════════════════════════════════════════════════════════════════
# CTCF-профиль (одиночный)
# ══════════════════════════════════════════════════════════════════════════════

def plot_ctcf_profile(bin_centers, density, algo, chrom, resolution, out_path, dpi=300):
    _ensure_parent_dir(out_path)
    fig, ax = plt.subplots(figsize=(8, 4))
    color = _algo_color(algo)
    ax.plot(bin_centers / 1000, density, color=color, lw=2, label=algo)
    ax.axvline(0, color="black", lw=1, ls="--", alpha=0.7)
    ax.fill_between(bin_centers / 1000, density, alpha=0.2, color=color)
    ax.set_xlabel("Distance from TAD boundary (kb)")
    ax.set_ylabel("CTCF density (peaks/boundary/kb)")
    ax.set_title(f"CTCF Enrichment | {algo} | {chrom} @ {resolution//1000} kb")
    ax.legend()
    plt.tight_layout()
    plt.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    logger.debug("CTCF profile: %s", out_path)


# ══════════════════════════════════════════════════════════════════════════════
# Jaccard heatmap
# ══════════════════════════════════════════════════════════════════════════════

def plot_jaccard_heatmap(jaccard_df, chrom, resolution, out_path, dpi=300):
    import seaborn as sns
    _ensure_parent_dir(out_path)
    fig, ax = plt.subplots(figsize=(6, 5))
    sns.heatmap(jaccard_df, annot=True, fmt=".3f", cmap="YlOrRd",
                vmin=0, vmax=1, square=True, ax=ax,
                cbar_kws={"label": "Jaccard Index"})
    ax.set_title(f"Pairwise Jaccard | {chrom} @ {resolution//1000} kb")
    plt.tight_layout()
    plt.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


# ══════════════════════════════════════════════════════════════════════════════
# Сводный CTCF-профиль
# ══════════════════════════════════════════════════════════════════════════════

def _load_arrowhead_for_chrom(arrowhead_path, chrom):
    try:
        df = pd.read_csv(arrowhead_path, sep="\t", comment="#", header=0, dtype=str)
        df = df.iloc[:, :3].copy()
        df.columns = ["chrom", "start", "end"]
        df["start"] = pd.to_numeric(df["start"], errors="coerce")
        df["end"]   = pd.to_numeric(df["end"],   errors="coerce")
        df = df.dropna(subset=["start", "end"])
        if not df["chrom"].iloc[0].startswith("chr"):
            df["chrom"] = "chr" + df["chrom"].astype(str)
        return df[df["chrom"] == chrom].reset_index(drop=True)
    except Exception as exc:
        logger.warning("Arrowhead load error (%s): %s", arrowhead_path, exc)
        return pd.DataFrame(columns=["chrom", "start", "end"])


def plot_ctcf_profile_all_algos(
    algo_results, ctcf_df, chrom, resolution, out_path,
    arrowhead_path=None, dpi=300,
    profile_range_bp=500_000, profile_bin_bp=10_000,
):
    from src.validation import compute_ctcf_profile as _ctcf_profile
    from scipy.ndimage import gaussian_filter1d
    _ensure_parent_dir(out_path)

    fig, ax = plt.subplots(figsize=(11, 5))
    any_plotted = False

    if arrowhead_path and os.path.exists(arrowhead_path):
        ah_df = _load_arrowhead_for_chrom(arrowhead_path, chrom)
        if not ah_df.empty:
            try:
                bins, dens = _ctcf_profile(ah_df, ctcf_df, chrom, resolution,
                                           profile_range_bp=profile_range_bp,
                                           profile_bin_bp=profile_bin_bp)
                dens_f = dens.astype(float)
                dens_norm = dens_f / (dens_f.mean() + 1e-9)
                dens_smooth = gaussian_filter1d(dens_norm, sigma=1.5)
                ax.plot(bins/1000, dens_smooth, color="black", lw=2.5, zorder=10,
                        label=f"Arrowhead (ref, n={len(ah_df)})")
                any_plotted = True
            except Exception as exc:
                logger.warning("Arrowhead profile error: %s", exc)

    # ── Слой слабого консенсуса: вертикальные линии на границах ──
    if weak_consensus_path is not None:
        _wdf = _load_consensus_bed(weak_consensus_path)
        if not _wdf.empty:
            _wdf = _wdf[
                (_wdf['chrom'] == chrom) &
                (_wdf['end']   > region_start) &
                (_wdf['start'] < region_end)
            ]
        for _, _row in _wdf.iterrows():
            _sup = min(int(_row['support']), 7)
            _col = CONSENSUS_COLORS.get(_sup, '#FFD700')
            for _bnd in [_row['start'], _row['end']]:
                ax_hic.axvline(
                    x=_bnd, color=_col, linewidth=1.0,
                    alpha=0.70, linestyle='--', zorder=3,
                )

    # ── Слой сильного консенсуса: заполненные прямоугольники ────────
    if strong_consensus_path is not None:
        _sdf = _load_consensus_bed(strong_consensus_path)
        if not _sdf.empty:
            _sdf = _sdf[
                (_sdf['chrom'] == chrom) &
                (_sdf['end']   > region_start) &
                (_sdf['start'] < region_end)
            ]
        for _, _row in _sdf.iterrows():
            _sup = min(int(_row['support']), 7)
            _col = CONSENSUS_COLORS.get(_sup, '#FFD700')
            ax_hic.axvspan(
                xmin=max(float(_row['start']), float(region_start)),
                xmax=min(float(_row['end']),   float(region_end)),
                alpha=0.18, color=_col, zorder=2, linewidth=0,
            )
            for _bnd in [_row['start'], _row['end']]:
                if region_start <= _bnd <= region_end:
                    ax_hic.axvline(
                        x=_bnd, color=_col, linewidth=2.0,
                        alpha=0.90, linestyle='-', zorder=4,
                    )

    for algo, df in algo_results.items():
        if df is None or df.empty:
            continue
        try:
            bins, dens = _ctcf_profile(df, ctcf_df, chrom, resolution,
                                       profile_range_bp=profile_range_bp,
                                       profile_bin_bp=profile_bin_bp)
            dens_f = dens.astype(float)
            dens_norm = dens_f / (dens_f.mean() + 1e-9)
            dens_smooth = gaussian_filter1d(dens_norm, sigma=1.5)
            ax.plot(bins/1000, dens_smooth, color=_algo_color(algo), lw=1.6,
                    alpha=0.85, label=f"{algo} (n={len(df)})")
            any_plotted = True
        except Exception as exc:
            logger.warning("CTCF profile %s skipped: %s", algo, exc)

    if not any_plotted:
        plt.close(fig)
        return

    ax.axvline(0, color="black", lw=1, ls="--", alpha=0.5, zorder=5)
    ax.axhline(1, color="gray", lw=0.8, ls=":", alpha=0.6, zorder=4)
    ax.set_xlabel("Distance from TAD boundary (kb)", fontsize=11)
    ax.set_ylabel("Relative CTCF density", fontsize=11)
    ax.set_title(f"CTCF profile \u00b1{profile_range_bp//1000} kb | "
                 f"{chrom} @ {resolution//1000} kb",
                 fontsize=12, fontweight="bold")
    ax.set_xlim(-profile_range_bp/1000, profile_range_bp/1000)
    ax.legend(fontsize=8, ncol=2 if len(ax.lines) > 4 else 1,
              loc="upper right", framealpha=0.85)
    ax.grid(axis="y", lw=0.4, alpha=0.5)
    plt.tight_layout()
    plt.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    logger.info("All-algo CTCF profile: %s", out_path)


# ══════════════════════════════════════════════════════════════════════════════
# Batch-визуализация (для интеграции с основным pipeline)
# ══════════════════════════════════════════════════════════════════════════════

def run_all_visualization(
    all_results, consensus_all, cfg,
    matrices_cache=None, all_loops=None,
):
    from src.data_prep import get_matrix
    from src.statistics import compute_pairwise_matrix
    from src.validation import load_ctcf_peaks

    figures_dir  = cfg["paths"]["figures_out"]
    hic_dir      = cfg["paths"].get("hic_figures_out",
                       os.path.join(figures_dir, "hic_maps"))
    chipseq_dir  = cfg["paths"].get("chipseq_figures_out",
                       os.path.join(figures_dir, "chipseq_profiles"))
    os.makedirs(hic_dir,    exist_ok=True)
    os.makedirs(chipseq_dir, exist_ok=True)
    resolutions = cfg["resolutions"]
    chromosomes = cfg["chromosomes"]["all"]
    dpi = int(cfg.get("visualization", {}).get("dpi", 300))
    styles = cfg.get("visualization", {}).get("styles", ["browser", "distance"])
    ctcf_df = load_ctcf_peaks(cfg["paths"]["ctcf_bed"])

    for res in resolutions:
        for chrom in chromosomes:
            algo_dfs = {}
            for algo in all_results:
                df = all_results[algo].get(chrom, {}).get(res)
                if df is not None and not df.empty:
                    algo_dfs[algo] = df
            if not algo_dfs:
                continue

            loop_dfs = {}
            if all_loops:
                for algo in all_loops:
                    df = all_loops[algo].get(chrom, {}).get(res)
                    if df is not None and not df.empty:
                        loop_dfs[algo] = df

            try:
                if matrices_cache and (chrom, res) in matrices_cache:
                    matrix = matrices_cache[(chrom, res)]
                else:
                    matrix = get_matrix(cfg, chrom, res)
            except Exception as exc:
                logger.error("Matrix load: %s @ %d: %s", chrom, res, exc)
                continue

            cons_df = None
            if chrom in consensus_all and res in consensus_all[chrom]:
                cons_df = consensus_all[chrom][res]

            if "browser" in styles:
                path = os.path.join(hic_dir, f"hic_browser_{chrom}_{res}bp.png")
                try:
                    plot_tad_browser_view(
                        matrix, algo_dfs, cons_df, chrom, res, path,
                        cfg=cfg, loop_results=loop_dfs or None, dpi=dpi,
                    weak_consensus_path=weak_path,
                    strong_consensus_path=strong_path)
                except Exception as exc:
                    logger.error("Browser-view %s@%d: %s", chrom, res, exc)

            if "distance" in styles:
                path = os.path.join(hic_dir, f"hic_distance_{chrom}_{res}bp.png")
                try:
                    plot_tad_distance_view(
                        matrix, algo_dfs, cons_df, chrom, res, path,
                        cfg=cfg, loop_results=loop_dfs or None, dpi=dpi)
                except Exception as exc:
                    logger.error("Distance-view %s@%d: %s", chrom, res, exc)

            for algo, df in algo_dfs.items():
                try:
                    from src.validation import compute_ctcf_profile as _cp
                    bins, dens = _cp(df, ctcf_df, chrom, res)
                    ppath = os.path.join(chipseq_dir,
                                         f"ctcf_profile_{algo}_{chrom}_{res}bp.png")
                    plot_ctcf_profile(bins, dens, algo, chrom, res, ppath, dpi)
                except Exception:
                    pass

            try:
                ah = cfg.get("paths", {}).get("arrowhead_ref")
                apath = os.path.join(chipseq_dir,
                                     f"ctcf_profile_all_{chrom}_{res}bp.png")
                plot_ctcf_profile_all_algos(algo_dfs, ctcf_df, chrom, res,
                                            apath, arrowhead_path=ah, dpi=dpi)
            except Exception:
                pass


# ══════════════════════════════════════════════════════════════════════════════
# Автономный запуск: python visualization.py --bed ... --raw ...
# ══════════════════════════════════════════════════════════════════════════════

def _demo():
    """
    Быстрый запуск из командной строки.

    Примеры:
        # С RAWobserved матрицей:
        python visualization.py \\
            --bed algo1_tads_armatus.bed algo2_tads_topdom.bed \\
            --raw chr1_25000bp.RAWobserved \\
            --chrom chr1 --resolution 25000

        # Без матрицы (синтетика):
        python visualization.py \\
            --bed algo1_tads_armatus.bed algo2_tads_topdom.bed \\
            --chrom chr1 --resolution 25000

        # С .cool файлом:
        python visualization.py \\
            --bed algo1_tads_armatus.bed algo2_tads_topdom.bed \\
            --cool data.mcool \\
            --chrom chr1 --resolution 25000
    """
    import argparse

    parser = argparse.ArgumentParser(
        description="Визуализация TAD-доменов из BED-файлов",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Примеры:
  # RAWobserved (Juicer формат):
  python visualization.py --bed algo1_tads_armatus.bed algo2_tads_topdom.bed \\
      --raw chr1_25000bp.RAWobserved --chrom chr1 --resolution 25000

  # .cool / .mcool:
  python visualization.py --bed *.bed --cool data.mcool --chrom chr1 --resolution 25000

  # Синтетическая матрица:
  python visualization.py --bed *.bed --chrom chr1 --resolution 25000
        """,
    )
    parser.add_argument(
        "--bed", nargs="+", required=True,
        help="BED-файлы с TAD-доменами",
    )
    parser.add_argument("--chrom", default=None,
                        help="Хромосома (если не указана — извлекается из имени --raw)")
    parser.add_argument("--resolution", type=int, default=None,
                        help="Разрешение bp (если не указано — извлекается из имени --raw)")
    parser.add_argument("--raw", default=None,
                        help="Hi-C матрица в Juicer RAWobserved формате")
    parser.add_argument("--cool", default=None,
                        help="Hi-C матрица в .cool/.mcool формате")
    parser.add_argument("--style", default="both",
                        choices=["browser", "distance", "both"])
    parser.add_argument("--max-distance", type=int, default=3_000_000,
                        help="Макс. расстояние от диагонали (bp)")
    parser.add_argument("--panel-size", type=int, default=10_000_000,
                        help="Ширина панели (bp)")
    parser.add_argument("--out-prefix", default="tad_vis",
                        help="Префикс выходных файлов")
    parser.add_argument("--dpi", type=int, default=200)

    args = parser.parse_args()

    # ── Определяем chrom и resolution ──
    chrom = args.chrom
    resolution = args.resolution

    # Пытаемся извлечь из имени RAW-файла
    if args.raw:
        if chrom is None:
            chrom = _parse_chrom_from_filename(args.raw)
            if chrom:
                print(f"  Хромосома из имени файла: {chrom}")
        if resolution is None:
            resolution = _parse_resolution_from_filename(args.raw)
            if resolution:
                print(f"  Разрешение из имени файла: {resolution} bp")

    # Дефолты
    if chrom is None:
        chrom = "chr1"
        print(f"  Хромосома не указана, используем: {chrom}")
    if resolution is None:
        resolution = 25000
        print(f"  Разрешение не указано, используем: {resolution} bp")

    # ── Загружаем BED-файлы ──
    algo_results: Dict[str, pd.DataFrame] = {}
    for path in args.bed:
        if not os.path.exists(path):
            print(f"ОШИБКА: файл не найден: {path}")
            return
        algo = _algo_name_from_filename(path)
        df = load_bed_domains(path, target_chrom=chrom)
        print(f"  {algo}: {len(df)} доменов на {chrom}  ← {path}")
        algo_results[algo] = df

    if not algo_results:
        print("Нет данных для визуализации")
        return

    # ── Размер матрицы ──
    all_ends = []
    for df in algo_results.values():
        if not df.empty:
            all_ends.append(df["end"].max())
    max_bp = max(all_ends) if all_ends else 250_000_000
    n_bins_estimate = int(np.ceil(max_bp / resolution))

    # ── Загружаем матрицу ──
    matrix = None

    # Вариант 1: RAWobserved
    if args.raw is not None:
        if not os.path.exists(args.raw):
            print(f"ОШИБКА: RAW файл не найден: {args.raw}")
            return
        try:
            matrix = load_raw_observed(args.raw, resolution, chrom_size_bp=max_bp)
            print(f"  Матрица из RAWobserved: {matrix.shape[0]}x{matrix.shape[1]}")
        except Exception as exc:
            print(f"  ОШИБКА загрузки RAWobserved: {exc}")
            return

    # Вариант 2: .cool/.mcool
    elif args.cool is not None:
        try:
            import cooler
            cool_uri = args.cool
            if args.cool.endswith(".mcool"):
                cool_uri = f"{args.cool}::resolutions/{resolution}"
            clr = cooler.Cooler(cool_uri)
            matrix = clr.matrix(balance=False).fetch(chrom)
            print(f"  Матрица из .cool: {matrix.shape[0]}x{matrix.shape[1]}")
        except Exception as exc:
            print(f"  ОШИБКА загрузки .cool: {exc}")
            print("  Используем синтетическую матрицу")

    # Вариант 3: Синтетика
    if matrix is None:
        print("  Генерируем синтетическую Hi-C матрицу...")
        n_bins = n_bins_estimate
        matrix = np.zeros((n_bins, n_bins), dtype=float)

        # Экспоненциальный decay
        max_diag = min(n_bins, args.max_distance // resolution + 50)
        for d in range(max_diag):
            val = 1000 * np.exp(-d / 50.0)
            idx = np.arange(0, n_bins - d)
            matrix[idx, idx + d] = val
            matrix[idx + d, idx] = val

        # Усиление внутри TAD
        first_df = list(algo_results.values())[0]
        for _, dom in first_df.iterrows():
            b0 = max(0, min(int(dom["start"] // resolution), n_bins - 1))
            b1 = max(0, min(int(dom["end"]   // resolution), n_bins))
            if b1 > b0:
                sz = b1 - b0
                boost = 200 * np.exp(
                    -np.abs(np.subtract.outer(np.arange(sz), np.arange(sz))) / 15.0
                )
                matrix[b0:b1, b0:b1] += boost

        noise = np.random.poisson(lam=5, size=matrix.shape).astype(float)
        matrix = (matrix + matrix.T) / 2.0 + noise
        print(f"  Синтетическая матрица: {matrix.shape[0]}x{matrix.shape[1]}")

    # ── Конфиг ──
    cfg = {
        "visualization": {
            "hic_colormap": "YlOrRd",
            "max_distance_bp": args.max_distance,
            "panel_bp": args.panel_size,
            "panel_overlap_bp": None,
            "dpi": args.dpi,
        }
    }

    # ── Рисуем ──
    styles = ["browser", "distance"] if args.style == "both" else [args.style]

    for style in styles:
        out_path = f"{args.out_prefix}_{style}_{chrom}_{resolution}bp.png"
        print(f"\n  Строим {style}-view → {out_path}")

        if style == "browser":
            plot_tad_browser_view(
                matrix=matrix, algo_results=algo_results,
                consensus_df=None, chrom=chrom,
                resolution=resolution, out_path=out_path,
                cfg=cfg, dpi=args.dpi,
                    weak_consensus_path=weak_path,
                    strong_consensus_path=strong_path)
        else:
            plot_tad_distance_view(
                matrix=matrix, algo_results=algo_results,
                consensus_df=None, chrom=chrom,
                resolution=resolution, out_path=out_path,
                cfg=cfg, dpi=args.dpi,
            )

        print(f"  \u2713 Сохранено: {out_path}")

    print("\nГотово!")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    _demo()
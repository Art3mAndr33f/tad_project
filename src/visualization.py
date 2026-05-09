"""
visualization.py
================
Визуализация результатов TAD-детекции:

1. Hi-C карта (log-scale, coolwarm)
2. TAD-домены каждого алгоритма (треугольники под картой)
3. Консенсусные границы (вертикальные линии)
4. CTCF-профиль обогащения
5. Статистические графики
6. Интерактивный HTML через plotly
"""

from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LogNorm
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from src.consensus import CONSENSUS_COLORS

logger = logging.getLogger(__name__)

# Цвета алгоритмов
ALGO_COLORS = {
    "armatus": "#1f77b4",
    "topdom":  "#ff7f0e",
    "scktld":  "#2ca02c",
    "coitad":  "#d62728",
}

# ──────────────────────────────────────────────────────────────────────────────
# Вспомогательные функции
# ──────────────────────────────────────────────────────────────────────────────

def _ensure_dir(path: str) -> None:
    Path(path).mkdir(parents=True, exist_ok=True)


def _log_matrix(matrix: np.ndarray) -> np.ndarray:
    """Логарифмическое преобразование для визуализации."""
    m = matrix.astype(np.float64)
    m = np.where(m > 0, m, np.nan)
    return np.log1p(m)


def _draw_tad_triangles(
    ax: plt.Axes,
    domains_df: pd.DataFrame,
    resolution: int,
    n_bins: int,
    row_offset: float,
    row_height: float,
    color: str,
    alpha: float = 0.7,
) -> None:
    """
    Нарисовать TAD-домены как треугольники под Hi-C картой.
    Треугольник: основание = [start_bin, end_bin], вершина = (start+end)/2.
    """
    for _, domain in domains_df.iterrows():
        i0 = int(domain["start"] // resolution)
        i1 = int(domain["end"]   // resolution)
        i0 = max(0, min(i0, n_bins - 1))
        i1 = max(0, min(i1, n_bins - 1))
        if i1 <= i0:
            continue

        mid = (i0 + i1) / 2.0
        triangle = plt.Polygon(
            [[i0, row_offset],
             [i1, row_offset],
             [mid, row_offset + row_height]],
            closed=True,
            facecolor=color,
            edgecolor=color,
            alpha=alpha,
            linewidth=0.5,
        )
        ax.add_patch(triangle)


# ──────────────────────────────────────────────────────────────────────────────
# Главная функция визуализации
# ──────────────────────────────────────────────────────────────────────────────

def plot_tad_comparison(
    matrix: np.ndarray,
    algo_results: Dict[str, pd.DataFrame],
    consensus_df: Optional[pd.DataFrame],
    chrom: str,
    resolution: int,
    out_path: str,
    cfg: Optional[dict] = None,
    dpi: int = 300,
    figsize: Tuple[int, int] = (18, 14),
) -> None:
    """
    Создать PNG-визуализацию:
      - Hi-C карта (верх)
      - Треугольники TAD для каждого алгоритма (4 строки)
      - Консенсусные границы (вертикальные линии)
      - Легенда

    Parameters
    ----------
    matrix       : dense numpy-матрица (n×n)
    algo_results : {algo_name: domains_DataFrame}
    consensus_df : DataFrame консенсусных границ
    chrom, resolution : метаданные
    out_path     : путь для сохранения PNG
    """
    _ensure_dir(os.path.dirname(out_path))

    n_bins = matrix.shape[0]
    # Высота: hic-карта + 4 строки треугольников
    n_algo_rows = len(algo_results)
    row_height  = n_bins * 0.12  # относительно размера матрицы

    fig = plt.figure(figsize=figsize, dpi=dpi)
    # GridSpec: верхняя часть — карта, нижняя — треугольники
    gs = fig.add_gridspec(
        2, 1,
        height_ratios=[4, 1],
        hspace=0.05,
    )
    ax_hic  = fig.add_subplot(gs[0])
    ax_tads = fig.add_subplot(gs[1], sharex=ax_hic)

    # ── Hi-C карта ──────────────────────────────────────────
    log_m = _log_matrix(matrix)
    vmin  = np.nanpercentile(log_m, 5)
    vmax  = np.nanpercentile(log_m, 99)

    im = ax_hic.imshow(
        log_m,
        aspect="auto",
        cmap=cfg["visualization"]["hic_colormap"] if cfg else "coolwarm",
        vmin=vmin, vmax=vmax,
        interpolation="nearest",
        origin="upper",
    )
    plt.colorbar(im, ax=ax_hic, fraction=0.02, pad=0.01, label="log(1+count)")
    ax_hic.set_title(
        f"Hi-C TAD Comparison | {chrom} @ {resolution//1000} kb",
        fontsize=13, fontweight="bold",
    )
    ax_hic.set_ylabel("Genomic bins")
    ax_hic.tick_params(labelbottom=False)

    # ── Консенсусные границы на карте ───────────────────────
    if consensus_df is not None and not consensus_df.empty:
        for _, row in consensus_df.iterrows():
            bin_pos = int(row["position"]) // resolution
            support = int(row["support"])
            color   = CONSENSUS_COLORS.get(support, "#999999")
            ax_hic.axvline(x=bin_pos, color=color, lw=0.8, alpha=0.9)
            ax_hic.axhline(y=bin_pos, color=color, lw=0.8, alpha=0.9)

    # ── Треугольники TAD ────────────────────────────────────
    ax_tads.set_xlim(0, n_bins)
    ax_tads.set_ylim(-0.5, n_algo_rows)
    ax_tads.set_yticks(np.arange(n_algo_rows) + 0.5)
    ax_tads.set_yticklabels(list(algo_results.keys()), fontsize=9)
    ax_tads.set_xlabel(f"Genomic bins (1 bin = {resolution//1000} kb)")

    for i, (algo, domains_df) in enumerate(algo_results.items()):
        color = ALGO_COLORS.get(algo, f"C{i}")
        if not domains_df.empty:
            _draw_tad_triangles(
                ax_tads, domains_df, resolution, n_bins,
                row_offset=i, row_height=0.9, color=color,
            )

    # Консенсусные линии на треугольниках тоже
    if consensus_df is not None and not consensus_df.empty:
        for _, row in consensus_df.iterrows():
            bin_pos = int(row["position"]) // resolution
            support = int(row["support"])
            color   = CONSENSUS_COLORS.get(support, "#999999")
            ax_tads.axvline(x=bin_pos, color=color, lw=0.7, alpha=0.8)

    # ── Легенда ─────────────────────────────────────────────
    algo_patches = [
        mpatches.Patch(color=ALGO_COLORS.get(a, f"C{i}"), label=a.capitalize())
        for i, a in enumerate(algo_results.keys())
    ]
    cons_patches = [
        mpatches.Patch(color=c, label=f"Consensus {n} algo")
        for n, c in CONSENSUS_COLORS.items()
    ]
    legend = ax_hic.legend(
        handles=algo_patches + cons_patches,
        loc="upper right",
        fontsize=7,
        framealpha=0.8,
        ncol=2,
    )

    plt.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    logger.info("Сохранено: %s", out_path)


# ──────────────────────────────────────────────────────────────────────────────
# Интерактивный plotly HTML
# ──────────────────────────────────────────────────────────────────────────────

def plot_tad_comparison_html(
    matrix: np.ndarray,
    algo_results: Dict[str, pd.DataFrame],
    consensus_df: Optional[pd.DataFrame],
    chrom: str,
    resolution: int,
    out_path: str,
) -> None:
    """
    Создать интерактивный HTML через plotly.
    """
    _ensure_dir(os.path.dirname(out_path))
    n_bins = matrix.shape[0]
    log_m  = _log_matrix(matrix)

    fig = make_subplots(
        rows=2, cols=1,
        row_heights=[0.75, 0.25],
        shared_xaxes=True,
        vertical_spacing=0.02,
        subplot_titles=[
            f"Hi-C Heatmap | {chrom} @ {resolution//1000} kb",
            "TAD Domains",
        ],
    )

    # Hi-C карта
    fig.add_trace(
        go.Heatmap(
            z=np.nan_to_num(log_m),
            colorscale="RdBu_r",
            showscale=True,
            colorbar=dict(title="log(1+count)", len=0.5, y=0.7),
            hovertemplate="bin_x=%{x}<br>bin_y=%{y}<br>log_count=%{z:.2f}<extra></extra>",
        ),
        row=1, col=1,
    )

    # TAD-границы как вертикальные линии
    for i, (algo, domains_df) in enumerate(algo_results.items()):
        color = ALGO_COLORS.get(algo, f"hsl({i*90},70%,50%)")
        if domains_df.empty:
            continue
        for _, d in domains_df.iterrows():
            b0 = int(d["start"]) // resolution
            b1 = int(d["end"])   // resolution
            for b in [b0, b1]:
                fig.add_vline(
                    x=b, line_width=0.5, line_color=color,
                    opacity=0.5, row=1, col=1,
                )

    # Консенсусные границы
    if consensus_df is not None and not consensus_df.empty:
        for _, row in consensus_df.iterrows():
            bin_pos = int(row["position"]) // resolution
            support = int(row["support"])
            color   = CONSENSUS_COLORS.get(support, "#888888")
            fig.add_vline(
                x=bin_pos, line_width=1.5, line_color=color,
                opacity=0.9, row=1, col=1,
            )

    # Горизонтальные полосы TAD в нижней части
    for i, (algo, domains_df) in enumerate(algo_results.items()):
        color = ALGO_COLORS.get(algo, f"hsl({i*90},70%,50%)")
        if domains_df.empty:
            continue
        for _, d in domains_df.iterrows():
            b0 = int(d["start"]) // resolution
            b1 = int(d["end"])   // resolution
            fig.add_shape(
                type="rect",
                x0=b0, x1=b1, y0=i, y1=i + 0.9,
                fillcolor=color, opacity=0.6,
                line_width=0,
                row=2, col=1,
            )

    fig.update_layout(
        height=900,
        title=f"TAD Comparison | {chrom} @ {resolution//1000} kb",
        showlegend=False,
        template="plotly_white",
    )
    fig.update_xaxes(title_text=f"Genomic bin ({resolution//1000} kb)", row=2, col=1)
    fig.update_yaxes(showticklabels=False, row=2, col=1)

    fig.write_html(out_path)
    logger.info("HTML сохранён: %s", out_path)


# ──────────────────────────────────────────────────────────────────────────────
# CTCF-профиль
# ──────────────────────────────────────────────────────────────────────────────

def plot_ctcf_profile(
    bin_centers: np.ndarray,
    density: np.ndarray,
    algo: str,
    chrom: str,
    resolution: int,
    out_path: str,
    dpi: int = 300,
) -> None:
    """Нарисовать профиль CTCF-обогащения вокруг TAD-границ."""
    _ensure_dir(os.path.dirname(out_path))

    fig, ax = plt.subplots(figsize=(8, 4))
    color = ALGO_COLORS.get(algo, "steelblue")

    ax.plot(bin_centers / 1000, density, color=color, lw=2, label=algo)
    ax.axvline(0, color="black", lw=1, linestyle="--", alpha=0.7)
    ax.fill_between(bin_centers / 1000, density, alpha=0.2, color=color)

    ax.set_xlabel("Distance from TAD boundary (kb)")
    ax.set_ylabel("CTCF density (peaks/boundary/kb)")
    ax.set_title(
        f"CTCF Enrichment Profile | {algo} | {chrom} @ {resolution//1000} kb"
    )
    ax.legend()
    plt.tight_layout()
    plt.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    logger.debug("CTCF-профиль: %s", out_path)


# ──────────────────────────────────────────────────────────────────────────────
# Матрица Jaccard (тепловая карта)
# ──────────────────────────────────────────────────────────────────────────────

def plot_jaccard_heatmap(
    jaccard_df: pd.DataFrame,
    chrom: str,
    resolution: int,
    out_path: str,
    dpi: int = 300,
) -> None:
    """Тепловая карта попарного Jaccard index."""
    import seaborn as sns
    _ensure_dir(os.path.dirname(out_path))

    fig, ax = plt.subplots(figsize=(6, 5))
    sns.heatmap(
        jaccard_df,
        annot=True, fmt=".3f",
        cmap="YlOrRd", vmin=0, vmax=1,
        square=True, ax=ax,
        cbar_kws={"label": "Jaccard Index"},
    )
    ax.set_title(f"Pairwise Jaccard | {chrom} @ {resolution//1000} kb")
    plt.tight_layout()
    plt.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


# ──────────────────────────────────────────────────────────────────────────────
# Batch-визуализация
# ──────────────────────────────────────────────────────────────────────────────

def run_all_visualization(
    all_results: Dict,
    consensus_all: Dict,
    cfg: dict,
    matrices_cache: Optional[Dict] = None,
) -> None:
    """
    Запустить полный цикл визуализации для всех хромосом и разрешений.
    """
    from src.data_prep import get_matrix
    from src.statistics import compute_pairwise_matrix
    from src.validation import compute_ctcf_profile, load_ctcf_peaks

    figures_dir = cfg["paths"]["figures_out"]
    resolutions  = cfg["resolutions"]
    chromosomes  = cfg["chromosomes"]["all"]
    dpi          = cfg["visualization"]["dpi"]
    figsize      = tuple(cfg["visualization"]["figsize"])
    gen_html     = cfg["visualization"]["generate_html"]
    ctcf_df      = load_ctcf_peaks(cfg["paths"]["ctcf_bed"])

    for res in resolutions:
        for chrom in chromosomes:
            # Собрать результаты
            algo_dfs: Dict[str, pd.DataFrame] = {}
            for algo in all_results:
                df = all_results[algo].get(chrom, {}).get(res)
                if df is not None and not df.empty:
                    algo_dfs[algo] = df

            if not algo_dfs:
                continue

            # Загрузить матрицу
            try:
                if matrices_cache and (chrom, res) in matrices_cache:
                    matrix = matrices_cache[(chrom, res)]
                else:
                    matrix = get_matrix(cfg, chrom, res)
            except Exception as exc:
                logger.error("Матрица не загружена: %s @ %d: %s", chrom, res, exc)
                continue

            # Консенсус
            cons_df = None
            if chrom in consensus_all and res in consensus_all[chrom]:
                cons_df = consensus_all[chrom][res]

            # PNG
            png_path = os.path.join(
                figures_dir, f"hic_tads_{chrom}_{res}bp.png"
            )
            try:
                plot_tad_comparison(
                    matrix, algo_dfs, cons_df, chrom, res, png_path, cfg, dpi, figsize
                )
            except Exception as exc:
                logger.error("Ошибка PNG %s @ %d: %s", chrom, res, exc)

            # HTML
            if gen_html:
                html_path = os.path.join(
                    figures_dir, f"hic_tads_{chrom}_{res}bp.html"
                )
                try:
                    plot_tad_comparison_html(
                        matrix, algo_dfs, cons_df, chrom, res, html_path
                    )
                except Exception as exc:
                    logger.error("Ошибка HTML %s @ %d: %s", chrom, res, exc)

            # CTCF-профили
            for algo, df in algo_dfs.items():
                try:
                    from src.validation import compute_ctcf_profile as _cp
                    bins, dens = _cp(df, ctcf_df, chrom, res)
                    profile_path = os.path.join(
                        figures_dir, f"ctcf_profile_{algo}_{chrom}_{res}bp.png"
                    )
                    plot_ctcf_profile(bins, dens, algo, chrom, res, profile_path, dpi)
                except Exception as exc:
                    logger.debug("CTCF-профиль пропущен: %s %s @ %d: %s",
                                 algo, chrom, res, exc)

            # Jaccard heatmap
            if len(algo_dfs) >= 2:
                try:
                    jac_df, _ = compute_pairwise_matrix(algo_dfs, res)
                    jac_path  = os.path.join(
                        figures_dir, f"jaccard_{chrom}_{res}bp.png"
                    )
                    plot_jaccard_heatmap(jac_df, chrom, res, jac_path, dpi)
                except Exception as exc:
                    logger.debug("Jaccard heatmap пропущен: %s @ %d: %s", chrom, res, exc)
#!/usr/bin/env python3
"""
run_chipseq_validation.py
=========================
Валидация TAD-границ по ChIP-seq трекам (RAD21, SMC3, H3K4me3, H3K27ac).
Переиспользует run_ctcf_profile_analysis() из ctcf_analysis.py.
"""
import logging, sys, os, glob
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import yaml
import pandas as pd
from src.ctcf_analysis import run_ctcf_profile_analysis

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)

# ── Конфиг ────────────────────────────────────────────────────────────────────
with open("config/config.yaml") as f:
    cfg = yaml.safe_load(f)

TRACKS = {
    "rad21":    cfg["paths"]["rad21_bed"],
    "smc3":     cfg["paths"]["smc3_bed"],
    "h3k4me3":  cfg["paths"]["h3k4me3_bed"],
    "h3k27ac":  cfg["paths"]["h3k27ac_bed"],
}

CHROM      = "chr1"
RESOLUTION = 100_000


# ── Загрузить результаты алгоритмов из BED-файлов ────────────────────────────
def load_algo_results(chrom: str, res: int) -> dict:
    """Найти и загрузить TAD-домены всех алгоритмов."""
    algo_results = {}

    # Ищем в нескольких возможных директориях
    search_dirs = ["results/tads", "results/domains", "results"]
    patterns = [
        f"*_{chrom}_{res}bp.bed",
        f"*_{chrom}_{res}.bed",
        f"{chrom}_{res}bp_*.bed",
    ]

    found_files = []
    for d in search_dirs:
        if not os.path.isdir(d):
            continue
        for pat in patterns:
            hits = glob.glob(os.path.join(d, pat))
            found_files.extend(hits)

    if not found_files:
        logger.error(
            "Нет TAD BED-файлов! Проверил директории: %s",
            search_dirs
        )
        logger.error(
            "Ожидаемый паттерн: *_%s_%dbp.bed", chrom, res
        )
        return {}

    for path in sorted(set(found_files)):
        basename = os.path.basename(path)
        # Извлечь имя алгоритма: убрать суффикс _chr1_100000bp.bed
        algo = basename.replace(f"_{chrom}_{res}bp.bed", "") \
                       .replace(f"_{chrom}_{res}.bed", "")
        # Пропустить консенсусные файлы
        if "consensus" in algo.lower():
            continue
        try:
            df = pd.read_csv(
                path, sep="\t", header=None, comment="#",
                names=["chrom", "start", "end"],
                usecols=[0, 1, 2],
                dtype={"chrom": str, "start": int, "end": int},
            )
            df = df[df["chrom"] == chrom].dropna().reset_index(drop=True)
            if df.empty:
                logger.warning("  %s: нет данных для %s", algo, chrom)
                continue
            algo_results[algo] = df
            logger.info("  ✓ %s: %d доменов", algo, len(df))
        except Exception as e:
            logger.warning("  Пропущен %s: %s", basename, e)

    return algo_results


# ── Загрузить ChIP-seq трек ───────────────────────────────────────────────────
def load_chipseq_bed(path: str, chrom: str) -> pd.DataFrame:
    """Загрузить BED/narrowPeak/broadPeak — только первые 3 колонки."""
    df = pd.read_csv(
        path, sep="\t", header=None, comment="#",
        usecols=[0, 1, 2],
        names=["chrom", "start", "end"],
        dtype={"chrom": str},
    )
    # Добавить chr-префикс если отсутствует
    if not df["chrom"].iloc[0].startswith("chr"):
        df["chrom"] = "chr" + df["chrom"].astype(str)

    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"]   = pd.to_numeric(df["end"],   errors="coerce")
    df = df.dropna().astype({"start": int, "end": int})
    df = df[df["chrom"] == chrom].reset_index(drop=True)
    logger.info(
        "  Загружен %s: %d пиков на %s",
        os.path.basename(path), len(df), chrom
    )
    return df


# ── Основной цикл ─────────────────────────────────────────────────────────────
def main():
    logger.info("=" * 60)
    logger.info("ChIP-seq валидация TAD-границ")
    logger.info("Хромосома: %s | Разрешение: %d bp", CHROM, RESOLUTION)
    logger.info("=" * 60)

    # Загрузить TAD-домены
    logger.info("\n[1/3] Загрузка TAD-доменов...")
    algo_results = load_algo_results(CHROM, RESOLUTION)
    if not algo_results:
        logger.error("Нет алгоритмов — прерываем.")
        sys.exit(1)
    logger.info("Алгоритмов: %d → %s", len(algo_results), list(algo_results.keys()))

    arrowhead_ref = cfg.get("paths", {}).get("arrowhead_ref")
    os.makedirs("results/figures", exist_ok=True)
    os.makedirs("results/stats",   exist_ok=True)

    all_metrics = []

    logger.info("\n[2/3] Анализ профилей...")
    for track_name, bed_path in TRACKS.items():
        if not os.path.exists(bed_path):
            logger.warning("⚠️  Файл не найден, пропускаю: %s", bed_path)
            continue

        logger.info("\n─── %s ───", track_name.upper())
        chipseq_df = load_chipseq_bed(bed_path, CHROM)
        if chipseq_df.empty:
            logger.warning("Нет пиков на %s для %s", CHROM, track_name)
            continue

        out_csv = f"results/stats/{track_name}_enrichment.csv"

        df_metrics = run_ctcf_profile_analysis(
            algo_results     = algo_results,
            ctcf_df          = chipseq_df,      # параметризовано — не только CTCF
            chrom            = CHROM,
            resolution       = RESOLUTION,
            arrowhead_path   = arrowhead_ref,
            profile_range_bp = 500_000,
            profile_bin_bp   =  10_000,
            near_window_bp   = 100_000,
            center_window_bp =  50_000,
            out_csv          = out_csv,
        )
        df_metrics["track"] = track_name
        all_metrics.append(df_metrics)

    # Сводная таблица
    logger.info("\n[3/3] Сводная таблица...")
    if all_metrics:
        df_summary = pd.concat(all_metrics, ignore_index=True)
        summary_path = "results/stats/chipseq_validation_summary.csv"
        df_summary.to_csv(summary_path, index=False)
        logger.info("Сводка сохранена: %s", summary_path)

        # Сводка вердиктов по трекам
        print("\n" + "═" * 60)
        print("  СВОДКА ВЕРДИКТОВ")
        print("═" * 60)
        pivot = df_summary.groupby(["track", "verdict"]).size().unstack(fill_value=0)
        print(pivot.to_string())

    # ── Генерация фигур профилей ────────────────────────────────────────────
    logger.info("\n[+] Генерация фигур профилей...")
    try:
        from src.validation import compute_ctcf_profile
        from src.consensus import extract_boundaries
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        for track_name, bed_path in TRACKS.items():
            if not os.path.exists(bed_path):
                continue
            chipseq_df = load_chipseq_bed(bed_path, CHROM)
            if chipseq_df.empty:
                continue

            fig, axes = plt.subplots(
                2, 4, figsize=(20, 8), sharey=False
            )
            axes = axes.flatten()

            algo_list = list(algo_results.keys())
            for idx, (algo, df) in enumerate(algo_results.items()):
                if idx >= 8:
                    break
                try:
                    bins, dens = compute_ctcf_profile(
                        df, chipseq_df, CHROM, RESOLUTION,
                        profile_range_bp=500_000, profile_bin_bp=10_000,
                    )
                    ax = axes[idx]
                    ax.plot(bins / 1000, dens, lw=1.5, color="#2196F3")
                    ax.axvline(0, color="red", lw=1, ls="--", alpha=0.7)
                    ax.axhline(dens.mean(), color="gray", lw=0.8, ls=":", alpha=0.6)
                    ax.set_title(algo, fontsize=9)
                    ax.set_xlabel("Distance to boundary (kb)", fontsize=7)
                    ax.set_ylabel("Peak density", fontsize=7)
                    ax.tick_params(labelsize=7)
                except Exception as e:
                    logger.warning("  Фигура %s/%s: %s", track_name, algo, e)

            # Скрыть пустые субплоты
            for idx in range(len(algo_results), 8):
                axes[idx].set_visible(False)

            fig.suptitle(
                f"{track_name.upper()} enrichment at TAD boundaries | {CHROM} @ {RESOLUTION//1000}kb",
                fontsize=11, fontweight="bold"
            )
            plt.tight_layout()
            chipseq_fig_dir = cfg.get("paths", {}).get(
                "chipseq_figures_out", "results/figures/chipseq_profiles")
            os.makedirs(chipseq_fig_dir, exist_ok=True)
            out_png = os.path.join(chipseq_fig_dir,
                f"{track_name}_profile_all_algos_{CHROM}_{RESOLUTION}bp.png")
            plt.savefig(out_png, dpi=150, bbox_inches="tight")
            plt.close()
            logger.info("  Сохранено: %s", out_png)

    except ImportError as e:
        logger.warning("plot_ctcf_profile недоступен: %s — пропускаю фигуры", e)

    logger.info("\n✅ ChIP-seq валидация завершена")


if __name__ == "__main__":
    main()

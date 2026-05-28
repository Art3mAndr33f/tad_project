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
import numpy as np
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

# ─────────────────────────────────────────────────────────────────────────────
# Консенсусная ChIP-seq валидация (strong vs weak)
# ─────────────────────────────────────────────────────────────────────────────

def _boundaries_from_domain_bed(bed_path: str, chrom: str, resolution: int) -> pd.DataFrame:
    """Извлечь уникальные граничные позиции из BED-файла доменов.

    Поддерживает оба формата:
      - BED без заголовка (стандарт: col0=chrom, col1=start, col2=end)
      - TSV с заголовком (chrom/start/end/support/...)

    Из каждого домена [start, end] берём start и end как отдельные границы.
    Возвращает DataFrame(chrom, start, end).
    """
    p = bed_path
    if not os.path.exists(p):
        logger.warning("_boundaries_from_domain_bed: не найден %s", p)
        return pd.DataFrame(columns=["chrom", "start", "end"])
    try:
        # Читаем первую строку чтобы понять есть ли заголовок
        with open(p) as _f:
            first_line = _f.readline().strip()
        has_header = (not first_line.startswith("chr")) and ("chrom" in first_line or "start" in first_line)

        if has_header:
            df = pd.read_csv(p, sep="\t", comment="#")
            # Переименовать на случай нестандартных имён
            if "chrom" not in df.columns and df.columns[0] not in ("chrom",):
                df.columns = ["chrom","start","end"] + list(df.columns[3:])
        else:
            # BED без заголовка — берём первые 3 колонки
            df = pd.read_csv(p, sep="\t", comment="#", header=None,
                             usecols=[0, 1, 2],
                             names=["chrom", "start", "end"],
                             dtype={"chrom": str, "start": int, "end": int})
    except Exception as exc:
        logger.warning("_boundaries_from_domain_bed: ошибка чтения %s: %s", p, exc)
        return pd.DataFrame(columns=["chrom", "start", "end"])

    if df.empty:
        return pd.DataFrame(columns=["chrom", "start", "end"])

    df["chrom"] = df["chrom"].astype(str)
    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"]   = pd.to_numeric(df["end"],   errors="coerce")
    df = df.dropna(subset=["start","end"])

    sub = df[df["chrom"] == chrom].copy()
    if sub.empty:
        logger.warning("_boundaries_from_domain_bed: нет строк для %s в %s", chrom, os.path.basename(p))
        return pd.DataFrame(columns=["chrom", "start", "end"])

    rows = []
    for _, row in sub.iterrows():
        for pos in [int(row["start"]), int(row["end"])]:
            rows.append({"chrom": chrom, "start": pos, "end": pos + resolution})
    result = pd.DataFrame(rows).drop_duplicates(subset=["start"]).reset_index(drop=True)
    logger.info(
        "  _boundaries_from_domain_bed: %s → %d границ на %s",
        os.path.basename(p), len(result), chrom,
    )
    return result


def run_consensus_chipseq_validation(
    algo_results: dict,
    chrom: str = "chr1",
    resolution: int = 100_000,
) -> pd.DataFrame:
    """Прогнать strong/weak консенсус через compute_ctcf_profile для всех ChIP-seq треков.

    Дополняет algo_results записями 'strong_consensus' и 'weak_consensus'.
    Возвращает DataFrame строк для дозаписи в chipseq_validation_summary.csv.
    """
    from src.validation import compute_ctcf_profile

    consensus_dir = cfg.get("paths", {}).get("consensus_out", "results/consensus")
    chipseq_fig_dir = cfg.get("paths", {}).get(
        "chipseq_figures_out", "results/figures/chipseq_profiles"
    )
    os.makedirs(chipseq_fig_dir, exist_ok=True)

    beds = {
        "weak_consensus":   os.path.join(consensus_dir, f"tad_consensus_{chrom}_{resolution}bp.bed"),
        "strong_consensus": os.path.join(consensus_dir, f"strong_tad_consensus_{chrom}_{resolution}bp.bed"),
    }

    # Стиль линий для сравнительного графика
    STYLES = {
        "strong_consensus": dict(lw=2.5, ls="-",  color="#8B008B", label="Strong consensus"),
        "weak_consensus":   dict(lw=1.4, ls="--", color="#FF8C00", label="Weak consensus (Jaccard≥0.5)"),
        "arrowhead":        dict(lw=1.5, ls=":",  color="#333333", label="Arrowhead (ref)"),
    }

    summary_rows = []

    for track_name, bed_path in TRACKS.items():
        if not os.path.exists(bed_path):
            continue
        chipseq_df = load_chipseq_bed(bed_path, chrom)
        if chipseq_df.empty:
            continue

        profiles: dict = {}   # label -> (bins, density)

        # Прогнать Arrowhead как эталон
        arrowhead_path = cfg.get("paths", {}).get("arrowhead_ref")
        if arrowhead_path and os.path.exists(arrowhead_path):
            try:
                # Arrowhead имеет заголовок: chr1 x1 x2 chr2 y1 y2 ...
                # Нужные колонки: chr1(0), x1(1), x2(2)
                arr_df = pd.read_csv(arrowhead_path, sep="\t", comment="#",
                                     header=0)
                # Переименовать первые три колонки
                cols = list(arr_df.columns)
                rename = {cols[0]: "chrom", cols[1]: "start", cols[2]: "end"}
                arr_df = arr_df.rename(columns=rename)
                arr_df["chrom"] = arr_df["chrom"].astype(str)
                arr_df["start"] = pd.to_numeric(arr_df["start"], errors="coerce")
                arr_df["end"]   = pd.to_numeric(arr_df["end"],   errors="coerce")
                arr_df = arr_df.dropna(subset=["start","end"])
                arr_df = arr_df[arr_df["chrom"] == chrom]
                if not arr_df.empty:
                    rows_a = []
                    for _, row in arr_df.iterrows():
                        for pos in [int(row["start"]), int(row["end"])]:
                            rows_a.append({"chrom": chrom,
                                           "start": pos,
                                           "end":   pos + resolution})
                    arr_bnd = pd.DataFrame(rows_a).drop_duplicates(
                        subset=["start"]).reset_index(drop=True)
                    logger.info("  Arrowhead границ на %s: %d", chrom, len(arr_bnd))
                    if not arr_bnd.empty:
                        bins_a, dens_a = compute_ctcf_profile(
                            arr_bnd, chipseq_df, chrom, resolution,
                            profile_range_bp=500_000, profile_bin_bp=10_000,
                        )
                        profiles["arrowhead"] = (bins_a, dens_a)
            except Exception as exc:
                logger.warning("  Arrowhead профиль (%s): %s", track_name, exc)

        # Прогнать weak и strong консенсус
        for label, bed in beds.items():
            bnd_df = _boundaries_from_domain_bed(bed, chrom, resolution)
            if bnd_df.empty:
                logger.warning("  %s: нет границ (%s)", label, os.path.basename(bed))
                continue
            if len(bnd_df) < 10:
                logger.warning("  %s: слишком мало границ (%d) — пропускаем", label, len(bnd_df))
                continue
            try:
                bins, dens = compute_ctcf_profile(
                    bnd_df, chipseq_df, chrom, resolution,
                    profile_range_bp=500_000, profile_bin_bp=10_000,
                )
                profiles[label] = (bins, dens)

                # Метрика: central_ratio = mean(центр ±50kb) / mean(фланги >200kb)
                mask_ctr  = (np.abs(bins) <= 50_000)
                mask_flnk = (np.abs(bins) >= 200_000)
                ctr_val   = float(np.mean(dens[mask_ctr]))  if mask_ctr.any()  else float("nan")
                flnk_val  = float(np.mean(dens[mask_flnk])) if mask_flnk.any() else float("nan")
                central_ratio = round(ctr_val / flnk_val, 3) if flnk_val > 0 else float("nan")

                summary_rows.append({
                    "algorithm":      label,
                    "track":          track_name,
                    "chrom":          chrom,
                    "resolution":     resolution,
                    "n_boundaries":   len(bnd_df),
                    "center_ratio":   central_ratio,
                    "verdict":        "STRONG" if central_ratio >= 1.3
                                      else "GOOD" if central_ratio >= 1.1
                                      else "WEAK",
                })
                logger.info(
                    "  %-20s %-8s ctr_ratio=%.3f", label, track_name, central_ratio
                )
            except Exception as exc:
                logger.error("  compute_ctcf_profile (%s, %s): %s", label, track_name, exc)

        # ── Сравнительный PNG ─────────────────────────────────────────────
        if len(profiles) < 2:
            logger.warning("  Нет профилей для сравнения (%s)", track_name)
            continue

        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(10, 4))
        for label, (bins_p, dens_p) in profiles.items():
            style = STYLES.get(label, dict(lw=1.0, ls="-", label=label))
            ax.plot(bins_p / 1000, dens_p, **style)
        ax.axvline(0, color="gray", ls=":", lw=0.8, alpha=0.6)
        ax.set_xlabel("Distance from TAD boundary (kb)", fontsize=10)
        ax.set_ylabel(f"{track_name.upper()} peak density", fontsize=10)
        ax.set_title(
            f"{track_name.upper()}: Strong vs Weak consensus | {chrom} @ {resolution // 1000}kb",
            fontsize=11, fontweight="bold",
        )
        ax.legend(fontsize=9, framealpha=0.85)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        plt.tight_layout()

        out_png = os.path.join(
            chipseq_fig_dir,
            f"{track_name}_strong_vs_weak_consensus_{chrom}_{resolution}bp.png",
        )
        plt.savefig(out_png, dpi=150, bbox_inches="tight")
        plt.close()
        logger.info("  Сохранён сравнительный PNG: %s", out_png)

    return pd.DataFrame(summary_rows) if summary_rows else pd.DataFrame()


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


    # ── Консенсусная валидация (strong vs weak) ────────────────────────────
    logger.info("\n[2.5/3] Консенсусная валидация (strong vs weak)...")
    df_consensus_val = run_consensus_chipseq_validation(
        algo_results=algo_results,
        chrom=CHROM,
        resolution=RESOLUTION,
    )
    if not df_consensus_val.empty:
        # Дозаписать strong/weak напрямую в summary CSV (схема совпадает с run_ctcf_profile_analysis)
        summary_path = "results/stats/chipseq_validation_summary.csv"
        if os.path.exists(summary_path):
            existing = pd.read_csv(summary_path)
            # Удалить старые строки strong/weak перед перезаписью
            if "algorithm" in existing.columns:
                existing = existing[~existing["algorithm"].isin(
                    ["strong_consensus", "weak_consensus"])]
            df_merged = pd.concat([existing, df_consensus_val], ignore_index=True)
        else:
            df_merged = df_consensus_val
        df_merged.to_csv(summary_path, index=False)
        logger.info("  Summary обновлён: %d строк (добавлено %d)",
                    len(df_merged), len(df_consensus_val))

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

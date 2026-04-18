"""
run_pipeline.py
===============
Главный скрипт пайплайна с CLI-интерфейсом.

Пример запуска:
  python pipeline/run_pipeline.py \\
    --resolution 25000 \\
    --chroms chr17 chr18 chr19 chr20 chr21 chr22 \\
    --algorithms armatus topdom scktld coitad

Флаги:
  --resolution INT         разрешение в bp (10000/25000/50000/100000)
  --chroms STR [STR ...]   список хромосом
  --algorithms STR [...]   алгоритмы для запуска
  --config PATH            путь к config.yaml (default: config/config.yaml)
  --prep-data              только подготовить матрицы
  --skip-stats             пропустить расчёт статистики
  --skip-validation        пропустить CTCF-валидацию
  --skip-report            не генерировать HTML-отчёт
  --only-viz               только визуализация (TAD-листы должны быть)
  --force                  перезаписать существующие результаты
  --log-level              DEBUG/INFO/WARNING (default: INFO)
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd
import yaml

# Добавляем корень проекта в sys.path,
# чтобы работали импорты вида `from src.algorithms import ...`
_PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

def setup_logging(level: str) -> None:
    logging.basicConfig(
        level=getattr(logging, level.upper(), logging.INFO),
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler("pipeline.log", mode="a"),
        ],
    )


def load_config(path: str) -> dict:
    with open(path) as f:
        return yaml.safe_load(f)


def ensure_dirs(cfg: dict) -> None:
    for key in ["processed", "tads_out", "consensus_out",
                "stats_out", "figures_out"]:
        Path(cfg["paths"][key]).mkdir(parents=True, exist_ok=True)


# ──────────────────────────────────────────────────────────────────────────────
# TAD-листы: сохранение / загрузка
# ──────────────────────────────────────────────────────────────────────────────

def save_tad_bed(df: pd.DataFrame, out_path: str) -> None:
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, sep="\t", header=False, index=False)


def load_tad_bed(path: str) -> pd.DataFrame:
    if not Path(path).exists():
        return pd.DataFrame(columns=["chrom", "start", "end"])
    df = pd.read_csv(path, sep="\t", header=None,
                     names=["chrom", "start", "end"])
    return df


def tad_path(cfg: dict, algo: str, chrom: str, resolution: int) -> str:
    return os.path.join(
        cfg["paths"]["tads_out"],
        f"{algo}_{chrom}_{resolution}bp.bed",
    )


# ──────────────────────────────────────────────────────────────────────────────
# Запуск алгоритмов
# ──────────────────────────────────────────────────────────────────────────────

def run_detection(
    cfg: dict,
    algorithms: list[str],
    chromosomes: list[str],
    resolutions: list[int],
    force: bool = False,
) -> dict:
    """
    Запустить все алгоритмы и вернуть словарь результатов.
    {algo: {chrom: {resolution: DataFrame}}}
    """
    from src.algorithms import ALGORITHM_REGISTRY
    logger = logging.getLogger("run_detection")

    results: dict = defaultdict(lambda: defaultdict(dict))

    for res in resolutions:
        logger.info("═══ Разрешение: %d bp (%d kb) ═══", res, res // 1000)
        for chrom in chromosomes:
            for algo in algorithms:
                bed = tad_path(cfg, algo, chrom, res)

                if Path(bed).exists() and not force:
                    logger.info("[%s] %s @ %d — загружаю из кэша", algo, chrom, res)
                    df = load_tad_bed(bed)
                    results[algo][chrom][res] = df
                    continue

                # scKTLD — проверка ограничений памяти
                if algo == "scktld":
                    limits = cfg["chromosomes"]["scktld_limits"].get(res, [])
                    if limits and chrom not in limits:
                        logger.info(
                            "[scKTLD] Пропуск %s @ %d (ограничение памяти)", chrom, res
                        )
                        results[algo][chrom][res] = pd.DataFrame(
                            columns=["chrom", "start", "end"]
                        )
                        continue

                func = ALGORITHM_REGISTRY.get(algo)
                if func is None:
                    logger.error("Неизвестный алгоритм: %s", algo)
                    continue

                logger.info("[%s] Запуск: %s @ %d bp", algo, chrom, res)
                try:
                    df = func(
                        chrom=chrom,
                        resolution=res,
                        data_path=cfg["paths"]["processed"],
                        cfg=cfg,
                    )
                    results[algo][chrom][res] = df
                    save_tad_bed(df, bed)
                    logger.info("[%s] %s @ %d → %d TADs", algo, chrom, res, len(df))
                except Exception as exc:
                    logger.error("[%s] Ошибка %s @ %d: %s", algo, chrom, res, exc,
                                 exc_info=True)
                    results[algo][chrom][res] = pd.DataFrame(
                        columns=["chrom", "start", "end"]
                    )

    return dict(results)


# ──────────────────────────────────────────────────────────────────────────────
# HTML-отчёт
# ──────────────────────────────────────────────────────────────────────────────

def generate_html_report(stats_dict: dict, cfg: dict) -> None:
    """Сгенерировать сводный HTML-отчёт с помощью pandas + jinja2."""
    from jinja2 import Environment, BaseLoader
    logger = logging.getLogger("report")

    TEMPLATE = """
<!DOCTYPE html>
<html lang="ru">
<head>
  <meta charset="UTF-8">
  <title>TAD Consensus Pipeline — Report</title>
  <style>
    body { font-family: Arial, sans-serif; margin: 20px; background: #f8f9fa; }
    h1   { color: #2c3e50; }
    h2   { color: #34495e; border-bottom: 2px solid #3498db; padding-bottom: 4px; }
    table.dataframe {
      border-collapse: collapse; width: 100%; font-size: 12px;
      background: white; box-shadow: 0 1px 3px rgba(0,0,0,0.1);
    }
    table.dataframe th {
      background: #3498db; color: white; padding: 8px;
      text-align: left; position: sticky; top: 0;
    }
    table.dataframe td { padding: 6px 8px; border-bottom: 1px solid #ecf0f1; }
    table.dataframe tr:hover td { background: #ebf5fb; }
    .section { margin-bottom: 40px; }
    .note { color: #7f8c8d; font-size: 11px; margin-top: 6px; }
  </style>
</head>
<body>
  <h1>🧬 TAD Consensus Pipeline — Отчёт</h1>
  <p class="note">Данные: GM12878 (hg19), Rao et al. 2014 | Генерировано автоматически</p>

  <div class="section">
    <h2>1. Базовая статистика алгоритмов</h2>
    {{ basic_html }}
  </div>

  <div class="section">
    <h2>2. Попарный Jaccard Index</h2>
    {{ jaccard_html }}
  </div>

  <div class="section">
    <h2>3. Попарный Boundary Overlap Rate</h2>
    {{ overlap_html }}
  </div>

  <div class="section">
    <h2>4. Сравнение с эталоном Rao 2014 (Arrowhead)</h2>
    {{ rao_html }}
  </div>

  <div class="section">
    <h2>5. CTCF Enrichment</h2>
    {{ ctcf_html }}
  </div>
</body>
</html>
"""

    def df_to_html(df: pd.DataFrame, max_rows: int = 200) -> str:
        if df is None or df.empty:
            return "<p><i>Нет данных</i></p>"
        return df.head(max_rows).to_html(
            classes="dataframe", index=False, na_rep="—",
            float_format=lambda x: f"{x:.4f}",
        )

    env      = Environment(loader=BaseLoader())
    template = env.from_string(TEMPLATE)

    html = template.render(
        basic_html   = df_to_html(stats_dict.get("basic")),
        jaccard_html = df_to_html(stats_dict.get("pairwise_jaccard")),
        overlap_html = df_to_html(stats_dict.get("pairwise_overlap")),
        rao_html     = df_to_html(stats_dict.get("rao_comparison")),
        ctcf_html    = df_to_html(stats_dict.get("ctcf_enrichment")),
    )

    out_path = cfg["paths"]["report"]
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(html)
    logger.info("HTML-отчёт: %s", out_path)


# ──────────────────────────────────────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="TAD Consensus Pipeline — сравнительный анализ Hi-C TAD-детекции",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--config", default="config/config.yaml",
        help="Путь к config.yaml (default: config/config.yaml)",
    )
    parser.add_argument(
        "--resolution", type=int, nargs="+",
        default=None,
        help="Разрешения в bp (default: из конфига)",
    )
    parser.add_argument(
        "--chroms", nargs="+", default=None,
        help="Список хромосом (default: из конфига)",
    )
    parser.add_argument(
        "--algorithms", nargs="+",
        default=["armatus", "topdom", "scktld", "coitad"],
        choices=["armatus", "topdom", "scktld", "coitad"],
        help="Алгоритмы для запуска",
    )
    parser.add_argument(
        "--prep-data", action="store_true",
        help="Только подготовить матрицы из .hic",
    )
    parser.add_argument(
        "--only-viz", action="store_true",
        help="Только визуализация (TAD-листы уже есть)",
    )
    parser.add_argument(
        "--skip-stats",       action="store_true", help="Пропустить статистику"
    )
    parser.add_argument(
        "--skip-validation",  action="store_true", help="Пропустить CTCF-валидацию"
    )
    parser.add_argument(
        "--skip-report",      action="store_true", help="Не генерировать HTML-отчёт"
    )
    parser.add_argument(
        "--force", action="store_true",
        help="Перезаписать существующие результаты",
    )
    parser.add_argument(
        "--log-level", default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
    )
    return parser.parse_args()


# ──────────────────────────────────────────────────────────────────────────────
# Главная функция
# ──────────────────────────────────────────────────────────────────────────────

def main() -> None:
    args = parse_args()
    setup_logging(args.log_level)
    logger = logging.getLogger("pipeline")

    logger.info("╔══════════════════════════════════════╗")
    logger.info("║    TAD Consensus Pipeline  v1.0      ║")
    logger.info("╚══════════════════════════════════════╝")

    # ── Конфигурация ────────────────────────────────────────────────────────
    cfg = load_config(args.config)
    ensure_dirs(cfg)

    resolutions = args.resolution or cfg["resolutions"]
    chromosomes = args.chroms    or cfg["chromosomes"]["all"]
    algorithms  = args.algorithms

    logger.info("Разрешения: %s", [f"{r//1000}kb" for r in resolutions])
    logger.info("Хромосомы:  %s", chromosomes)
    logger.info("Алгоритмы:  %s", algorithms)

    # ── Подготовка данных ────────────────────────────────────────────────────
    if args.prep_data:
        from src.data_prep import prepare_all_matrices
        logger.info("Подготовка матриц из .hic ...")
        prepare_all_matrices(cfg, resolutions, chromosomes, force=args.force)
        logger.info("Подготовка завершена.")
        return

    # ── Только визуализация ──────────────────────────────────────────────────
    if args.only_viz:
        from src.visualization import run_all_visualization
        from src.consensus import compute_all_consensus

        logger.info("Режим только визуализации")
        # Загружаем готовые TAD-листы
        all_results: dict = defaultdict(lambda: defaultdict(dict))
        for algo in algorithms:
            for res in resolutions:
                for chrom in chromosomes:
                    bed = tad_path(cfg, algo, chrom, res)
                    all_results[algo][chrom][res] = load_tad_bed(bed)

        consensus_all = compute_all_consensus(
            all_results, cfg, out_dir=cfg["paths"]["consensus_out"]
        )
        run_all_visualization(all_results, consensus_all, cfg)
        return

    # ── Детекция TAD ─────────────────────────────────────────────────────────
    logger.info("Запуск детекции TAD ...")
    all_results = run_detection(cfg, algorithms, chromosomes, resolutions, args.force)

    # ── Консенсус ────────────────────────────────────────────────────────────
    logger.info("Расчёт консенсусных границ ...")
    from src.consensus import compute_all_consensus
    consensus_all = compute_all_consensus(
        all_results, cfg, out_dir=cfg["paths"]["consensus_out"]
    )

    stats_dict = {}

    # ── Статистика ───────────────────────────────────────────────────────────
    if not args.skip_stats:
        logger.info("Расчёт статистики ...")
        from src.statistics import aggregate_statistics, save_statistics
        stats = aggregate_statistics(all_results, cfg, consensus_all)
        save_statistics(stats, cfg["paths"]["stats_out"])
        stats_dict.update(stats)

    # ── Валидация CTCF ───────────────────────────────────────────────────────
    if not args.skip_validation:
        logger.info("CTCF-валидация ...")
        try:
            from src.validation import run_all_validation
            ctcf_stats = run_all_validation(all_results, consensus_all, cfg)
            ctcf_path  = os.path.join(cfg["paths"]["stats_out"], "ctcf_enrichment.csv")
            ctcf_stats.to_csv(ctcf_path, index=False)
            stats_dict["ctcf_enrichment"] = ctcf_stats
            logger.info("CTCF-валидация завершена: %s", ctcf_path)
        except Exception as exc:
            logger.error("Ошибка CTCF-валидации: %s", exc, exc_info=True)

    # ── Визуализация ─────────────────────────────────────────────────────────
    logger.info("Визуализация ...")
    try:
        from src.visualization import run_all_visualization
        run_all_visualization(all_results, consensus_all, cfg)
    except Exception as exc:
        logger.error("Ошибка визуализации: %s", exc, exc_info=True)

    # ── HTML-отчёт ───────────────────────────────────────────────────────────
    if not args.skip_report and stats_dict:
        logger.info("Генерация HTML-отчёта ...")
        generate_html_report(stats_dict, cfg)

    logger.info("╔══════════════════════════════════════╗")
    logger.info("║       Пайплайн завершён успешно      ║")
    logger.info("╚══════════════════════════════════════╝")


if __name__ == "__main__":
    main()
# TAD Consensus Pipeline — Полная реализация

## Структура проекта

```
tad_consensus/
├── README.md
├── .gitignore
├── requirements.txt
├── environment.yml
├── config/
│   └── config.yaml
├── src/
│   ├── __init__.py
│   ├── data_prep.py
│   ├── algorithms/
│   │   ├── __init__.py
│   │   ├── run_armatus.py
│   │   ├── run_topdom.py
│   │   ├── run_scktld.py
│   │   └── run_coitad.py
│   ├── consensus.py
│   ├── statistics.py
│   ├── validation.py
│   └── visualization.py
├── pipeline/
│   └── run_pipeline.py
├── notebooks/
│   └── exploration.ipynb
└── tests/
    ├── test_consensus.py
    └── test_statistics.py
```

---

# TAD Consensus Pipeline

Сравнительный анализ детекции топологически ассоциированных доменов (TAD)
на bulk Hi-C данных GM12878 (Rao et al. 2014, hg19). Запускает четыре
алгоритма (Armatus, TopDom, scKTLD, coiTAD), строит консенсус границ и
оценивает качество через обогащение CTCF и сравнение с эталоном Arrowhead.

---

## Консенсусная логика и цветовая схема

Граница считается консенсусной, если её поддерживают ≥2 алгоритмов
(допуск совпадения ±1 bin):

| Поддержка        | Цвет                        | Hex       |
|------------------|-----------------------------|-----------|
| 1 алгоритм       | —  (не консенсус)           | —         |
| 2 алгоритма      | 🟡 Жёлтый   (слабый)        | `#FFD700` |
| 3 алгоритма      | 🟠 Оранжевый (умеренный)    | `#FF8C00` |
| 4 алгоритма      | 🟢 Зелёный  (сильный)       | `#00C800` |

Консенсус вычисляется независимо для каждого разрешения и хромосомы.

---

## Требования

- **OS**: Ubuntu 22.04
- **Python**: 3.10+
- **Armatus**: бинарник в `tools/armatus/armatus` (собран из исходников)
- **RAM**: ≥32 GB рекомендуется для 10 kb

---

## Установка

```bash
# 1. Клонировать репозиторий
git clone <repo_url> && cd tad_project

# 2. Создать conda-окружение
conda env create -f environment.yml
conda activate tad_pipeline

# 3. Доустановить pip-зависимости (если нужно)
pip install -r requirements.txt
```

---

## Подготовка данных

Положить файлы в `data/raw/`:

```
data/raw/
├── GSE63525_GM12878_insitu_primary_30.hic
├── chr1_10kb.RAWobserved      # и аналогичные файлы для других хромосом/разрешений
├── ...
data/reference/
├── GSE63525_GM12878_primary+replicate_Arrowhead_domainlist.txt
└── GM12878_CTCF_peaks_hg19.bed
```

RAWobserved-файлы можно извлечь автоматически:
```bash
python pipeline/run_pipeline.py --prep-data --resolution 25000
```

---

## Запуск

```bash
# Полный пайплайн на chr17-chr22, разрешение 25kb
python pipeline/run_pipeline.py \
    --resolution 25000 \
    --chroms chr17 chr18 chr19 chr20 chr21 chr22 \
    --algorithms armatus topdom scktld coitad

# Только детекция + консенсус, без отчёта
python pipeline/run_pipeline.py \
    --resolution 100000 --chroms chr1 chr2 \
    --algorithms armatus topdom coitad \
    --skip-report

# Только визуализация (TAD-листы уже есть)
python pipeline/run_pipeline.py \
    --resolution 50000 --chroms chr21 \
    --only-viz
```

---

## Структура результатов

```
results/
├── tads/          # <algorithm>_<chrom>_<resolution>.bed
├── consensus/     # consensus_<chrom>_<resolution>.bed  (score = кол-во алгоритмов)
├── stats/         # per_algorithm_stats.csv, pairwise_jaccard.csv,
│                  # rao_comparison.csv, ctcf_enrichment.csv
├── figures/       # hic_tads_<chrom>_<resolution>.png
└── report.html    # сводный интерактивный отчёт
```
<<<<<<< HEAD

## Обновление правил проекта
```bash
python scripts/update_rules.py check        # проверить
python scripts/update_rules.py add-algorithm <Name>  # новый алгоритм
python scripts/update_rules.py changelog    # зафиксировать изменения
```
=======
>>>>>>> e516dc2480df2b5d9bd89745443958a8bcb5fb74

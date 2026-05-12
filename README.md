# TAD Consensus Pipeline

Сравнительный анализ детекции топологически ассоциированных доменов (TAD)
на bulk Hi-C данных GM12878 (Rao et al. 2014, hg19). Запускает **семь**
алгоритмов (Armatus, TopDom, scKTLD, coiTAD, DI+HMM, OnTAD, ModularityTAD),
строит консенсус границ и оценивает качество через обогащение CTCF
и сравнение с эталоном Arrowhead.

---

## Структура проекта

```
tad_consensus/
├── README.md
├── .gitignore
├── requirements.txt
├── environment.yml
├── download_data.sh
├── config/
│   └── config.yaml
├── data/
│   ├── raw/           # RAWobserved файлы — НЕ в git
│   ├── processed/     # .npy матрицы — НЕ в git
│   └── reference/     # Arrowhead + CTCF BED — НЕ в git
├── tools/
│   ├── armatus/       # бинарник armatus
│   ├── coiTAD/        # Python-версия coiTAD
│   ├── ontad/         # бинарник OnTAD
│   └── ontad_src/     # исходники OnTAD (git clone)
├── src/
│   ├── __init__.py
│   ├── data_prep.py
│   ├── algorithms/
│   │   ├── __init__.py            # ALGORITHM_REGISTRY
│   │   ├── run_armatus.py
│   │   ├── run_topdom.py
│   │   ├── run_scktld.py
│   │   ├── run_coitad.py
│   │   ├── run_dihmm.py           # DI + Hidden Markov Model
│   │   ├── run_ontad.py           # OnTAD (иерархический)
│   │   └── run_modularity_tad.py  # Graph modularity DP
│   ├── consensus.py
│   ├── statistics.py
│   ├── validation.py
│   └── visualization.py
├── pipeline/
│   └── run_pipeline.py
├── tests/
│   ├── test_consensus.py
│   └── test_statistics.py
└── notebooks/
    └── exploration.ipynb
```

---

## Алгоритмы

| Алгоритм | Метод | Статус |
|----------|-------|--------|
| Armatus | Dynamic programming (gamma) | ✅ |
| TopDom | Insulation Score | ✅ |
| scKTLD | Spectral clustering (GPU) | ✅ |
| coiTAD | OI matrix + IS fallback | ✅ |
| DI+HMM | Directionality Index + HMM | ✅ |
| OnTAD | Hierarchical sliding average | ⚠️ калибровка |
| ModularityTAD | Graph modularity + DP | ⚠️ калибровка |

---

## Консенсусная логика и цветовая схема

| Поддержка | Цвет | Hex |
|-----------|------|-----|
| 1 алгоритм | — (не консенсус) | — |
| 2 алгоритма | 🟡 Жёлтый | #FFD700 |
| 3 алгоритма | 🟠 Оранжевый | #FF8C00 |
| ≥4 алгоритма | 🟢 Зелёный | #00C800 |

---

## Требования

- OS: Ubuntu 22.04
- Python: 3.10+
- GPU: NVIDIA A100 (опционально, для scKTLD)
- RAM: ≥32 GB рекомендуется

---

## Установка

```bash
# 1. Клонировать репозиторий
git clone <repo_url> && cd tad_project

# 2. Создать conda-окружение
conda env create -f environment.yml
conda activate tad_pipeline

# 3. Доустановить hmmlearn
pip install hmmlearn

# 4. Собрать Armatus (детали в rules.md §4.2)
cd tools/armatus && mkdir -p build && cd build
cmake .. -DCMAKE_POLICY_VERSION_MINIMUM=3.5 ...
make -j2

# 5. Собрать OnTAD (детали в rules.md §4.2)
cd tools/
git clone https://github.com/anlin00007/OnTAD.git ontad_src
cd ontad_src/src
# Применить патч straw.cpp (см. rules.md §4.2)
g++ -std=c++11 -I${CONDA_PREFIX}/include -L${CONDA_PREFIX}/lib \
    -Wl,-rpath,${CONDA_PREFIX}/lib \
    main.cpp step1.cpp step2.cpp step3.cpp step4.cpp common.cpp straw.cpp \
    -lm -lcurl -lz -o OnTAD
mkdir -p ~/tad_project/tools/ontad
cp OnTAD ~/tad_project/tools/ontad/OnTAD
```

---

## Подготовка данных

```bash
# Скачать данные GEO + CTCF
bash download_data.sh

# Структура после скачивания:
# data/raw/GM12878_primary/<res>kb_resolution_intrachromosomal/...
# data/reference/GSE63525_GM12878_primary+replicate_Arrowhead_domainlist.txt
# data/reference/GM12878_CTCF_peaks_hg19.bed

# Подготовить матрицы
python pipeline/run_pipeline.py --prep-data --resolution 25000
```

---

## Запуск

```bash
# Переменные окружения (обязательно на сервере brain-lab)
export CUDA_VISIBLE_DEVICES=3   # GPU 3: ~40 GB свободно
export OMP_NUM_THREADS=2
export MKL_NUM_THREADS=2

# Полный пайплайн, 4 классических алгоритма
python pipeline/run_pipeline.py \
    --resolution 25000 \
    --chroms chr17 chr18 chr19 chr20 chr21 chr22 \
    --algorithms armatus topdom scktld coitad

# Новые алгоритмы (DI+HMM готов; OnTAD и ModularityTAD в калибровке)
python pipeline/run_pipeline.py \
    --resolution 25000 \
    --chroms chr17 chr18 chr19 chr20 chr21 chr22 \
    --algorithms dihmm ontad modularity_tad \
    --force

# Все 7 алгоритмов
python pipeline/run_pipeline.py \
    --resolution 25000 \
    --chroms chr17 chr18 chr19 chr20 chr21 chr22 \
    --algorithms armatus topdom scktld coitad dihmm ontad modularity_tad \
    --force

# Только визуализация (TAD-листы уже есть)
python pipeline/run_pipeline.py \
    --resolution 50000 --chroms chr21 --only-viz

# Без CTCF и отчёта (быстро)
python pipeline/run_pipeline.py \
    --resolution 100000 --chroms chr1 chr2 \
    --algorithms armatus topdom coitad \
    --skip-validation --skip-report
```

---

## Структура результатов

```
results/
├── tads/       # <algo>_<chrom>_<res>bp.bed
├── consensus/  # consensus_<chrom>_<res>bp.bed
├── stats/      # basic.csv, pairwise_jaccard.csv,
│               # rao_comparison.csv, ctcf_enrichment.csv
├── figures/    # hic_tads_<chrom>_<res>bp.png
│               # ctcf_profile_<algo>_<chrom>_<res>bp.png
│               # jaccard_<chrom>_<res>bp.png
└── report.html
```

---

## Известные проблемы

| Проблема | Статус | Описание |
|----------|--------|----------|
| OnTAD мало TADs | ⚠️ | depth=1 даёт 10–16 TADs вместо 50–100; нужен max-nonoverlapping набор |
| ModularityTAD много TADs | ⚠️ | 300+ TADs вместо 50–100; нужна переработка objective |
| CTCF нет пика | 🔍 | Нет чёткого обогащения у 4 классических алгоритмов; нужна диагностика |

Подробности и гипотезы — в `rules.md` §5.7, §5.8, §7.4.

---

## Тесты

```bash
pytest tests/ -v
pytest tests/ --cov=src --cov-report=html
```
```

---

## Следующие шаги (приоритизированы)

```
🔴 P1 — CTCF диагностика:
    head data/reference/GM12878_CTCF_peaks_hg19.bed
    cut -f1 data/reference/GM12878_CTCF_peaks_hg19.bed | sort -u
    # Ожидаем: chr1, chr2, ... (с префиксом chr)
    # Если без префикса → добавить chr-fix в validation.py

🔴 P2 — OnTAD: заменить depth=1 на max-nonoverlapping (листовые TAD)

🔴 P3 — ModularityTAD: заменить mean(B) на intra/inter ratio objective

🟡 P4 — Запустить DI+HMM на chr17–chr22 @ 25kb/50kb/100kb

🟢 P5 — Полный прогон всех 7 алгоритмов после калибровки P2–P3
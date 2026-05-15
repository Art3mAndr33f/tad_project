# TAD Consensus Pipeline

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Ubuntu 22.04](https://img.shields.io/badge/OS-Ubuntu%2022.04-orange.svg)](https://ubuntu.com/)

Сравнительный анализ алгоритмов детекции топологически ассоциированных доменов (TAD)
на bulk Hi-C данных GM12878 (Rao et al. 2014, hg19).

**Научная задача:** показать, что существующий консенсусный подход к разметке
данных для обучения TAD-нейросетей биологически неинформативен, и предложить
пайплайн на основе статистических методов с валидацией через CTCF ChIP-seq.

> **Полный контекст проекта:** `rules.md`
> **Документация по темам:** `agent_docs/` (список файлов → `agent_docs/00_index.md`)
> **Новый участник:** сразу читай раздел [Для новых участников](#для-новых-участников)

---

## Содержание

- [Научный контекст](#научный-контекст)
- [Алгоритмы](#алгоритмы)
- [Данные](#данные)
- [Установка](#установка)
- [Запуск](#запуск)
- [Структура проекта](#структура-проекта)
- [Результаты](#структура-результатов)
- [Параллельная работа](#параллельная-работа--workstreams)
- [Для новых участников](#для-новых-участников)
- [Известные проблемы](#известные-проблемы)
- [Текущий статус](#текущий-статус-исследования)
- [Тесты](#тесты)

---

## Научный контекст

**TAD (Topologically Associating Domains)** — структурные единицы хроматина,
в пределах которых геномные элементы (гены, энхансеры) взаимодействуют
преимущественно друг с другом. Нарушения границ TAD связаны с онкогенезом
и врождёнными заболеваниями.

**CTCF** — ключевой архитектурный белок хроматина. Совместно с когезином
формирует петли ДНК, заякоривающие границы TAD. CTCF ChIP-seq — прямой
биологический эксперимент на тех же клетках (GM12878), что делает его
**биологически информированным ground truth** для оценки качества границ TAD.

**Проблема текущих методов:** Нейросети для детекции TAD (deepTAD и аналоги)
обучаются на консенсусе нескольких алгоритмов. Мы показываем, что этот консенсус
биологически неинформативен: детектированные границы не имеют обогащения CTCF.

**Наш вклад:** Пайплайн статистического отбора границ с CTCF-валидацией
и публикация размеченных Hi-C карт для обучения нейросетей.

---

## Алгоритмы

Семь алгоритмов детекции TAD, покрывающих все основные классы методов:

| Алгоритм | Класс метода | Ключевой механизм | Статья | Статус |
|----------|-------------|-------------------|--------|--------|
| **Armatus** | DP-оптимизация | gamma-регуляризация плотности TAD | Filippova 2014 | ✅ |
| **TopDom** | Локальный | Insulation Score, скользящее окно | Shin 2016 | ✅ |
| **scKTLD** | Спектральный | kNN-граф + eigen-декомпозиция | Zheng 2024 | ✅ |
| **coiTAD** | Локальный | Observed/Insulation ratio | — | ✅ |
| **DI+HMM** | Вероятностный | Directionality Index + GaussianHMM | Dixon 2012 | ✅ |
| **OnTAD** | DP-иерархия | Иерархическая детекция, скольз. среднее | An 2019 | ✅ |
| **ModularityTAD** | Граф | 1D graph modularity + DP | — | ✅ |

**Почему такой набор:** Каждый алгоритм оптимизирует **разное свойство**
контактной матрицы. Это позволяет показать, что разные методы детектируют
принципиально разные границы — и простой консенсус усредняет биологический сигнал.

Подробная механика, параметры и результаты → `agent_docs/01_algorithms.md`

---

## Данные

| Датасет | Источник | Описание |
|---------|---------|----------|
| Hi-C контактные матрицы | GEO [GSE63525](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63525) | GM12878 primary, RAWobserved, 25/50/100kb |
| Эталон границ | Arrowhead (Rao et al. 2014) | `data/reference/` |
| CTCF ChIP-seq | ENCODE [ENCFF796WRU](https://www.encodeproject.org/files/ENCFF796WRU/) | GM12878, hg19 — **биол. ground truth** |

**Геном:** hg19 (GRCh37) — координатная система везде одна.
**Формат матриц:** 3-колоночный TSV (bin_i, bin_j, count), только верхний треугольник.

---

## Установка

### Требования

- OS: Ubuntu 22.04
- Python: 3.10+
- RAM: ≥ 32 GB рекомендуется
- GPU: NVIDIA с CUDA 12.1 (опционально, для scKTLD)

### Шаг 1: conda-окружение

```bash
git clone <repo_url> && cd tad_project
conda env create -f environment.yml
conda activate tad_pipeline
pip install hmmlearn
```

### Шаг 2: torch (опционально, для GPU-ускорения scKTLD)

```bash
pip install torch --index-url https://download.pytorch.org/whl/cu121
python -c "import torch; print(torch.cuda.is_available())"  # → True
```

### Шаг 3: Armatus

```bash
cd tools/armatus && mkdir -p build && cd build
cmake .. -DCMAKE_POLICY_VERSION_MINIMUM=3.5
make -j2
```

### Шаг 4: OnTAD

```bash
conda install -c conda-forge libcurl gcc gxx make -y
cd tools/ && git clone https://github.com/anlin00007/OnTAD.git ontad_src
cd ontad_src/src

# Патч для GCC 15
python3 -c "
text = open('straw.cpp').read()
old = '#include \"straw.h\"'
new = '#include \"straw.h\"\n#include <cstdint>'
if '<cstdint>' not in text:
    text = text.replace(old, new, 1)
    open('straw.cpp', 'w').write(text)
    print('Patch applied')
else:
    print('Already patched')
"

g++ -std=c++11 \
    -I${CONDA_PREFIX}/include \
    -L${CONDA_PREFIX}/lib \
    -Wl,-rpath,${CONDA_PREFIX}/lib \
    main.cpp step1.cpp step2.cpp step3.cpp step4.cpp common.cpp straw.cpp \
    -lm -lcurl -lz -o OnTAD

mkdir -p ~/tad_project/tools/ontad
cp OnTAD ~/tad_project/tools/ontad/OnTAD && chmod +x ~/tad_project/tools/ontad/OnTAD
```

Детальная инструкция с разбором ошибок → `agent_docs/04_build_ontad.md`

---

## Данные — скачивание

```bash
bash download_data.sh

# После скачивания:
# data/raw/GM12878_primary/<res>kb_resolution_intrachromosomal/...
# data/reference/GSE63525_GM12878_primary+replicate_Arrowhead_domainlist.txt
# data/reference/GM12878_CTCF_peaks_hg19.bed

# Подготовить матрицы (.npy кеш)
python pipeline/run_pipeline.py --prep-data --resolution 25000
```

---

## Запуск

```bash
# Обязательно на сервере brain-lab
export CUDA_VISIBLE_DEVICES=3   # GPU 3: ~40 GB свободно
export OMP_NUM_THREADS=2
export MKL_NUM_THREADS=2

# ── Основные режимы ─────────────────────────────────────────────────────────

# Все 7 алгоритмов @ 25kb (chr17–chr22, ~6 часов)
python pipeline/run_pipeline.py \
    --resolution 25000 \
    --chroms chr17 chr18 chr19 chr20 chr21 chr22 \
    --algorithms armatus topdom scktld coitad dihmm ontad modularity_tad \
    --force

# Только новые алгоритмы
python pipeline/run_pipeline.py \
    --resolution 25000 \
    --chroms chr17 chr18 chr19 chr20 chr21 chr22 \
    --algorithms dihmm ontad modularity_tad

# Более широкий прогон @ 100kb (все хромосомы)
python pipeline/run_pipeline.py \
    --resolution 100000 \
    --chroms chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 \
             chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 \
             chr17 chr18 chr19 chr20 chr21 chr22 chrX \
    --algorithms armatus topdom coitad dihmm ontad modularity_tad

# ── Специальные режимы ───────────────────────────────────────────────────────

# Только визуализация (TAD-листы уже посчитаны)
python pipeline/run_pipeline.py \
    --resolution 50000 --chroms chr21 --only-viz

# Без CTCF и отчёта (быстрый прогон)
python pipeline/run_pipeline.py \
    --resolution 100000 --chroms chr1 chr2 \
    --algorithms armatus topdom coitad \
    --skip-validation --skip-report

# ── Запуск в фоне с логом ────────────────────────────────────────────────────

mkdir -p logs
nohup python pipeline/run_pipeline.py \
    --resolution 25000 \
    --chroms chr17 chr18 chr19 chr20 chr21 chr22 \
    --algorithms armatus topdom scktld coitad dihmm ontad modularity_tad \
    --force \
    > logs/pipeline_25kb_$(date +%Y%m%d_%H%M%S).log 2>&1 &
echo "PID: $!"
```

---

## Структура проекта

```
tad_consensus/
├── rules.md                  ← главный контекст (читать первым)
├── README.md
├── CONTRIBUTING.md           ← правила для новых участников
├── create_agent_docs.py      ← скрипт создания agent_docs/
├── .gitignore
├── requirements.txt
├── environment.yml
├── download_data.sh
│
├── agent_docs/               ← подробная дока по темам
│   ├── 00_index.md           ← индекс: что и когда читать
│   ├── 01_algorithms.md      ← механика и параметры алгоритмов
│   ├── 02_data.md            ← форматы данных, пути
│   ├── 03_ctcf_validation.md ← CTCF как биол. ground truth
│   ├── 04_build_ontad.md     ← сборка OnTAD (C++)
│   ├── 05_known_issues.md    ← баги и решения
│   ├── 06_thesis_structure.md← план диплома по главам
│   ├── 07_visualization.md   ← задачи и API визуализации
│   ├── 08_gpu_server.md      ← brain-lab настройка
│   ├── 09_workstreams.md     ← параллельная работа
│   └── 10_statistical_methods.md ← Фурье, спектр, new методы
│
├── config/
│   └── config.yaml           ← ВСЕ параметры
│
├── data/                     ← НЕ в git
│   ├── raw/                  ← RAWobserved файлы
│   ├── processed/            ← .npy матрицы (кеш)
│   └── reference/            ← Arrowhead + CTCF BED
│
├── tools/
│   ├── armatus/              ← бинарник armatus
│   ├── coiTAD/               ← Python-версия coiTAD
│   ├── ontad/                ← бинарник OnTAD
│   └── ontad_src/            ← исходники OnTAD (git clone)
│
├── src/
│   ├── __init__.py
│   ├── data_prep.py          ← загрузка матриц (4-уровневый fallback)
│   ├── algorithms/
│   │   ├── __init__.py       ← ALGORITHM_REGISTRY (7 алгоритмов)
│   │   ├── run_armatus.py
│   │   ├── run_topdom.py
│   │   ├── run_scktld.py
│   │   ├── run_coitad.py
│   │   ├── run_dihmm.py
│   │   ├── run_ontad.py
│   │   └── run_modularity_tad.py
│   ├── consensus.py          ← жадная кластеризация границ
│   ├── statistics.py         ← Jaccard, boundary overlap, Arrowhead
│   ├── validation.py         ← CTCF ChIP-seq обогащение
│   └── visualization.py      ← heatmap, профили, dashboards
│
├── pipeline/
│   └── run_pipeline.py       ← оркестратор
│
├── tests/
│   ├── test_consensus.py
│   └── test_statistics.py
│
├── notebooks/
│   └── exploration.ipynb
│
├── thesis/                   ← текст диплома
│   ├── chapters/
│   └── figures/              ← симлинки из results/figures/
│
├── logs/                     ← логи запусков (НЕ в git)
│
└── results/                  ← НЕ в git
    ├── tads/                 ← <algo>_<chrom>_<res>bp.bed
    ├── consensus/            ← consensus_<chrom>_<res>bp.bed
    ├── stats/                ← basic.csv, pairwise_jaccard.csv, ...
    ├── figures/              ← PNG графики
    └── report.html
```

---

## Структура результатов

```
results/
├── tads/
│   └── <algo>_<chrom>_<res>bp.bed
│       # Пример: armatus_chr21_25000bp.bed
│       # Формат BED3: chrom start end
│
├── consensus/
│   └── consensus_<chrom>_<res>bp.bed
│       # BED3 + колонка support (2/3/4+)
│
├── stats/
│   ├── basic.csv                  # n_tads, median_size_kb, coverage_pct, ...
│   ├── pairwise_jaccard.csv       # попарный Jaccard между алгоритмами
│   ├── rao_comparison.csv         # Recall/Precision/F1 vs Arrowhead
│   └── ctcf_enrichment.csv        # CTCF обогащение по алгоритмам
│
└── figures/
    ├── hic_tads_<chrom>_<res>bp.png           # Hi-C + overlay границ
    ├── ctcf_profile_wide_<chrom>_<res>bp.png  # CTCF профиль ±500kb
    ├── ctcf_profile_<algo>_<chrom>_<res>bp.png
    ├── jaccard_<chrom>_<res>bp.png            # попарный Jaccard heatmap
    └── dashboard_<chrom>_<res>bp.png          # genomic dashboard
```

---

## Консенсусная логика

Жадная кластеризация: границы двух алгоритмов объединяются в кластер
если расстояние между ними ≤ 1 бин (tolerance_bins=1).

| Support | Цвет | Смысл |
|---------|------|-------|
| 1 алгоритм | — | Не консенсус, не включается |
| 2 алгоритма | 🟡 `#FFD700` | Минимальный консенсус |
| 3 алгоритма | 🟠 `#FF8C00` | Средняя уверенность |
| ≥4 алгоритма | 🟢 `#00C800` | Высокая уверенность |

---

## Параллельная работа — Workstreams

Проект разбит на изолированные рабочие потоки с чёткими границами ответственности:

| Workstream | Git-ветка | Зона ответственности | Статус |
|------------|-----------|---------------------|--------|
| **WS-1: Algorithms** | `ws/algorithms` | `src/algorithms/`, `tools/` | 🟡 Active |
| **WS-2: Visualization** | `ws/visualization` | `src/visualization.py`, `results/figures/` | 🔲 Открыт |
| **WS-3: Validation** | `ws/validation` | `src/validation.py`, `src/statistics.py`, `src/consensus.py` | 🔴 Blocked |
| **WS-4: Thesis** | `ws/thesis` | `thesis/`, `notebooks/` | 🔲 Открыт |

**Золотые правила:**
- Работай только в своей зоне ответственности
- `config/config.yaml` — только чтение, изменения через PR
- Интерфейс `run_<algo>()` — неприкосновенный контракт
- Хочешь изменить чужой файл → создай cross-ws issue

**Полная документация** → `agent_docs/09_workstreams.md`

---

## Для новых участников

```bash
# 1. Прочти контекст (5–10 минут)
cat rules.md

# 2. Узнай воркстримы, найди свой
cat agent_docs/09_workstreams.md

# 3. Загрузи детали по своей задаче
cat agent_docs/07_visualization.md   # если работаешь с графиками
cat agent_docs/06_thesis_structure.md  # если пишешь диплом
cat agent_docs/01_algorithms.md      # если работаешь с алгоритмами
cat agent_docs/03_ctcf_validation.md # если работаешь с CTCF

# 4. Создай свою ветку
git checkout main && git pull
git checkout -b ws/<workstream_name>

# 5. Убедись что тесты проходят
pytest tests/ -v
```

**Не знаешь с чего начать?** → `agent_docs/09_workstreams.md` §Текущие задачи

---

## Текущий статус исследования

### 🔴 В работе (блокирует основные результаты)

| Задача | Документация |
|--------|-------------|
| Диагностика CTCF-валидации (coord, chr-prefix, окно) | `agent_docs/03_ctcf_validation.md` |
| Широкий CTCF-профиль ±500kb для всех алгоритмов | `agent_docs/07_visualization.md` |
| Воспроизведение результатов алгоритмов из статей | `agent_docs/01_algorithms.md` |
| Обоснование гиперпараметров по хромосомам | `agent_docs/01_algorithms.md` |

### 🟡 Запланировано (формирует аргументы диплома)

| Задача | Документация |
|--------|-------------|
| Улучшенный overlay границ на Hi-C heatmap | `agent_docs/07_visualization.md` |
| Консенсус по TAD-пересечениям (не по границам) | `agent_docs/10_statistical_methods.md` |
| CTCF как единственный ground truth | `agent_docs/03_ctcf_validation.md` |
| Прогон @ 50kb и 100kb (chr1–chr22) | — |
| Начать написание диплома (главы 3 и 6) | `agent_docs/06_thesis_structure.md` |

### 🟢 Будущее (основной научный вклад)

| Задача | Документация |
|--------|-------------|
| Фурье-анализ и спектральные методы для поиска границ | `agent_docs/10_statistical_methods.md` |
| Доказательство превосходства нового датасета через CTCF | `agent_docs/03_ctcf_validation.md` |
| GNN на консенсусных границах (proof-of-concept) | — |
| Публикация размеченных Hi-C карт | — |

---

## Известные проблемы

| Проблема | Статус | Документация |
|----------|--------|-------------|
| CTCF: нет чёткого пика обогащения | 🔴 Диагностируется | `agent_docs/03_ctcf_validation.md` |
| Гиперпараметры: единый TAD/Mb диапазон для всех хромосом | 🟡 Открыт | `agent_docs/01_algorithms.md` |
| coiTAD: частота fallback на IS неизвестна | 🟡 Открыт | `agent_docs/05_known_issues.md` |
| OnTAD: depth=1 → мало TADs | ✅ Исправлено v2.1 | `agent_docs/05_known_issues.md` |
| ModularityTAD: 300+ TADs | ✅ Исправлено v2.1 | `agent_docs/05_known_issues.md` |
| OnTAD: компиляция (uint64_t, libcurl) | ✅ Исправлено v2.0 | `agent_docs/04_build_ontad.md` |

---

## Тесты

```bash
# Запуск всех тестов
pytest tests/ -v

# С coverage-репортом
pytest tests/ --cov=src --cov-report=html
open htmlcov/index.html

# Только быстрые тесты (без интеграционных)
pytest tests/ -v -k "not Integration"
```

**Тесты обязательны** перед каждым PR — см. `CONTRIBUTING.md`.

---

## Цитирование данных

Если используете данные или результаты из этого репозитория:

```
Rao et al. (2014). A 3D Map of the Human Genome at Kilobase Resolution
Reveals Principles of Chromatin Looping. Cell, 159(7), 1665–1680.
GEO: GSE63525
```

---

## Требования к системе

```
OS:     Ubuntu 22.04
Python: 3.10+
RAM:    ≥ 32 GB (рекомендуется ≥ 64 GB для chr1 @ 25kb)
GPU:    NVIDIA CUDA 12.1+ (опционально, для scKTLD)
Disk:   ≥ 100 GB (данные + результаты)
```

---

*TAD Consensus Pipeline | hg19 | GSE63525 GM12878 | v3.1 | 2026-05-14*
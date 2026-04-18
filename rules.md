# TAD Consensus Pipeline — Project Rules for LLM Context

> Этот файл — единственный источник правды о проекте.
> Кидай его в начало любого нового чата, чтобы модель сразу знала контекст.

---

## 1. Что это за проект

Сравнительный анализ детекции топологически ассоциированных доменов (TAD)
на bulk Hi-C данных клеточной линии **GM12878** (human lymphoblastoid,
эталонный датасет Rao et al. 2014, GEO: **GSE63525**, координатная система
**hg19**).

Пайплайн запускает четыре алгоритма детекции TAD, строит консенсус границ,
оценивает качество через сравнение с эталоном Arrowhead и обогащение CTCF.

**OS**: Ubuntu 22.04  
**Python**: 3.10+  
**Conda env**: `tad_pipeline`

---

## 2. Структура проекта (каноническая)

```
tad_consensus/
├── README.md
├── .gitignore
├── requirements.txt
├── environment.yml
├── download_data.sh          # скрипт скачивания данных (см. раздел 7)
│
├── config/
│   └── config.yaml           # ВСЕ пути и параметры — только здесь
│
├── data/
│   ├── raw/                  # архив + распакованные RAWobserved — НЕ в git
│   ├── processed/            # .npy матрицы + .RAWobserved — НЕ в git
│   └── reference/            # Arrowhead domainlist + CTCF BED — НЕ в git
│
├── tools/
│   ├── armatus/              # бинарник armatus (собран из исходников)
│   └── coiTAD/               # кастомная Python-версия coiTAD
│       ├── __init__.py
│       ├── coitad_core.py    # ядро: build_oi_matrix, detect_tads
│       └── detector.py       # класс CoiTADDetector
│
├── src/
│   ├── __init__.py
│   ├── data_prep.py
│   ├── algorithms/
│   │   ├── __init__.py       # ALGORITHM_REGISTRY
│   │   ├── run_armatus.py
│   │   ├── run_topdom.py
│   │   ├── run_scktld.py
│   │   └── run_coitad.py
│   ├── consensus.py
│   ├── statistics.py
│   ├── validation.py
│   └── visualization.py
│
├── pipeline/
│   └── run_pipeline.py       # CLI точка входа
│
├── tests/
│   ├── test_consensus.py
│   └── test_statistics.py
│
├── notebooks/
│   └── exploration.ipynb
│
└── results/                  # генерируется — НЕ в git
    ├── tads/                 # <algo>_<chrom>_<resolution>bp.bed
    ├── consensus/            # consensus_<chrom>_<resolution>bp.bed
    ├── stats/                # *.csv
    ├── figures/              # *.png / *.html
    └── report.html
```

---

## 3. Данные

### 3.1 Источники и файлы

| Файл | Путь | Описание |
|------|------|----------|
| Интрахромосомные матрицы | `data/raw/` | RAWobserved по хромосомам и разрешениям |
| Arrowhead domain list | `data/reference/GSE63525_GM12878_primary+replicate_Arrowhead_domainlist.txt` | Эталон Rao 2014 |
| CTCF ChIP-seq | `data/reference/GM12878_CTCF_peaks_hg19.bed` | ENCODE ENCFF796WRU, hg19 |

⚠️ .hic файл НЕ используется и НЕ скачивается.
Источник данных — RAWobserved файлы из data/raw/:
data/raw/GM12878_primary/<res>kb_resolution_intrachromosomal/<chrom>/MAPQGE30/<chrom>_<res>kb.RAWobserved

### 3.2 Скрипт скачивания (`download_data.sh`)

```bash
#!/bin/bash
set -e
mkdir -p data/raw data/reference

BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl"

# Интрахромосомные матрицы (primary, ~7.1 Gb)
wget -c -P data/raw/ \
    "${BASE}/GSE63525_GM12878_primary_intrachromosomal_contact_matrices.tar.gz"
tar -xzf data/raw/GSE63525_GM12878_primary_intrachromosomal_contact_matrices.tar.gz \
    -C data/raw/

# Arrowhead domain list (эталон, 242 Kb)
wget -c -P data/reference/ \
    "${BASE}/GSE63525_GM12878_primary+replicate_Arrowhead_domainlist.txt.gz"
gunzip data/reference/GSE63525_GM12878_primary+replicate_Arrowhead_domainlist.txt.gz

# CTCF ChIP-seq (ENCODE)
wget -c -O data/reference/GM12878_CTCF_peaks_hg19.bed.gz \
    "https://www.encodeproject.org/files/ENCFF796WRU/@@download/ENCFF796WRU.bed.gz"
gunzip data/reference/GM12878_CTCF_peaks_hg19.bed.gz
```

После скачивания структура `data/raw/` содержит директории вида:
`GM12878/primary/10kb_resolution_intrachromosomal/chr1/MAPQGE30/chr1_10kb.RAWobserved`

### 3.3 Разрешения и хромосомы

**Разрешения**: `[25000, 50000, 100000]` (в bp)

**Хромосомы по умолчанию**: `chr1–chr22, chrX`

**Ограничения памяти для scKTLD** (dense-матрицы):

| Разрешение | Допустимые хромосомы |
|-----------|----------------------|
| 10 kb | chr21, chr22 |
| 25 kb | chr17–chr22 |
| 50 kb | chr1–chr22 |
| 100 kb | chr1–chr22 |

### 3.4 Формат RAWobserved

3-колоночный TSV без заголовка: `bin_i<TAB>bin_j<TAB>count`  
Координаты в bp (не индексы бинов).  
Только верхний треугольник (включая диагональ), нули пропущены.

---

## 4. Технический стек

### 4.1 Зависимости

```
# Hi-C / геномика
hic-straw>=1.3.1      # чтение .hic файлов
cooler>=0.9.3
pybedtools>=0.9.1
pyranges>=0.0.129

# Вычисления
numpy>=1.24
scipy>=1.11
pandas>=2.0
scikit-learn>=1.3

# Визуализация
matplotlib>=3.7
seaborn>=0.12
plotly>=5.18

# Утилиты
PyYAML>=6.0
tqdm>=4.66
jinja2>=3.1
click>=8.1

# Тесты
pytest>=7.4
pytest-cov>=4.1
```

### 4.2 Внешние инструменты

- **Armatus**: C++-бинарник в `tools/armatus/armatus`.
  Сборка через CMake (не Make):
    sudo apt-get install libboost-dev libboost-program-options-dev \
                         libboost-iostreams-dev zlib1g-dev
    mkdir build && cd build
    cmake .. -DCMAKE_EXE_LINKER_FLAGS="-lboost_iostreams -lz"
    make -j$(nproc)
  Бинарник: `build/src/armatus`
- **coiTAD**: кастомная Python-версия в `tools/coiTAD/`.
  Импортируется динамически через `sys.path`.

---

## 5. Алгоритмы

### 5.1 Единый интерфейс (ОБЯЗАТЕЛЬНО)

Каждый алгоритм реализует функцию:

```python
def run_<algorithm>(
    chrom: str,           # 'chr1', 'chrX', ...
    resolution: int,      # в bp: 10000 / 25000 / 50000 / 100000
    data_path: str,       # путь к data/processed/
    cfg: Optional[dict],  # конфиг-словарь (может быть None)
    **kwargs,
) -> pd.DataFrame:
    # Возвращает DataFrame с колонками: chrom, start, end
    # start и end в базовых парах (hg19)
    # При ошибке или пропуске: pd.DataFrame(columns=['chrom','start','end'])
```

### 5.2 Armatus (`src/algorithms/run_armatus.py`)

- Вызов через `subprocess.run()`
- Флаги: -S (sparse 3-column TSV), -N (без нормализации), -c <chrom>
- Флаг -R — для директории Rao, НЕ использовать
- Флаг -z — не существует, НЕ использовать
- Флаг -j — меняет логику на "только gamma_max", НЕ использовать при переборе gamma

- Выходной файл: <prefix>.consensus.txt  (не .txt и не _level_0.txt)
candidates = [
    f"{prefix}.consensus.txt",   # ← первый приоритет
    f"{prefix}.txt",
    f"{prefix}_level_0.txt",
    f"{prefix}_domain.txt",
]
- Перебор gamma: `[0.1, 0.5, 1.0, 2.0, 5.0]`
- Выбор лучшей gamma: наибольшая средняя стабильность
  (средний Jaccard на уровне доменов со всеми остальными gamma)
- Парсинг выхода: файл `<prefix>.txt` или `<prefix>_level_0.txt`
- Timeout: 600 секунд на один запуск

### 5.3 TopDom (`src/algorithms/run_topdom.py`)

- Чистый Python, без внешних зависимостей
- Входные данные: dense numpy-матрица (загружается из `.npy`)
- Перебор window: `[3, 5, 10]`
- Выбор лучшего окна: максимальный intra/inter ratio
- Bin-score[i] = mean(upstream_triangle) - mean(downstream_triangle)
- Границы = локальные минимумы bin-score ниже mean−std

### 5.4 scKTLD (`src/algorithms/run_scktld.py`)

- **КРИТИЧНО**: принимает dense numpy n×n матрицу **без нормализации**
  (RAW counts, `balance=False`) — требование авторов алгоритма
- Алгоритм: RBF-ядро → нормализованный Лапласиан → eigsh →
  спектральное вложение → DP-сегментация
- `dimension=32`
- `knn_k=20`
- `penalty`: автоподбор через elbow на логарифмической сетке
- **Перед запуском**: обязательная проверка `scktld_limits[resolution]`
- При нарушении лимита — возвращать пустой DataFrame, НЕ падать

### 5.5 coiTAD (`src/algorithms/run_coitad.py`)

- Импорт из `tools/coiTAD/` через динамический `sys.path`
- Точки входа (в порядке приоритета):
  1. `CoiTADDetector().detect(matrix, resolution=resolution)`
  2. `detect_tads(matrix, resolution=resolution)`
  3. `run(matrix, resolution=resolution)`
- **Fallback**: если `tools/coiTAD/` недоступен → встроенная реализация
  на основе OI-матрицы + Directional Index
- Авторские модификации из `tools/coiTAD/` — СОХРАНЯТЬ

### 5.6 ALGORITHM_REGISTRY (`src/algorithms/__init__.py`)

```python
ALGORITHM_REGISTRY = {
    "armatus": run_armatus,
    "topdom":  run_topdom,
    "scktld":  run_scktld,
    "coitad":  run_coitad,
}
```

---

## 6. Консенсус границ

### 6.1 Алгоритм

1. Для каждой хромосомы × разрешение собрать все границы TAD
   (start и end каждого домена) от всех алгоритмов
2. Округлить позиции до ближайшего бина: `round(pos / resolution) * resolution`
3. Жадная кластеризация: граница входит в кластер если расстояние
   до текущего центра кластера ≤ `tolerance_bins * resolution`
4. Центр кластера = медиана всех позиций в кластере
5. Support = число алгоритмов, у которых хотя бы одна граница
   попадает в кластер

### 6.2 Параметры (из конфига)

```yaml
consensus:
  tolerance_bins: 1    # ±1 bin
  min_support: 2       # минимум для консенсуса
```

### 6.3 Цветовая схема (ФИКСИРОВАНА, не менять)

| Support | Цвет | Hex | Значение |
|---------|------|-----|----------|
| 1 | — | — | не консенсус |
| 2 | 🟡 Жёлтый | `#FFD700` | слабый |
| 3 | 🟠 Оранжевый | `#FF8C00` | умеренный |
| 4 | 🟢 Зелёный | `#00C800` | сильный |

```python
# Канонический словарь — не изменять
CONSENSUS_COLORS = {
    2: "#FFD700",
    3: "#FF8C00",
    4: "#00C800",
}
```

### 6.4 Формат BED-файла консенсуса

```
chrom  start  end  name                  score  strand
chr17  500000 525000 consensus_support3  3      .
```

`start = position`, `end = position + resolution`, `score = support`

---

## 7. Метрики и статистика

### 7.1 Базовая статистика (на каждый ал��оритм × хромосому × разрешение)

| Метрика | Описание |
|---------|----------|
| `n_tads` | Число TAD |
| `median_size_kb` | Медианный размер TAD (kb) |
| `mean_size_kb` | Средний размер TAD (kb) |
| `std_size_kb` | Стандартное отклонение размера |
| `coverage_pct` | % покрытия хромосомы |
| `n_consensus_bnd` | Число границ, поддержанных ≥2 алгоритмами |
| `pct_consensus_bnd` | % консенсусных границ от всех границ алгоритма |

### 7.2 Попарное сравнение

- **Jaccard index** на уровне доменов (точное совпадение start + end)
- **Boundary overlap rate**: % границ алгоритма A, найденных алгоритмом B
  в пределах ±1 bin

### 7.3 Сравнение с эталоном Rao 2014

Эталон: `data/reference/GSE63525_GM12878_primary+replicate_Arrowhead_domainlist.txt`

Метрики вычисляются при допусках **±1 bin** и **±2 bin**:

| Метрика | Формула |
|---------|---------|
| Boundary Recall | # ref-границ найденных алгоритмом / # всех ref-границ |
| Boundary Precision | # pred-границ совпавших с ref / # всех pred-границ |
| F1-score | 2 × P × R / (P + R) |
| Jaccard (домены) | \|A ∩ B\| / \|A ∪ B\| |

### 7.4 CTCF-валидация

- Источник: `data/reference/GM12878_CTCF_peaks_hg19.bed`
- Окно вокруг границы: `±1 bin * resolution`
- **1000 пермутаций** случайных окон (seed=42)
- Метрики: `enrichment_score = observed / mean(random)`, `z_score`, `p_value`
- Профиль обогащения: ±500 kb от границы, бин 10 kb

---

## 8. Конфигурация

**Всё** хранится в `config/config.yaml`. Пути и параметры **не хардкодить** в коде.

Ключевые поля:

```yaml
project:
  seed: 42             # везде где есть случайность

paths:
  raw_data:       "data/raw"
  processed:      "data/processed"
  reference:      "data/reference"
  results:        "results"
  tads_out:       "results/tads"
  consensus_out:  "results/consensus"
  stats_out:      "results/stats"
  figures_out:    "results/figures"
  report:         "results/report.html"
  armatus_bin:    "tools/armatus/armatus"
  coitad_dir:     "tools/coiTAD"
  hic_file:       ""
  arrowhead_file: "data/reference/GSE63525_GM12878_primary+replicate_Arrowhead_domainlist.txt"
  ctcf_bed:       "data/reference/GM12878_CTCF_peaks_hg19.bed"

resolutions:      [25000, 50000, 100000]

chromosomes:
  all: ["chr1".."chr22", "chrX"]
  scktld_limits:
    10000:  ["chr21","chr22"]
    25000:  ["chr17","chr18","chr19","chr20","chr21","chr22"]
    50000:  ["chr1".."chr22"]
    100000: ["chr1".."chr22"]

algorithms:
  armatus:
    gamma_values:    [0.1, 0.5, 1.0, 2.0, 5.0]
    stability_top_n: 3
  topdom:
    window_sizes: [3, 5, 10]
    default_window: 5
  scktld:
    dimension: 32
    knn_k: 20
    balance: false
    auto_penalty: true
  coitad:
    default_params: {}

consensus:
  tolerance_bins: 1
  min_support:    2
  colors:
    2: "#FFD700"
    3: "#FF8C00"
    4: "#00C800"

validation:
  ctcf_window_bp:              1
  n_permutations:              1000
  enrichment_profile_range_bp: 500000
  enrichment_profile_bin_bp:   10000

rao_comparison:
  tolerances_bins: [1, 2]

visualization:
  hic_colormap:  "coolwarm"
  hic_log_scale: true
  dpi:           300
  format:        "png"
  figsize:       [18, 14]
  generate_html: true
  algorithm_colors:
    armatus: "#1f77b4"
    topdom:  "#ff7f0e"
    scktld:  "#2ca02c"
    coitad:  "#d62728"
```

---

## 9. Соглашения по именованию

### 9.1 Файлы результатов

| Тип | Шаблон | Пример |
|-----|--------|--------|
| TAD-лист | `results/tads/<algo>_<chrom>_<res>bp.bed` | `armatus_chr17_25000bp.bed` |
| Консенсус | `results/consensus/consensus_<chrom>_<res>bp.bed` | `consensus_chr17_25000bp.bed` |
| Статистика | `results/stats/<name>.csv` | `basic.csv`, `pairwise_jaccard.csv` |
| Фигура PNG | `results/figures/hic_tads_<chrom>_<res>bp.png` | `hic_tads_chr17_25000bp.png` |
| Фигура HTML | `results/figures/hic_tads_<chrom>_<res>bp.html` | |
| CTCF-профиль | `results/figures/ctcf_profile_<algo>_<chrom>_<res>bp.png` | |
| Jaccard heatmap | `results/figures/jaccard_<chrom>_<res>bp.png` | |

### 9.2 Матрицы в `data/processed/`

| Тип | Шаблон | Пример |
|-----|--------|--------|
| Dense numpy | `<chrom>_<res>bp.npy` | `chr17_25000bp.npy` |
| RAWobserved | `<chrom>_<res>bp.RAWobserved` | `chr17_25000bp.RAWobserved` |

### 9.3 Переменные и функции

```python
# Хромосомы — всегда с префиксом 'chr'
chrom = "chr17"          # ✅
chrom = "17"             # ❌

# Разрешение — всегда в bp (int)
resolution = 25000       # ✅
resolution = "25kb"      # ❌ (только в CLI-аргументах и именах файлов)

# Позиции — всегда в bp (int или np.int64)
start = 500_000          # ✅

# Индексы бинов — только внутри функций, не возвращать наружу
bin_i = start // resolution   # только внутри

# DataFrame колонки TAD
df.columns == ["chrom", "start", "end"]   # ✅ строго этот порядок

# Пустой DataFrame при ошибке/пропуске
pd.DataFrame(columns=["chrom", "start", "end"])   # ✅
```

### 9.4 Логирование

```python
import logging
logger = logging.getLogger(__name__)   # ✅ всегда __name__

# Уровни:
logger.debug(...)    # детали вычислений, пути файлов
logger.info(...)     # ключевые шаги, результаты (число TAD, метрики)
logger.warning(...)  # пропуск из-за ограничений, fallback
logger.error(...)    # ошибки, не останавливающие пайплайн
# raise — только для критических ошибок конфигурации
```

---

## 10. Правила и запреты

### 10.1 ЗАПРЕЩЕНО

```python
# ❌ Хардкодить пути
matrix = np.load("data/processed/chr17_25000bp.npy")

# ❌ Хардкодить параметры алгоритмов
gamma_values = [0.1, 0.5, 1.0, 2.0, 5.0]

# ❌ Нормализовать матрицу для scKTLD
matrix = matrix / matrix.max()   # scKTLD требует RAW counts

# ❌ Запускать scKTLD без проверки памяти
run_scktld("chr1", 10000, ...)   # OOM на ~3GB матрице

# ❌ Менять CONSENSUS_COLORS
CONSENSUS_COLORS[2] = "#FF0000"  # цветовая схема фиксирована

# ❌ Возвращать None из run_<algorithm>
return None   # всегда pd.DataFrame(columns=["chrom","start","end"])

# ❌ print() в src/ и pipeline/
print("debug")   # только logger.*

# ❌ Менять seed
np.random.seed(123)   # везде seed=42 из конфига

# ❌ Использовать абсолютные пути
"/home/user/tad_consensus/data/..."   # только относительные от корня проекта

# ❌ Коммитить данные и результаты
git add data/      # ❌
git add results/   # ❌
```

### 10.2 ОБЯЗАТЕЛЬНО

```python
# ✅ Загружать конфиг из файла
cfg = load_config("config/config.yaml")

# ✅ Читать параметры из конфига
gamma_values = cfg["algorithms"]["armatus"]["gamma_values"]

# ✅ Проверять лимиты памяти scKTLD до запуска
limits = cfg["chromosomes"]["scktld_limits"].get(resolution, [])
if limits and chrom not in limits:
    return pd.DataFrame(columns=["chrom", "start", "end"])

# ✅ Создавать директории перед записью
Path(out_path).parent.mkdir(parents=True, exist_ok=True)

# ✅ seed=42 для всех генераторов случайных чисел
rng = np.random.default_rng(42)

# ✅ balance=False для scKTLD (RAW counts, без нормализации)
matrix = get_matrix(cfg, chrom, resolution)   # уже RAW из hicstraw NONE

# ✅ Все хромосомы с префиксом 'chr' в DataFrame
df["chrom"] = "chr17"

# ✅ Консенсус считается отдельно для каждого разрешения
for res in resolutions:
    consensus = compute_consensus(algo_results, chrom, res)

# ✅ Сравнение с Rao 2014 — строго в координатах hg19
```

### 10.3 Обработка ошибок

```python
# ✅ Паттерн для алгоритмов: ловить все, логировать, возвращать пустой DF
try:
    df = run_algorithm(...)
except Exception as exc:
    logger.error("[algo] Ошибка %s @ %d: %s", chrom, resolution, exc, exc_info=True)
    df = pd.DataFrame(columns=["chrom", "start", "end"])

# ✅ coiTAD: try original → except → fallback
try:
    _core = _import_coitad(cfg)
    ...
except Exception as exc:
    logger.warning("[coiTAD] Использую fallback: %s", exc)
    df = _coitad_fallback(matrix, chrom, resolution)

# ✅ Armatus subprocess: проверять returncode, обрабатывать TimeoutExpired
proc = subprocess.run(cmd, timeout=600, ...)
if proc.returncode != 0:
    logger.warning("Armatus ненулевой код: %d", proc.returncode)
    return _empty_df()
```

---

## 11. data_prep.py — ключевые функции

```python
# Загрузить конфиг
cfg = load_config("config/config.yaml") -> dict

# Извлечь dense-матрицу из .hic (через hicstraw, normalization="NONE")
matrix = extract_dense_matrix(hic_path, chrom, resolution, normalization="NONE")
# → np.ndarray shape (n, n), dtype float32, симметричная

# Получить матрицу из кэша (или извлечь если нет)
matrix = get_matrix(cfg, chrom, resolution) -> np.ndarray

# Получить путь к RAWobserved (или создать если нет)
raw_path = get_rawobserved_path(cfg, chrom, resolution) -> Path

# Подготовить все матрицы (batch)
prepare_all_matrices(cfg, resolutions, chromosomes, force=False)

get_matrix теперь имеет 4-уровневый fallback:
  1. data/processed/<chrom>_<res>bp.npy           (кэш)
  2. data/processed/<chrom>_<res>bp.RAWobserved   (кэш)
  3. data/raw/GM12878_primary/<res>kb_.../MAPQGE30/<chrom>_<res>kb.RAWobserved
  4. .hic через hicstraw (только если файл существует)

При сохранении кэша сохраняются ОБА файла: .npy и .RAWobserved

Новая функция:
get_rawobserved_path_raw(cfg, chrom, resolution) -> Path
  Возвращает путь к RAWobserved прямо в data/raw/ (без processed/)
```

**Кэш**: `data/processed/<chrom>_<res>bp.npy` и `<chrom>_<res>bp.RAWobserved`

---

## 12. CLI пайплайна

```bash
# Полный пайплайн
python pipeline/run_pipeline.py \
    --resolution 25000 \
    --chroms chr17 chr18 chr19 chr20 chr21 chr22 \
    --algorithms armatus topdom scktld coitad

# Только визуализация (TAD-листы уже есть в results/tads/)
python pipeline/run_pipeline.py --resolution 50000 --chroms chr17 --only-viz

# Пропустить тяжёлые шаги
python pipeline/run_pipeline.py \
    --resolution 100000 --chroms chr1 chr2 \
    --skip-validation --skip-report

# Перезаписать существующие результаты
python pipeline/run_pipeline.py --resolution 50000 --force

# Доступные флаги:
# --config PATH         путь к config.yaml (default: config/config.yaml)
# --resolution INT+     разрешения в bp
# --chroms STR+         хромосомы
# --algorithms STR+     armatus topdom scktld coitad
# --prep-data           только матрицы
# --only-viz            только визуализация
# --skip-stats          пропустить статистику
# --skip-validation     пропустить CTCF
# --skip-report         не генерировать HTML
# --force               перезаписать
# --log-level           DEBUG/INFO/WARNING/ERROR
```

---

## 13. Тестирование

```bash
# Запустить все тесты
pytest tests/ -v

# С покрытием
pytest tests/ --cov=src --cov-report=html

# Только быстрые тесты (без интеграционных)
pytest tests/ -v -k "not Integration"
```

**Что покрывают тесты**:

- `tests/test_consensus.py`: `extract_boundaries`, `cluster_boundaries`,
  `compute_consensus`, цветовая схема, tolerance matching, edge cases
- `tests/test_statistics.py`: `compute_basic_stats`, `jaccard_domains`,
  `boundary_overlap_rate`, `compute_pairwise_matrix`, `compare_with_reference`,
  `HG19_CHROM_SIZES`

**При добавлении новой функции** — добавить unit-тест.

---

## 14. Размеры хромосом hg19 (справочник)

```python
HG19_CHROM_SIZES = {
    "chr1":  249_250_621, "chr2":  243_199_373, "chr3":  198_022_430,
    "chr4":  191_154_276, "chr5":  180_915_260, "chr6":  171_115_067,
    "chr7":  159_138_663, "chr8":  146_364_022, "chr9":  141_213_431,
    "chr10": 135_534_747, "chr11": 135_006_516, "chr12": 133_851_895,
    "chr13": 115_169_878, "chr14": 107_349_540, "chr15": 102_531_392,
    "chr16":  90_354_753, "chr17":  81_195_210, "chr18":  78_077_248,
    "chr19":  59_128_983, "chr20":  63_025_520, "chr21":  48_129_895,
    "chr22":  51_304_566, "chrX":  155_270_560,
}
```

---

## 15. Визуализация

### Что генерируется для каждой хромосомы × разрешение

1. **PNG** (`hic_tads_<chrom>_<res>bp.png`, 300 dpi):
   - Hi-C карта: log1p(count), colormap=`coolwarm`
   - TAD-домены: цветные треугольники (4 строки = 4 алгоритма)
   - Консенсусные границы: вертикальные линии с цветом по схеме
   - Легенда: алгоритм → цвет, support → цвет

2. **HTML** (`hic_tads_<chrom>_<res>bp.html`):
   - Интерактивная plotly-версия (при `generate_html: true`)

3. **CTCF-профиль** (`ctcf_profile_<algo>_<chrom>_<res>bp.png`):
   - Линейный график ±500 kb от границы, бин 10 kb

4. **Jaccard heatmap** (`jaccard_<chrom>_<res>bp.png`):
   - Seaborn heatmap попарного Jaccard

### Цвета алгоритмов (фиксированы)

```python
ALGO_COLORS = {
    "armatus": "#1f77b4",   # синий
    "topdom":  "#ff7f0e",   # оранжевый
    "scktld":  "#2ca02c",   # зелёный
    "coitad":  "#d62728",   # красный
}
```

---

## 16. Типичные задачи и куда смотреть

| Задача | Файл |
|--------|------|
| Добавить новый алгоритм | `src/algorithms/run_<new>.py` + `ALGORITHM_REGISTRY` |
| Изменить параметры | `config/config.yaml` |
| Изменить логику консенсуса | `src/consensus.py` → `compute_consensus()` |
| Добавить метрику | `src/statistics.py` → `compute_basic_stats()` / `compare_with_reference()` |
| Изменить визуализацию | `src/visualization.py` → `plot_tad_comparison()` |
| Изменить CTCF-анализ | `src/validation.py` → `compute_ctcf_enrichment()` |
| Добавить шаг в пайплайн | `pipeline/run_pipeline.py` → функция + вызов в `main()` |
| Скачать данные заново | `bash download_data.sh` |
| Извлечь матрицы из .hic | `python pipeline/run_pipeline.py --prep-data` |

---

## 17. Частые подводные камни

| Проблема | Причина | Решение |
|----------|---------|---------|
| OOM при scKTLD | Запуск на крупных хромосомах | Проверить `scktld_limits` в конфиге |
| Armatus не находит output | Разные версии пишут разные имена файлов | Проверять несколько кандидатов: `.txt`, `_level_0.txt`, `_domain.txt` |
| hicstraw не видит хромосому | Разные форматы: `chr1` vs `1` | `_chrom_strip()` убирает префикс `chr` для hicstraw |
| NaN в матрице | Нули заменены NaN при log | `np.nan_to_num(matrix)` перед вычислениями, кроме scKTLD |
| Несовпадение координат | Смешивание hg19/hg38 | Все данные строго hg19, проверять при новых источниках |
| coiTAD не импортируется | `tools/coiTAD/` отсутствует или битый | Автоматический fallback на OI+DI реализацию |
| Разные `chr`-форматы в BED | ENCODE иногда даёт `1` вместо `chr1` | `if not df["chrom"].str.startswith("chr"): df["chrom"] = "chr" + df["chrom"]` |
| `No module named 'algorithms'` | `pipeline/` в sys.path вместо корня | Добавить `_PROJECT_ROOT` в sys.path в начале `run_pipeline.py` |
| Armatus 0 TADs, нет ошибки сборки | Флаг `-R` ожидает директорию, не файл | Использовать `-S -N -c <chrom>` |
| Armatus output не найден | Файл пишется как `.consensus.txt` | Добавить в candidates первым |
| scKTLD зависает на часы | `M@M.T` на dense матрице O(n³) + eigsh k=128 | KNN-sparse + dimension=32 |
| scKTLD 0–14 TADs | penalty_grid слишком широкий, elbow выбирает край | Уменьшить диапазон, добавить min_tads фильтр |
| coiTAD 300–700 TADs | `tools/coiTAD/` отсутствует, fallback пересегментирует | Исправить fallback или создать tools/coiTAD/ |
| Arrowhead / CTCF не найдены | `data/reference/` пуста | Скачать отдельно (не входит в основной архив) |
| `No module named 'src.visualization'` | `src/visualization.py` не создан | Создать файл |

---

## 18. Как обновлять этот файл

> rules.md — живой документ. При любом архитектурном изменении проекта
> он обновляется **до** изменения кода (docs-first подход).

### 18.1 Версионирование

Строка версии в последней строке файла:

```
*Версия rules.md: <MAJOR>.<MINOR> | ...*
```

| Тип изменения | Что менять |
|---------------|------------|
| Новый алгоритм, новый источник данных, новая секция | MAJOR + 1 |
| Новый параметр, новая метрика, правка запретов | MINOR + 1 |
| Опечатки, форматирование | версия не меняется |

---

### 18.2 Чеклист: добавить новый алгоритм

**Шаг 1 — Обновить rules.md** (этот файл):

- [ ] Секция **5** («Алгоритмы»): добавить подраздел `5.N <AlgoName>`
  с описанием входа, выхода, параметров, особенностей
- [ ] Секция **5.6** (`ALGORITHM_REGISTRY`): добавить строку
  `"<algo>": run_<algo>`
- [ ] Секция **8** (`config.yaml`): добавить блок под `algorithms:`
- [ ] Секция **9.1** (Именование файлов): убедиться что шаблон
  `<algo>_<chrom>_<res>bp.bed` подходит
- [ ] Секция **15** (Цвета алгоритмов `ALGO_COLORS`): добавить цвет
- [ ] Секция **19** (Changelog): добавить запись

**Шаг 2 — Обновить код**:

- [ ] Создать `src/algorithms/run_<algo>.py` с функцией
  `run_<algo>(chrom, resolution, data_path, cfg, **kwargs) -> pd.DataFrame`
- [ ] Зарегистрировать в `src/algorithms/__init__.py`
- [ ] Добавить параметры в `config/config.yaml`
- [ ] Добавить цвет в `ALGO_COLORS` (`src/visualization.py`)
- [ ] Добавить в `requirements.txt` / `environment.yml` если нужны
  новые зависимости
- [ ] Написать тесты в `tests/test_<algo>.py`

**Шаблон подраздела для секции 5**:

```markdown
### 5.N <AlgoName> (`src/algorithms/run_<algo>.py`)

- **Входные данные**: [dense numpy matrix / RAWobserved / другое]
- **Нормализация**: [RAW / KR / VC — указать явно]
- **Параметры** (из конфига):
  - `param1`: описание, диапазон значений
  - `param2`: описание, диапазон значений
- **Выбор параметров**: [как подбирается оптимальное значение]
- **Ограничения памяти**: [если есть — добавить в `scktld_limits` или
  аналогичную структуру в конфиге]
- **Внешние зависимости**: [subprocess / pip-пакет / нет]
- **Fallback при ошибке**: [описание]
- **Особые требования**: [если есть]
```

**Шаблон блока в `config.yaml`**:

```yaml
algorithms:
  <algo>:
    param1: default_value
    param2: default_value
    # ограничения по хромосомам (если нужны)
    chrom_limits:
      10000:  ["chr21", "chr22"]
      25000:  ["chr17", "chr18", "chr19", "chr20", "chr21", "chr22"]
      50000:  []   # пусто = все хромосомы
      100000: []
```

---

### 18.3 Чеклист: добавить новое разрешение

- [ ] Секция **3.3**: добавить в таблицу «Разрешения»
- [ ] Секция **3.3**: обновить `scktld_limits` (и аналоги для новых
  алгоритмов с ограничениями памяти)
- [ ] Секция **8**: добавить в `resolutions:` и в `scktld_limits`
- [ ] Секция **19**: добавить запись в Changelog
- [ ] Обновить `config/config.yaml`
- [ ] Проверить `prepare_all_matrices()` — новое разрешение будет
  обработано автоматически если оно в `cfg["resolutions"]`

---

### 18.4 Чеклист: добавить новый источник данных

- [ ] Секция **3.1** (таблица файлов): добавить строку
- [ ] Секция **3.2** (`download_data.sh`): добавить команду wget/curl
- [ ] Секция **8** (`paths:`): добавить ключ пути
- [ ] Секция **9.2** (именование в `data/processed/`): если новый
  формат — описать шаблон
- [ ] Секция **19**: добавить запись
- [ ] Обновить `download_data.sh`
- [ ] Обновить `config/config.yaml` (ключ под `paths:`)
- [ ] Обновить `.gitignore` если нужно исключить новый формат файлов

---

### 18.5 Чеклист: добавить новую метрику

- [ ] Секция **7** («Метрики и статистика»): добавить в
  соответствующую таблицу
- [ ] Секция **19**: добавить запись
- [ ] Добавить вычисление в `src/statistics.py` или `src/validation.py`
- [ ] Добавить колонку в соответствующий CSV-отчёт
- [ ] Добавить в HTML-шаблон отчёта (`pipeline/run_pipeline.py` →
  `generate_html_report()`)
- [ ] Написать unit-тест в `tests/test_statistics.py`

---

### 18.6 Чеклист: изменить структуру проекта

- [ ] Секция **2** (дерево структуры): обновить дерево
- [ ] Секция **16** (таблица «Куда смотреть»): обновить строки
- [ ] Секция **9.1** (именование файлов): обновить шаблоны если изменились
- [ ] Секция **19**: добавить запись
- [ ] Обновить `.gitignore` если появились новые директории с данными
- [ ] Обновить `config/config.yaml` (ключи под `paths:`)

---

### 18.7 Правила оформления секций

```markdown
## N. Название секции          ← ## для верхнего уровня

### N.M Подраздел              ← ### для подразделов

| Колонка 1 | Колонка 2 |      ← таблицы для справочной информации
|-----------|-----------|

```python                       ← блоки кода с указанием языка
# Комментарий
```

- [ ] Пункт чеклиста           ← чеклисты для процедур обновления

> Заметка / предупреждение     ← blockquote для важных замечаний
```

---

### 18.8 Формат записи в Changelog (секция 19)

```markdown
### vMAJOR.MINOR — YYYY-MM-DD

**Добавлено**
- Краткое описание что добавлено

**Изменено**
- Краткое описание что изменено

**Удалено**
- Краткое описание что удалено

**Затронутые секции**: 5.N, 8, 15
**Затронутые файлы**: `src/algorithms/run_<algo>.py`, `config/config.yaml`
```

---

## 19. Changelog

### v1.2 — (дата следующего изменения)

> Шаблон для следующей записи. Скопировать, заполнить, удалить эту строку.

**Добавлено**
- ...

**Изменено**
- ...

**Удалено**
- ...

**Затронутые секции**: ...
**Затронутые файлы**: ...

---

### v1.1 — 2026-04-11

**Изменено**
- Разрешения: убрано 10kb, теперь [25000, 50000, 100000]
- download_data.sh: потоковая распаковка без сохранения архива,
  только 25kb/50kb/100kb
- .hic файл не используется — данные читаются из RAWobserved напрямую
- get_matrix: 4-уровневый fallback, при кэшировании сохраняет .npy и .RAWobserved
- Armatus: флаги исправлены (-S -N -c вместо -R -z), 
  output файл .consensus.txt
- Armatus: сборка через CMake + libboost-iostreams-dev
- scKTLD: dimension 128→32, добавлен knn_k=20,
  KNN-sparse RBF вместо полного M@M.T,
  _dp_segmentation_fast через cumsum,
  _auto_penalty с фильтром min/max TADs

**Добавлено**
- get_rawobserved_path_raw() в data_prep.py
- sys.path fix в pipeline/run_pipeline.py

**Затронутые секции**: 3, 4.2, 5.2, 5.4, 8, 11, 12, 17
**Затронутые файлы**: `download_data.sh`, `config/config.yaml`,
  `src/data_prep.py`, `src/algorithms/run_armatus.py`,
  `src/algorithms/run_scktld.py`, `pipeline/run_pipeline.py`

### v1.0 — 2025-01-01

**Добавлено**
- Начальная версия rules.md
- Секции 1–17: стек, архитектура, данные, алгоритмы (Armatus / TopDom /
  scKTLD / coiTAD), консенсус, метрики, конфигурация, соглашения,
  запреты, CLI, тесты, визуализация
- `download_data.sh` с источниками GEO GSE63525 и ENCODE CTCF

**Затронутые секции**: все (первая версия)
**Затронутые файлы**: `rules.md` (создан)

---

*Версия rules.md: 1.1 | Проект: TAD Consensus Pipeline | Геном: hg19 | Данные: GSE63525 GM12878*
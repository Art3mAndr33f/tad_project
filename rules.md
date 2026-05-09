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

  **Структура в git**: только `.gitkeep` (сам бинарник и `build/` в `.gitignore`)

    **Сборка без sudo (через conda, CMake 4.x):**
  ```bash
  conda install -c conda-forge boost boost-cpp cmake make zlib -y

  cd tools/armatus
  git clone https://github.com/kingsfordgroup/armatus.git .
  mkdir -p build && cd build

  cmake .. \
      -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
      -DCMAKE_C_COMPILER="${CONDA_PREFIX}/bin/gcc" \
      -DCMAKE_CXX_COMPILER="${CONDA_PREFIX}/bin/g++" \
      -DCMAKE_EXE_LINKER_FLAGS="-L${CONDA_PREFIX}/lib -Wl,-rpath,${CONDA_PREFIX}/lib -lboost_iostreams -lz" \
      -DBOOST_ROOT="${CONDA_PREFIX}" \
      -DBOOST_INCLUDEDIR="${CONDA_PREFIX}/include" \
      -DBOOST_LIBRARYDIR="${CONDA_PREFIX}/lib" \
      -DCMAKE_PREFIX_PATH="${CONDA_PREFIX}" \
      -DBoost_NO_SYSTEM_PATHS=ON \
      -DCMAKE_BUILD_TYPE=Release

  make -j2   # ← не j$(nproc): правило сервера — не нагружать CPU
  cp src/armatus ../../armatus
  chmod +x ../../armatus

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
- Флаги -R, -z, -j — НЕ использовать (см. §17)
- Выходной файл: `<prefix>.consensus.txt` (первый приоритет)
- Перебор gamma: `[0.1, 0.3, 0.5, 0.6, 0.65, 0.7, 0.75, 0.8, 0.9, 1.0, 2.0]`
- Выбор лучшей gamma: **valid zone + stability**
  1. Вычислить n_tads и средний Jaccard для каждой gamma
  2. Valid zone: `chrom_mb × 0.8 ≤ n_tads ≤ chrom_mb × 2.5`
     (эмпирическая плотность TAD ~0.8–2.5 TAD/Mb, GM12878 @ 25kb)
  3. Среди валидных — выбрать gamma с max(stability)
  4. Если valid zone пуста — выбрать ближайший к `target = (min+max)/2`
- Переменные в `run_armatus`: результаты хранятся в `results` (не `gamma_results`)
- `stability_top_n` читается из `cfg["algorithms"]["armatus"]["stability_top_n"]`
- Формат лога gamma: `%.2f` (не `%.1f` — иначе 0.65 и 0.7 неотличимы)
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
- Алгоритм: RBF-ядро → нормализованный Лапласиан → eigendecomposition →
  спектральное вложение → DP-сегментация
- `dimension=32`
- `knn_k=20`
- `penalty`: автоподбор через elbow на логарифмической сетке
- **Перед запуском**: обязательная проверка `scktld_limits[resolution]`
- При нарушении лимита — возвращать пустой DataFrame, НЕ падать
- **GPU-детектирование**: `DEVICE = _get_device()` на уровне модуля
  (использует `print`, не `logger` — вызывается до инициализации logging)
- **`_GPU_DENSE_THRESHOLD = 5000`**: n < 5000 бинов → GPU-ветка
  (chr17–chr22 @ 25kb: 1900–3250 бинов → всегда GPU)
- **`_spectral_embedding` ветки**:
  - GPU (n < 5000): `torch.linalg.eigh` на CUDA, плотная матрица ~100 MB
  - GPU (n ≥ 5000): fallback на CPU `eigsh` (OOM риск)
  - CPU: `scipy.sparse.linalg.eigsh` (текущая реализация без GPU)
- **torch**: установлен через `pip install torch --index-url https://download.pytorch.org/whl/cu121`

### 5.5 coiTAD (`src/algorithms/run_coitad.py`)

- Импорт из `tools/coiTAD/` через динамический `sys.path`
- Точки входа (в порядке приоритета):
  1. `CoiTADDetector().detect(matrix, resolution=resolution)`
  2. `detect_tads(matrix, resolution=resolution)`
  3. `run(matrix, resolution=resolution)`
- **Fallback**: Insulation Score (Crane et al. 2015) если `tools/coiTAD/` недоступен

#### Параметры fallback (откалиброваны на GM12878 @ 25kb):
| Параметр | Значение | Обоснование |
|----------|----------|-------------|
| `window_bins` | `max(5, 125_000 // resolution)` | 5 бинов @ 25kb = 125kb |
| `sigma` | 1.5 | умеренное сглаживание insulation score |
| `prominence_factor` | 0.20 | порог глубины минимума = mean - 0.20×std |
| `min_tad_kb` | 100 | минимальный размер TAD |

#### Адаптивная логика fallback:
- Если получено < 10 TAD → ослабить sigma на -0.5, prominence × 0.3
- Если получено > 200 TAD → ужесточить sigma на +1.0, prominence × 1.5

#### Ожидаемый результат:
- 30–80 TAD на хромосому при 25kb
- median_size_kb: 500–2000 kb

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
    # Расширенная сетка для покрытия переходной зоны при 25kb
    gamma_values: [0.1, 0.3, 0.5, 0.6, 0.65, 0.7, 0.75, 0.8, 0.9, 1.0, 2.0]
    stability_top_n: 3
    # Плотность TAD для valid zone (TAD/Mb):
    tad_density_min: 0.8
    tad_density_max: 2.5
  coitad:
    default_params:
      window_bins: null     # null = auto (max(5, 125000 // resolution))
      sigma: 1.5
      prominence_factor: 0.20
      min_tad_kb: 100
  topdom:
    window_sizes: [3, 5, 10]
    default_window: 5
  scktld:
    dimension: 32
    knn_k: 20
    balance: false
    auto_penalty: true

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
  generate_html: false  # отключено: Plotly рендер chr1@25kb занимает 39+ мин
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
   - ⚠️ `plt.Polygon`: использовать `facecolor=` + `edgecolor=` (не `color=` — вызывает дубль kwargs)
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
| `gamma_results` NameError | опечатка в имени переменной | использовать `results` (dict объявлен как `results: dict[...]`) |
| `logger` NameError в `_get_device()` | функция вызывается на уровне модуля до `logger = logging.getLogger(__name__)` | использовать `print()` внутри `_get_device()` |
| HTML-визуализация висит часами | Plotly рендерит матрицу chr1@25kb (~9970 бинов) | установить `generate_html: false` в config.yaml |
| `... chrX` в --chroms | `...` воспринимается буквально как имя хромосомы | перечислять хромосомы явно: `chr1 chr2 chr3 ... chr22 chrX` |
| stability_top_n NameError | не читается из cfg | добавить `stability_top_n = cfg[...].get("stability_top_n", 3)` в тело `run_armatus` |
| coiTAD 2–11 TADs | prominence=0.60 слишком строго | использовать 0.20 + window=125kb |
| coiTAD нет CTCF-пика | TADs слишком крупные (>3Mb) | исправить prominence_factor |
| Нагрузка CPU на сервере | скрипты используют numpy многопоток | `OMP_NUM_THREADS=2 MKL_NUM_THREADS=2` |
| GPU не используется | нет CUDA_VISIBLE_DEVICES | установить перед запуском |

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

### v1.4 — (дата следующего изменения)

> Шаблон для следующей записи. Скопировать, заполнить, удалить эту строку.

**Добавлено**
- ...

**Изменено**
- ...

**Удалено**
- ...

**Затронутые секции**: ...
**Затронутые файлы**: ...

### v1.3 — 2026-05-03

**Добавлено**
- `_get_device()` в run_scktld.py: print (не logger) т.к. module-level до logging
- `DEVICE = _get_device()` — module-level константа
- `_GPU_DENSE_THRESHOLD = 5000` — порог бинов для GPU-ветки
- GPU-ветка `_spectral_embedding`: `torch.linalg.eigh` на CUDA (n < 5000)
- torch установлен: `pip install torch --index-url .../cu121`
- Секция 17: три новых подводных камня (logger в _get_device, HTML медленно, `...` в --chroms)

**Изменено**
- `_spectral_embedding`: CPU eigsh → GPU/CPU двухветочная реализация
- `src/visualization.py` ~90: `color=` → `facecolor=`, убран дублирующийся `edgecolor=`
- `config/config.yaml`: `generate_html: true` → `false`
- Секции 5.4, 8, 15, 17, 20 обновлены

**Результаты прогона**
- chr1–chr22 + chrX @ 25kb / 50kb / 100kb: все TAD-листы готовы
- scKTLD пропускает chr1–chr16 @ 25kb (лимит памяти) — штатно
- Консенсус: 27–119 границ support≥2 на хромосому

**Затронутые секции**: 5.4, 8, 15, 17, 20
**Затронутые файлы**: `src/algorithms/run_scktld.py`, `src/visualization.py`, `config/config.yaml`, `rules.md`

---

### v1.2 — 2026-04-14

**Изменено**
- Armatus `_select_best_gamma`: добавлен valid zone фильтр по числу TAD
  (плотность 0.8–2.5 TAD/Mb × размер хромосомы в Mb)
- Armatus gamma_values: [0.1, 0.5, 1.0, 2.0, 5.0] → [0.1, 0.3, 0.5, 0.6, 0.65, 0.7, 0.75, 0.8, 0.9, 1.0, 2.0]
- Armatus лог формат: %.1f → %.2f (0.65 и 0.70 были неотличимы)
- Исправлен NameError: `gamma_results` → `results`, `stability_top_n` добавлен в scope
- coiTAD fallback: prominence_factor 0.60→0.20, sigma 2.5→1.5, window 250kb→125kb
- coiTAD адаптивное ослабление: factor×0.5→×0.3

**Добавлено**
- Секция 20: GPU-адаптация для сервера с правилами использования
- `_HG19_CHROM_SIZES_MB` в run_armatus.py
- `_TAD_DENSITY_MIN/MAX` константы
- `_get_device()` утилита для детектирования GPU
- `CUDA_VISIBLE_DEVICES` инструкции в README и rules.md
- GPU-версия scKTLD через cupy/torch (src/algorithms/run_scktld.py)

**Исправлено**
- Arrowhead warning spam: 90+ одинаковых предупреждений → однократное

**Затронутые секции**: 5.2, 5.5, 8, 17, 20 (новая)
**Затронутые файлы**: `src/algorithms/run_armatus.py`, `src/algorithms/run_coitad.py`,
  `src/algorithms/run_scktld.py`, `config/config.yaml`, `README.md`

*Версия rules.md: 1.3 | Проект: TAD Consensus Pipeline | Геном: hg19 | Данные: GSE63525 GM12878*

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

## 20. GPU-адаптация (сервер)

### 20.1 Правила сервера

⚠️ **На сервере ЗАПРЕЩЕНО нагружать CPU.** Все вычислительно тяжёлые операции
должны выполняться на GPU.

```bash
# Перед запуском — проверить свободный GPU:
nvtop
nvidia-smi

# Запуск с явным GPU:
CUDA_VISIBLE_DEVICES=0 python pipeline/run_pipeline.py ...
```

### 20.2 Какие части пайплайна используют GPU

| Компонент | CPU (локально) | GPU (сервер) | Ускорение |
|-----------|---------------|--------------|-----------|
| scKTLD: RBF-ядро | numpy/scipy | cupy | ~10–30× |
| scKTLD: eigsh | scipy.sparse.linalg | torch.lobpcg / cupy.sparse | ~5–20× |
| scKTLD: DP | numpy | numpy (O(n²) — не критично) | — |
| TopDom | numpy | numpy (быстро и так) | — |
| Armatus | subprocess C++ | subprocess C++ (без GPU) | — |
| coiTAD fallback | numpy/scipy | numpy (быстро и так) | — |

### 20.3 Требования для GPU-режима

```yaml
# environment.yml — добавить:
  - cupy-cuda12x>=13.0    # или cupy-cuda11x в зависимости от CUDA версии
  - pytorch>=2.0          # для torch.linalg альтернативы
```

**Фактически установлено на brain-lab (CUDA 13.0, Driver 580.82.07):**
```bash
pip install torch --index-url https://download.pytorch.org/whl/cu121
# Проверка: CUDA available: True | GPU 0: NVIDIA A100 80GB PCIe | Free: ~62 GB
```

### 20.4 Детектирование GPU в коде

```python
# Каноническая проверка — использовать везде где нужен GPU
def _get_device():
    """Определить устройство: 'cuda' если доступно, иначе 'cpu'."""
    try:
        import cupy as cp
        cp.cuda.runtime.getDeviceCount()   # проверка без создания массива
        return "cuda"
    except Exception:
        return "cpu"

DEVICE = _get_device()
```

### 20.5 Переменные окружения

```bash
# Обязательно перед запуском пайплайна на сервере:
export CUDA_VISIBLE_DEVICES=0   # или 1 — выбрать свободный GPU
export OMP_NUM_THREADS=2        # ограничить CPU-потоки
export MKL_NUM_THREADS=2
```

*Версия rules.md: 1.3 | Проект: TAD Consensus Pipeline | Геном: hg19 | Данные: GSE63525 GM12878*
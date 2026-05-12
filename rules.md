# Обновлённый `rules.md` (v2.0) и исправленный `README.md`

---

## `rules.md` — полная версия

```markdown
# TAD Consensus Pipeline — Project Rules for LLM Context

> Этот файл — единственный источник правды о проекте.
> Кидай его в начало любого нового чата, чтобы модель сразу знала контекст.

---

## 1. Что это за проект

Сравнительный анализ детекции топологически ассоциированных доменов (TAD)
на bulk Hi-C данных клеточной линии **GM12878** (human lymphoblastoid,
эталонный датасет Rao et al. 2014, GEO: **GSE63525**, координатная система
**hg19**).

Пайплайн запускает **семь** алгоритмов детекции TAD, строит консенсус границ,
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
├── download_data.sh
│
├── config/
│   └── config.yaml
│
├── data/
│   ├── raw/
│   ├── processed/
│   └── reference/
│
├── tools/
│   ├── armatus/              # бинарник armatus
│   ├── coiTAD/               # кастомная Python-версия coiTAD
│   │   ├── __init__.py
│   │   ├── coitad_core.py
│   │   └── detector.py
│   ├── ontad/                # бинарник OnTAD (собран из ontad_src)
│   │   └── OnTAD
│   └── ontad_src/            # исходники OnTAD (git clone)
│       └── src/
│           └── (исходники C++)
│
├── src/
│   ├── __init__.py
│   ├── data_prep.py
│   ├── algorithms/
│   │   ├── __init__.py            # ALGORITHM_REGISTRY (7 алгоритмов)
│   │   ├── run_armatus.py
│   │   ├── run_topdom.py
│   │   ├── run_scktld.py
│   │   ├── run_coitad.py
│   │   ├── run_dihmm.py           # NEW: DI + HMM
│   │   ├── run_ontad.py           # NEW: OnTAD (subprocess + Python fallback)
│   │   └── run_modularity_tad.py  # NEW: Graph modularity DP
│   ├── consensus.py
│   ├── statistics.py
│   ├── validation.py
│   └── visualization.py
│
├── pipeline/
│   └── run_pipeline.py
│
├── tests/
│   ├── test_consensus.py
│   └── test_statistics.py
│
├── notebooks/
│   └── exploration.ipynb
│
└── results/
    ├── tads/
    ├── consensus/
    ├── stats/
    ├── figures/
    └── report.html
```

---

## 3. Данные

### 3.1 Источники и файлы

| Файл | Путь | Описание |
|------|------|----------|
| Интрахромосомные матрицы | `data/raw/` | RAWobserved по хромосомам |
| Arrowhead domain list | `data/reference/GSE63525_GM12878_primary+replicate_Arrowhead_domainlist.txt` | Эталон Rao 2014 |
| CTCF ChIP-seq | `data/reference/GM12878_CTCF_peaks_hg19.bed` | ENCODE ENCFF796WRU, hg19 |

⚠️ .hic файл НЕ используется. Источник — RAWobserved файлы.

### 3.2 Скрипт скачивания (`download_data.sh`)

```bash
#!/bin/bash
set -e
mkdir -p data/raw data/reference
BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl"
wget -c -P data/raw/ \
    "${BASE}/GSE63525_GM12878_primary_intrachromosomal_contact_matrices.tar.gz"
tar -xzf data/raw/GSE63525_GM12878_primary_intrachromosomal_contact_matrices.tar.gz \
    -C data/raw/
wget -c -P data/reference/ \
    "${BASE}/GSE63525_GM12878_primary+replicate_Arrowhead_domainlist.txt.gz"
gunzip data/reference/*.gz
wget -c -O data/reference/GM12878_CTCF_peaks_hg19.bed.gz \
    "https://www.encodeproject.org/files/ENCFF796WRU/@@download/ENCFF796WRU.bed.gz"
gunzip data/reference/GM12878_CTCF_peaks_hg19.bed.gz
```

### 3.3 Разрешения и хромосомы

**Разрешения**: `[25000, 50000, 100000]`
**Хромосомы**: `chr1–chr22, chrX`

**Ограничения памяти для scKTLD**:

| Разрешение | Допустимые хромосомы |
|-----------|----------------------|
| 25 kb | chr17–chr22 |
| 50 kb | chr1–chr22 |
| 100 kb | chr1–chr22 |

### 3.4 Формат RAWobserved

3-колоночный TSV: `bin_i<TAB>bin_j<TAB>count`. Только верхний треугольник.

---

## 4. Технический стек

### 4.1 Зависимости

```
numpy>=1.24, scipy>=1.11, pandas>=2.0, scikit-learn>=1.3
matplotlib>=3.7, seaborn>=0.12, plotly>=5.18
PyYAML>=6.0, tqdm>=4.66, jinja2>=3.1, click>=8.1
hmmlearn>=0.3.3        # для run_dihmm.py
cooler>=0.9.3, pybedtools>=0.9.1, pyranges>=0.0.129
pytest>=7.4, pytest-cov>=4.1
```

### 4.2 Внешние инструменты

- **Armatus**: C++ бинарник в `tools/armatus/armatus`
- **OnTAD**: C++ бинарник в `tools/ontad/OnTAD`

**Сборка OnTAD (без sudo, с conda GCC и libcurl):**
```bash
conda install -c conda-forge libcurl gcc gxx make -y

cd tools/
git clone https://github.com/anlin00007/OnTAD.git ontad_src
cd ontad_src/src

# Патч: добавить #include <cstdint> (GCC 15 строже)
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

# Сборка с явными conda-путями
g++ -std=c++11 \
    -I${CONDA_PREFIX}/include \
    -L${CONDA_PREFIX}/lib \
    -Wl,-rpath,${CONDA_PREFIX}/lib \
    main.cpp step1.cpp step2.cpp step3.cpp step4.cpp common.cpp straw.cpp \
    -lm -lcurl -lz -o OnTAD

mkdir -p ~/tad_project/tools/ontad
cp OnTAD ~/tad_project/tools/ontad/OnTAD
chmod +x ~/tad_project/tools/ontad/OnTAD
```

⚠️ **Не использовать** `make` напрямую — он не передаёт `-I${CONDA_PREFIX}/include`.
⚠️ **Не отключать** `#include <curl/curl.h>` в straw.cpp — это сломает CURL-типы.

---

## 5. Алгоритмы

### 5.1 Единый интерфейс (ОБЯЗАТЕЛЬНО)

```python
def run_<algorithm>(
    chrom: str,
    resolution: int,
    data_path: str,
    cfg: Optional[dict],
    **kwargs,
) -> pd.DataFrame:
    # columns: ["chrom", "start", "end"]
    # При ошибке: pd.DataFrame(columns=["chrom","start","end"])
```

### 5.2 Armatus (`src/algorithms/run_armatus.py`)

- Subprocess, флаги: `-S -N -c <chrom>`
- Gamma sweep: `[0.1, 0.3, 0.5, 0.6, 0.65, 0.7, 0.75, 0.8, 0.9, 1.0, 2.0]`
- Valid zone: `chrom_mb × 0.8 ≤ n_tads ≤ chrom_mb × 2.5`
- Output: `.consensus.txt`
- Timeout: 600 с

### 5.3 TopDom (`src/algorithms/run_topdom.py`)

- Чистый Python
- Window sweep: `[3, 5, 10]`
- Выбор: max intra/inter ratio

### 5.4 scKTLD (`src/algorithms/run_scktld.py`)

- RAW counts, `balance=False`
- `dimension=32`, `knn_k=20`
- GPU-ветка (n < 5000): `torch.linalg.eigh` на CUDA
- Ограничения по хромосомам: см. §3.3

### 5.5 coiTAD (`src/algorithms/run_coitad.py`)

- Импорт из `tools/coiTAD/`
- Fallback: Insulation Score

### 5.6 DI+HMM (`src/algorithms/run_dihmm.py`) ← NEW

- **Метод**: Directionality Index (Dixon et al. 2012) + GaussianHMM (3 состояния)
- **Входные данные**: dense numpy matrix (RAW counts)
- **Нормализация**: не требуется
- **Зависимости**: `hmmlearn>=0.3.3` (pip)
- **Параметры** (из конфига):
  - `window_bins: 10` — окно для DI (10 бинов × 25kb = 250kb)
  - `n_states: 3` — upstream / unbiased / downstream
  - `n_iter: 100` — итерации EM для HMM
  - `min_tad_kb: 100.0`
- **Стратегия извлечения TAD** (трёхуровневая):
  1. Dixon паттерн: блоки `state=2 → state=0`
  2. Boundary-mode: переходы состояний HMM
  3. DI-minima fallback: нулевые пересечения сглаженного DI (σ=2)
- **Результаты** (GM12878 @ 25kb):
  - chr22: 54 TADs, median=488kb, max=2575kb ✅
  - chr21: 68 TADs, median=450kb, max=9650kb ✅
- **Статус**: ✅ Готов к полному прогону

### 5.7 OnTAD (`src/algorithms/run_ontad.py`) ← NEW

- **Метод**: иерархическая детекция через скользящее среднее (An et al. 2019)
- **Входные данные**: N×N матрица TAB-separated (записывается во временный файл)
- **Нормализация**: `-log2` флаг (log2(x+1)), RAW counts
- **Внешние зависимости**: C++ бинарник `tools/ontad/OnTAD`
- **Fallback**: Python-реализация OnTAD через 2D prefix-суммы если subprocess не даёт TAD
- **Параметры** (из конфига):
  - `penalty: null` — null = auto-sweep
  - `minsz: 3` — минимальный TAD в бинах
  - `maxsz: 200` — максимальный TAD в бинах
  - `log2: true`
- **Формат `.tad` файла** (5 колонок, TAB):
  ```
  start_bin  end_bin  depth  score1  score2
  1          2050     0      ...     ...    ← root (depth=0), пропускаем
  685        747      1      ...     ...    ← top-level TAD (depth=1) ← используем
  685        731      2      ...     ...    ← sub-TAD, пропускаем
  ```
  Конвертация: `start_bp = (start_bin - 1) * resolution`, `end_bp = end_bin * resolution`
- **⚠️ ТЕКУЩАЯ ПРОБЛЕМА** (нерешена):
  - Penalty sweep `[1.0, 0.5, 0.3, 0.2, 0.1, 0.05, 0.02, 0.01]` даёт только 16 TADs chr22 (target 51–128)
  - `depth=1` возвращает только крупные домены верхнего уровня независимо от penalty
  - **Гипотеза**: нужно использовать не `depth=1` а **максимальный непересекающийся набор TAD** (листья дерева), но только те что >= min_tad_kb
  - **Следующий шаг**: реализовать `_get_max_nonoverlapping(df_all, min_tad_bins)` — жадный выбор листовых TAD снизу вверх по иерархии
- **Результаты** (GM12878 @ 25kb):
  - chr22: 16 TADs @ depth=1 ⚠️ (мало, target 51–128)
  - chr21: 10 TADs @ depth=1 ⚠️
- **Статус**: ⚠️ Требует исправления уровня иерархии

### 5.8 ModularityTAD (`src/algorithms/run_modularity_tad.py`) ← NEW

- **Метод**: 1D-адаптация Newman graph modularity Q для Hi-C, DP-сегментация
- **Входные данные**: dense numpy matrix (RAW counts)
- **Нормализация**: не требуется
- **Нет внешних зависимостей** (только numpy)
- **Параметры** (из конфига):
  - `max_dist_mb: 5.0` — максимальное расстояние для контактов
  - `penalty: null` — null = auto-sweep (рекомендуется)
  - `min_tad_kb: 100.0`
  - `max_tad_kb: 3000.0`
- **Auto-penalty**: sweep `[2.0, 1.0, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.001]` после нормировки score на max
- **⚠️ ТЕКУЩАЯ ПРОБЛЕМА** (нерешена):
  - 331/328 TADs (target 51–128), median=100kb = ровно min_tad_bins
  - Auto-penalty sweep не срабатывает: в логе видно только одно значение `penalty=0.0500` → значит `penalty: 0.05` захардкожен в config.yaml
  - Score после нормировки: median=0.023, max=1.0 → penalty=0.05 > median, но DP всё равно дробит на минимальные куски
  - **Гипотеза 1**: нормировка score/max неверна — нужно нормировать score на `(l * l)` (площадь TAD), а не на глобальный max
  - **Гипотеза 2**: DP objective неверен — нужна не сумма mean_B, а разница intra vs inter (как в OnTAD)
  - **Следующий шаг**: переключиться на intra/inter ratio objective (аналог TopDom score), это даёт интуитивно более правильную функцию: `score(i,j) = mean_A_intra(i,j) - mean_A_flanks(i,j)`, без модулярного лапласиана
- **Результаты** (GM12878 @ 25kb):
  - chr22: 331 TADs ⚠️ (мало, target 51–128)
  - chr21: 328 TADs ⚠️
- **Статус**: ⚠️ Требует переработки objective function

### 5.9 ALGORITHM_REGISTRY (`src/algorithms/__init__.py`)

```python
from .run_armatus        import run_armatus
from .run_topdom         import run_topdom
from .run_scktld         import run_scktld
from .run_coitad         import run_coitad
from .run_dihmm          import run_dihmm
from .run_ontad          import run_ontad
from .run_modularity_tad import run_modularity_tad

ALGORITHM_REGISTRY = {
    "armatus":        run_armatus,
    "topdom":         run_topdom,
    "scktld":         run_scktld,
    "coitad":         run_coitad,
    "dihmm":          run_dihmm,
    "ontad":          run_ontad,
    "modularity_tad": run_modularity_tad,
}
```

---

## 6. Консенсус границ

### 6.1 Алгоритм

1. Собрать все границы TAD (start и end) от всех алгоритмов
2. Округлить до ближайшего бина
3. Жадная кластеризация: tolerance ≤ `tolerance_bins * resolution`
4. Support = число алгоритмов с границей в кластере

### 6.2 Параметры

```yaml
consensus:
  tolerance_bins: 1
  min_support: 2
```

### 6.3 Цветовая схема (ФИКСИРОВАНА)

| Support | Цвет | Hex |
|---------|------|-----|
| 2 | 🟡 Жёлтый | `#FFD700` |
| 3 | 🟠 Оранжевый | `#FF8C00` |
| 4+ | 🟢 Зелёный | `#00C800` |

```python
CONSENSUS_COLORS = {2: "#FFD700", 3: "#FF8C00", 4: "#00C800"}
```

---

## 7. Метрики и статистика

### 7.1 Базовая статистика

| Метрика | Описание |
|---------|----------|
| `n_tads` | Число TAD |
| `median_size_kb` | Медианный размер TAD (kb) |
| `mean_size_kb` | Средний размер TAD (kb) |
| `coverage_pct` | % покрытия хромосомы |
| `n_consensus_bnd` | Число консенсусных границ (support≥2) |
| `pct_consensus_bnd` | % консенсусных границ от всех |

### 7.2 Попарное сравнение

- **Jaccard index** на уровне доменов
- **Boundary overlap rate** ±1 bin

### 7.3 Сравнение с эталоном Rao 2014

Метрики при допусках ±1 bin и ±2 bin: Recall, Precision, F1, Jaccard

### 7.4 CTCF-валидация

- Окно: `±1 bin * resolution`
- 1000 пермутаций (seed=42)
- **⚠️ ИЗВЕСТНАЯ ПРОБЛЕМА**: CTCF-обогащение у первых 4 алгоритмов (armatus, topdom, scktld, coitad) не показывает чёткого пика на графиках
- **Гипотезы**:
  1. Слишком широкое окно вокруг границы → сигнал размыт
  2. RAW матрица без нормализации → границы смещены
  3. Координатное несоответствие CTCF BED и TAD (hg19 vs hg38, или chr-prefix)
  4. Слишком мало хромосом в анализе (chr17–chr22 @ 25kb)
- **Задача**: диагностировать CTCF-валидацию:
  1. Проверить координатный формат CTCF BED: `head data/reference/GM12878_CTCF_peaks_hg19.bed`
  2. Проверить chr-prefix: `cut -f1 data/reference/GM12878_CTCF_peaks_hg19.bed | sort -u`
  3. Запустить валидацию на chr1 @ 100kb (много TAD, лучше статистика)
  4. Построить профиль обогащения ±500kb и проверить форму кривой
  5. Сравнить с эталоном Rao 2014 Arrowhead: если у Arrowhead тоже нет пика → проблема в данных CTCF

---

## 8. Конфигурация

`config/config.yaml`:

```yaml
project:
  seed: 42

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
  ontad_bin:      "tools/ontad/OnTAD"
  coitad_dir:     "tools/coiTAD"
  hic_file:       ""
  arrowhead_file: "data/reference/GSE63525_GM12878_primary+replicate_Arrowhead_domainlist.txt"
  ctcf_bed:       "data/reference/GM12878_CTCF_peaks_hg19.bed"

resolutions: [25000, 50000, 100000]

chromosomes:
  all: ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9",
        "chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17",
        "chr18","chr19","chr20","chr21","chr22","chrX"]
  scktld_limits:
    25000:  ["chr17","chr18","chr19","chr20","chr21","chr22"]
    50000:  ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9",
             "chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17",
             "chr18","chr19","chr20","chr21","chr22"]
    100000: ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9",
             "chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17",
             "chr18","chr19","chr20","chr21","chr22","chrX"]

algorithms:
  armatus:
    gamma_values: [0.1, 0.3, 0.5, 0.6, 0.65, 0.7, 0.75, 0.8, 0.9, 1.0, 2.0]
    stability_top_n: 3
    tad_density_min: 0.8
    tad_density_max: 2.5
  coitad:
    default_params:
      window_bins: null
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
  dihmm:
    window_bins: 10
    n_states: 3
    n_iter: 100
    min_tad_kb: 100.0
  ontad:
    penalty: null          # null = auto-sweep [1.0..0.01]
    minsz: 3
    maxsz: 200
    log2: true
  modularity_tad:
    max_dist_mb: 5.0
    penalty: null          # null = auto-sweep [2.0..0.001]
    min_tad_kb: 100.0
    max_tad_kb: 3000.0

consensus:
  tolerance_bins: 1
  min_support: 2
  colors:
    2: "#FFD700"
    3: "#FF8C00"
    4: "#00C800"

validation:
  ctcf_window_bp: 1
  n_permutations: 1000
  enrichment_profile_range_bp: 500000
  enrichment_profile_bin_bp: 10000

rao_comparison:
  tolerances_bins: [1, 2]

visualization:
  hic_colormap: "coolwarm"
  hic_log_scale: true
  dpi: 300
  format: "png"
  figsize: [18, 14]
  generate_html: false
  algorithm_colors:
    armatus:        "#1f77b4"
    topdom:         "#ff7f0e"
    scktld:         "#2ca02c"
    coitad:         "#d62728"
    dihmm:          "#9467bd"
    ontad:          "#e377c2"
    modularity_tad: "#17becf"
```

---

## 9. Соглашения по именованию

### 9.1 Файлы результатов

| Тип | Шаблон | Пример |
|-----|--------|--------|
| TAD-лист | `results/tads/<algo>_<chrom>_<res>bp.bed` | `dihmm_chr17_25000bp.bed` |
| Консенсус | `results/consensus/consensus_<chrom>_<res>bp.bed` | |
| Статистика | `results/stats/<name>.csv` | |
| Фигура | `results/figures/hic_tads_<chrom>_<res>bp.png` | |

### 9.2 Матрицы в `data/processed/`

| Тип | Шаблон |
|-----|--------|
| Dense numpy | `<chrom>_<res>bp.npy` |
| RAWobserved | `<chrom>_<res>bp.RAWobserved` |

### 9.3 Переменные и функции

```python
chrom = "chr17"        # ✅ всегда с префиксом chr
resolution = 25000     # ✅ в bp (int)
df.columns == ["chrom", "start", "end"]  # ✅ строго
pd.DataFrame(columns=["chrom", "start", "end"])  # ✅ при ошибке
```

### 9.4 Логирование

```python
logger = logging.getLogger(__name__)
# print() разрешён ТОЛЬКО в _get_device() — вызывается до init logging
```

---

## 10. Правила и запреты

### 10.1 ЗАПРЕЩЕНО

```python
# ❌ Хардкодить пути, параметры, seed
# ❌ print() в src/ и pipeline/ (кроме _get_device())
# ❌ Нормализовать матрицу для scKTLD
# ❌ Менять CONSENSUS_COLORS
# ❌ Возвращать None из run_<algorithm>
# ❌ Менять seed (везде 42)
# ❌ Коммитить data/ и results/
# ❌ Отключать #include <curl/curl.h> в straw.cpp OnTAD
# ❌ Использовать make без явных -I${CONDA_PREFIX}/include (OnTAD)
```

### 10.2 ОБЯЗАТЕЛЬНО

```python
# ✅ cfg = load_config("config/config.yaml")
# ✅ Все параметры из конфига
# ✅ Проверять scktld_limits до запуска scKTLD
# ✅ Path(out).parent.mkdir(parents=True, exist_ok=True)
# ✅ rng = np.random.default_rng(42)
# ✅ balance=False для scKTLD
```

### 10.3 Обработка ошибок

```python
try:
    df = run_algorithm(...)
except Exception as exc:
    logger.error("[algo] %s @ %d: %s", chrom, resolution, exc, exc_info=True)
    df = pd.DataFrame(columns=["chrom", "start", "end"])
```

---

## 11. data_prep.py — ключевые функции

```python
cfg = load_config("config/config.yaml")
matrix = get_matrix(cfg, chrom, resolution)  # 4-уровневый fallback
raw_path = get_rawobserved_path(cfg, chrom, resolution)
prepare_all_matrices(cfg, resolutions, chromosomes, force=False)
```

---

## 12. CLI пайплайна

```bash
# Полный пайплайн (4 классических алгоритма)
python pipeline/run_pipeline.py \
    --resolution 25000 \
    --chroms chr17 chr18 chr19 chr20 chr21 chr22 \
    --algorithms armatus topdom scktld coitad

# Новые алгоритмы
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

# Переменные окружения (обязательно на сервере)
export CUDA_VISIBLE_DEVICES=3   # GPU 3: ~40GB свободно
export OMP_NUM_THREADS=2
export MKL_NUM_THREADS=2
```

---

## 13. Тестирование

```bash
pytest tests/ -v
pytest tests/ --cov=src --cov-report=html
pytest tests/ -v -k "not Integration"
```

---

## 14. Размеры хромосом hg19

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

### Цвета алгоритмов (ФИКСИРОВАНЫ)

```python
ALGO_COLORS = {
    "armatus":        "#1f77b4",   # синий
    "topdom":         "#ff7f0e",   # оранжевый
    "scktld":         "#2ca02c",   # зелёный
    "coitad":         "#d62728",   # красный
    "dihmm":          "#9467bd",   # фиолетовый
    "ontad":          "#e377c2",   # розовый
    "modularity_tad": "#17becf",   # циан
}
```

---

## 16. Типичные задачи и куда смотреть

| Задача | Файл |
|--------|------|
| Добавить новый алгоритм | `src/algorithms/run_<new>.py` + `ALGORITHM_REGISTRY` |
| Изменить параметры | `config/config.yaml` |
| Изменить консенсус | `src/consensus.py` |
| Добавить метрику | `src/statistics.py` |
| Изменить CTCF-анализ | `src/validation.py` |
| Добавить шаг в пайплайн | `pipeline/run_pipeline.py` |

---

## 17. Частые подводные камни

| Проблема | Причина | Решение |
|----------|---------|---------|
| OOM при scKTLD | Крупные хромосомы | Проверять `scktld_limits` |
| Armatus 0 TADs | Флаг `-R` вместо `-S -N` | Использовать `-S -N -c <chrom>` |
| OnTAD не компилируется: `curl/curl.h` | libcurl не в PATH GCC | `conda install libcurl`, явный `-I${CONDA_PREFIX}/include` |
| OnTAD не компилируется: `uint64_t` | GCC 15 строже | Добавить `#include <cstdint>` после `#include "straw.h"` |
| OnTAD: `make` не находит libcurl | make не передаёт conda пути | Использовать явный `g++ -I${CONDA_PREFIX}/include ...` |
| OnTAD: двойной `#include <cstdint>` | Патч применён дважды | Проверить `grep -c cstdint straw.cpp` перед патчем |
| OnTAD: 0 TADs в subprocess | Парсер `parts[3]` вместо `parts[2]` для depth | depth = `parts[2]`, не `parts[3]` |
| OnTAD: 16 TADs @ depth=1 | depth=1 = только крупные домены | Нужен max-nonoverlapping набор листовых TAD |
| ModularityTAD: 300+ TADs | penalty из config переопределяет auto | Удалить `penalty:` из config или поставить `null` |
| ModularityTAD: DP предпочитает min_tad_bins | Score не зависит от размера TAD | Нормировать intra/inter, а не mean(B) |
| logger NameError в `_get_device()` | Функция вызывается до init logging | Использовать `print()` |
| HTML-визуализация зависает | Plotly chr1@25kb | `generate_html: false` |
| CTCF нет чёткого пика | Возможна проблема координат или данных | Диагностика: §7.4 |
| pytadbit не устанавливается | Нет wheel для Python 3.10 | Не использовать; заменён на DI+HMM |

---

## 18. Как обновлять этот файл

> docs-first подход: rules.md обновляется ДО изменения кода.

### 18.1 Версионирование

| Тип изменения | Что менять |
|---------------|------------|
| Новый алгоритм, новый источник данных | MAJOR + 1 |
| Новый параметр, новая метрика, правка запретов | MINOR + 1 |
| Опечатки, форматирование | не меняется |

### 18.2 Чеклист: добавить новый алгоритм

**Шаг 1 — rules.md:**
- [ ] §5: добавить подраздел `5.N`
- [ ] §5.9: добавить в ALGORITHM_REGISTRY
- [ ] §8: добавить блок в `algorithms:`
- [ ] §15: добавить цвет в `ALGO_COLORS`
- [ ] §19: Changelog

**Шаг 2 — код:**
- [ ] `src/algorithms/run_<algo>.py`
- [ ] `src/algorithms/__init__.py`
- [ ] `config/config.yaml`
- [ ] `src/visualization.py` (ALGO_COLORS)
- [ ] `requirements.txt` / `environment.yml`
- [ ] `tests/test_<algo>.py`

---

## 19. Changelog

### v2.0 — 2026-05-12

**Добавлено**
- Три новых алгоритма: DI+HMM (§5.6), OnTAD (§5.7), ModularityTAD (§5.8)
- `hmmlearn>=0.3.3` в зависимости
- `tools/ontad/` — бинарник OnTAD (сборка с conda gcc + libcurl)
- Секция §20 GPU-адаптация обновлена (GPU 3 = 40GB свободно)
- Инструкция сборки OnTAD в §4.2 (патч straw.cpp + явный g++)
- Документация текущих проблем OnTAD и ModularityTAD в §5.7 и §5.8
- Задача диагностики CTCF-валидации в §7.4
- Удалена ссылка на несуществующий `scripts/update_rules.py`

**Изменено**
- §5.9 ALGORITHM_REGISTRY: 4 → 7 алгоритмов
- §8 config.yaml: добавлены блоки `dihmm`, `ontad`, `modularity_tad`, `ontad_bin`
- §12 CLI: добавлены команды для новых алгоритмов
- §15 ALGO_COLORS: добавлены 3 новых цвета
- §17: добавлены 7 новых подводных камней (OnTAD сборка и парсинг, ModularityTAD, CTCF)
- README.md: структура проекта обновлена (новые файлы), удалён `scripts/update_rules.py`

**Результаты smoke-test** (chr21, chr22 @ 25kb):
- DI+HMM: 54/68 TADs ✅
- OnTAD: 16/10 TADs ⚠️ (требует исправления уровня иерархии)
- ModularityTAD: 331/328 TADs ⚠️ (требует переработки objective)

**Затронутые секции**: 2, 4.2, 5 (5.6–5.9), 7.4, 8, 12, 15, 17, 19, 20
**Затронутые файлы**: `src/algorithms/run_dihmm.py`, `src/algorithms/run_ontad.py`,
  `src/algorithms/run_modularity_tad.py`, `src/algorithms/__init__.py`,
  `config/config.yaml`, `src/visualization.py`, `README.md`, `rules.md`

### v1.3 — 2026-05-03

*(см. предыдущую версию)*

---

## 20. GPU-адаптация (сервер brain-lab)

### 20.1 Правила сервера

⚠️ Не нагружать CPU. Тяжёлые вычисления — на GPU.

```bash
# Статус GPU (актуально на 2026-05-12):
# GPU 0: 73989/81154 MiB занято — НЕ использовать
# GPU 1: 55443/81154 MiB занято — НЕ использовать
# GPU 2: 76549/81154 MiB занято — НЕ использовать
# GPU 3: 41114/81154 MiB занято, ~40 GB свободно ← ИСПОЛЬЗОВАТЬ

export CUDA_VISIBLE_DEVICES=3
export OMP_NUM_THREADS=2
export MKL_NUM_THREADS=2
```

### 20.2 Использование GPU по компонентам

| Компонент | CPU | GPU | Ускорение |
|-----------|-----|-----|-----------|
| scKTLD (n<5000) | scipy eigsh | torch.linalg.eigh CUDA | ~10–20× |
| TopDom, Armatus, coiTAD | numpy/subprocess | — | — |
| DI+HMM | numpy + hmmlearn | — (достаточно CPU) | — |
| OnTAD | subprocess C++ | — | — |
| ModularityTAD | numpy | — (достаточно CPU) | — |

### 20.3 Установка torch (brain-lab, CUDA 12.1)

```bash
pip install torch --index-url https://download.pytorch.org/whl/cu121
# Проверка: python -c "import torch; print(torch.cuda.is_available())"
```

*Версия rules.md: 2.0 | Проект: TAD Consensus Pipeline | Геном: hg19 | Данные: GSE63525 GM12878*
```
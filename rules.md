# TAD Consensus Pipeline — Project Rules for LLM Context
# Version: 3.4 | 2026-05-20 | Genome: hg19 | Data: GSE63525 GM12878

> **docs-first:** rules.md обновляется ДО изменения кода.
> Кидай этот файл в начало любого нового чата — это единственный источник правды.
> Подробная документация по темам: см. `agent_docs/` (загружай только нужный файл).

---

## 1. Суть проекта (WHY + WHAT)

**Научная задача диплома:**
Существующий подход к разметке данных для обучения TAD-детекторов
(deepTAD-style консенсус нескольких алгоритмов) биологически неинформативен.
Задача — доказать это и предложить пайплайн на основе статистических методов
с биологически информированной валидацией через CTCF ChIP-seq.

**Итоговый импакт:**
Общий пайплайн сбора данных для обучения Hi-C TAD-экстракторов, основанный
на статистическом отборе границ регуляторных элементов, сопрягаемый с
биологически информированной валидацией через CTCF ChIP-seq, с публикацией
размеченных Hi-C карт для обучения нейросетей.

**Пайплайн (текущая реализация):**
1. Запускает 7 алгоритмов детекции TAD (GM12878, hg19, RAWobserved)
2. Строит консенсус границ — **baseline для критики**
3. Валидирует через CTCF ChIP-seq — **биологический ground truth**
4. Предлагает статистические методы (Фурье, спектр, TAD-overlap) — **новизна**

**OS:** Ubuntu 22.04 | **Python:** 3.10+ | **Conda env:** `tad_pipeline`
**GPU:** `CUDA_VISIBLE_DEVICES=2` (brain-lab, ~20 GB свободно)

---

## 2. Быстрый старт

```bash
conda activate tad_pipeline
export CUDA_VISIBLE_DEVICES=2 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2

# Все 7 алгоритмов
python pipeline/run_pipeline.py \
    --resolution 25000 \
    --chroms chr17 chr18 chr19 chr20 chr21 chr22 \
    --algorithms armatus topdom scktld coitad dihmm ontad modularity_tad \
    --force

# Тесты
pytest tests/ -v

# Новый участник — контекст за 10 минут
cat rules.md
cat agent_docs/09_workstreams.md
cat agent_docs/<твой_файл>.md
git checkout -b ws/<workstream>
```

---

## 3. Архитектура (карта кода)

```
src/algorithms/          ← 7 детекторов TAD (единый интерфейс → §5)
src/consensus.py         ← жадная кластеризация границ + TAD-консенсус (Jaccard)
src/validation.py        ← CTCF ChIP-seq валидация (биол. ground truth)
src/ctcf_analysis.py     ← профильный анализ ChIP-seq треков (параметризован)
src/statistics.py        ← Jaccard, boundary overlap, сравнение с Arrowhead
src/visualization.py     ← Hi-C heatmap, CTCF профили, Jaccard matrix
src/data_prep.py         ← загрузка матриц (4-уровневый fallback)
pipeline/run_pipeline.py ← оркестратор
scripts/run_chipseq_validation.py ← валидация по RAD21/SMC3/H3K4me3/H3K27ac
config/config.yaml       ← ВСЕ параметры (не хардкодить!)
agent_docs/              ← подробная дока по темам (загружай только нужное)
thesis/                  ← текст диплома (главы в Markdown/LaTeX)
logs/                    ← логи запусков (не в git)
```

**Карта agent_docs/ — загружай только нужный файл:**

| Файл | Загружай когда... |
|------|-------------------|
| `00_index.md` | Всегда — это карта |
| `01_algorithms.md` | Работаешь с алгоритмами, параметрами |
| `02_data.md` | Работаешь с данными, форматами |
| `03_ctcf_validation.md` | Работаешь с CTCF, биол. интерпретацией |
| `04_build_ontad.md` | Собираешь OnTAD из исходников |
| `05_known_issues.md` | Что-то сломалось |
| `06_thesis_structure.md` | Пишешь диплом |
| `07_visualization.md` | Работаешь с графиками |
| `08_gpu_server.md` | Запускаешь на brain-lab |
| `09_workstreams.md` | Параллельная работа |
| `10_statistical_methods.md` | Фурье, спектр, новые методы |

---

## 4. Данные

| Источник | Путь | Описание |
|----------|------|----------|
| Hi-C матрицы | `data/raw/` | RAWobserved, GM12878 primary |
| Arrowhead эталон | `data/reference/GSE63525_..._Arrowhead_domainlist.txt` | Rao 2014 |
| CTCF ChIP-seq | `data/reference/GM12878_CTCF_peaks_hg19.bed` | ENCODE ENCFF796WRU, 44 217 пиков |
| RAD21 ChIP-seq | `data/reference/GM12878_RAD21_peaks_hg19.bed` | ENCODE ENCFF001VFE (narrowPeak), 23 947 пиков ✅ 2026-05-19 |
| SMC3 ChIP-seq | `data/reference/GM12878_SMC3_peaks_hg19.bed` | ENCODE ENCFF001VFH (narrowPeak), 64 597 пиков ✅ 2026-05-19 |
| H3K4me3 ChIP-seq | `data/reference/GM12878_H3K4me3_peaks_hg19.bed` | UCSC wgEncodeBroadHistone (broadPeak), 57 476 пиков ✅ 2026-05-19 |
| H3K27ac ChIP-seq | `data/reference/GM12878_H3K27ac_peaks_hg19.bed` | UCSC wgEncodeBroadHistone (broadPeak), 56 069 пиков ✅ 2026-05-19 |

**Пути ChIP-seq треков в config/config.yaml (секция paths):**
```yaml
ctcf_bed:            "data/reference/GM12878_CTCF_peaks_hg19.bed"
rad21_bed:           "data/reference/GM12878_RAD21_peaks_hg19.bed"
smc3_bed:            "data/reference/GM12878_SMC3_peaks_hg19.bed"
h3k4me3_bed:         "data/reference/GM12878_H3K4me3_peaks_hg19.bed"
h3k27ac_bed:         "data/reference/GM12878_H3K27ac_peaks_hg19.bed"
hic_figures_out:     "results/figures/hic_maps"
chipseq_figures_out: "results/figures/chipseq_profiles"
arrowhead_ref:       "data/reference/GSE63525_GM12878_primary+replicate_Arrowhead_domainlist.txt"
```

**Разрешения:** `[25000, 50000, 100000]`
**Хромосомы:** `chr1–chr22, chrX`
**⚠️ .hic файл НЕ используется** — только RAWobserved (3-кол. TSV: bin_i, bin_j, count)
**⚠️ Координатная система везде hg19** — без исключений

**Ограничения памяти scKTLD:**
| Разрешение | Хромосомы |
|-----------|-----------|
| 25kb | chr17–chr22 |
| 50kb | chr1–chr22 |
| 100kb | chr1–chr22, chrX |

---

## 5. Алгоритмы — 7 штук, единый интерфейс

### Обязательный контракт (нарушение ломает pipeline)

```python
def run_<algorithm>(
    chrom: str,          # "chr17" — всегда с префиксом chr
    resolution: int,     # 25000 / 50000 / 100000 (в bp)
    data_path: str,
    cfg: Optional[dict],
    **kwargs,
) -> pd.DataFrame:
    # Обязательные колонки: ["chrom", "start", "end"]
    # При ЛЮБОЙ ошибке: return pd.DataFrame(columns=["chrom","start","end"])
```

### Реестр алгоритмов

| Алгоритм | Метод | Класс | Статья | Статус |
|----------|-------|-------|--------|--------|
| `armatus` | DP gamma-регуляризация | DP | Filippova 2014 | ✅ |
| `topdom` | Insulation Score | Локальный | Shin 2016 | ✅ |
| `scktld` | Spectral kNN clustering | Спектральный | Zheng 2024 | ✅ |
| `coitad` | OI matrix + IS fallback | Локальный | — | ✅ |
| `dihmm` | Directionality Index + HMM | Вероятностный | Dixon 2012 | ✅ |
| `ontad` | Иерархическая детекция | DP-иерархия | An 2019 | ✅ |
| `modularity_tad` | Graph modularity + DP | Граф | — | ✅ |

**Подробные параметры, механика, результаты → `agent_docs/01_algorithms.md`**

### ALGORITHM_REGISTRY (`src/algorithms/__init__.py`)

```python
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

**Алгоритм:** Собрать все границы (start/end) → округлить до бина →
жадная кластеризация с `tolerance_bins=1` → support = число алгоритмов в кластере.

**TAD-консенсус (Jaccard):** `compute_tad_consensus()` — Union-Find по парам TAD
из разных алгоритмов с Jaccard ≥ jaccard_threshold. Интегрирован в
`compute_all_consensus()`. Выход: `results/consensus/tad_consensus_{chrom}_{res}bp.bed`.

```yaml
consensus:
  tolerance_bins: 1
  min_support: 2
  jaccard_threshold: 0.5   # порог Jaccard для TAD-консенсуса
```

**Цветовая схема CONSENSUS_COLORS (расширена до 7 уровней — не менять):**

| Support | Цвет | Hex |
|---------|------|-----|
| 2 | 🟡 Жёлтый | `#FFD700` |
| 3 | 🟠 Оранжевый | `#FF8C00` |
| 4 | 🟢 Зелёный | `#00C800` |
| 5 | 🌲 Тёмно-зелёный | `#008000` |
| 6 | 🔵 Синий | `#0000CD` |
| 7 | 🟣 Фиолетовый | `#8B008B` |

---

## 7. CTCF-валидация — КРИТИЧЕСКИ ВАЖНО

**CTCF ChIP-seq — биологический ground truth.** CTCF совместно с когезином
формирует петли ДНК, заякоривающие границы TAD. Хороший детектор TAD ОБЯЗАН
давать обогащение CTCF на границах.

**Текущий статус:** нет чёткого пика у 4 классических алгоритмов →
**это центральный аргумент диплома**.

### Статус биологических валидаций (chr1 @ 100kb)

| Приоритет | Валидация | CSV | PNG |
|-----------|-----------|-----|-----|
| ✅ P1 | CTCF профиль ±500kb | `results/stats/ctcf_enrichment.csv` | ✅ `ctcf_profile_all_algos_chr1_*` |
| ✅ P2 | RAD21 профиль | `results/stats/rad21_enrichment.csv` | ✅ `chipseq_profiles/rad21_profile_all_algos_chr1_100000bp.png` |
| ✅ P2 | SMC3 профиль | `results/stats/smc3_enrichment.csv` | ✅ `chipseq_profiles/smc3_profile_all_algos_chr1_100000bp.png` |
| ✅ P4 | H3K4me3 профиль | `results/stats/h3k4me3_enrichment.csv` | ✅ `chipseq_profiles/h3k4me3_profile_all_algos_chr1_100000bp.png` |
| ✅ P4 | H3K27ac профиль | `results/stats/h3k27ac_enrichment.csv` | ✅ `chipseq_profiles/h3k27ac_profile_all_algos_chr1_100000bp.png` |
| 🟡 P3 | TSS / housekeeping genes | — | — |
| 🟢 P5 | Alu/SINE обогащение | — | — |

**Сводная таблица:** `results/stats/chipseq_validation_summary.csv`

### Ключевые результаты валидации (chr1 @ 100kb)

**RAD21 / SMC3 (когезин — специфичный маркер петлевых якорей):**
| Алгоритм | RAD21 | SMC3 |
|----------|-------|------|
| Arrowhead | STRONG (ctr=1.46) | STRONG (ctr=1.25) |
| ontad | STRONG (ctr=1.42) | STRONG (ctr=1.19) |
| modularity_tad | STRONG (ctr=1.31) | STRONG (ctr=1.10) |
| dihmm | GOOD (ctr=1.03) | SHIFTED (+475kb) |
| armatus | WEAK (ctr=1.15) | WEAK (ctr=1.04) |
| scktld | WEAK (ctr=1.02) | WEAK (ctr=1.06) |
| topdom | SHIFTED (+355kb) | SHIFTED (+355kb) |
| coitad | NO_DATA (15 бнд) | NO_DATA (15 бнд) |

**H3K4me3 / H3K27ac (активационные марки):**
- Сигнал слабее RAD21/SMC3 — ожидаемо (энхансеры/промоторы внутри TAD, не на границах)
- STRONG: только Arrowhead
- GOOD: modularity_tad, armatus (H3K4me3)
- topdom: SHIFTED +155–355kb на всех треках (систематический артефакт)

### Скрипт валидации
`scripts/run_chipseq_validation.py` — переиспользует `run_ctcf_profile_analysis()`
из `src/ctcf_analysis.py`. Параметр `ctcf_df` принимает любой ChIP-seq трек.

### Ключевые функции src/validation.py
- `compute_ctcf_profile(df, chipseq_df, chrom, resolution, ...)` → `(bins, density)`
- `compute_ctcf_enrichment(...)` → статистика обогащения
- ⚠️ `plot_ctcf_profile` **НЕ СУЩЕСТВУЕТ** — строить фигуры через matplotlib напрямую

**Подробности → `agent_docs/03_ctcf_validation.md`**

---

## 8. Правила кода

```python
# ❌ ЗАПРЕЩЕНО
# Хардкодить пути / параметры / seed
# print() в src/ и pipeline/ (кроме _get_device())
# Нормализовать матрицу для scKTLD (balance: false ОБЯЗАТЕЛЬНО)
# Менять CONSENSUS_COLORS или ALGO_COLORS
# Возвращать None из run_<algorithm>
# Менять seed (везде 42)
# Коммитить data/ и results/
# Отключать #include <curl/curl.h> в straw.cpp OnTAD
# Использовать make без явных -I${CONDA_PREFIX}/include (OnTAD)

# ✅ ОБЯЗАТЕЛЬНО
cfg = load_config("config/config.yaml")
rng = np.random.default_rng(42)
Path(out).parent.mkdir(parents=True, exist_ok=True)
logger = logging.getLogger(__name__)
```

**Цвета алгоритмов (ФИКСИРОВАНЫ):**
```python
ALGO_COLORS = {
    "armatus":        "#1f77b4",
    "topdom":         "#ff7f0e",
    "scktld":         "#2ca02c",
    "coitad":         "#d62728",
    "dihmm":          "#9467bd",
    "ontad":          "#e377c2",
    "modularity_tad": "#17becf",
}
```

---

## 9. Соглашения по именованию

```
results/tads/<algo>_<chrom>_<res>bp.bed                    # TAD-лист
results/consensus/consensus_<chrom>_<res>bp.bed            # консенсус (границы)
results/consensus/tad_consensus_<chrom>_<res>bp.bed        # консенсус (домены, Jaccard)
results/stats/<name>.csv                                   # метрики
results/figures/hic_maps/hic_browser_<chrom>_<res>bp.png      # Hi-C browser-style
results/figures/hic_maps/hic_distance_<chrom>_<res>bp.png     # Hi-C distance-style
results/figures/chipseq_profiles/<track>_profile_all_algos_<chrom>_<res>bp.png  # ChIP-seq профили

chrom = "chr17"    # всегда с префиксом chr
resolution = 25000 # в bp (int)
df.columns == ["chrom", "start", "end"]  # строго
```

---

## 10. Параллельная работа — Workstreams

| Workstream | Ветка | Зона ответственности |
|------------|-------|----------------------|
| WS-1: Algorithms | `ws/algorithms` | `src/algorithms/`, `tools/` |
| WS-2: Visualization | `ws/visualization` | `src/visualization.py`, `results/figures/` |
| WS-3: Validation | `ws/validation` | `src/validation.py`, `src/ctcf_analysis.py`, `src/consensus.py` |
| WS-4: Thesis | `ws/thesis` | `thesis/` |

**Подробности → `agent_docs/09_workstreams.md`**

---

## 11. Research Backlog

### 🔴 P0 — Критично

- [x] ~~CTCF как ground truth: диагностика~~ ✅ v3.2
- [x] ~~Широкий CTCF-профиль ±500kb~~ ✅ v3.2
- [x] ~~compute_tad_consensus() — Jaccard ≥ 0.5~~ ✅ v3.3

### 🟡 P1 — Важно

- [x] ~~**PNG-фигуры профилей ChIP-seq (RAD21/SMC3/H3K4me3/H3K27ac)** ✅ v3.4~~
  В `results/figures/chipseq_profiles/`

- [x] ~~**Изменить визуализацию Hi-C на visualization_sveta.py** ✅ v3.4~~
  browser-style + distance-style. Jaccard удалён из pipeline.

- [x] ~~**tad_consensus по всем хромосомам** ✅ v3.4~~
  22 файла @ 100kb в `results/consensus/`. chr1: 142 домена.

- [x] ~~**Прогон @ 50kb и 100kb, chr1–chr22** ✅ v3.4~~
  TAD-файлы уже существовали; консенсус пересчитан напрямую.
  ```bash
  python pipeline/run_pipeline.py --resolution 50000 \
      --chroms chr1 ... chr22 \
      --algorithms armatus topdom scktld coitad dihmm ontad modularity_tad --force
  ```
  Исключения: coitad @ 100kb (unusable), scKTLD только chr17–chr22 @ 25kb.

- [x] ~~**Наглядная overlay-визуализация границ** ✅ v3.5~~
  weak consensus (dashed lines) + strong consensus (filled axvspan)
  поверх Hi-C ленты в `plot_tad_browser_view`. → `src/visualization.py`

### 🟢 P2 — Развитие

- [ ] **Написание диплома — глава по биол. валидации:**
  Интерпретировать ChIP-seq результаты. Ключевые тезисы:
  ontad + modularity_tad STRONG по когезину; topdom систематич. смещение +355kb;
  H3K4me3/H3K27ac слабее RAD21/SMC3 (ожидаемо).
  → `agent_docs/06_thesis_structure.md`

- [ ] **TSS / housekeeping genes (P3):** GENCODE hg19 GTF
- [ ] **Alu/SINE обогащение (P5):** UCSC RepeatMasker
- [ ] **Статистические методы:** Фурье, спектр, change-point → `agent_docs/10_statistical_methods.md`
- [ ] **GNN на консенсусных границах:** proof-of-concept

---

## 12. Структура диплома

| # | Глава | Суть | Статус |
|---|-------|------|--------|
| 1 | Введение + мотивация | Hi-C, TAD, зачем нужен качественный датасет | 🔲 |
| 2 | Обзор техник | 7 алгоритмов + deepTAD | 🔲 |
| 3 | Критика консенсуса | Почему deepTAD-style плохой | 🔲 |
| 4 | Механика алгоритмов | Почему разные алгоритмы → разные результаты | 🔲 |
| 5 | Новый пайплайн | Статистические методы + CTCF-валидация | 🔲 |
| 6 | Результаты | ChIP-seq валидация, сравнение алгоритмов | 🔲 |
| 7 | Заключение | Вклад, публикация датасета | 🔲 |

**Подробный план → `agent_docs/06_thesis_structure.md`**

---

## 13. Зависимости

```
numpy>=1.24, scipy>=1.11, pandas>=2.0, scikit-learn>=1.3
matplotlib>=3.7, seaborn>=0.12, plotly>=5.18
PyYAML>=6.0, tqdm>=4.66, jinja2>=3.1, click>=8.1
hmmlearn>=0.3.3
ruptures>=1.1.7
cooler>=0.9.3, pybedtools>=0.9.1, pyranges>=0.0.129
pytest>=7.4, pytest-cov>=4.1
torch (CUDA 12.1)
```

**Внешние бинарники:**
- Armatus: `tools/armatus/armatus`
- OnTAD: `tools/ontad/OnTAD` (сборка → `agent_docs/04_build_ontad.md`)

---

## 14. Типичные задачи → куда смотреть

| Задача | Файл |
|--------|------|
| Добавить новый алгоритм | `src/algorithms/run_<new>.py` + `ALGORITHM_REGISTRY` + §5 |
| Изменить параметры | `config/config.yaml` |
| Изменить консенсус границ | `src/consensus.py` |
| Изменить TAD-консенсус (Jaccard) | `src/consensus.py::compute_tad_consensus()` |
| Добавить метрику | `src/statistics.py` |
| Изменить CTCF/ChIP-seq анализ | `src/ctcf_analysis.py` |
| Запустить ChIP-seq валидацию | `scripts/run_chipseq_validation.py` |
| Добавить шаг в пайплайн | `pipeline/run_pipeline.py` |
| Добавить/изменить график | `src/visualization.py` → `agent_docs/07_visualization.md` |
| Проблема с компиляцией OnTAD | `agent_docs/04_build_ontad.md` |
| Что-то сломалось | `agent_docs/05_known_issues.md` |
| Параллельная работа | `agent_docs/09_workstreams.md` |

---

## 15. Частые подводные камни

| Проблема | Причина | Решение |
|----------|---------|---------|
| OOM при scKTLD | Крупные хромосомы | Проверять `scktld_limits` в config |
| Armatus 0 TADs | Флаг `-R` вместо `-S -N` | `-S -N -c <chrom>` |
| CTCF нет пика | Координаты / окно / chr-prefix | Диагностика → `agent_docs/03_ctcf_validation.md` |
| OnTAD не компилируется: `curl/curl.h` | libcurl не в PATH | `conda install libcurl` + явный `-I${CONDA_PREFIX}/include` |
| OnTAD не компилируется: `uint64_t` | GCC 15 строже | `#include <cstdint>` после `#include "straw.h"` |
| OnTAD: 0 TADs | Парсер depth: `parts[3]` | `depth = parts[2]` |
| ModularityTAD: 300+ TADs | penalty из config переопределяет auto | `penalty: null` в config |
| logger NameError в `_get_device()` | До init logging | `print()` в этой функции |
| HTML-визуализация зависает | Plotly chr1@25kb | `generate_html: false` |
| scKTLD: balance=True | Нормализация ломает структуру | `balance: false` ВСЕГДА |
| coitad @ 100kb: 16 TADs | IS-fallback огрубляет при 100kb | Исключить из консенсуса @ 100kb |
| armatus @ 100kb: 200+ TADs | gamma sweep выбирает малый gamma | Известная проблема, зафиксирована |
| chr22 CTCF-профиль: пик смещён | Активный регион только 35Mb, CTCF неравномерно | Не использовать chr22. Эталон: chr1@100kb |
| CTCF-профиль на малых хромосомах | chr21/22 дают ложные сдвиги | Минимум: chr1–chr6 @ ≥50kb |
| CTCF-валидация зависает (50kb) | Python for-loop по 44k пиков | `_count_ctcf_overlaps` векторизован в v3.2 |
| `from __future__` не первой строкой | `import os` вставился выше | Держать `from __future__ import annotations` строкой №1 |
| Arrowhead: `invalid literal 'x1'` | Заголовок chr1/x1/x2 не числа | Читать с `header=0`, переименовывать колонки |
| `coitad_..._Xbp.bed` в results/tads/ | Битое имя файла от артефакта glob-паттерна | Игнорировать, в консенсус не попадает |
| topdom: систематический сдвиг +355kb | Артефакт алгоритма (воспроизводится на всех 4 треках) | Отразить в дипломе как ограничение topdom |
| coitad ChIP-seq: NO_DATA | 15 границ < min_boundaries=30 на chr1@100kb | Ограничение разрешения, отразить в дипломе |
| SMC3 из wgEncodeRegTfbsClustered | col4 = имя TF только для основных TF, SMC3 — в других колонках | Использовать ENCODE Portal: ENCFF001VFH |

---

## 16. Как обновлять этот файл

> docs-first: rules.md обновляется ДО изменения кода.
> Держать файл ≤ 300 строк — детали выносить в agent_docs/.

**Чеклист: добавить новый алгоритм:**
- [ ] §5: строка в таблице реестра
- [ ] §8: цвет в ALGO_COLORS
- [ ] §16: Changelog
- [ ] `agent_docs/01_algorithms.md`
- [ ] `src/algorithms/run_<algo>.py` + `ALGORITHM_REGISTRY`
- [ ] `config/config.yaml`
- [ ] `src/visualization.py` (ALGO_COLORS)
- [ ] `tests/test_<algo>.py`

---

## 17. Changelog

### v3.5 — 2026-05-28

**Добавлено:**
- §6: `compute_strong_tad_consensus()` — строгий TAD-консенсус:
  TAD включается только если |start_i−start_j| ≤ tol·res AND |end_i−end_j| ≤ tol·res
  у ≥ min_support алгоритмов. Файлы: `results/consensus/strong_tad_consensus_{chrom}_{res}bp.bed`
- §6: `config/consensus.strong_boundary_tolerance_bins: 1` добавлен в config/config.yaml
- §7 (визуализация): `plot_tad_browser_view` расширен двумя слоями:
  weak consensus → вертикальные dashed-линии; strong consensus → заполненные axvspan
- §7: `_load_consensus_bed()` — хелпер чтения BED с/без заголовка (поддержка обоих форматов)
- §11: Backlog P1 "наглядная overlay-визуализация" ✅ закрыта
- §15: новые known issues (Arrowhead заголовок, `from __future__` позиция)
- `scripts/run_chipseq_validation.py`: `run_consensus_chipseq_validation()`,
  `_boundaries_from_domain_bed()` — прогон strong/weak консенсусов через ChIP-seq профили
- Summary обновлён: 8 строк (strong_consensus + weak_consensus × 4 трека)
- 4 сравнительных PNG: `{track}_strong_vs_weak_consensus_chr1_100000bp.png`

**Ключевые результаты (chr1 @ 100kb):**
| Трек    | strong_consensus | weak_consensus | Arrowhead |
|---------|-----------------|----------------|-----------|
| RAD21   | STRONG 1.382    | STRONG 1.317   | 1.457     |
| SMC3    | GOOD 1.187      | GOOD 1.230     | 1.252     |
| H3K4me3 | GOOD 1.218      | GOOD 1.210     | 1.142     |
| H3K27ac | GOOD 1.107      | GOOD 1.140     | 1.116     |

Strong > weak по RAD21 (cohesin-специфичный маркер) — строгий критерий
границ захватывает биологически реальные loop-якоря точнее Jaccard-перекрытия.

**Затронутые файлы:**
`src/consensus.py`, `src/visualization.py`, `config/config.yaml`,
`scripts/run_chipseq_validation.py`, `rules.md`

### v3.4 — 2026-05-20

**Добавлено:**
- §3: `src/visualization.py` заменён на `visualization_sveta.py` (browser-style + distance-style)
- §4: `hic_figures_out`, `chipseq_figures_out`, `arrowhead_ref` добавлены в config/config.yaml
- §6: `compute_tad_consensus()` реально интегрирован в `compute_all_consensus()` (был задекларирован в v3.3, но не вызывался)
- §7: PNG-профили RAD21/SMC3/H3K4me3/H3K27ac сгенерированы → `results/figures/chipseq_profiles/`
- §9: структура `results/figures/` разделена на `hic_maps/` + `chipseq_profiles/`
- §11: Backlog P1-задачи закрыты
- §15: добавлен known issue `coitad_..._Xbp.bed`
- `config/config.yaml`: `styles: [browser, distance]`, `max_distance_bp`, `panel_bp`, все 7 цветов

**tad_consensus @ 100kb:** 22 файла, chr1–chr22 (coitad исключён)
**tad_consensus @ 50kb:** 22 файла, chr1–chr22 (coitad исключён); chr1: 177 доменов
**Структура figures/:** `hic_maps/` (69 PNG), `chipseq_profiles/` (417 PNG + 4 strong_vs_weak)
**Jaccard heatmap:** удалён из `run_all_visualization` pipeline

**Затронутые файлы:**
`src/consensus.py`, `src/visualization.py`, `config/config.yaml`,
`scripts/run_chipseq_validation.py`, `rules.md`

### v3.3 — 2026-05-19

**Добавлено:**
- §4: ChIP-seq треки RAD21 (ENCFF001VFE, 23 947 пиков), SMC3 (ENCFF001VFH, 64 597),
  H3K4me3 (57 476), H3K27ac (56 069) — скачаны в `data/reference/`
- §4: пути rad21/smc3/h3k4me3/h3k27ac_bed добавлены в config/config.yaml
- §6: CONSENSUS_COLORS расширен до 7 уровней (был до 4); исправлен `min(support,4)→min(support,7)`
- §6: `consensus.jaccard_threshold: 0.5` добавлен в config/config.yaml
- §6: `compute_tad_consensus()` (Union-Find + Jaccard) интегрирован в `compute_all_consensus()`
- §7: таблица валидаций обновлена; добавлены результаты по всем 4 трекам
- §7: сводная таблица `results/stats/chipseq_validation_summary.csv`
- §14: добавлена строка `scripts/run_chipseq_validation.py`
- §15: новые known issues — PNG-фигуры, topdom +355kb, coitad NO_DATA, SMC3 clustered
- Новый скрипт: `scripts/run_chipseq_validation.py`
- Новый модуль: `src/ctcf_analysis.py` (параметризован под любой ChIP-seq трек)

**Smoke-тест consensus.py пройден:** `compute_consensus` + `compute_tad_consensus` ✅

**Затронутые файлы:**
`src/consensus.py`, `src/ctcf_analysis.py`, `config/config.yaml`,
`scripts/run_chipseq_validation.py`, `rules.md`

### v3.2 — 2026-05-15

**Исправлено:**
- `run_modularity_tad.py` v2.2: O/E нормализация вместо log1p
- `run_ontad.py` v2.3: recursive top-down subdivision
- `pipeline/run_pipeline.py`: `--algorithms` динамический из ALGORITHM_REGISTRY
- `src/validation.py`: `_count_ctcf_overlaps` векторизован (~100x ускорение)
- `config/config.yaml`: `modularity_tad.penalty: null`

**Результаты полного прогона (покрытие):**
| Алгоритм | 25kb | 50kb | 100kb |
|----------|------|------|-------|
| armatus | 23 chr, avg 140 | 23 chr, avg 151 ⚠️ | 23 chr, avg 202 ⚠️ |
| topdom | 23 chr, avg 180 | 23 chr, avg 96 | 23 chr, avg 63 |
| scktld | 6 chr, avg 46 | 22 chr, avg 52 | 22 chr, avg 56 |
| coitad | 23 chr, avg 63 | 23 chr, avg 32 | 23 chr, avg 16 ⚠️ |
| dihmm | 6 chr, avg 73 | 22 chr, avg 69 | 22 chr, avg 53 |
| ontad | 6 chr, avg 46 | 22 chr, avg 65 | 22 chr, avg 33 |
| modularity_tad | 6 chr, avg 55 | 22 chr, avg 113 | 22 chr, avg 116 |

### v3.1 — 2026-05-14
agent_docs/ создана (10 файлов), workstreams, research backlog, структура диплома.

### v2.0–v1.3
*(см. историю git)*

---

*rules.md v3.4 | TAD Consensus Pipeline | hg19 | GSE63525 GM12878*
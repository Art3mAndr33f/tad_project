# TAD Consensus Pipeline — Project Rules for LLM Context
# Version: 3.1 | 2026-05-14 | Genome: hg19 | Data: GSE63525 GM12878

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
**GPU:** `CUDA_VISIBLE_DEVICES=3` (brain-lab, ~40 GB свободно на GPU 3)

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
cat agent_docs/09_workstreams.md   # найди свой workstream
cat agent_docs/<твой_файл>.md      # загрузи детали
git checkout -b ws/<workstream>    # создай свою ветку
```

---

## 3. Архитектура (карта кода)

```
src/algorithms/          ← 7 детекторов TAD (единый интерфейс → §5)
src/consensus.py         ← жадная кластеризация границ
src/validation.py        ← CTCF ChIP-seq валидация (биол. ground truth)
src/statistics.py        ← Jaccard, boundary overlap, сравнение с Arrowhead
src/visualization.py     ← Hi-C heatmap, CTCF профили, Jaccard matrix
src/data_prep.py         ← загрузка матриц (4-уровневый fallback)
pipeline/run_pipeline.py ← оркестратор
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
| CTCF ChIP-seq | `data/reference/GM12878_CTCF_peaks_hg19.bed` | ENCODE ENCFF796WRU, hg19 |

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

**Алгоритм:** Собрать все границы (start/end) → округлить до бина →
жадная кластеризация с `tolerance_bins=1` → support = число алгоритмов в кластере.

```yaml
consensus:
  tolerance_bins: 1
  min_support: 2
```

**Цветовая схема (ФИКСИРОВАНА — не менять):**

| Support | Цвет | Hex |
|---------|------|-----|
| 2 | 🟡 Жёлтый | `#FFD700` |
| 3 | 🟠 Оранжевый | `#FF8C00` |
| ≥4 | 🟢 Зелёный | `#00C800` |

---

## 7. CTCF-валидация — КРИТИЧЕСКИ ВАЖНО

**CTCF ChIP-seq — биологический ground truth.** Это не просто метрика.
CTCF совместно с когезином формирует петли ДНК, заякоривающие границы TAD.
Хороший детектор TAD ОБЯЗАН давать обогащение CTCF на границах.
Если алгоритм даёт границы без пика CTCF → его разметка биологически нерелевантна
→ такой консенсус не годится для обучения нейросети.

**Текущий статус:** нет чёткого пика у 4 классических алгоритмов →
**это центральный аргумент диплома**.

**⚠️ ACTIVE ISSUE — диагностировать в первую очередь:**
```bash
head -5 data/reference/GM12878_CTCF_peaks_hg19.bed  # формат
cut -f1 data/reference/GM12878_CTCF_peaks_hg19.bed | sort -u  # chr-prefix?
# Ожидаем: chr1, chr2, ... — если без chr → fix в src/validation.py
```

**Метрика успеха (для нового пайплайна):**
Распределение смещений CTCF от границ → ближе к δ(0) (узкий высокий пик).
Текущий консенсус → размытое распределение → плохой датасет.

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

# ✅ При ошибке алгоритма
try:
    df = run_algorithm(...)
except Exception as exc:
    logger.error("[algo] %s @ %d: %s", chrom, resolution, exc, exc_info=True)
    df = pd.DataFrame(columns=["chrom", "start", "end"])
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
results/tads/<algo>_<chrom>_<res>bp.bed          # TAD-лист
results/consensus/consensus_<chrom>_<res>bp.bed  # консенсус
results/stats/<name>.csv                         # метрики
results/figures/hic_tads_<chrom>_<res>bp.png     # фигуры
data/processed/<chrom>_<res>bp.npy               # матрица

chrom = "chr17"    # всегда с префиксом chr
resolution = 25000 # в bp (int)
df.columns == ["chrom", "start", "end"]  # строго
```

---

## 10. Параллельная работа — Workstreams

Проект разбит на **изолированные workstreams** с чёткими границами:

| Workstream | Ветка | Зона ответственности |
|------------|-------|----------------------|
| WS-1: Algorithms | `ws/algorithms` | `src/algorithms/`, `tools/` |
| WS-2: Visualization | `ws/visualization` | `src/visualization.py`, `results/figures/` |
| WS-3: Validation | `ws/validation` | `src/validation.py`, `src/statistics.py`, `src/consensus.py` |
| WS-4: Thesis | `ws/thesis` | `thesis/`, `notebooks/` |

**Читают ВСЕ (только чтение):** `config/config.yaml`, `data/`, `src/data_prep.py`
**Меняют по договорённости:** интерфейсы между WS (API-контракты)

**Подробности → `agent_docs/09_workstreams.md`**

---

## 11. Research Backlog

> Полная документация по задачам → соответствующие `agent_docs/` файлы.

### 🔴 P0 — Критично (блокирует основной результат)

- [ ] **CTCF как ground truth:** Подтвердить биологическую информированность.
  Диагностировать отсутствие пика (координаты, chr-prefix, ширина окна).
  → `agent_docs/03_ctcf_validation.md`

- [ ] **Широкий CTCF-профиль ±500kb:** Заменить текущий ±1bin на профиль
  плотности CTCF в окне ±500kb с binning 10kb. Показать на одном графике
  все алгоритмы + Arrowhead эталон.
  → `agent_docs/07_visualization.md` §P1

- [ ] **Воспроизвести результаты алгоритмов из оригинальных статей:**
  Armatus (Filippova 2014), TopDom (Shin 2016), OnTAD (An 2019),
  scKTLD (Zheng 2024). Проверить на разных разрешениях.
  → `agent_docs/01_algorithms.md`

- [ ] **Гиперпараметры — обоснование:** Перепроверить, почему выбирается
  конкретное число TAD для конкретной хромосомы. Правильно ли использовать
  единый диапазон TAD/Mb для всех хромосом и разрешений?
  → `agent_docs/01_algorithms.md` §OPEN QUESTIONS

### 🟡 P1 — Важно (формирует основной аргумент)

- [ ] **Наглядная визуализация границ:** Overlay алгоритмов и консенсуса
  на Hi-C heatmap. Два режима: per_algo (tracks) и consensus (weighted lines).
  Zoom на конкретные локусы.
  → `agent_docs/07_visualization.md` §P2

- [ ] **Консенсус по TAD, а не по границам:** Добавить поиск консенсусных
  границ не по отдельным границам, а по пересечениям целых TAD.
  Проверить реализацию в deepTAD.
  → `agent_docs/10_statistical_methods.md` §Метод 3

- [ ] **Биологическая валидация как ground truth:** Добавить проверку
  консенсусных границ на CTCF. Показать, что текущий консенсус не проходит.
  Ответить: можно ли использовать CTCF как единственный критерий датасета?
  → `agent_docs/03_ctcf_validation.md` §P3

- [ ] **Запустить прогон @ 50kb и 100kb:** Полный прогон chr1–chr22 для сбора
  данных к главам 3 и 6 диплома.

### 🟢 P2 — Развитие (основной научный вклад)

- [ ] **Статистические методы поиска границ:** Фурье-анализ insulation score,
  спектральные методы (Fiedler vector), change-point detection.
  → `agent_docs/10_statistical_methods.md`

- [ ] **Доказать превосходство нового датасета:** Количественно показать,
  что новый пайплайн даёт более узкий/высокий пик в CTCF профиле vs консенсус.
  → `agent_docs/03_ctcf_validation.md` §P2

- [ ] **GNN на консенсусных границах:** Небольшая графовая нейросеть для
  детекции TAD на полученных консенсусных границах как proof-of-concept.

- [ ] **Написание диплома:** Начать с глав 3 (критика консенсуса) и 6
  (результаты) — данные уже есть.
  → `agent_docs/06_thesis_structure.md`

---

## 12. Структура диплома

| # | Глава | Суть | Статус |
|---|-------|------|--------|
| 1 | Введение + мотивация | Hi-C, TAD, зачем нужен качественный датасет | 🔲 |
| 2 | Обзор техник | 7 алгоритмов + deepTAD, метод сбора датасета | 🔲 |
| 3 | Критика консенсуса | Почему deepTAD-style плохой (≥10 стр, много визуализаций) | 🔲 |
| 4 | Механика алгоритмов | Почему разные алгоритмы → разные результаты | 🔲 |
| 5 | Новый пайплайн | Статистические методы + CTCF-валидация | 🔲 |
| 6 | Результаты | Новый датасет проходит CTCF-валидацию — main result | 🔲 |
| 7 | Заключение | Вклад, публикация датасета, следующий шаг (нейросеть) | 🔲 |

**Подробный план каждой главы → `agent_docs/06_thesis_structure.md`**

---

## 13. Зависимости

```
numpy>=1.24, scipy>=1.11, pandas>=2.0, scikit-learn>=1.3
matplotlib>=3.7, seaborn>=0.12, plotly>=5.18
PyYAML>=6.0, tqdm>=4.66, jinja2>=3.1, click>=8.1
hmmlearn>=0.3.3
cooler>=0.9.3, pybedtools>=0.9.1, pyranges>=0.0.129
pytest>=7.4, pytest-cov>=4.1
torch (CUDA 12.1)  # pip install torch --index-url https://download.pytorch.org/whl/cu121
```

**Внешние бинарники:**
- Armatus: `tools/armatus/armatus`
- OnTAD: `tools/ontad/OnTAD` (сборка → `agent_docs/04_build_ontad.md`)

---

## 14. Типичные задачи → куда смотреть

| Задача | Файл |
|--------|------|
| Добавить новый алгоритм | `src/algorithms/run_<new>.py` + `ALGORITHM_REGISTRY` + §5 этого файла |
| Изменить параметры | `config/config.yaml` |
| Изменить консенсус | `src/consensus.py` |
| Добавить метрику | `src/statistics.py` |
| Изменить CTCF-анализ | `src/validation.py` |
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
| coitad @ 100kb: 16 TADs, median 12Mb | IS-fallback огрубляет при 100kb | Исключить из консенсуса @ 100kb |
| armatus @ 100kb: 200+ TADs | gamma sweep выбирает слишком малый gamma | Известная проблема, зафиксирована |
| CTCF-валидация зависает (50kb) | Python for-loop по 44k пиков × 1000 перм | `_count_ctcf_overlaps` векторизован в v3.2 |

---

## 16. Как обновлять этот файл

> docs-first: rules.md обновляется ДО изменения кода.
> Держать файл ≤ 250 строк — детали выносить в agent_docs/.

**Чеклист: добавить новый алгоритм:**
- [ ] §5: добавить строку в таблицу реестра
- [ ] §8: добавить цвет в ALGO_COLORS
- [ ] §16: Changelog
- [ ] `agent_docs/01_algorithms.md`: подраздел с механикой и параметрами
- [ ] `src/algorithms/run_<algo>.py`
- [ ] `src/algorithms/__init__.py` (ALGORITHM_REGISTRY)
- [ ] `config/config.yaml` (блок параметров)
- [ ] `src/visualization.py` (ALGO_COLORS)
- [ ] `requirements.txt` / `environment.yml`
- [ ] `tests/test_<algo>.py`

---

## 17. Changelog

### v3.2 — 2026-05-15

**Исправлено:**
- `run_modularity_tad.py` v2.2: O/E нормализация (observed/expected по диагоналям)
  вместо log1p — полностью убирает distance-decay без bias к размеру TAD;
  переписан `_dp_segment` (убран pass-through механизм, теперь только реальные сегменты)
- `run_ontad.py` v2.3: recursive top-down subdivision вместо score-based greedy —
  рекурсивно заменяет крупные TAD детьми до достижения target; исправлена опечатка
  `best_dist = best_rows` → `best_dist = dist`
- `pipeline/run_pipeline.py`: `--algorithms` теперь принимает все 7 алгоритмов
  (убран `choices=[...]` захардкоженный список, заменён на динамический из ALGORITHM_REGISTRY)
- `src/validation.py`: `_count_ctcf_overlaps` векторизован через numpy broadcasting
  (~100x ускорение; 30 мин → 2–5 мин для 50kb chr1–chr22)
- `config/config.yaml`: `modularity_tad.penalty: null` (был 0.05 — блокировал auto-sweep)

**Добавлено:**
- Документация known issues: coitad @ 100kb (16 TADs, median 12Mb — unusable),
  armatus @ 100kb (avg 202 TADs — завышено); оба зафиксированы в §15

**GPU:** GPU2 стал предпочтительным (20GB свободно); GPU3 — резервный

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

**Затронутые файлы:**
`src/algorithms/run_ontad.py`, `src/algorithms/run_modularity_tad.py`,
`src/validation.py`, `pipeline/run_pipeline.py`, `config/config.yaml`,
`rules.md`, `README.md`

### v3.1 — 2026-05-14

**Добавлено:**
- `agent_docs/` — 10 файлов подробной документации по темам
- `create_agent_docs.py` — скрипт автоматического создания agent_docs/
- §10: Workstreams (карта параллельной работы с ветками git)
- §11: Research Backlog (все задачи из комментариев тимлида 12.05.2026)
- §12: Структура диплома (7 глав с описанием)
- `thesis/` директория для текста диплома
- `logs/` директория для логов запусков

**Изменено:**
- rules.md реструктурирован по принципу Progressive Disclosure (~250 строк)
- Детали алгоритмов, данных, CTCF, OnTAD-сборки вынесены в agent_docs/
- §7 CTCF: расширен биологический контекст и диагностика
- README.md: добавлены разделы для новых участников, workstreams, статус

**Задачи из комментариев тимлида (12.05.2026) — добавлены в Backlog:**
- Широкий CTCF-профиль ±500kb
- Наглядная визуализация границ (overlay + zoom)
- Консенсус по TAD-пересечениям (deepTAD review)
- Биол. валидация как ground truth
- Проверка консенсуса на CTCF
- Воспроизведение результатов алгоритмов из статей
- GNN на консенсусных границах
- Статистические методы (Фурье, спектральный, change-point)
- Обоснование гиперпараметров по хромосомам
- Доказательство превосходства нового датасета через CTCF

**Затронутые файлы:**
`rules.md`, `README.md`, `agent_docs/` (создана), `create_agent_docs.py` (создан),
`thesis/` (создана), `logs/` (создана)

### v2.1 — 2026-05-13

**Исправлено:**
- OnTAD: `depth==1` → `_get_max_nonoverlapping()` (листовые TAD из иерархии)
- ModularityTAD: `mean(B)` → `intra/inter ratio` objective
- Auto-penalty: нормировка на `percentile_95(scores)` вместо `max`

**Затронутые файлы:**
`src/algorithms/run_ontad.py`, `src/algorithms/run_modularity_tad.py`, `rules.md`

### v2.0 — 2026-05-12

**Добавлено:** DI+HMM, OnTAD, ModularityTAD (3 новых алгоритма).
ALGORITHM_REGISTRY: 4 → 7 алгоритмов. Инструкция сборки OnTAD.

### v1.3 — 2026-05-03

*(см. историю git)*

---

*rules.md v3.1 | TAD Consensus Pipeline | hg19 | GSE63525 GM12878*
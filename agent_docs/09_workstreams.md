# Workstreams — параллельная работа

> Читать ВСЕМ участникам перед началом работы.
> Последнее обновление: 2026-05-14

---

## Принцип параллельной работы

Проект разбит на **изолированные workstreams** с чёткими границами ответственности.
Каждый workstream работает в своей git-ветке и трогает только свои файлы.

**Правило:** если тебе нужно изменить файл из чужого workstream → сначала
обсудить с владельцем, создать issue, договориться об интерфейсе.

---

## Карта workstreams

```
┌─────────────────────────────────────────────────────────────────────┐
│                         TAD Pipeline Project                        │
├──────────────┬──────────────┬───────────────┬───────────────────────┤
│  WS-1: ALGO  │  WS-2: VIZ   │  WS-3: VALID  │  WS-4: THESIS         │
│  Алгоритмы   │ Визуализация │  Валидация    │  Диплом               │
│              │              │  + Статистика │                       │
├──────────────┼──────────────┼───────────────┼───────────────────────┤
│ src/algo/    │src/visuali-  │src/validation │notebooks/             │
│ tools/       │zation.py     │.py            │thesis/                │
│              │results/fig/  │src/statistics │                       │
│              │              │.py            │                       │
└──────────────┴──────────────┴───────────────┴───────────────────────┘
             ↓ ЧИТАЮТ ВСЕ ↓
       config/config.yaml (только чтение)
       data/ (только чтение)
       src/data_prep.py (только чтение)
```

---

## WS-1: Алгоритмы (Algorithms)

**Git-ветка:** `ws/algorithms`
**Документация:** `agent_docs/01_algorithms.md`

**Зона ответственности:**
- `src/algorithms/run_*.py`
- `tools/armatus/`, `tools/ontad/`, `tools/coiTAD/`
- `tests/test_*.py` (для алгоритмов)

**Текущие задачи:**
- [ ] Воспроизвести результаты из оригинальных статей (Armatus, TopDom, OnTAD, scKTLD)
- [ ] Проверить качество алгоритмов на разных разрешениях (25/50/100kb)
- [ ] Обосновать выбор гиперпараметров для каждой хромосомы
- [ ] Запустить полный прогон chr1–chr22 @ 50kb, 100kb

**НЕ трогать:**
- `src/consensus.py` — принадлежит WS-3
- `src/visualization.py` — принадлежит WS-2
- `config/config.yaml` — только через PR с описанием

**Интерфейс (контракт):**
```python
# Любой run_<algo>() обязан возвращать:
pd.DataFrame(columns=["chrom", "start", "end"])
# Любое другое поведение — баг, блокирующий всех
```

---

## WS-2: Визуализация (Visualization)

**Git-ветка:** `ws/visualization`
**Документация:** `agent_docs/07_visualization.md`

**Зона ответственности:**
- `src/visualization.py`
- `results/figures/`
- `notebooks/exploration.ipynb`

**Текущие задачи:**
- [ ] Широкий CTCF-профиль ±500kb (P1, детали в 07_visualization.md)
- [ ] Улучшенный overlay границ на Hi-C heatmap (P2)
- [ ] Genomic dashboard (P3)
- [ ] Jaccard heatmap с дендрограммой (P4)

**НЕ трогать:**
- `src/validation.py` — данные для CTCF профиля берутся через её API
- `src/algorithms/` — TAD-листы берутся из `results/tads/` (файлы)
- Менять ALGO_COLORS или CONSENSUS_COLORS — только через rules.md

**Входные данные (читать, не писать):**
```python
# TAD-листы (уже посчитаны WS-1)
pd.read_csv("results/tads/<algo>_<chrom>_<res>bp.bed",
            sep="\t", names=["chrom","start","end"])

# CTCF данные
"data/reference/GM12878_CTCF_peaks_hg19.bed"

# Hi-C матрицы
get_matrix(cfg, chrom, resolution)  # из src/data_prep.py
```

---

## WS-3: Валидация и статистика (Validation + Statistics)

**Git-ветка:** `ws/validation`
**Документация:** `agent_docs/03_ctcf_validation.md`, `agent_docs/10_statistical_methods.md`

**Зона ответственности:**
- `src/validation.py`
- `src/statistics.py`
- `src/consensus.py`

**Текущие задачи:**
- [ ] Диагностика CTCF (chr-prefix, окно, координаты)
- [ ] Консенсус по пересечениям TAD (не по границам)
- [ ] Статистические методы поиска границ (Фурье, спектральный)
- [ ] Метрика качества датасета на основе CTCF профиля

**НЕ трогать:**
- `src/algorithms/` — только читать результаты
- `src/visualization.py` — только заказывать новые функции у WS-2

**API для WS-2 (предоставляет WS-3):**
```python
# WS-3 обязуется предоставить:
def get_ctcf_density_profile(
    boundaries: pd.DataFrame,
    ctcf_bed: str,
    window_bp: int = 500_000,
    bin_size_bp: int = 10_000,
) -> np.ndarray:
    # shape: (n_bins,) — плотность CTCF пиков по окну
    # WS-2 использует это для построения графика
```

---

## WS-4: Диплом (Thesis)

**Git-ветка:** `ws/thesis`
**Документация:** `agent_docs/06_thesis_structure.md`

**Зона ответственности:**
- `thesis/` (создать директорию)
- `notebooks/` (analysis notebooks)

**Текущие задачи:**
- [ ] Начать главу 3 (критика консенсуса) — есть данные
- [ ] Начать главу 2 (обзор алгоритмов) — есть agent_docs/01_algorithms.md
- [ ] Собрать все готовые рисунки из results/figures/ с подписями

**НЕ трогать:**
- Любой src/ файл
- config/config.yaml

**Источники данных (только читать):**
```
results/figures/   — готовые рисунки
results/stats/     — csv с метриками
results/tads/      — TAD-листы для примеров
agent_docs/        — документация для понимания методов
```

---

## Workflow: как вести параллельную работу

### Git workflow

```bash
# 1. Создать свою ветку от main
git checkout main && git pull
git checkout -b ws/visualization

# 2. Работать в своей ветке
git add src/visualization.py
git commit -m "viz: add wide CTCF profile ±500kb"

# 3. Перед merge — обновить rules.md (Changelog)
# 4. PR → code review → merge в main
```

### Как заказать изменение в чужом workstream

1. Создать issue с тегом `cross-ws`
2. Описать: что нужно, зачем, какой интерфейс ожидается
3. Дождаться ОК от владельца workstream
4. После реализации — обновить `09_workstreams.md`

### Как добавить новую зависимость

1. Обсудить в команде
2. Добавить в `requirements.txt` И `environment.yml`
3. Обновить `rules.md` §4.1
4. Commit с тегом `deps:`

---

## Статус workstreams (2026-05-14)

| Workstream | Ветка | Ответственный | Статус |
|------------|-------|---------------|--------|
| WS-1: Algorithms | `ws/algorithms` | — | 🟡 Active |
| WS-2: Visualization | `ws/visualization` | — | 🔲 Not started |
| WS-3: Validation | `ws/validation` | — | 🔴 Blocked (CTCF debug) |
| WS-4: Thesis | `ws/thesis` | — | 🔲 Not started |

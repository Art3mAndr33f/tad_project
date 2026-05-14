# Известные проблемы и решения

> Перед созданием нового issue — убедись, что его нет здесь.

---

## Активные проблемы

### 🔴 CTCF: нет чёткого пика обогащения

**Симптом:** График CTCF enrichment плоский или без выраженного пика у 4 алгоритмов.
**Гипотезы и диагностика:** см. `agent_docs/03_ctcf_validation.md`
**Статус:** В работе

### 🟡 Гиперпараметры: единый диапазон TAD/Mb для всех хромосом

**Симптом:** Armatus использует 0.8–2.5 TAD/Mb для всех хромосом и разрешений.
**Проблема:** chr1 и chr22 биологически разные; одинаковая плотность TAD неправомерна.
**Задача:** Изучить распределение размеров TAD в Arrowhead по хромосомам.
**Статус:** Открыт

### 🟡 coiTAD: частота срабатывания fallback неизвестна

**Симптом:** При падении OI-метода срабатывает IS fallback (= TopDom).
**Задача:** Добавить счётчик fallback в run_coitad.py, логировать в stats.
**Статус:** Открыт

---

## Решённые проблемы

### ✅ OnTAD: depth=1 → мало TADs (10–16 вместо 50+)

**Решение:** `_get_max_nonoverlapping()` — жадный выбор листовых TAD из иерархии.
**Файл:** `src/algorithms/run_ontad.py`
**Исправлено:** v2.1

### ✅ ModularityTAD: 300+ TADs вместо 50–100

**Решение:** Заменить `mean(B)` objective на `intra/inter ratio` + нормировка penalty
на `percentile_95(scores)`.
**Файл:** `src/algorithms/run_modularity_tad.py`
**Исправлено:** v2.1

### ✅ OnTAD: компиляция (uint64_t, libcurl)

**Решение:** Патч `#include <cstdint>` + явный `-I${CONDA_PREFIX}/include`.
**Детали:** `agent_docs/04_build_ontad.md`
**Исправлено:** v2.0

---

## Подводные камни (справочник)

| Проблема | Причина | Решение |
|----------|---------|---------|
| OOM при scKTLD | Крупные хромосомы | Проверять scktld_limits в config |
| Armatus 0 TADs | Флаг `-R` вместо `-S -N` | Использовать `-S -N -c <chrom>` |
| logger NameError в _get_device() | До init logging | Использовать print() |
| HTML-визуализация зависает | Plotly chr1@25kb (огромно) | generate_html: false |
| scKTLD: balance=True сломает | Нормализация меняет структуру | balance: false ОБЯЗАТЕЛЬНО |
| Seed не 42 | Случайная инициализация | Везде rng = np.random.default_rng(42) |

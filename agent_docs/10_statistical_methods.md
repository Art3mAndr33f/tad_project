# Статистические методы поиска границ TAD

> Workstream: WS-3 (Validation)
> Статус: ПЛАНИРОВАНИЕ — методы ещё не реализованы

---

## Контекст: зачем статистические методы?

**Проблема консенсуса:** жадная кластеризация границ не учитывает форму сигнала
в матрице. Алгоритм с bad boundaries "голосует" наравне с хорошим.

**Идея:** вместо консенсуса "кто из алгоритмов согласен" — найти границы напрямую
из матрицы через анализ сигнала. Это биологически более обоснованно.

---

## Метод 1: Фурье-анализ insulation score

**Идея:**
1. Вычислить insulation score IS(i) для всей хромосомы
2. Применить FFT к IS-профилю
3. Отфильтровать шумовые высокочастотные компоненты
4. Найти границы как минимумы восстановленного сигнала

**Преимущество:** автоматическая фильтрация шума без подбора window.

```python
from scipy.fft import fft, ifft, fftfreq

def find_boundaries_fourier(
    matrix: np.ndarray,
    resolution: int,
    cutoff_freq: float = 0.1,   # в единицах 1/bin
    min_tad_kb: float = 100.0,
) -> pd.DataFrame:
    # 1. IS = mean(matrix[i-w:i, i:i+w]) для каждого i
    # 2. FFT(IS) → фильтр → IFFT → IS_smooth
    # 3. Минимумы IS_smooth → кандидаты в границы
    # 4. Фильтр по min_tad_kb
    ...
```

---

## Метод 2: Спектральный анализ матрицы

**Идея:**
1. Вычислить eigenvalues/eigenvectors матрицы контактов (или её Лапласиана)
2. Первый ненулевой eigenvector (Fiedler vector) кодирует глобальное разбиение
3. Нулевые пересечения Fiedler vector → кандидаты в границы TAD

**Связь с scKTLD:** scKTLD использует похожую идею но для кластеризации.
Здесь мы используем eigenvectors напрямую для поиска границ.

```python
from scipy.linalg import eigh

def find_boundaries_spectral(
    matrix: np.ndarray,
    resolution: int,
    n_eigvec: int = 3,          # сколько eigenvectors использовать
    min_tad_kb: float = 100.0,
) -> pd.DataFrame:
    # 1. Нормализовать матрицу → Лапласиан L = D - A
    # 2. n_eigvec наименьших ненулевых eigenvectors
    # 3. Zero crossings → кандидаты в границы
    # 4. Объединить сигналы нескольких eigenvectors
    ...
```

---

## Метод 3: Консенсус по пересечениям TAD (не по границам)

**Идея (из deepTAD review):** вместо кластеризации отдельных границ —
найти геномные регионы, которые *всегда попадают внутрь TAD* у всех алгоритмов,
и *всегда разделяют TAD* у всех алгоритмов.

```python
def consensus_by_tad_overlap(
    tads_by_algo: dict,          # {algo: pd.DataFrame[chrom,start,end]}
    chrom_size: int,
    resolution: int,
    min_support: int = 4,        # минимум алгоритмов согласных на разрыв
) -> pd.DataFrame:
    # 1. Для каждого бина → считаем сколько алгоритмов имеют TAD-границу рядом
    # 2. Бины с support >= min_support → кандидаты в финальные границы
    # 3. В отличие от жадной кластеризации — учитываем СТРУКТУРУ TAD
    ...
```

**Ссылка для изучения:** deepTAD GitHub — проверить как реализован их data collection.

---

## Метод 4: Change-point detection

**Идея:** insulation score = временной ряд. Границы TAD = точки резкого изменения.

```python
# Библиотека: ruptures (pip install ruptures)
import ruptures as rpt

def find_boundaries_changepoint(
    matrix: np.ndarray,
    resolution: int,
    n_bkps: int = None,         # None = автоподбор через BIC/AIC
    min_tad_kb: float = 100.0,
) -> pd.DataFrame:
    is_profile = compute_insulation_score(matrix)
    model = rpt.Pelt(model="rbf").fit(is_profile)
    breakpoints = model.predict(pen=...)
    ...
```

---

## План сравнения методов

| Метод | Биол. интерпретируем | Без гиперпараметров | CTCF пик? |
|-------|---------------------|--------------------|----|
| Жадный консенсус (текущий) | ❌ | ❌ | ❌ плохой |
| Фурье IS | ✅ | ✅ (cutoff_freq) | ? |
| Спектральный | ✅ | ✅ (n_eigvec) | ? |
| TAD-overlap консенсус | ✅ | ❌ (min_support) | ? |
| Change-point | ✅ | ✅ (pen) | ? |

**Метрика сравнения:** ширина пика в CTCF профиле на ±500kb графике.
Чем уже и выше пик — тем лучше метод.

---

## Порядок реализации

1. Реализовать `get_ctcf_density_profile()` (нужно WS-2 для графиков)
2. Реализовать TAD-overlap консенсус (простейший из новых методов)
3. Проверить CTCF профиль → если лучше консенсуса → аргумент для диплома
4. Реализовать Фурье и спектральный методы
5. Сравнить все методы через CTCF профиль

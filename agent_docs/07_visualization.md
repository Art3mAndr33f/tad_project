# Визуализация — руководство для разработчика

> Workstream: Visualization
> Зона ответственности: src/visualization.py, results/figures/
> НЕ трогать: src/algorithms/, src/consensus.py, src/validation.py

---

## Текущее состояние

Реализовано в `src/visualization.py`:
- Hi-C heatmap с overlay границ TAD
- CTCF enrichment profile (базовый)
- Jaccard heatmap между алгоритмами

---

## Фиксированные константы (НЕ МЕНЯТЬ)

```python
# Цвета алгоритмов — ФИКСИРОВАНЫ, менять только через rules.md
ALGO_COLORS = {
    "armatus":        "#1f77b4",
    "topdom":         "#ff7f0e",
    "scktld":         "#2ca02c",
    "coitad":         "#d62728",
    "dihmm":          "#9467bd",
    "ontad":          "#e377c2",
    "modularity_tad": "#17becf",
}

# Цвета консенсуса — ФИКСИРОВАНЫ
CONSENSUS_COLORS = {2: "#FFD700", 3: "#FF8C00", 4: "#00C800"}

# Параметры рендеринга из config.yaml
visualization:
  hic_colormap: "coolwarm"
  hic_log_scale: true
  dpi: 300
  format: "png"
  figsize: [18, 14]
  generate_html: false   # Plotly для chr1@25kb зависает — держать false
```

---

## Задачи (приоритизированы)

### 🔴 P1: Широкий CTCF-профиль (±500kb)

**Цель:** Заменить текущий ±1bin профиль на ±500kb с binning 10kb.

```python
def plot_ctcf_boundary_profile(
    boundaries: pd.DataFrame,        # chrom, start, end
    ctcf_bed_path: str,
    window_bp: int = 500_000,
    bin_size_bp: int = 10_000,
    ax=None,
    label: str = "",
    color: str = "blue",
) -> plt.Axes:
    """
    Строит профиль плотности CTCF пиков относительно границ TAD.
    Ось X: -500kb .. +500kb от границы
    Ось Y: нормализованная плотность (пиков / kb / границу)
    Пик в 0 → алгоритм биологически информирован.
    """
    ...
```

**Где вызывать:** После compute_ctcf_enrichment() в pipeline/run_pipeline.py
**Выходной файл:** `results/figures/ctcf_profile_wide_<chrom>_<res>bp.png`

**Нужно показать на одном графике:**
- Каждый алгоритм отдельной линией (ALGO_COLORS)
- Arrowhead эталон — чёрной пунктирной линией
- Легенда, подпись осей, title с chrom и resolution

---

### 🔴 P2: Улучшенная визуализация границ на Hi-C heatmap

**Цель:** Сделать overlay границ информативным и читаемым.

**Текущая проблема:** При 7 алгоритмах + консенсус = каша линий.

**Решение — два режима:**

**Режим A: По алгоритмам (для сравнения)**
- Каждый алгоритм = отдельная строка треков под heatmap (как genome browser)
- Цвет линии = ALGO_COLORS[algo]
- Высота трека = 0.3 × resolution

**Режим B: Консенсус (для финального результата)**
- На heatmap: вертикальные линии с толщиной пропорциональной support
- Цвет = CONSENSUS_COLORS[support]
- Support=1 → не показывать

```python
def plot_hic_with_boundaries(
    matrix: np.ndarray,
    boundaries_by_algo: dict,        # {algo_name: pd.DataFrame}
    consensus: pd.DataFrame,
    chrom: str,
    resolution: int,
    mode: str = "consensus",         # "per_algo" | "consensus"
    region: tuple = None,            # (start_bp, end_bp) для zoom
    output_path: str = None,
) -> plt.Figure:
```

**Дополнительно — zoom режим:**
- Аргумент `region=(start_bp, end_bp)` для отображения конкретного локуса
- Полезно для Главы 3 диплома (показать конкретный пример расхождения)

---

### 🟡 P3: Сравнительный dashboard

Один PNG с несколькими subplots для одной хромосомы:
```
┌──────────────────────────────┐
│   Hi-C heatmap (log scale)   │
├──────────────────────────────┤
│   Tracks: границы алгоритмов │
├──────────────────────────────┤
│   Insulation Score (TopDom)  │
├──────────────────────────────┤
│   CTCF peaks (BED)           │
└──────────────────────────────┘
```

**Функция:**
```python
def plot_genomic_dashboard(
    matrix, boundaries_by_algo, ctcf_bed,
    chrom, resolution, region=None, output_path=None
) -> plt.Figure:
```

---

### 🟡 P4: Jaccard heatmap улучшения

Текущий Jaccard heatmap — базовый. Добавить:
- Аннотации значений в ячейках
- Дендрограмму (иерархическая кластеризация алгоритмов)
- Цветовую шкалу с явными метками (0.0 = нет перекрытия, 1.0 = идентично)

---

## Технические правила для этого workstream

```python
# ✅ Обязательно
import matplotlib.pyplot as plt
import seaborn as sns
cfg = load_config("config/config.yaml")
dpi = cfg["visualization"]["dpi"]           # 300
fmt = cfg["visualization"]["format"]        # "png"
figsize = cfg["visualization"]["figsize"]   # [18, 14]

# ✅ Сохранение
fig.savefig(output_path, dpi=dpi, bbox_inches="tight", format=fmt)
plt.close(fig)  # ОБЯЗАТЕЛЬНО — иначе memory leak при 23 хромосомах

# ❌ Не использовать plt.show() в src/ — только в notebooks/
# ❌ Не хардкодить figsize, dpi, colors
# ❌ generate_html: false — не включать без согласования
```

---

## Папка с результатами

```
results/figures/
├── hic_tads_<chrom>_<res>bp.png           # основной heatmap
├── ctcf_profile_wide_<chrom>_<res>bp.png  # новый широкий профиль
├── ctcf_profile_<algo>_<chrom>_<res>bp.png
├── jaccard_<chrom>_<res>bp.png
└── dashboard_<chrom>_<res>bp.png          # новый dashboard
```

# Алгоритмы детекции TAD — подробная документация

> Связанные файлы: src/algorithms/, config/config.yaml §algorithms

---

## Единый интерфейс (ОБЯЗАТЕЛЬНО соблюдать)

```python
def run_<algorithm>(
    chrom: str,          # "chr17" — всегда с префиксом
    resolution: int,     # 25000, 50000, 100000 (в bp)
    data_path: str,      # путь к processed данным
    cfg: Optional[dict], # секция алгоритма из config.yaml
    **kwargs,
) -> pd.DataFrame:
    # Обязательные колонки: ["chrom", "start", "end"]
    # При любой ошибке: return pd.DataFrame(columns=["chrom","start","end"])
```

**Нарушение интерфейса ломает pipeline/run_pipeline.py и тесты.**

---

## 1. Armatus — DP с gamma-регуляризацией

**Статья:** Filippova et al. 2014, Genome Research
**Механика:** Минимизирует функцию `sum(density(TAD)) - gamma * n_TAD` через DP.
Параметр `gamma` штрафует за количество доменов → маленький gamma = крупные TAD.

```yaml
armatus:
  gamma_values: [0.1, 0.3, 0.5, 0.6, 0.65, 0.7, 0.75, 0.8, 0.9, 1.0, 2.0]
  tad_density_min: 0.8   # TAD / Mb — нижняя граница валидного gamma
  tad_density_max: 2.5   # TAD / Mb — верхняя граница
  stability_top_n: 3
```

**Subprocess:** `armatus -S -N -c <chrom> -g <gamma> -r <resolution> -i <matrix>`
**Timeout:** 600с

**⚠️ ОТКРЫТЫЙ ВОПРОС:**
Диапазон 0.8–2.5 TAD/Mb используется одинаково для всех хромосом и разрешений.
Это биологически не обосновано — chr1 и chr22 имеют разную архитектуру.
Задача: проверить, нужен ли отдельный диапазон для разных хромосом.

**Воспроизведение из статьи:** TODO — запустить на chr1 @ 40kb (как в Filippova 2014)

---

## 2. TopDom — Insulation Score

**Статья:** Shin et al. 2016, PLOS Computational Biology
**Механика:** Строит insulation score IS(i) = mean(contacts в окне WxW вокруг i).
Границы TAD = локальные минимумы IS. Выбирается window с max intra/inter ratio.

```yaml
topdom:
  window_sizes: [3, 5, 10]  # в бинах
  default_window: 5
```

**⚠️ Известная проблема:**
При window=3 на 25kb разрешении чувствительность = 75kb — грубее масштаба CTCF.
Max intra/inter ratio может давать false positives на RAW матрицах без нормализации.

**Воспроизведение из статьи:** TODO — запустить на IMR90 @ 40kb (как в Shin 2016)

---

## 3. scKTLD — спектральная кластеризация

**Статья:** Zheng et al. 2024
**Механика:** Hi-C матрица → kNN-граф → спектральная декомпозиция (eigenvectors) →
k-means кластеризация. TAD = кластеры геномных бинов.

```yaml
scktld:
  dimension: 32   # число eigenvectors
  knn_k: 20       # k для kNN-графа
  balance: false  # ОБЯЗАТЕЛЬНО false — RAW counts, не нормировать
  auto_penalty: true
```

**GPU ветка:** если n_bins < 5000 → `torch.linalg.eigh` на CUDA (10–20x быстрее)

**Ограничения памяти:**
| Разрешение | Хромосомы |
|-----------|-----------|
| 25kb | chr17–chr22 |
| 50kb | chr1–chr22 |
| 100kb | chr1–chr22, chrX |

**⚠️ ОТКРЫТЫЙ ВОПРОС:** Число кластеров (TAD) нужно задавать явно.
Как автоматически подобрать k для каждой хромосомы? Gap statistic? Eigengap?

---

## 4. coiTAD — OI matrix

**Механика:** Observed/Insulation matrix → пики = границы TAD.
Fallback на IS если OI не даёт результат.

```yaml
coitad:
  default_params:
    window_bins: null
    sigma: 1.5
    prominence_factor: 0.20
    min_tad_kb: 100
```

**⚠️ Проблема:** Fallback на IS делает результаты похожими на TopDom.
Нужно проверить: в каком % случаев срабатывает fallback?

**Воспроизведение из статьи:** TODO — найти оригинальную статью coiTAD

---

## 5. DI+HMM — Directionality Index + Hidden Markov Model

**Статья:** Dixon et al. 2012, Nature (оригинальная работа по TAD)
**Механика:**
1. DI(i) = (D - U) / (D + U) где D = downstream contacts, U = upstream contacts
2. GaussianHMM (3 состояния): upstream-biased / unbiased / downstream-biased
3. TAD = блок `upstream → unbiased → downstream`

```yaml
dihmm:
  window_bins: 10   # окно для DI (10 × 25kb = 250kb)
  n_states: 3
  n_iter: 100
  min_tad_kb: 100.0
```

**Трёхуровневый fallback:**
1. Dixon паттерн: `state=2 → state=0`
2. Boundary-mode: любые переходы состояний HMM
3. DI-minima: нулевые пересечения сглаженного DI (σ=2)

**Результаты (25kb):** chr22=54 TAD, chr21=68 TAD ✅

---

## 6. OnTAD — иерархическая детекция

**Статья:** An et al. 2019, Genome Biology
**Механика:** Скользящее среднее по диагонали матрицы → иерархия TAD (дерево).
Глубина 0 = весь геном, глубина 1 = крупные домены, глубина N = субTAD.

```yaml
ontad:
  penalty: null   # null = auto-sweep [1.0..0.01]
  minsz: 3        # минимальный TAD в бинах
  maxsz: 200      # максимальный TAD в бинах
  log2: true      # -log2(x+1) нормализация входа
```

**Формат .tad файла:**
```
start_bin  end_bin  depth  score1  score2
1          2050     0      ...     ...    ← root, пропускаем
685        747      1      ...     ...    ← top-level TAD ← берём
685        731      2      ...     ...    ← sub-TAD, опционально
```

**Стратегия выбора:** `_get_max_nonoverlapping(df_all, min_tad_bins)` —
жадный выбор листовых TAD, непересекающийся набор.

**Конвертация координат:**
`start_bp = (start_bin - 1) * resolution`
`end_bp = end_bin * resolution`

**Результаты (25kb):** chr22≈60 TAD, chr21≈55 TAD ✅

---

## 7. ModularityTAD — граф-модулярность

**Механика:** 1D-адаптация модулярности Ньюмана для Hi-C.
`score(i,j) = mean_intra(i,j) - mean_flanks(i,j)` через 2D prefix-суммы + DP.

```yaml
modularity_tad:
  max_dist_mb: 5.0
  penalty: null       # null = auto-sweep [2.0..0.001]
  min_tad_kb: 100.0
  max_tad_kb: 3000.0
```

**Auto-penalty:** sweep по `pen_norm × percentile_95(scores)`

**Результаты (25kb):** chr22≈70 TAD, chr21≈65 TAD ✅

---

## Сравнение алгоритмов (для диплома)

| Алгоритм | Что оптимизирует | Нормализация нужна? | Иерархия? | Чувствителен к... |
|----------|-----------------|--------------------|-----------|--------------------|
| Armatus | Плотность инт-й внутри TAD | Нет | Нет | gamma, размер хром |
| TopDom | Локальная изоляция | Желательна | Нет | window size |
| scKTLD | Глобальная структура матрицы | Нет (RAW) | Нет | k (число TAD) |
| coiTAD | Observed/Insulation ratio | Нет | Нет | prominence |
| DI+HMM | Направленность контактов | Нет | Нет | window_bins |
| OnTAD | Локальное скольз. среднее | log2 | ДА | penalty, уровень |
| ModularityTAD | intra/inter ratio | Нет | Нет | penalty |

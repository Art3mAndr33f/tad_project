# Данные — форматы, пути, скачивание

---

## Источники

| Датасет | GEO/URL | Файл |
|---------|---------|------|
| Hi-C GM12878 primary | GSE63525 | RAWobserved по хромосомам |
| Arrowhead boundaries | GSE63525 suppl | `GSE63525_GM12878_primary+replicate_Arrowhead_domainlist.txt` |
| CTCF ChIP-seq | ENCODE ENCFF796WRU | `GM12878_CTCF_peaks_hg19.bed` |

**Координатная система:** hg19 (GRCh37) — ВЕЗДЕ, без исключений.
**.hic файл НЕ используется** — только RAWobserved.

---

## Скачивание

```bash
bash download_data.sh
# Результат:
# data/raw/GM12878_primary/<res>kb_resolution_intrachromosomal/<chrom>/...
# data/reference/GSE63525_GM12878_primary+replicate_Arrowhead_domainlist.txt
# data/reference/GM12878_CTCF_peaks_hg19.bed
```

---

## Формат RAWobserved

3-колоночный TSV, только верхний треугольник матрицы:
```
bin_i <TAB> bin_j <TAB> count
0           25000       142
0           50000       87
25000       25000       891
```
Координаты в bp (не в номерах бинов).

---

## Разрешения и хромосомы

```
Разрешения: [25000, 50000, 100000]
Хромосомы:  chr1–chr22, chrX

scKTLD ограничения памяти:
  25kb  → chr17–chr22
  50kb  → chr1–chr22
  100kb → chr1–chr22, chrX
```

---

## Ключевые функции data_prep.py

```python
cfg = load_config("config/config.yaml")

# Загрузить матрицу (4-уровневый fallback: .npy → RAWobserved → rebin → error)
matrix = get_matrix(cfg, chrom="chr17", resolution=25000)
# → numpy array shape (N, N), dtype float64, симметричная

# Путь к RAWobserved
raw_path = get_rawobserved_path(cfg, chrom="chr17", resolution=25000)

# Подготовить все матрицы (кеш в .npy)
prepare_all_matrices(cfg, resolutions=[25000], chromosomes=["chr17"], force=False)
```

---

## Форматы результатов

| Тип | Путь | Формат |
|-----|------|--------|
| TAD-лист | `results/tads/<algo>_<chrom>_<res>bp.bed` | BED3: chrom start end |
| Консенсус | `results/consensus/consensus_<chrom>_<res>bp.bed` | BED3 + support колонка |
| Dense матрица | `data/processed/<chrom>_<res>bp.npy` | numpy float64, (N,N) |

---

## Размеры хромосом hg19

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

## Формат Arrowhead (эталон Rao 2014)

TSV с заголовком:
```
chr1  start1  end1  chr2  start2  end2  color  score  ...
```
Только хромосома 1 в обоих столбцах (intrachromosomal). Координаты в bp.

---

## Формат CTCF BED

Standard BED6+:
```
chr1  903000  904000  peak_1  100  .  ...
```
**Важно:** проверить chr-prefix перед использованием:
```bash
cut -f1 data/reference/GM12878_CTCF_peaks_hg19.bed | sort -u | head
# Ожидаем: chr1, chr2, ... (с chr)
# Если без chr → нужен fix в src/validation.py
```

#!/bin/bash
# download_data.sh — корректные имена файлов GSE63525

set -e  # выход при ошибке

mkdir -p data/raw data/reference

BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl"

echo "=== Скачиваем интрахромосомные матрицы (primary, 7.1 Gb) ==="
wget -c -P data/raw/ \
    "${BASE}/GSE63525_GM12878_primary_intrachromosomal_contact_matrices.tar.gz"

echo "=== Распаковываем ==="
tar -xzf data/raw/GSE63525_GM12878_primary_intrachromosomal_contact_matrices.tar.gz \
    -C data/raw/

echo "=== Скачиваем Arrowhead domain list (242 Kb) — эталон ==="
wget -c -P data/reference/ \
    "${BASE}/GSE63525_GM12878_primary+replicate_Arrowhead_domainlist.txt.gz"
gunzip data/reference/GSE63525_GM12878_primary+replicate_Arrowhead_domainlist.txt.gz

echo "=== Скачиваем CTCF ChIP-seq GM12878 hg19 (ENCODE) ==="
wget -c -O data/reference/GM12878_CTCF_peaks_hg19.bed.gz \
    "https://www.encodeproject.org/files/ENCFF796WRU/@@download/ENCFF796WRU.bed.gz"
gunzip data/reference/GM12878_CTCF_peaks_hg19.bed.gz

echo "=== Готово ==="

#!/bin/bash
# download_data.sh

set -e

mkdir -p data/raw data/reference

BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl"

echo "=== Скачиваем и распаковываем матрицы (только 25kb / 50kb / 100kb) ==="
echo "    Архив не сохраняется на диск — потоковая распаковка"
echo "    Ожидаемый объём после распаковки: ~11.5 GiB"

wget -q --show-progress -O - \
    "${BASE}/GSE63525_GM12878_primary_intrachromosomal_contact_matrices.tar.gz" \
| tar -xzf - \
    -C data/raw/ \
    --wildcards \
    "*/25kb_resolution_intrachromosomal/*" \
    "*/50kb_resolution_intrachromosomal/*" \
    "*/100kb_resolution_intrachromosomal/*"

echo "=== Скачиваем Arrowhead domain list (эталон) ==="
wget -c -P data/reference/ \
    "${BASE}/GSE63525_GM12878_primary+replicate_Arrowhead_domainlist.txt.gz"
gunzip -f data/reference/GSE63525_GM12878_primary+replicate_Arrowhead_domainlist.txt.gz

echo "=== Скачиваем CTCF ChIP-seq GM12878 hg19 (ENCODE) ==="
wget -c -O data/reference/GM12878_CTCF_peaks_hg19.bed.gz \
    "https://www.encodeproject.org/files/ENCFF796WRU/@@download/ENCFF796WRU.bed.gz"
gunzip -f data/reference/GM12878_CTCF_peaks_hg19.bed.gz

echo "=== Готово. Итоговый размер данных: ==="
du -sh data/raw/GM12878_primary/*/
du -sh data/reference/
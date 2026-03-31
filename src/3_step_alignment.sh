#!/usr/bin/env bash

# Configuración del match de patrones
set -euo pipefail
shopt -s nullglob

# Creación de carpetas
DATA="./data"
LOGS="./logs"
RESULTS="./results"
CLEAN="${DATA}/fastp"
INDEXES="${DATA}/indexes"
STAR_PE="${RESULTS}/star_pe"
STAR_SE="${RESULTS}/star_se"
HISAT2_PE="${RESULTS}/hisat2_pe"
HISAT2_SE="${RESULTS}/hisat2_se"

mkdir -p "$STAR_PE" "$HISAT2_PE" "$HISAT2_SE" "$LOGS" "$STAR_SE"

# Índices de alineamiento
STAR_IDX="${INDEXES}/mm39.gencode.M36.star"
HISAT2_IDX="${INDEXES}/mm39.gencode.M36.hisat/mm39.gencode.M36.hisat"

[[ ! -f "${HISAT2_IDX}.1.ht2" ]] && echo "[ERROR] Índice HISAT2 no existente en $HISAT2_IDX" && exit 1
[[ ! -d "$STAR_IDX" ]] && echo "[ERROR] Directorio de índice STAR_PE no existente en $STAR_IDX" && exit 1

# Configuración de hilos
threads_star=8
threads_hisat2=8

echo "========== HISAT2 SE =========="

files=("$CLEAN"/SRR*_clean_1.fastq)

for f in "${files[@]}"; do
    base=$(basename "$f" _clean_1.fastq)

    # Alineamiento SE con HISAT2
    hisat2 \
        -p "$threads_hisat2" \
        -x "$HISAT2_IDX" \
        -U "$f" \
        -S "${HISAT2_SE}/${base}_hisat2_se.sam" \
        -t \
        --no-unal \
        --summary-file "${HISAT2_SE}/${base}_hisat2_se_summary.txt" \
        2> ${HISAT2_SE}/${base}_hisat2_se_time.txt

done

echo "========== HISAT2 PE =========="

for f in "${files[@]}"; do
    base=$(basename "$f" _clean_1.fastq)

    # Alineamiento PE con HISAT2
    hisat2 \
        -p "$threads_hisat2" \
        -x "$HISAT2_IDX" \
        -1 "$f" \
        -2 "$CLEAN/${base}_clean_2.fastq" \
        -S "$HISAT2_PE/${base}_hisat2_pe.sam" \
        -t \
        --no-unal \
        --summary-file "$HISAT2_PE/${base}_hisat2_pe_summary.txt" \
        2> ${HISAT2_PE}/${base}_hisat2_pe_time.txt
        
done

echo "========== STAR PE =========="

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate star

for f in "${files[@]}"; do
    base=$(basename "$f" _clean_1.fastq)

    # Alineamiento PE con STAR
    /usr/bin/time -f "%e" STAR \
        --runThreadN "$threads_star" \
        --genomeDir "$STAR_IDX" \
        --readFilesIn "$f" "${CLEAN}/${base}_clean_2.fastq" \
        --outFileNamePrefix "${STAR_PE}/${base}_star_pe_" \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped None \
        2> ${STAR_PE}/${base}_star_pe_time.txt

done

echo "========== STAR SE =========="

for f in "${files[@]}"; do
    base=$(basename "$f" _clean_1.fastq)

    # Alineamiento PE con STAR
    /usr/bin/time -f "%e" STAR \
        --runThreadN "$threads_star" \
        --genomeDir "$STAR_IDX" \
        --readFilesIn "$f" \
        --outFileNamePrefix "${STAR_SE}/${base}_star_se_" \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped None \
        2> ${STAR_SE}/${base}_star_se_time.txt

done

echo "Proceso terminado."
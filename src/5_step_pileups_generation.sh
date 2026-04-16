#!/usr/bin/env bash

set -euo pipefail
shopt -s nullglob

# Activación del entorno deeptools
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate deeptools

# Directorios y archivos
RESULTS="./results"
BLACKLIST="${RESULTS}/reference/mm38-blacklist.bed"

declare -A rutas
rutas[HISAT2_PE]="${RESULTS}/hisat2_pe"
rutas[HISAT2_SE]="${RESULTS}/hisat2_se"
rutas[STAR_PE]="${RESULTS}/star_pe"
rutas[STAR_SE]="${RESULTS}/star_se"

# Variables útiles
wks=8
BIN_SIZE=20
BAM_FILTERED="BAM_FILTERED_SORTED"

echo "================================================================================================="
current_date_time="$(date "+%Y-%m-%d %H:%M:%S")"
echo "Iniciando la generación de pileups - $current_date_time"

for key in "${!rutas[@]}"; do
    ruta="${rutas[$key]}/$BAM_FILTERED"
    OUT_DIR="${rutas[$key]}/PILEUPS"

    mkdir -p "$OUT_DIR"

    for f in "$ruta"/*.bam; do
        base="$(basename "$f" .bam)"

        echo "-------------------------------------------------------------------------------------------------"

        if [[ "$key" == *_PE ]]; then
            FLAG_INCLUDE=(--samFlagInclude 64)
            EXTEND_OPTION=(--extendReads)
        else
            FLAG_INCLUDE=()
            EXTEND_OPTION=()
        fi

        if [[ -f "$BLACKLIST" ]]; then
            echo "Procesando $f con la blacklist $BLACKLIST - [$(date "+%Y-%m-%d %H:%M:%S")]"

            bamCoverage \
                -p "$wks" \
                -b "$f" \
                -o "${OUT_DIR}/${base}.bw" \
                --binSize "$BIN_SIZE" \
                --blackListFileName "$BLACKLIST" \
                --normalizeUsing BPM \
                --skipNAs \
                --ignoreDuplicates \
                "${FLAG_INCLUDE[@]}" \
                "${EXTEND_OPTION[@]}"
        else
            echo "Procesando $f sin blacklist - [$(date "+%Y-%m-%d %H:%M:%S")]"

            bamCoverage \
                -p "$wks" \
                -b "$f" \
                -o "${OUT_DIR}/${base}.bw" \
                --binSize "$BIN_SIZE" \
                --normalizeUsing BPM \
                --skipNAs \
                --ignoreDuplicates \
                "${FLAG_INCLUDE[@]}" \
                "${EXTEND_OPTION[@]}"
        fi

        echo "Procesamiento finalizado - [$(date "+%Y-%m-%d %H:%M:%S")]"
    done
done

echo "-------------------------------------------------------------------------------------------------"
echo "Eliminando archivos SAM de los alineadores HISAT2"
rm -fr "${rutas[HISAT2_PE]}/SAM/" "${rutas[HISAT2_SE]}/SAM/"

echo "-------------------------------------------------------------------------------------------------"
current_date_time="$(date "+%Y-%m-%d %H:%M:%S")"
echo "Script finalizado con éxito - $current_date_time"
echo "================================================================================================="

#!/usr/bin/env bash

set -euo pipefail
shopt -s nullglob

# Directorios de trabajo
DATA="./data"
RESULTS="./results"
FASTQ_CLEANED="${DATA}/fastp"
SALMON_DIR="${RESULTS}/SALMON"

mkdir -p "$SALMON_DIR"

# Variables útiles
wks=8
READ_TYPE="paired"
SALMON_IDX="${DATA}/indexes/mm39.gencode.M36.salmon"

# Entorno de salmon
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate salmon

echo "================================================================================================="
current_date_time="$(date "+%Y-%m-%d %H:%M:%S")"
echo "Iniciando pseudoalineamientos con SALMON - $current_date_time"
echo ""
echo "-------------------------------------------------------------------------------------------------"

# Obtención de nombres de muestra
fastq_array=("${FASTQ_CLEANED}"/*.fastq)
mapfile -t nombres < <(
    for f in "${fastq_array[@]}"; do
        base=$(basename "$f")
        echo "${base%%_clean*}"
    done | sort -u
)

for sample in "${nombres[@]}"; do
    echo "-------------------------------------------------------------------------------------------------"
    echo "Procesando muestra ${sample}... - [$(date "+%Y-%m-%d %H:%M:%S")]"
    echo "Tipo de lectura: $READ_TYPE"
    echo ""

    OUT="${SALMON_DIR}/$sample"

    if [[ "$READ_TYPE" == "paired" ]]; then

        R1="${FASTQ_CLEANED}/${sample}_clean_1.fastq"
        R2="${FASTQ_CLEANED}/${sample}_clean_2.fastq"

        if [[ ! -f "$R1" || ! -f "$R2" ]]; then
            echo "No está uno de los pares de la muestra $sample"
            continue
        fi

        salmon quant \
            -i "$SALMON_IDX" \
            -l A \
            -1 "$R1" \
            -2 "$R2" \
            -p "$wks" \
            --validateMappings \
            -o "$OUT"

    else

        R="${FASTQ_CLEANED}/${sample}_clean.fastq"

        if [[ ! -f "$R" ]]; then
            echo "Sin archivo para la muestra $sample"
            continue
        fi

        salmon quant \
            -i "$SALMON_IDX" \
            -l A \
            -r "$R" \
            -p "$wks" \
            --validateMappings \
            -o "$OUT"
    fi

    echo "Fin del procesamiento - [$(date "+%Y-%m-%d %H:%M:%S")]"
    echo "-------------------------------------------------------------------------------------------------"
done

current_date_time="$(date "+%Y-%m-%d %H:%M:%S")"
echo "Script finalizado con éxito - $current_date_time"
echo "================================================================================================="
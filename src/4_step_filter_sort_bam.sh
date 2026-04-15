#!/usr/bin/env bash

set -euo pipefail
shopt -s nullglob

# Directorios de trabajo
RESULTS="./results"

declare -A rutas
rutas[HISAT2_PE]="${RESULTS}/hisat2_pe"
rutas[HISAT2_SE]="${RESULTS}/hisat2_se"
rutas[STAR_PE]="${RESULTS}/star_pe"
rutas[STAR_SE]="${RESULTS}/star_se"

# Variables útiles
SAM="sam"
BAM="bam"
SORTED="bam_filtered_sorted"
wks=8

echo "================================================================================================="
current_date_time="$(date "+%Y-%m-%d %H:%M:%S")"
echo "Iniciando conversión de SAM a BAM, filtrado y ordenamiento - $current_date_time"

# Iterando sobre las rutas
for key in "${!rutas[@]}"; do
    echo "-------------------------------------------------------------------------------------------------"
    ruta="${rutas[$key]}"

    RAW_BAM_DIR="${ruta}/${BAM^^}"
    SORTED_DIR="${ruta}/${SORTED^^}"

    mkdir -p "$RAW_BAM_DIR" "$SORTED_DIR"

    # Archivos de entrada y salida dependiendo del alineador
    if [[ "$key" == HISAT2* ]]; then
        INPUT_DIR="${ruta}/${SAM^^}"
        ARCHIVOS=("$INPUT_DIR"/*.sam)

        echo "Procesando directorio ${INPUT_DIR}/"
        echo "Alineamientos de HISAT2"

    else
        INPUT_DIR="${ruta}/${BAM^^}"
        ARCHIVOS=("$INPUT_DIR"/*.bam)

        echo "Procesando directorio ${INPUT_DIR}/"
        echo "Alineamientos de STAR"
    fi

    # Iterando sobre el arreglo con los archivos
    for f in "${ARCHIVOS[@]}"; do
        if [[ "$key" == HISAT2* ]]; then
            base="$(basename "$f" .sam)"
            RAW_BAM="${RAW_BAM_DIR}/${base}.bam"
        else
            base="$(basename "$f" .bam)"
            RAW_BAM="$f"
        fi

        SORTED_BAM="${SORTED_DIR}/${base}_filtered_sorted.bam"

        echo "Procesando archivo: $f"

        tmp_dir=$(mktemp -d)

        # Sólo HISAT2_(SE|PE) no tiene BAM crudo
        if [[ "$key" == HISAT2* ]]; then
            echo "Convirtiendo SAM a BAM crudo: $RAW_BAM"
            if ! samtools view -bS "$f" > "$RAW_BAM"; then
                echo "[ERROR] Fallo al convertir SAM a BAM con samtools view en $f" >&2
                rm -rf "$tmp_dir"
                exit 1
            fi
        fi

        echo "-------------------------------------------------------------------------------------------------"

        echo "Filtrando y ordenando hacia: $SORTED_BAM"
        if samtools view -F 4 -q 10 -u "$RAW_BAM" | \
           samtools sort -@ $wks -T "${tmp_dir}/sort" -o "$SORTED_BAM"; then
            echo "Procesado con éxito: $f"
        else
            echo "[ERROR] Falló al filtrar/ordenar con samtools en $f" >&2
            rm -rf "$tmp_dir"
            exit 1
        fi

        echo "-------------------------------------------------------------------------------------------------"

        echo "Indexando BAM: $SORTED_BAM"
        if samtools index -@ $wks "$SORTED_BAM"; then
            echo "Índice generado: ${SORTED_BAM}.bai"
        else
            echo "[ERROR] Falló el indexado de $SORTED_BAM" >&2
            rm -rf "$tmp_dir"
            exit 1
        fi

        rm -rf "$tmp_dir"
    done
done

echo "-------------------------------------------------------------------------------------------------"
current_date_time="$(date "+%Y-%m-%d %H:%M:%S")"
echo "Script finalizado con éxito - $current_date_time"
echo "================================================================================================="

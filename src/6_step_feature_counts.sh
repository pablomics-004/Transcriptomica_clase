#!/usr/bin/env bash

set -euo pipefail
shopt -s nullglob

# Descarga de la anotación
URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M38/gencode.vM38.annotation.gtf.gz"
ANNOTATION_DIR="./data/annotation"

mkdir -p "$ANNOTATION_DIR"

# Ejecución de Feature Counts
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate subread

# Variables útiles
wks=8

echo "================================================================================================="
current_date_time="$(date "+%Y-%m-%d %H:%M:%S")"
echo "Iniciando descarga de GTF - $current_date_time"
echo ""

# Descarga de datos
ANNOTATION_GZ="${ANNOTATION_DIR}/gencode.vM38.annotation.gtf.gz"
ANNOTATION="${ANNOTATION_DIR}/gencode.vM38.annotation.gtf"

# Evita volver a descargar si ya estaba
if [[ ! -f "$ANNOTATION" ]]; then
    wget -P "$ANNOTATION_DIR" "$URL" && gunzip -f "$ANNOTATION_GZ"
fi

echo "-------------------------------------------------------------------------------------------------"
current_date_time="$(date "+%Y-%m-%d %H:%M:%S")"
echo "Descarga exitosa - $current_date_time"
echo "================================================================================================="

# Directorios de trabajo
RESULTS="./results"
BAM_FILTERED="BAM_FILTERED_SORTED"

declare -A rutas
rutas[HISAT2_PE]="${RESULTS}/hisat2_pe"
rutas[HISAT2_SE]="${RESULTS}/hisat2_se"
rutas[STAR_PE]="${RESULTS}/star_pe"
rutas[STAR_SE]="${RESULTS}/star_se"

# Renombrado de archivos
declare -A edades
edades[SRR9126754]="18"
edades[SRR9127244]="3"
edades[SRR9127382]="18"
edades[SRR9127568]="18"
edades[SRR9126950]="3"
edades[SRR9127139]="3"
edades[SRR9126924]="18"

echo "================================================================================================="
current_date_time="$(date "+%Y-%m-%d %H:%M:%S")"
echo "Iniciando renombrado de archivos BAM y BAI - $current_date_time"
echo ""

# Iterando sobre las rutas de cada alineador
for key in "${!rutas[@]}"; do
    dir="${rutas[$key]}/${BAM_FILTERED}"

    # Evalúa si se trata de un directorio o no
    [[ -d "$dir" ]] || continue

    archivos_bam=()
    archivos=("$dir"/*.bam "$dir"/*.bai)

    for file in "${archivos[@]}"; do
        base=$(basename "$file")
        ext="${base##*.}"
        name="${base%.*}"

        # Separa el nombre por "_"
        IFS='_' read -r srr align tipo _ <<< "$name"

        edad="${edades[$srr]:-}"

        if [[ -z "$edad" ]]; then
            echo "No encontré edad para $srr, se omite: $base"
            continue
        fi

        new_name="${srr}_age${edad}_${align}_${tipo}_filtered.${ext}"

        # Desplazando archivos BAI a su respectivo directorio
        if [[ "${ext}" == "bai" ]]; then
            new_dir="${rutas[$key]}/BAI"
            mkdir -p "$new_dir"
            new_path="${new_dir}/${new_name}"
        else
            new_path="${dir}/${new_name}"
        fi

        # Obtención de los nombres de los archivos BAM
        if [[ "${ext}" == "bam" ]]; then
            archivos_bam+=("$new_path")
        fi

        echo "$base -> $new_name"
        echo "Destino: $new_path"
        mv -n -- "$file" "$new_path" # No sobreescribe archivos existentes
        echo "-------------------------------------------------------------------------------------------------"
    done

    # Ordenando los nombres de los archivos para featureCounts dentro de un arreglo
    mapfile -t bams_ordenados < <(printf "%s\n" "${archivos_bam[@]}" | sort)

    # Generación de listas ordenadas
    printf "%s\n" "${bams_ordenados[@]}" > "${rutas[$key]}/bam_${key}_sorted_names.txt"

    echo "================================================================================================="
    current_date_time="$(date "+%Y-%m-%d %H:%M:%S")"
    echo "Iniciando FeatureCounts - $current_date_time"
    echo ""

    if [[ ${#bams_ordenados[@]} -eq 0 ]]; then
        echo "No hay BAMs para $key, se omite featureCounts"
        continue
    fi

    # Revisión del tipo de lectura
    if [[ "$key" == *_PE ]]; then
        FEATURECOUNTS_PE=(-p -B)
    else
        FEATURECOUNTS_PE=()
    fi

    # Directorio de salida
    OUT="${rutas[$key]}/FeatureCounts"
    mkdir -p "$OUT"

    featureCounts \
        -T $wks \
        -a "$ANNOTATION" \
        -o "$OUT/counts_${key}_table.txt" \
        --largestOverlap \
        -s 0 \
        "${FEATURECOUNTS_PE[@]}" \
        "${bams_ordenados[@]}"

    current_date_time="$(date "+%Y-%m-%d %H:%M:%S")"
    echo "FeatureCounts finalizado, tabla guardada en "$OUT" - $current_date_time"
    echo "================================================================================================="
done

current_date_time="$(date "+%Y-%m-%d %H:%M:%S")"
echo "Renombrado exitoso - $current_date_time"
echo "================================================================================================="


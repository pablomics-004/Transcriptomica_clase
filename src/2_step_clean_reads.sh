#!/usr/bin/env bash

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate fastp

set -euo pipefail # Si un comando falla, el script se detiene
shopt -s nullglob # Evita errores si no se encuentran archivos coincidentes

# Creación de carpetas
DATA="./data"
FASTQ="$DATA/fastq"
CLEAN="$DATA/fastp"
FASTQC_CLEAN_1="$DATA/fastqc_1_clean"
FASTQC_CLEAN_2="$DATA/fastqc_2_clean"
MULTIQC_DIR="$DATA/multiqc"

mkdir -p "$FASTQ" "$CLEAN" "$MULTIQC_DIR" "$FASTQC_CLEAN_1" "$FASTQC_CLEAN_2"

# Configuración de hilos y control de procesos
wks=3
threads_fastp=2
threads_fastqc=2

# Manejo de procesos en paralelo
pids=()
batch_processing() {
    if [[ ${#pids[@]} -eq $wks ]]; then
        wait "${pids[@]}"
        current_date_time="$(date "+%Y-%m-%d %H:%M:%S")"
        echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
        echo "$current_date_time"
        echo "Lote de $wks procesos completados. Continuando con el script..."
        pids=()
    fi
}

echo "========== FASTP =========="

files=("$FASTQ"/SRR*_1.fastq)
erase_bp=17

for f in "${files[@]}"; do
    base=$(basename "$f" _1.fastq)

    # Limpieza PE por plataforma Ilumina NovaSeq
    fastp \
        -i "$FASTQ/${base}_1.fastq" \
        -I "$FASTQ/${base}_2.fastq" \
        -o "$CLEAN/${base}_clean_1.fastq" \
        -O "$CLEAN/${base}_clean_2.fastq" \
        -w "$threads_fastp" \
        --trim_poly_g \
        --trim_front1 $erase_bp \
        --trim_front2 $erase_bp \
        --detect_adapter_for_pe &

    pids+=("$!")
    batch_processing
done

# Manejo de procesos restantes
if [[ ${#pids[@]} -gt 0 ]]; then
    wait "${pids[@]}"
    pids=()
fi

echo "========== FASTQC SOBRE READS LIMPIOS =========="

conda activate multiqc

clean_files_1=("$CLEAN"/*_clean_1.fastq)

# Stats de la calidad resultante en las reads
for f1 in "${clean_files_1[@]}"; do
    base=$(basename "$f1" _clean_1.fastq)
    f2="$CLEAN/${base}_clean_2.fastq"

    fastqc \
        -o "$FASTQC_CLEAN_1" \
        -t "$threads_fastqc" \
        "$f1" &
    pids+=("$!")
    batch_processing

    if [[ -f "$f2" ]]; then
        fastqc \
            -o "$FASTQC_CLEAN_2" \
            -t "$threads_fastqc" \
            "$f2" &
        pids+=("$!")
        batch_processing
    else
        echo "Falta archivo par: $f2"
    fi
done

echo "========== MULTIQC =========="

# Comparación global de las stats antes y después de la limpieza
multiqc "$DATA" -o "$MULTIQC_DIR"

echo "Proceso terminado."
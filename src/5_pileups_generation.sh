#!/usr/bin/env bash

set -euo pipefail
shopt -s nullglob

# Activación del entorno deeptools
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate deeptools

# Directorios y archivos
RESULTS="./results"
BLACKLIST="${RESULTS}/reference/mm39-blacklist.bed"

declare -A rutas
rutas[HISAT2_PE]="${RESULTS}/hisat2_pe"
rutas[HISAT2_SE]="${RESULTS}/hisat2_se"
rutas[STAR_PE]="${RESULTS}/star_pe"
rutas[STAR_SE]="${RESULTS}/star_se"

# Variables útiles
wks=8
BAM_FILTERED="bam_filtered_sorted"




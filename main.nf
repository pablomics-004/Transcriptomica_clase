#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Parámetros
params.fastq_dir = "./data/fastq"
params.output_dir = "./data/"
params.fastqc_tag = ""
params.trim_front = 15
params.is_pe = null

// ================= PROCESOS =================

// ----------------- FASTQC -----------------
process FASTQC_RAW {

  publishDir { "${params.output_dir}/${subdir}" }, mode: "copy"

  input:
  tuple path(read), val(subdir)

  output:
  path "*_fastqc.{html,zip}"

  script:
  """
  fastqc ${read}
  """
}

process FASTQC_CLEAN {

  publishDir { "${params.output_dir}/${subdir}" }, mode: "copy"

  input:
  tuple path(read), val(subdir)

  output:
  path "*_fastqc.{html,zip}"

  script:
  """
  fastqc ${read}
  """
}

// ----------------- FASTP -----------------
process FASTP_PE {

  publishDir "${params.output_dir}/fastp", mode: "copy"

  input:
  tuple val(sample), path(reads)

  output:
  tuple val(sample), path("${sample}_clean_1.fastq"), path("${sample}_clean_2.fastq")

  script:
  """
  fastp \
    -i ${reads[0]} \
    -I ${reads[1]} \
    -o ${sample}_clean_1.fastq \
    -O ${sample}_clean_2.fastq \
    -w ${task.cpus} \
    --trim_poly_g \
    --trim_front1 ${params.trim_front} \
    --trim_front2 ${params.trim_front} \
    --detect_adapter_for_pe
  """
}

process FASTP_SE {

  publishDir "${params.output_dir}/fastp", mode: "copy"

  input:
  tuple val(sample), path(read)

  output:
  tuple val(sample), path("${sample}_clean.fastq")

  script:
  """
  fastp \
    -i ${read} \
    -o ${sample}_clean.fastq \
    -w ${task.cpus} \
    --trim_poly_g \
    --trim_front1 ${params.trim_front}
  """
}

// ----------------- MULTIQC -----------------
process MULTIQC {

  publishDir "${params.output_dir}/multiqc", mode: "copy"

  input:
  path qc_files

  output:
  path "multiqc_report.html"
  path "multiqc_data"

  script:
  """
  multiqc . -o .
  """
}

// ================= MAIN =================
workflow {
  // Validando parámetros
  if (params.is_pe == null) {
    error "Debe indicarse si se trata de lecturas PE o SE con el parámetro --is_pe (true|false)"
  }
  if (!file(params.fastq_dir).exists()) {
    error "El directorio no existe: ${params.fastq_dir}"
  }

  is_pe = params.is_pe.toString().toBoolean()
  tag = params.fastqc_tag ?: ""

  // ----------------- PE -----------------
  if (is_pe) {

    reads = Channel.fromFilePairs(
      "${params.fastq_dir}/*_{1,2}.fastq*",
      checkIfExists: true
    )

    // Visualización inicial de la calidad
    reads_for_fastqc = reads.flatMap { sample, pair ->
      [
        tuple(pair[0], "fastqc_1${tag}"),
        tuple(pair[1], "fastqc_2${tag}")
      ]
    }

    FASTQC_RAW(reads_for_fastqc)

    // Limpieza de los archivos
    clean = FASTP_PE(reads)

    clean_for_fastqc = clean.flatMap { sample, r1, r2 ->
      [
        tuple(r1, "fastqc_1_clean"),
        tuple(r2, "fastqc_2_clean")
      ]
    }

  // ----------------- SE -----------------
  } else {

    reads = Channel.fromPath(
      "${params.fastq_dir}/*.fastq*",
      checkIfExists: true
    ).map { read ->
      def sample = read.name.replaceAll(/\.fastq(\.gz)?$/, "")
      tuple(sample, read)
    }

    // Visualización inicial de la calidad
    reads_for_fastqc = reads.map { sample, read ->
      tuple(read, "fastqc${tag}")
    }

    FASTQC_RAW(reads_for_fastqc)

    // Limpieza de los archivos
    clean = FASTP_SE(reads)

    clean_for_fastqc = clean.map { sample, read ->
      tuple(read, "fastqc_clean")
    }
  }

  // ----------------- PARA AMBOS -----------------

  // Visualización final
  FASTQC_CLEAN(clean_for_fastqc)

  // Visualización global
  multiqc_input = FASTQC_RAW.out.mix(FASTQC_CLEAN.out).collect()
  MULTIQC(multiqc_input)
}
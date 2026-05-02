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
  tuple val(sample), path(read), val(subdir), val(name_suffix), val(ext)

  output:
  path "${sample}${name_suffix}_fastqc.{html,zip}"

  script:
  """
  ln -s ${read} ${sample}${name_suffix}${ext}
  fastqc ${sample}${name_suffix}${ext}
  """
}

process FASTQC_CLEAN {

  publishDir { "${params.output_dir}/${subdir}" }, mode: "copy"

  input:
  tuple val(sample), path(read), val(subdir), val(name_suffix), val(ext)

  output:
  path "${sample}${name_suffix}_fastqc.{html,zip}"

  script:
  """
  ln -s ${read} ${sample}${name_suffix}${ext}
  fastqc ${sample}${name_suffix}${ext}
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

// ================= FUNCIONES =================
def get_ext(read) {
    read.name.endsWith(".fastq.gz") ? ".fastq.gz" : ".fastq"
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
    
    // Visualización inicial de la calidad
    reads = Channel.fromFilePairs(
      "${params.fastq_dir}/*_{1,2}.fastq*",
      checkIfExists: true
    )

    reads_for_fastqc = reads.flatMap { sample, pair ->
      [
        tuple(sample, pair[0], "fastqc_1${tag}", "${tag}_1", get_ext(pair[0])),
        tuple(sample, pair[1], "fastqc_2${tag}", "${tag}_2", get_ext(pair[1]))
      ]
    }

    FASTQC_RAW(reads_for_fastqc)

    // Limpieza de los archivos
    clean = FASTP_PE(reads)

    clean_for_fastqc = clean.flatMap { sample, r1, r2 ->
      [
        tuple(sample, r1, "fastqc_1_clean", "_clean_1", get_ext(r1)),
        tuple(sample, r2, "fastqc_2_clean", "_clean_2", get_ext(r2))
      ]
    }
  
  // ----------------- SE -----------------
  } else {
    
    // Visualización inicial de la calidad
    reads = Channel.fromPath(
      "${params.fastq_dir}/*.fastq*",
      checkIfExists: true
    ).map { read ->
      def sample = read.name.replaceAll(/\.fastq(\.gz)?$/, "")
      tuple(sample, read) 
    }

    reads_for_fastqc = reads.map { sample, read ->
      tuple(sample, read, "fastqc${tag}", tag, get_ext(read))
    }

    FASTQC_RAW(reads_for_fastqc)

    // Limpieza de los archivos
    clean = FASTP_SE(reads)

    clean_for_fastqc = clean.map { sample, read ->
      tuple(sample, read, "fastqc_clean", "_clean", get_ext(read))
    }

  // ----------------- PARA AMBOS -----------------
  
  // Visualización final
  FASTQC_CLEAN(clean_for_fastqc)

  // Visualización global 
  multiqc_input = FASTQC_RAW.out.mix(FASTQC_CLEAN.out).collect()
  MULTIQC(multiqc_input)
  }
}
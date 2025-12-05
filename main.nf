#!/usr/bin/env nextflow
nextflow.enable.dsl=2

def idr_output = params.idr_output ?: "idr_output"

process idr_call {
  tag "${pair}"
  stageInMode 'symlink'
  stageOutMode 'move'

  publishDir "${params.project_folder}/${idr_output}", mode: 'copy'

  input:
    tuple val(pair), path(peaks1), path(peaks2)
  
  output:
    path("${pair}_idr.narrowPeak")
    path("${pair}_idr.txt")
    path("${pair}_idr.log")
    path("${pair}_idr.txt.png")

  script:
  """
  set -eux

  mkdir -p ${params.project_folder}/${idr_output}

  idr \\
    --samples ${peaks1} ${peaks2} \\
    --input-file-type narrowPeak \\
    --rank signal.value \\
    --output-file ${pair}_idr.txt \\
    --log-output-file ${pair}_idr.log \\
    --plot \\
    --idr-threshold ${params.idr_threshold}

  awk 'BEGIN{OFS="\\t"} !/^#/ {print \$1,\$2,\$3,"'${pair}'_IDR_peak_"NR,1000,".",\$7,\$8,\$9,\$10}' \\
      ${pair}_idr.txt > ${pair}_idr.narrowPeak
  """
}

process make_pseudo_reps {
  tag "${rep_name}"
  stageInMode 'symlink'
  stageOutMode 'move'

  publishDir "${params.project_folder}/${idr_output}/pseudo_bams", mode: 'copy'

  input:
    tuple val(rep_name), path(bam)

  output:
    tuple val(rep_name), path("${rep_name}.pseudo1.bam"), path("${rep_name}.pseudo2.bam")

  script:
  """
  set -eux

  samtools view ${bam} | awk '{print \$1}' | sort -u | shuf > ${rep_name}.all_readnames.txt

  total=\$(wc -l < ${rep_name}.all_readnames.txt)
  half=\$(( (total + 1) / 2 ))

  head -n \$half ${rep_name}.all_readnames.txt > ${rep_name}.names1.txt
  tail -n \$((half + 1)) ${rep_name}.all_readnames.txt > ${rep_name}.names2.txt

  samtools view -b -N ${rep_name}.names1.txt ${bam} -o ${rep_name}.pseudo1.bam
  samtools view -b -N ${rep_name}.names2.txt ${bam} -o ${rep_name}.pseudo2.bam

  samtools index ${rep_name}.pseudo1.bam
  samtools index ${rep_name}.pseudo2.bam
  """
}

process pseudo_idr {
  tag "${rep_name}"
  stageInMode 'symlink'
  stageOutMode 'move'

  publishDir "${params.project_folder}/${idr_output}", mode: 'copy'

  input:
    tuple val(rep_name), path(pseudo1_bam), path(pseudo2_bam)

  output:
    path("${rep_name}_pseudo_idr.narrowPeak")
    path("${rep_name}_pseudo_idr.txt")
    path("${rep_name}_pseudo_idr.log")
    path("${rep_name}_pseudo_idr.txt.png")

  script:
  """
  set -eux

  macs2 callpeak \\
    -t ${pseudo1_bam} \\
    -n ${rep_name}_pseudo1 \\
    -f BAM \\
    -g ${params.macs3_genome} \\
    --outdir . \\
    --keep-dup all \\
    -q 0.01

  macs2 callpeak \\
    -t ${pseudo2_bam} \\
    -n ${rep_name}_pseudo2 \\
    -f BAM \\
    -g ${params.macs3_genome} \\
    --outdir . \\
    --keep-dup all \\
    -q 0.01

  idr \\
    --samples ${rep_name}_pseudo1_peaks.narrowPeak ${rep_name}_pseudo2_peaks.narrowPeak \\
    --input-file-type narrowPeak \\
    --rank signal.value \\
    --output-file ${rep_name}_pseudo_idr.txt \\
    --log-output-file ${rep_name}_pseudo_idr.log \\
    --plot \\
    --idr-threshold ${params.idr_threshold}

  awk 'BEGIN{OFS="\\t"} !/^#/ {print \$1,\$2,\$3,"'${rep_name}'_PSEUDO_IDR_peak_"NR,1000,".",\$7,\$8,\$9,\$10}' \\
      ${rep_name}_pseudo_idr.txt > ${rep_name}_pseudo_idr.narrowPeak
  """
}


workflow {

  def idr_csv = params.idr_pairs_csv ?: "idr_pairs.csv"

  rows = Channel
    .fromPath(idr_csv, checkIfExists: true)
    .splitCsv(header: true)

  rows = rows.filter { row ->
    ! file("${params.project_folder}/${idr_output}/${row.pair_name}_idr.narrowPeak").exists()
  }

  pairs_ch = rows.map { row ->
    tuple( row.pair_name, file(row.rep1_peaks), file(row.rep2_peaks) )
  }

  idr_call(pairs_ch)

  if( params.do_pseudo_idr ) {

    def pseudo_bam_csv = params.pseudo_idr_bam_csv ?: "pseudo_idr_bam.csv"

    pseudo_rows = Channel
      .fromPath(pseudo_bam_csv, checkIfExists: true)
      .splitCsv(header: true)

    pseudo_rows = pseudo_rows.filter { row ->
      ! file("${params.project_folder}/${idr_output}/${row.rep_name}_pseudo_idr.narrowPeak").exists()
    }

    pseudo_input_ch = pseudo_rows.map { row ->
      tuple( row.rep_name, file(row.bam) )
    }

    pseudo_input_ch | make_pseudo_reps | pseudo_idr
  }
}

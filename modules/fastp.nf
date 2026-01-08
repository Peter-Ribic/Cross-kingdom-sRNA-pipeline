#!/usr/bin/env nextflow

process FASTP_TRIM {

  tag "$sample_id"

  conda 'bioconda::fastp'
  publishDir "results/trimming/${sample_id}", mode: 'symlink'

  input:
  tuple val(sample_id), path(reads)

  output:
  tuple val(sample_id), path("*_trimmed.fq.gz"), emit: trimmed_reads
  path("*_fastp.html"), emit: qc_html
  path("*_fastp.json"), emit: qc_json

  script:
  """
  set -euo pipefail

  for r in ${reads.join(' ')}; do
      base=\$(basename "\$r")
      prefix="\${base%%.*}"

      fastp \\
        -i "\$r" \\
        -o "\${prefix}_trimmed.fq.gz" \\
        --disable_adapter_trimming \\
        --qualified_quality_phred 20 \\
        --length_required 18 \\
        --html "\${prefix}_fastp.html" \\
        --json "\${prefix}_fastp.json"
  done
  """
}

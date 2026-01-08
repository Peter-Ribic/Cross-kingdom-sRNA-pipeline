#!/usr/bin/env nextflow

process TRIM_GALORE {
    tag "$sample_id"
    container "community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18"
    publishDir "results/trimming/${sample_id}", mode: 'symlink'

    input:
    tuple val(sample_id), path(reads)
    path(adapter_table)

    output:
    tuple val(sample_id),
          path("*_trimmed.fq.gz"), emit: trimmed_reads
    path("*_trimming_report.txt"), emit: trimming_reports
    path("*_trimmed_fastqc.{zip,html}"), emit: fastqc_reports
    path "${task.process}_${sample_id}.tsv", emit: log_info

    script:
    """
    set -euo pipefail

    # Iterate each reads file and run trim_galore with its matching adapter
    for r in ${reads.join(' ')}; do
      reads_base=\$(basename "\$r")
      reads_key="\${reads_base%%.*}"   # prefix before first dot

      adapter=\$(
        awk -v key="\$reads_key" '
          BEGIN { FS="[\t, ]+" }
          \$1 == key { print \$2; found=1; exit }
          END { if (!found) exit 2 }
        ' "${adapter_table}"
      ) || {
        echo "ERROR: No adapter found in ${adapter_table} for reads key: \$reads_key (from file: \$reads_base)" >&2
        exit 1
      }

      echo "Running trim_galore for \$reads_base using adapter \$adapter" >&2

      trim_galore \\
        --fastqc --gzip --length 19 --stringency 10\\
        --adapter "\$adapter" \\
        "\$r"
    done

    num_outputed_reads=\$(zcat *_trimmed.fq.gz | awk 'END{print NR/4}')
    echo -e "${task.process}\\t${sample_id}\\t\$num_outputed_reads" >> ${task.process}_${sample_id}.tsv
    """
}

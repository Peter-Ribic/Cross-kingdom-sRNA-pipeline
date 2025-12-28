#!/usr/bin/env nextflow

process TRIM_GALORE {
    tag "$sample_id"
    container "community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18"
    publishDir "results/trimming/${sample_id}", mode: 'symlink'

    input:
    tuple val(sample_id), path(reads)

    // output:
    // path "${reads.simpleName}_trimmed.fq.gz", emit: trimmed_reads
    // path "${reads}_trimming_report.txt", emit: trimming_reports
    // path "${reads.simpleName}_trimmed_fastqc.{zip,html}", emit: fastqc_reports

    output:
        tuple val(sample_id), 
        path("*_trimmed.fq.gz"), emit: trimmed_reads
        path("*_trimming_report.txt"), emit: trimming_reports
        path("*_trimmed_fastqc.{zip,html}"), emit: fastqc_reports
    script:
    """
    trim_galore --fastqc --gzip --small_rna --length 19 ${reads.join(' ')}
    """
}

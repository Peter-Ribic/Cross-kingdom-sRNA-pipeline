#!/usr/bin/env nextflow

process MERGE_READS {
    input:
    tuple val(sample_id), path(reads)
    output:
    tuple val(sample_id), path("${sample_id}_merged.fq.gz"), emit: merged_reads

    script:
    """
    cat ${reads.join(' ')} > ${sample_id}_merged.fq.gz
    """
}

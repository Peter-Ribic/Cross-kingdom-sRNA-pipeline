process FILTER_PATHOGEN_READS {
    tag "$sample_id"
    container "quay.io/biocontainers/seqkit:2.10.1--he881be0_0"
    publishDir "results/${sample_id}/pathogen_specific_reads", mode: 'symlink'

    input:
    tuple val(sample_id), path(trimmed_reads), path (pathogen_unique_ids)

    output:
    tuple val(sample_id), path("${sample_id}_pathogen_specific.fq.gz")

    script:
    """
    zcat ${trimmed_reads.join(' ')} | \
        seqkit grep -f ${pathogen_unique_ids} -o ${sample_id}_pathogen_specific.fq
    gzip ${sample_id}_pathogen_specific.fq
    """
}

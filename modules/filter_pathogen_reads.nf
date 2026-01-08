process FILTER_PATHOGEN_READS {
    tag "$sample_id"
    container "quay.io/biocontainers/seqkit:2.10.1--he881be0_0"
    publishDir "results/pathogen_specific_reads/${sample_id}", mode: 'symlink'

    input:
    tuple val(sample_id), path(trimmed_reads), path (pathogen_unique_ids)

    output:
    tuple val(sample_id), path("${sample_id}_pathogen_specific.fq.gz"), emit: pathogen_specific_reads
    path "${task.process}_${sample_id}.tsv", emit: log_info

    script:
    """
    zcat ${trimmed_reads.join(' ')} | \
        seqkit grep -f ${pathogen_unique_ids} -o ${sample_id}_pathogen_specific.fq

    num_sequences=\$(awk 'END{print NR/4}' ${sample_id}_pathogen_specific.fq)
    echo -e "${task.process}\\t${sample_id}\\t\$num_sequences" > ${task.process}_${sample_id}.tsv
    
    gzip ${sample_id}_pathogen_specific.fq
    """
}

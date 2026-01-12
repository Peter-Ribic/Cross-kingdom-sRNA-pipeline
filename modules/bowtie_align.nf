process BOWTIE_ALIGN {
    tag "$sample_id"
    container "biocontainers/bowtie2:v2.4.1_cv1"

    input:
    tuple val(sample_id), path(reads)
    path index_files
    val(index_prefix)
    
    output:
    tuple val(sample_id), path(reads), path("${sample_id}_${index_prefix}.sam"), emit: list_input
    path("${sample_id}_${index_prefix}.log"), emit: log

    script:
    """
    bowtie2 \
        --local \
        -L 19 -N 0 \
        --no-unal \
        -p ${task.cpus} \
        -x ${index_prefix} \
        -U ${reads} \
        -S ${sample_id}_${index_prefix}.sam 2> ${sample_id}_${index_prefix}.log   
    """
}
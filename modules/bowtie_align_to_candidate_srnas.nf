process BOWTIE_ALIGN_TO_CANDIDATE_SRNAS {
    tag "$sample_id"
    container "biocontainers/bowtie2:v2.4.1_cv1"
    publishDir "results/candidate_srnas_alignments/${sample_id}", mode: 'symlink'
    input:
    tuple val(sample_id), path(reads)
    path index_files

    output:
    tuple val(sample_id), path(reads), path("${sample_id}_candidate_srnas_0mm.sam"), emit: list_input
    path("${sample_id}_candidate_srnas_0mm.log"), emit: log

    script:
    """
    bowtie2 \
        --end-to-end \
        --score-min L,0,-0.1 \
        --no-unal \
        -p ${task.cpus} \
        -x candidate_srnas \
        -U ${reads} \
        -S ${sample_id}_candidate_srnas_0mm.sam 2> ${sample_id}_candidate_srnas_0mm.log   
    """
}
process BOWTIE_ALIGN_TO_HOST {
    tag "$sample_id"
    memory '100 GB'
    container "biocontainers/bowtie2:v2.4.1_cv1"
    publishDir "results/${sample_id}/host_alignments", mode: 'symlink'

    input:
    tuple val(sample_id), path(reads), path(sam_file), path(pathogen_ids)
    path index_prefix
    
    output:
    tuple val(sample_id), path(reads), path("${sample_id}_host_0mm.sam"), path(pathogen_ids), emit: list_input

    script:
    """
    bowtie2 \
        --very-sensitive \
        --no-unal \
        -p ${task.cpus} \
        -x host_index \
        -U ${reads.join(',')} \
        -S ${sample_id}_host_0mm.sam
    """
}
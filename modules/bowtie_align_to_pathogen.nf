process BOWTIE_ALIGN_TO_PATHOGEN {
    tag "$sample_id"
    container "biocontainers/bowtie2:v2.4.1_cv1"
    publishDir "results/${sample_id}/pathogen_alignments", mode: 'symlink'
    input:
    tuple val(sample_id), path(reads)
    path index_files

    output:
    tuple val(sample_id), path(reads), path("${sample_id}_pathogen_0mm.sam"), path("${sample_id}_pathogen_ids.txt"), emit: results

    script:
    """
    bowtie2 \
        --end-to-end \
        --score-min L,0,-0.1 \
        --no-unal \
        -p ${task.cpus} \
        -x pathogen_index \
        -U ${reads} \
        -S ${sample_id}_pathogen_0mm.sam
    

    grep -v "^@" ${sample_id}_pathogen_0mm.sam | cut -f1 | sort | uniq > ${sample_id}_pathogen_ids.txt
    """
}
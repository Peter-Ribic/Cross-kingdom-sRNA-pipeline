process BOWTIE_ALIGN_TO_HOST {
    tag "$sample_id"
    memory '100 GB'
    container "biocontainers/bowtie2:v2.4.1_cv1"
    publishDir "results/host_alignments/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(reads), path(sam_file), path(pathogen_ids)
    path index_prefix
    
    output:
    tuple val(sample_id), path(reads), path("${sample_id}_host_0mm.sam"), path(pathogen_ids), emit: list_input
    path("${sample_id}_host_0mm.log"), emit: log
    path "${task.process}_${sample_id}.tsv", emit: log_info


    script:
    """
    bowtie2 \
        --local \
        -L 19 -N 0 \
        --no-unal \
        -p ${task.cpus} \
        -x host_index \
        -U ${reads} \
        -S ${sample_id}_host_0mm.sam 2> ${sample_id}_host_0mm.log   

    total_aligned_records=\$(grep -vc "^@" ${sample_id}_host_0mm.sam)
    echo -e "${task.process}\\t${sample_id}\\t\$total_aligned_records" > ${task.process}_${sample_id}.tsv
    """
}
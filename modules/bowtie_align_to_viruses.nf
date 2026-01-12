process BOWTIE_ALIGN_TO_VIRUSES {
    tag "$sample_id"
    container "biocontainers/bowtie2:v2.4.1_cv1"
    publishDir "results/viruses/viruses_alignments/${sample_id}", mode: 'symlink'

    input:
    tuple val(sample_id), path(reads)
    path index_files
    
    output:
    tuple val(sample_id), path(reads), path("${sample_id}_viruses.sam"), emit: results
    path("${sample_id}_viruses.log"), emit: log
    path "${task.process}_${sample_id}.tsv", emit: log_info

    script:
    """
    bowtie2 \
        --local \
        -L 19 -N 0 \
        --no-unal \
        -p ${task.cpus} \
        -x viruses_index \
        -U ${reads} \
        -S ${sample_id}_viruses.sam 2> ${sample_id}_viruses.log   

    total_aligned_records=\$(grep -vc "^@" ${sample_id}_viruses.sam)
    echo -e "${task.process}\\t${sample_id}\\t\$total_aligned_records" > ${task.process}_${sample_id}.tsv
 """
}
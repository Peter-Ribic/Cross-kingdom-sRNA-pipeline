process BOWTIE_ALIGN_TO_VIRUSES {
    tag "$sample_id"
    container "biocontainers/bowtie2:v2.4.1_cv1"
    publishDir "results/viruses_alignments/${sample_id}", mode: 'symlink'

    input:
    tuple val(sample_id), path(reads)
    path index_files
    
    output:
    tuple val(sample_id), path(reads), path("${sample_id}_viruses_0mm.sam"), emit: results
    path("${sample_id}_viruses_0mm.log"), emit: log
    path "${task.process}_${sample_id}.tsv", emit: log_info

    script:
    """
    bowtie2 \
        --end-to-end \
        --score-min L,0,0 \
        --no-unal \
        -p ${task.cpus} \
        -x viruses_index \
        -U ${reads} \
        -S ${sample_id}_viruses_0mm.sam 2> ${sample_id}_viruses_0mm.log   

    total_aligned_records=\$(grep -vc "^@" ${sample_id}_viruses_0mm.sam)
    echo -e "${task.process}\\t${sample_id}\\t\$total_aligned_records" > ${task.process}_${sample_id}.tsv
 """
}
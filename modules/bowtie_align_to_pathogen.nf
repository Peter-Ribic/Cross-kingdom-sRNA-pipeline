process BOWTIE_ALIGN_TO_PATHOGEN {
    tag "$sample_id"
    container "biocontainers/bowtie2:v2.4.1_cv1"
    publishDir "results/pathogen_alignments/${sample_id}", mode: 'copy'
    input:
    tuple val(sample_id), path(reads)
    path index_files

    output:
    tuple val(sample_id), path(reads), path("${sample_id}_pathogen_0mm.sam"), path("${sample_id}_pathogen_ids.txt"), emit: results
    path("${sample_id}_pathogen_0mm.log"), emit: log
    path "${task.process}_${sample_id}.tsv", emit: log_info
    
    script:
    """
    bowtie2 \
        --local \
        -L 19 -N 0 \
        --no-unal \
        -p ${task.cpus} \
        -x pathogen_index \
        -U ${reads} \
        -S ${sample_id}_pathogen_0mm.sam 2> ${sample_id}_pathogen_0mm.log   

    total_aligned_records=\$(grep -vc "^@" ${sample_id}_pathogen_0mm.sam)
    echo -e "${task.process}\\t${sample_id}\\t\$total_aligned_records" > ${task.process}_${sample_id}.tsv

    grep -v "^@" ${sample_id}_pathogen_0mm.sam | cut -f1 | sort | uniq > ${sample_id}_pathogen_ids.txt
    """
}
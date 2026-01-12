process LIST_PATHOGEN_READS {
    tag "$sample_id"
    memory '120 GB'
    cpus 40
    container "quay.io/biocontainers/samtools:1.22.1--h96c455f_0"
    publishDir "results/main_filtering/pathogen_reads_list/${sample_id}", mode: 'symlink'

    input:
    tuple val(sample_id), path(reads), path(sam_file), path(pathogen_ids)
    
    output:
    tuple val(sample_id), path(reads), path("${sample_id}_fungal_unique_ids.txt"), emit: filter_input
    path("${sample_id}_host_ids.txt"), emit: host_ids

    script:
    """
     samtools view -@ 8 "${sam_file}" \
    | awk '{print \$1}' \
    | awk '!seen[\$1]++' > "${sample_id}_host_ids.txt"

    sort -u "${sample_id}_host_ids.txt" -o "${sample_id}_host_ids.txt"

    comm -23 ${pathogen_ids} ${sample_id}_host_ids.txt > ${sample_id}_fungal_unique_ids.txt
    """
}
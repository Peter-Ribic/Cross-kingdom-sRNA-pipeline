process BOWTIE_ALIGN_TO_HOST {
    tag "$sample_id"
    container "quay.io/biocontainers/bowtie:1.3.1--py39h9046dc2_10"
    publishDir "results/${sample_id}/host_alignments", mode: 'symlink'

    input:
    tuple val(sample_id), path(reads), path(sam_file), path(pathogen_ids)
    path index_prefix
    
    output:
    tuple val(sample_id), path("${sample_id}_host_0mm.sam"), emit: sam
    tuple val(sample_id), path("${sample_id}_host_ids.txt"), emit: host_ids
    tuple val(sample_id), path(reads), path("${sample_id}_fungal_unique_ids.txt"), emit: filter_input

    script:
    """
    bowtie -v 0 -a --best --strata -p ${task.cpus} -x host_index \
        -q ${reads.join(' ')} \
        -S ${sample_id}_host_0mm.sam

    grep -v "^@" ${sample_id}_host_0mm.sam | cut -f1 | sort | uniq > ${sample_id}_host_ids.txt
    comm -23 ${pathogen_ids} ${sample_id}_host_ids.txt > ${sample_id}_fungal_unique_ids.txt
    """
}
process BOWTIE_ALIGN_TO_PATHOGEN {
    tag "$sample_id"
    container "quay.io/biocontainers/bowtie:1.3.1--py39h9046dc2_10"
    publishDir "results/${sample_id}/pathogen_alignments", mode: 'symlink'
    input:
    tuple val(sample_id), path(reads)
    path index_files

    output:
    tuple val(sample_id), path(reads), path("${sample_id}_pathogen_0mm.sam"), path("${sample_id}_pathogen_ids.txt"), emit: results

    script:
    """
    bowtie -v 0 -a --best --strata -p ${task.cpus} -x pathogen_index \
        -q ${reads.join(' ')} \
        -S ${sample_id}_pathogen_0mm.sam

    grep -v "^@" ${sample_id}_pathogen_0mm.sam | cut -f1 | sort | uniq > ${sample_id}_pathogen_ids.txt
    """
}
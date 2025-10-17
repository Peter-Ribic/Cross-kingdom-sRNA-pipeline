process BOWTIE_BUILD_PATHOGEN {
    container "quay.io/biocontainers/bowtie:1.3.1--py39h9046dc2_10"
    publishDir "results/pathogen_index", mode: 'copy'

    input:
    path fasta

    output:
    path "pathogen_index*", emit: index_files

    script:
    """
    bowtie-build $fasta pathogen_index
    """
}


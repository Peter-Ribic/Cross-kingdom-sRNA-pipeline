process BOWTIE_BUILD_PATHOGEN {
    container "biocontainers/bowtie2:v2.4.1_cv1"
    publishDir "results/pathogen_index", mode: 'copy'

    input:
    path fasta

    output:
    path "pathogen_index*", emit: index_files

    script:
    """
    bowtie2-build $fasta pathogen_index
    """
}


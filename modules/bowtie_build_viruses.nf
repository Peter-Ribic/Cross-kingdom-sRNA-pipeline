process BOWTIE_BUILD_VIRUSES {
    container "biocontainers/bowtie2:v2.4.1_cv1"

    input:
    path fasta

    output:
    path "viruses_index*", emit: index_files

    script:
    """
    bowtie2-build $fasta viruses_index
    """
}


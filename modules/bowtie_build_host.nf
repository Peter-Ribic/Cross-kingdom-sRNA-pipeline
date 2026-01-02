process BOWTIE_BUILD_HOST {
    memory '20 GB'
    cpus 10
    container "biocontainers/bowtie2:v2.4.1_cv1"

    input:
    path fasta

    output:
    path "host_index*", emit: index_prefix

    script:
    """
    bowtie2-build $fasta host_index
    """
}

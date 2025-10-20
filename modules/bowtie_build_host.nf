process BOWTIE_BUILD_HOST {
    container "biocontainers/bowtie2:v2.4.1_cv1"
    publishDir "results/host_index", mode: 'copy'

    input:
    path fasta

    output:
    path "host_index*", emit: index_prefix

    script:
    """
    bowtie2-build $fasta host_index
    """
}

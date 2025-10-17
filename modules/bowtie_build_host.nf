process BOWTIE_BUILD_HOST {
    container "quay.io/biocontainers/bowtie:1.3.1--py39h9046dc2_10"
    publishDir "results/host_index", mode: 'copy'

    input:
    path fasta

    output:
    path "host_index*", emit: index_prefix

    script:
    """
    bowtie-build $fasta host_index
    """
}

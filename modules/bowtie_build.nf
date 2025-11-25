process BOWTIE_BUILD {
    container "biocontainers/bowtie2:v2.4.1_cv1"
    publishDir "results/viruses_index", mode: 'copy'

    input:
    path fasta
    val(index_prefix)

    output:
    path "${index_prefix}*", emit: index_files

    script:
    """
    bowtie2-build $fasta ${index_prefix}
    """
}


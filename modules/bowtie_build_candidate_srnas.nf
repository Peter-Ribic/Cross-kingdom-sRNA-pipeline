process BOWTIE_BUILD_CANDIDATE_SRNAS {
    container "biocontainers/bowtie2:v2.4.1_cv1"
    publishDir "results/candidate_srnas", mode: 'copy'

    input:
    path fasta

    output:
    path "candidate_srnas*", emit: index_files

    script:
    """
    bowtie2-build $fasta candidate_srnas
    """
}


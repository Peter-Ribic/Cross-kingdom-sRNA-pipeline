process HISAT2_BUILD {
    container "community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e"
    publishDir "results/index", mode: 'copy'

    input:
    path fasta

    output:
    path "genome_index.*.ht2", emit: index

    script:
    """
    hisat2-build -p ${task.cpus} $fasta genome_index
    """
}

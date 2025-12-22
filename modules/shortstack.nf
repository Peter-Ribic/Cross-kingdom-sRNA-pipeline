process SHORTSTACK {

    tag "$sample_id"

    container "quay.io/biocontainers/shortstack:4.1.2--hdfd78af_0"

    publishDir "results/shortstack/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(trimmed_reads)
    path genome_fasta

    output:
    tuple val(sample_id), path("ShortStack_out"), emit: shortstack_out
    tuple val(sample_id), path("ShortStack_out/*condensed.bam"), emit: bam
    tuple val(sample_id), path("ShortStack_out/*condensed.fa"), emit: fasta
    tuple val(sample_id), path("ShortStack_out/Results.txt"), emit: results
    tuple val(sample_id), path("${sample_id}_shortstack_majorRNAs.fasta"), emit: major_rnas



    script:
    """
    ShortStack \
        --readfile ${trimmed_reads} \
        --genomefile ${genome_fasta} \
        --outdir ShortStack_out \
        --threads 10 \
        --mincov 100 \
        --dicermin 18

    # Extract major RNAs from ShortStack Results.txt and create FASTA file
    awk -F'\\t' 'NR>1 && \$12 > 0 {print ">"\$2"\\n"\$11}' ShortStack_out/Results.txt > ${sample_id}_shortstack_majorRNAs.fasta
    """
}

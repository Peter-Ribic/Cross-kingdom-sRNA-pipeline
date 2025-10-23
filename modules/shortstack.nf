process SHORTSTACK {

    tag "$sample_id"

    container "quay.io/biocontainers/shortstack:4.1.2--hdfd78af_0"

    publishDir "results/${sample_id}/shortstack", mode: 'symlink'

    input:
    tuple val(sample_id), path(trimmed_reads)
    path genome_fasta

    output:
    tuple val(sample_id), path("ShortStack_out"), emit: shortstack_out

    script:
    """
    ShortStack --readfile ${trimmed_reads.join(' ')} \
               --genomefile ${genome_fasta} \
               --outdir ShortStack_out \
               --threads 10 \
               --mmap u \
               --mincov 1
    """
}

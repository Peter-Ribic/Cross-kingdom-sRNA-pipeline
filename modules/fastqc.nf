process FASTQC {
    tag "$sample_id"
    container "community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18"
    publishDir "results/fastqc/${sample_id}", mode: 'symlink'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*.zip", emit: zip
    path "*.html", emit: html

    script:
    """
    fastqc ${reads.join(' ')} --outdir .
    """
}

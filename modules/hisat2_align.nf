process HISAT2_ALIGN {
   // container "community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e"
   // container "quay.io/biocontainers/hisat2:2.2.1--h503566f_8"
    tag "$sample_id"
    publishDir "results/${sample_id}/hisat2", mode: 'symlink'

    input:
    tuple val(sample_id), path(reads)
    path index_files

    output:
    tuple val(sample_id), path("${sample_id}.bam"), emit: bam
    path("${sample_id}.hisat2.log"), emit: log

    script:
    """
    hisat2 -x genome_index \
        -U $reads \
        --new-summary --summary-file ${sample_id}.hisat2.log \
        -S ${sample_id}.sam

    samtools sort -o ${sample_id}.bam ${sample_id}.sam
    samtools index ${sample_id}.bam
    """
}

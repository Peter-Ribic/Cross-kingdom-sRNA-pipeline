process FASTQ_TO_FASTA {
    tag "$sample_id"
    container "quay.io/biocontainers/mulled-v2-0fd299cadb7a80e2cc704b5d903ccc54893c512d:377fcafe3b6e7ef703094f8f47a64d081622ee09-2"
    publishDir "results/fasta_conversion/${sample_id}", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("${sample_id}.fasta"), emit: fasta
    
    script:
    """
    seqtk seq -A ${reads} > ${sample_id}.fasta
    """
}
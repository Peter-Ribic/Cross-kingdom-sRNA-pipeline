process KEEP_ONLY_PATHOGEN_READS {
    tag "$sample_id"
    memory '200 GB'
    cpus 20
    container "quay.io/biocontainers/seqkit:2.10.1--he881be0_0"
    publishDir "results/${sample_id}/keep_pathogen_reads_only", mode: 'copy'

    input:
    tuple val(sample_id), path(reads), path(pathogen_reads)

    output:
    tuple val(sample_id), path("${sample_id}_filtered.fq.gz"), emit: filtered_reads

    script:
    """
    # Extract unique pathogen sequences (just the sequence strings)
    seqkit seq -s $pathogen_reads | sort | uniq > pathogen_unique.txt

    # Keep only reads from hop sample that exactly match pathogen sequences
    seqkit grep -s -f pathogen_unique.txt --threads 20 -I $reads | gzip > ${sample_id}_filtered.fq.gz
    """
}
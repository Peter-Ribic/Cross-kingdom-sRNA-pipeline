process FILTER_SRNA_LENGTH {

    tag "$sample_id"
    memory '120 GB'
    cpus 40
    container "quay.io/biocontainers/seqkit:2.10.1--he881be0_0"

    publishDir "results/filtered_srna/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_20_22nt.fq.gz"), emit: filtered_reads
    path("${sample_id}_length_distribution.tsv"), emit: length_distribution

    script:
    """
    echo "Processing sample: ${sample_id}"
    echo "Input file: ${reads}"

    # Generate read length distribution table
    echo "Calculating read length distribution..."
    seqkit fx2tab -l -i -n -H ${reads} | cut -f2 | sort | uniq -c | awk '{print \$2"\\t"\$1}' > ${sample_id}_length_distribution.tsv

    # Filter reads between 20 and 22 nt
    echo "Filtering reads between 20 and 22 nt..."
    seqkit seq -m 20 -M 24 ${reads} | gzip > ${sample_id}_20_22nt.fq.gz

    echo "Done. Outputs:"
    echo "- ${sample_id}_20_22nt.fq.gz"
    echo "- ${sample_id}_length_distribution.tsv"
    """
}

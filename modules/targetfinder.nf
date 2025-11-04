process TARGETFINDER {

    tag "$sample_id"
    container "quay.io/biocontainers/targetfinder:1.7--0"

    publishDir "results/${sample_id}/target_prediction", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    path host_transcriptome_fasta

    output:
    path "${sample_id}_TargetFinder_results.txt", emit: targetfinder_results

    script:
    """
     # Convert FASTQ to FASTA
    echo "Converting FASTQ to FASTA..."
    zcat ${reads} | awk 'NR%4==1 {gsub(/^@/,">"); sub(/ .*/,""); print} NR%4==2 {print}' > ${sample_id}_predicted_sRNAs.fa

    echo "Running TargetFinder directly against the transcriptome genome..."
    targetfinder.pl \
        -s ${sample_id}_predicted_sRNAs.fa \
        -d ${host_transcriptome_fasta} \
        -c 1000 \
        > ${sample_id}_TargetFinder_results.txt

    echo "Done. Results saved to ${sample_id}_TargetFinder_results.txt"
    """
}
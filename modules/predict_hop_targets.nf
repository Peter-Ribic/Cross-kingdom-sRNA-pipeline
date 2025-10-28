process PREDICT_HOP_TARGETS {

    tag "$sample_id"
    container "quay.io/biocontainers/targetfinder:1.7--0"

    publishDir "results/${sample_id}/target_prediction", mode: 'copy'

    input:
    tuple val(sample_id), path(shortstack_out)
    path host_genome_fasta

    output:
    path "${sample_id}_TargetFinder_results.txt", emit: targetfinder_results

    script:
    """
    echo "Extracting small RNA sequences from ShortStack results..."
    awk 'NR>1 {print ">"\$1"\\n"\$11}' ${shortstack_out}/Results.txt > ${sample_id}_predicted_sRNAs.fa
    sed 's/T/U/g' ${host_genome_fasta} > hop_genome_rna.fa

    echo "Running TargetFinder directly against the host genome..."
    targetfinder.pl \
        -s ${sample_id}_predicted_sRNAs.fa \
        -d hop_genome_rna.fa \
        -q 7 \
        -r \
        > ${sample_id}_TargetFinder_results.txt

    echo "Done. Results saved to ${sample_id}_TargetFinder_results.txt"
    """
}

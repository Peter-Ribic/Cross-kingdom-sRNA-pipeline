process PREDICT_HOP_TARGETS {

    tag "$sample_id"
    container "./mirnatarget.sif"

    publishDir "results/${sample_id}/target_prediction", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    path host_transcriptome_fasta

    output:
    path "${sample_id}_MiRNATarget_results.txt", emit: mirnatarget_results

    script:
    """
    echo "Working directory: \$PWD"

    # Convert FASTQ to FASTA
    echo "Converting FASTQ to FASTA..."
    zcat ${reads} | awk 'NR%4==1 {gsub(/^@/,\">\" ); print} NR%4==2 {print}' > ${sample_id}_predicted_sRNAs.fa

    # Copy host genome to working directory
    cp ${host_transcriptome_fasta} ./ref.fasta

    # Run ssearch36 alignment to a file
    echo "Running ssearch36 alignment..."
    ssearch36 -f -8 -g -3 -E 10000 -T 8 -b 200 -r +4/-3 -n -U -W 10 -N 20000 \
        ${sample_id}_predicted_sRNAs.fa ./ref.fasta > ${sample_id}_ssearch.out

    # Parse ssearch output
    echo "Parsing ssearch..."
    /MiRNATarget/parse_ssearch.py -i ${sample_id}_ssearch.out > ${sample_id}_parsed_ssearch.tsv

    # Predict miRNA targets
    echo "Parsing miRNA targets..."
    /MiRNATarget/parse_mirna_targets.py -i ${sample_id}_parsed_ssearch.tsv -o ${sample_id}_MiRNATarget_results.txt

    # If output does not exist, create empty placeholder
    if [ ! -s ${sample_id}_MiRNATarget_results.txt ]; then
        echo "# No targets found" > ${sample_id}_MiRNATarget_results.txt
    fi

    echo "MiRNATarget finished. Results saved to ${sample_id}_MiRNATarget_results.txt"
    """
}

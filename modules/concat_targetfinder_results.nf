process CONCAT_TARGETFINDER_RESULTS {
    tag "$sample_id"
    container "quay.io/biocontainers/mulled-v2-0fd299cadb7a80e2cc704b5d903ccc54893c512d:377fcafe3b6e7ef703094f8f47a64d081622ee09-2"
    publishDir "results/concat_targetfinder_results", mode: 'copy'
    
    input:
    tuple val(sample_id), path(log_files)
    
    output:
    path("${sample_id}_combined_targetfinder_results.txt"), emit: combined_results
    path("${sample_id}_targetfinder_summary.txt"), emit: summary
    
    script:
    """
    # Concatenate all TargetFinder log files
    cat $log_files > ${sample_id}_combined_targetfinder_results.txt
    
    # Create a summary file with basic statistics
    echo "TargetFinder Results Summary for $sample_id" > ${sample_id}_targetfinder_summary.txt
    echo "==============================================" >> ${sample_id}_targetfinder_summary.txt
    echo "Total input files processed: \$(echo $log_files | wc -w)" >> ${sample_id}_targetfinder_summary.txt
    echo "Total lines in combined results: \$(wc -l < ${sample_id}_combined_targetfinder_results.txt)" >> ${sample_id}_targetfinder_summary.txt
    echo "Combined results file: ${sample_id}_combined_targetfinder_results.txt" >> ${sample_id}_targetfinder_summary.txt
    """
}
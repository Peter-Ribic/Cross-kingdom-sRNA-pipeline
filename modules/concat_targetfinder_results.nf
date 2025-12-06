process CONCAT_TARGETFINDER_RESULTS {
    tag "$sample_id"
    container "quay.io/biocontainers/mulled-v2-0fd299cadb7a80e2cc704b5d903ccc54893c512d:377fcafe3b6e7ef703094f8f47a64d081622ee09-2"
    publishDir "results/concat_targetfinder_results", mode: 'copy'
    
    input:
    tuple val(sample_id), path(log_files)
    
    output:
    tuple val(sample_id), path("${sample_id}_combined_targetfinder_results.txt"), emit: combined_results
    path("${sample_id}_targetfinder_summary.txt"), emit: summary
    
    script:
    """
    # Initialize combined results file
    touch ${sample_id}_combined_targetfinder_results.txt
    
    # Counter for files with actual results
    files_with_results=0
    total_files=0
    
    # Process each file - only concatenate if it doesn't start with "No results for"
    for file in $log_files; do
        total_files=\$((total_files + 1))
        # Check if file exists and has content, and first line doesn't start with "No results for"
        if [[ -s "\$file" ]] && ! head -n 1 "\$file" 2>/dev/null | grep -q "^No results for"; then
            cat "\$file" >> ${sample_id}_combined_targetfinder_results.txt
            files_with_results=\$((files_with_results + 1))
        fi
    done
    
    # Create a summary file with basic statistics
    echo "TargetFinder Results Summary for $sample_id" > ${sample_id}_targetfinder_summary.txt
    echo "==============================================" >> ${sample_id}_targetfinder_summary.txt
    echo "Total input files processed: \$total_files" >> ${sample_id}_targetfinder_summary.txt
    echo "Files with actual results: \$files_with_results" >> ${sample_id}_targetfinder_summary.txt
    echo "Files with no results: \$((total_files - files_with_results))" >> ${sample_id}_targetfinder_summary.txt
    
    # Check if combined file has content
    if [[ -s ${sample_id}_combined_targetfinder_results.txt ]]; then
        echo "Total lines in combined results: \$(wc -l < ${sample_id}_combined_targetfinder_results.txt)" >> ${sample_id}_targetfinder_summary.txt
    else
        echo "Total lines in combined results: 0" >> ${sample_id}_targetfinder_summary.txt
        echo "No target predictions found for any sRNA sequences." > ${sample_id}_combined_targetfinder_results.txt
    fi
    
    echo "Combined results file: ${sample_id}_combined_targetfinder_results.txt" >> ${sample_id}_targetfinder_summary.txt
    """
}
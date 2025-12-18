process MULTIQC {
    tag "$output_name"
    container "quay.io/biocontainers/multiqc:1.32--pyhdfd78af_0"
    publishDir "results/multiqc/end", mode: 'symlink'

    input:
    val output_name
    path qc_files

    output:
    path "${output_name}.html", emit: report
    path "${output_name}_data", emit: data
    path "${output_name}_merged_pathogen_filter.tsv", emit: merged_pathogen
    path "${output_name}_merged_treated_minus_control.tsv", emit: merged_treated_minus_control

    script:
    """
    # Run MultiQC on QC files
    multiqc ${qc_files.join(' ')} -n ${output_name}.html

    # Merge all *_pathogen_filter.tsv files
    pathogen_files=( *_pathogen_filter.tsv )

    merged_file="${output_name}_merged_pathogen_filter.tsv"

    if [ \${#pathogen_files[@]} -gt 0 ]; then
        # Take header from first non-comment line
        head -n 20 "\${pathogen_files[0]}" | grep -v '^#' | head -n 1 > "\$merged_file"

        # Append data rows from all files (skip comments and header)
        for f in "\${pathogen_files[@]}"; do
            grep -v '^#' "\$f" | tail -n +2 >> "\$merged_file"
        done
    else
        # No pathogen files found → create empty table
        echo -e "sample\ttotal_reads\tpathogen_unique\tmatched_reads\tpercent_matched" > "\$merged_file"
    fi
   
    # Merge all *_treated_minus_control.tsv files
    treated_minus_control_files=( *_treated_minus_control.tsv )

    merged_file_treated_minus_control="${output_name}_merged_treated_minus_control.tsv"

    if [ \${#treated_minus_control_files[@]} -gt 0 ]; then
        # Take header from first non-comment line
        head -n 20 "\${treated_minus_control_files[0]}" | grep -v '^#' | head -n 1 > "\$merged_file_treated_minus_control"

        # Append data rows from all files (skip comments and header)
        for f in "\${treated_minus_control_files[@]}"; do
            grep -v '^#' "\$f" | tail -n +2 >> "\$merged_file_treated_minus_control"
        done
    else
        # No pathogen files found → create empty table
        echo -e "sample\ttotal_reads\tpathogen_unique\tmatched_reads\tpercent_matched" > "\$merged_file_treated_minus_control"
    fi
    """
}

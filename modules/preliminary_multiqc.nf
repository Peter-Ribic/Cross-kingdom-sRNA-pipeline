process PRELIMINARY_MULTIQC {
    tag "$output_name"
    container "quay.io/biocontainers/multiqc:1.32--pyhdfd78af_0"
    publishDir "results/multiqc/preliminary", mode: 'symlink'

    input:
    val output_name
    path qc_files

    output:
    path "${output_name}.html", emit: report
    path "${output_name}_data", emit: data

    script:
    """
    multiqc ${qc_files.join(' ')} -n ${output_name}.html
    """
}

#!/usr/bin/env nextflow
process MULTIQC {
    tag "$output_name"
    container "community.wave.seqera.io/library/pip_multiqc:ad8f247edb55897c"
    publishDir "results/multiqc", mode: 'symlink'

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
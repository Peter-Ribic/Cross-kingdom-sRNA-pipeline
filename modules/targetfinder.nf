process TARGETFINDER {

    tag "$sample_id"
    container "quay.io/biocontainers/targetfinder:1.7--hdfd78af_4"

    publishDir "results/targetfinder/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    path host_transcriptome_fasta

    output:
    path("${sample_id}.log"), emit: log

    script:
    """
    targetfinder.pl \
        -s TTCTTGATTAGGTCTTGGAATA \
        -t ${task.cpus} \
        -d ${host_transcriptome_fasta} \
        -c 4 > ${sample_id}.log
    """
}
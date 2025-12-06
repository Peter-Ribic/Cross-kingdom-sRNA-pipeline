process TARGETFINDER {
    cpus 2
    memory '6 GB'
    tag "$sample_id"
    container "quay.io/biocontainers/targetfinder:1.7--hdfd78af_4"
    
    input:
    tuple val(sample_id), val(fasta_id), val(sequence)
    path host_transcriptome_fasta
    
    output:
    tuple val(sample_id), path("*.log"), emit: log

    
    script:
    """
    targetfinder.pl \
        -s ${sequence} \
        -q ${fasta_id} \
        -t 3 \
        -d ${host_transcriptome_fasta} \
        -c 0.5 > ${fasta_id}.log
    """
}
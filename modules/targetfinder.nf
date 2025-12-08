process TARGETFINDER {
    cpus '4'
    memory '10 GB'
    container "quay.io/biocontainers/targetfinder:1.7--hdfd78af_4"
    
    input:
    val tuples_array
    path host_transcriptome_fasta
    
    output:
    path "*.log", emit: log

    
    script:
 // Build a multi-line bash script
    def commands = new StringBuilder()
    tuples_array.each { t ->
        def (sample_id, fasta_id, sequence) = t
        commands << """
        echo "Running targetfinder for ${fasta_id}..."
        targetfinder.pl \\
            -s ${sequence} \\
            -q ${fasta_id} \\
            -t 4 \\
            -d ${host_transcriptome_fasta} \\
            -c 0.5 > ${fasta_id}.log
        """
    }

    // Pass the commands to bash
    """
    ${commands.toString()}
    """
}
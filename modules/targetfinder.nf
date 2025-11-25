process TARGETFINDER {

    tag "$sample_id"
    container "quay.io/biocontainers/targetfinder:1.7--hdfd78af_4"

    publishDir "results/targetfinder/${sample_id}", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    path host_transcriptome_fasta
    
    output:
    path("*.log"), emit: log
    
    script:
    """
    # Parse FASTA file and extract sequence IDs
    grep '^>' $reads | sed 's/^>//' | while read fasta_id; do
        # Extract the sequence for this ID
        # Use awk to get the sequence corresponding to this header
        sequence=\$(awk -v id="\$fasta_id" '
            BEGIN { found=0 }
            /^>/ { 
                if (\$0 ~ ">" id) { found=1; next } 
                else { found=0 }
            }
            found && !/^>/ { print; exit }
        ' $reads)
        
        # Run targetfinder for each sequence
        targetfinder.pl \\
            -s "\$sequence" \\
            -t ${task.cpus} \\
            -d ${host_transcriptome_fasta} \\
            -q "\$fasta_id" \\
            -c 0.5 >> "target_results.log"
        echo "targetfinder.pl \\
            -s "\$sequence" \\
            -t ${task.cpus} \\
            -d ${host_transcriptome_fasta} \\
            -q "\$fasta_id" \\
            -c 0.5" >> "target_results.log"
    done
    """
}
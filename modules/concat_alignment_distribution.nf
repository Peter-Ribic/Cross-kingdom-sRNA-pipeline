process CONCAT_ALIGNMENT_DISTRIBUTION {
    tag "concat_alignment_distribution"
    container "quay.io/biocontainers/samtools:1.22--h96c455f_0"
    publishDir "results/viruses/virus_alignment_dist", mode: 'symlink'

    input:
    path stats_files 

    output:
    path("virus_alignment_dist.txt"), emit: results

    script:
    """
    header_written=0

    for file in ${stats_files.join(' ')}
    do
        if [ \$header_written -eq 0 ]; then
            head -n 1 \$file > virus_alignment_dist.txt
            header_written=1
        fi
        tail -n +2 \$file >> virus_alignment_dist.txt
    done
    """
}

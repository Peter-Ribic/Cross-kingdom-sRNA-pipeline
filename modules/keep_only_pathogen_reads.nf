process KEEP_ONLY_PATHOGEN_READS {
    tag "$sample_id"
    memory '200 GB'
    cpus 20
    container "quay.io/biocontainers/bioawk:1.0--h577a1d6_13"
    publishDir "results/${sample_id}/keep_pathogen_reads_only", mode: 'copy'

    input:
    tuple val(sample_id), path(reads), val(pathogen_sample_id), path(pathogen_reads)

    output:
    tuple val(sample_id), path("${sample_id}_filtered.fq.gz"), emit: filtered_reads

    // """
    // # Extract unique pathogen sequences (just the sequence strings)
    // seqkit seq -s $pathogen_reads | sort | uniq > pathogen_unique.txt

    // # Keep only reads from hop sample that exactly match pathogen sequences
    // seqkit grep -s -f pathogen_unique.txt --threads 30 -I $reads | gzip > ${sample_id}_filtered.fq.gz
    // """
    script:
    """
    echo "=== DEBUG: Extracting pathogen sequences ==="
    zcat $pathogen_reads | bioawk 'NR % 4 == 2 {print toupper(\$0)}' | tr -d '\\r' | sort | uniq > pathogen_unique.txt
    pathogen_count=\$(wc -l < pathogen_unique.txt)
    echo "Number of unique pathogen sequences: \$pathogen_count"
    echo "Sample pathogen sequences (first 5):"
    head -n 5 pathogen_unique.txt

    echo "=== DEBUG: Filtering host reads ==="
    zcat $reads | bioawk '
        BEGIN {
            host_count = 0
            matched_count = 0
            while ((getline line < "pathogen_unique.txt") > 0) {
                keep[line] = 1
            }
            pathogen_count = length(keep)
        }
        NR % 4 == 1 { name = substr(\$0, 2) }
        NR % 4 == 2 { seq = toupper(\$0); host_count++ }
        NR % 4 == 0 {
            qual = \$0
            if (seq in keep) {
                if (matched_count < 20) {
                    print "@"name"\\n"seq"\\n+\\n"qual
                }
                matched_count++
                print "@"name"\\n"seq"\\n+\\n"qual > "'${sample_id}_filtered.fq.gz.tmp'"
            }
        }
       END {
        pct = (host_count > 0) ? matched_count / host_count * 100 : 0
        printf "Matched %d host reads against %d pathogen sequences, %d reads matched (%.2f%%)\\n", host_count, pathogen_count, matched_count, pct > "/dev/stderr"
    }   
    '

    gzip -c '${sample_id}_filtered.fq.gz.tmp' > ${sample_id}_filtered.fq.gz
    rm '${sample_id}_filtered.fq.gz.tmp'

    echo "=== DEBUG: Finished filtering ==="
    """
}

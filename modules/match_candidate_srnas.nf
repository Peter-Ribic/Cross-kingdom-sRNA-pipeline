process MATCH_WITH_CANDIDATE_SRNAs {
    tag "$sample_id"
    publishDir "results/candidate_reads", mode: 'copy'

    input:
    tuple val(sample_id), path(sample_fq)
    path(candidate_srna_fasta)

    output:
    path "${sample_id}_candidate_reads_stats.txt"

    script:
    """
    # Extract sequences from FASTQ
    zcat ${sample_fq} | awk 'NR%4==2{print}' > ${sample_id}_reads.txt

    # Create associative array of candidate sequences
    awk '/^>/ {id=\$0; next} {seq=id"\t"\$0; print seq}' ${candidate_srna_fasta} > candidates.tsv

    # Count exact matches using join + sort
    # Sort read seqs and candidate seqs for join
    sort -k2,2 candidates.tsv > candidates_sorted.tsv
    sort ${sample_id}_reads.txt | uniq -c | awk '{print ">"NR"\\t"\$2"\\t"\$1}' > reads_sorted.tsv

    awk '
    NR==FNR {cand[\$2]=\$1; next}
    {seq=\$2; count[\$2]=\$3}
    END {
        print "candidate\tmatched_reads"
        for (s in count) {
        if (s in cand)
            print cand[s]"\t"count[s]
        }
    }
    ' candidates_sorted.tsv reads_sorted.tsv > ${sample_id}_candidate_reads_stats.txt

    """
}

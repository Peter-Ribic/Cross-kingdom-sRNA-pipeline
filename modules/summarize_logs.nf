process SUMMARIZE_LOGS {
    tag "summarize_logs"
    publishDir "results/summarized_logs", mode: 'copy'

    input:
    path log_files

    output:
    path "all_logs.tsv", emit: concatenated
    path "summed_by_process.tsv", emit: summed

    script:
    """
    cat ${log_files.join(' ')} > all_logs.tsv

    awk -F '\\t' '{sum[\$1] += \$NF} END {for (p in sum) print p"\\t"sum[p]}' all_logs.tsv \
      > summed_by_process.tsv
    """
}

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
    # 1) Concatenate all one-line TSV logs
    cat ${log_files.join(' ')} > all_logs.tsv

    # 2) Sum last column grouped by first column (process name)
    awk -F '\\t' '{sum[\$1] += \$NF} END {for (p in sum) print p"\\t"sum[p]}' all_logs.tsv \
      > summed_by_process.tsv
    """
}

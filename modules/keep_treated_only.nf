process KEEP_TREATED_ONLY {
    tag "$sample_id"
    memory '200 GB'
    cpus 20
    container "quay.io/biocontainers/bioawk:1.0--h577a1d6_13"
    publishDir "results/main_filtering/treated_minus_control/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(treated_reads), path(control_reads)

    output:
    tuple val(sample_id), path("${sample_id}_treated_only.fq.gz"), emit: filtered_reads
    path("${sample_id}_treated_minus_control.tsv"), emit: log
    path "${task.process}_${sample_id}.tsv", emit: log_info

    script:
    """
    # Extract unique control sequences
    zcat $control_reads | bioawk 'NR % 4 == 2 {print toupper(\$0)}' | tr -d '\\r' | sort | uniq > control_unique.txt
    control_count=\$(wc -l < control_unique.txt)

    # Filter treated reads by removing sequences present in control
    zcat $treated_reads | bioawk -v sid="${sample_id}" '
        BEGIN {
            treated_count = 0
            kept_count = 0
            while ((getline line < "control_unique.txt") > 0) {
                discard[line] = 1
            }
            control_count = length(discard)
        }
        NR % 4 == 1 { header = \$0; name = substr(\$0, 2) }
        NR % 4 == 2 { seq = toupper(\$0); treated_count++ }
        NR % 4 == 3 { plus = \$0 }
        NR % 4 == 0 {
            qual = \$0
            if (!(seq in discard)) {
                kept_count++
                print header"\\n"seq"\\n"plus"\\n"qual > "'${sample_id}_treated_only.fq'"
            }
        }
        END {
            pct = (treated_count > 0) ? kept_count / treated_count * 100 : 0
            
            print "# plot_type: general_stats" > "'${sample_id}_treated_minus_control.tsv'"
            print "sample\\ttotal_reads\\tcontrol_unique\\tkept_reads\\tpercent_kept" >> "'${sample_id}_treated_minus_control.tsv'"
            printf "%s\\t%d\\t%d\\t%d\\t%.2f\\n", sid, treated_count, control_count, kept_count, pct >> "'${sample_id}_treated_minus_control.tsv'"
        }
    '

    # Compress filtered reads
    gzip -c '${sample_id}_treated_only.fq' > ${sample_id}_treated_only.fq.gz

    num_sequences=\$(awk 'END{print NR/4}' ${sample_id}_treated_only.fq)
    echo -e "${task.process}\\t${sample_id}\\t\$num_sequences" > ${task.process}_${sample_id}.tsv
    """
}

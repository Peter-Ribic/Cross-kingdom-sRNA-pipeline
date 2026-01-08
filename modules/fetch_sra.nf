process FETCH_SRA {
    tag "$sample_id"
    conda 'bioconda::sra-tools=3.0.3'
    publishDir "results/fetch_sra/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), val(accessions)

    output:
    tuple val(sample_id), path("*.fastq"), emit: reads
    path "${task.process}_${sample_id}.tsv", emit: log_info

    script:
    """
    for acc in ${accessions.join(' ')}; do
        prefetch \$acc
        fasterq-dump \$acc --split-files -O .

        for f in \${acc}*.fastq; do
            mv "\$f" "${sample_id}_\${f}"
        done
    done

    num_sequences=\$(awk 'END{print NR/4}' *.fastq)
        echo -e "${task.process}\\t${sample_id}\\t\$num_sequences" > ${task.process}_${sample_id}.tsv
    """
}

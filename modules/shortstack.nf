process SHORTSTACK {

    tag "$sample_id"

    container "quay.io/biocontainers/shortstack:4.1.2--hdfd78af_0"

    publishDir "results/shortstack/shortstack_results/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(trimmed_reads)
    path genome_fasta

    output:
    tuple val(sample_id), path("ShortStack_out"), emit: shortstack_out
    tuple val(sample_id), path("ShortStack_out/*condensed.bam"), emit: bam
    tuple val(sample_id), path("ShortStack_out/*condensed.fa"), emit: fasta
    tuple val(sample_id), path("ShortStack_out/Results.txt"), emit: results
    tuple val(sample_id), path("ShortStack_out/${sample_id}_MajorRNA.fa"), emit: majorrna_fasta
    path "${task.process}_${sample_id}.tsv", emit: log_info

    script:
    """
    ShortStack \
        --readfile ${trimmed_reads} \
        --genomefile ${genome_fasta} \
        --outdir ShortStack_out \
        --threads 10 \
        --mincov 400 \
        --dicermin 20

    # MajorRNA FASTA (one per cluster): column 2 = Name, column 11 = MajorRNA
    awk -F \$'\\t' 'NR>1 && \$11 != "" {
        printf(">%s\\n%s\\n", \$2, \$11)
    }' ShortStack_out/Results.txt > ShortStack_out/${sample_id}_MajorRNA.fa

    num_clusters=\$(awk 'END{print NR}' ShortStack_out/Results.gff3)
    echo -e "${task.process}\\t${sample_id}\\t\$num_clusters" > ${task.process}_${sample_id}.tsv
    """
}

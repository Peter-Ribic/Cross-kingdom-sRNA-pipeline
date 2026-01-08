process MIRNA_TARGET {

    tag "$sample_id"

    conda 'bioconda::fasta3 conda-forge::python=3.11'

    publishDir "results/mirna_target/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(shortstack_results)
    path transcriptome_fasta
    path mirna_target_repo

    output:
    tuple val(sample_id), path("${sample_id}_identified_targets.tsv"), emit: results
    tuple val(sample_id), path("${sample_id}_shortstack_majorRNAs.fasta"), emit: extracted_fasta
    path "${task.process}_${sample_id}.tsv", emit: log_info

    script:
    """
    # Extract MajorRNAs from ShortStack Results.txt
    # Use short numeric-safe IDs (c1, c2, ...) to avoid FASTA36/parser truncation
    awk -F '\\t' '
      NR>1 && \$11 != "" {
        id=\$2
        sub(/^Cluster_/, "", id)
        print ">c"id"\\n"\$11
      }
    ' ${shortstack_results} > ${sample_id}_shortstack_majorRNAs.fasta
    
    ssearch36 -i -f -8 -g -3 -E 10000 -T 8 -b 200 -r +4/-3 -n -U -W 10 -N 20000 \
    ${sample_id}_shortstack_majorRNAs.fasta ${transcriptome_fasta} \
    | ${mirna_target_repo}/parse_ssearch.py \
    | ${mirna_target_repo}/parse_mirna_targets.py \
        --E_cutoff 3 \
        --num_mismatch_seed 2 \
        --hsp_cutoff 18 \
        --maximum_alignment_length 24 \
    > ${sample_id}_identified_targets.tsv

    num_targets=\$(awk 'END{print NR-1}' ${sample_id}_identified_targets.tsv)
    echo -e "${task.process}\\t${sample_id}\\t\$num_targets" > ${task.process}_${sample_id}.tsv
    """
}
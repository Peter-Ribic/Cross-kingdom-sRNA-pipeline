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
    path "${task.process}_${sample_id}.tsv", emit: log_info

    script:
    """
    
    # 2) Run ssearch36 and SAVE raw output
    ssearch36 -i -f -8 -g -3 -E 10000 -T 8 -b 200 -r +4/-3 -n -U -W 10 -N 20000 \\
      ${shortstack_results} ${transcriptome_fasta} \\
      > ${sample_id}_ssearch36.raw.txt

    # 3) Parse ssearch output and SAVE parsed intermediate
    ${mirna_target_repo}/parse_ssearch.py \\
      < ${sample_id}_ssearch36.raw.txt \\
      > ${sample_id}_parsed_ssearch.tsv

    # 4) Final miRNA target parsing (unchanged parameters)
    ${mirna_target_repo}/parse_mirna_targets.py \\
        --E_cutoff 3 \\
        --num_mismatch_seed 2 \\
        --hsp_cutoff 18 \\
        --maximum_alignment_length 24 \\
      < ${sample_id}_parsed_ssearch.tsv \\
      > ${sample_id}_identified_targets.tsv

    num_targets=\$(awk 'END{print NR-1}' ${sample_id}_identified_targets.tsv)
    echo -e "${task.process}\\t${sample_id}\\t\$num_targets" > ${task.process}_${sample_id}.tsv
    """
}
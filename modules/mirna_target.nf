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

    script:
    """
     > ${sample_id}_shortstack_majorRNAs.fasta
    ssearch36 -f -8 -g -3 -E 10000 -T 8 -b 200 -r +4/-3 -n -U -W 10 -N 20000 \
    ${sample_id}_shortstack_majorRNAs.fasta ${transcriptome_fasta} \
    | ${mirna_target_repo}/parse_ssearch.py \
    | ${mirna_target_repo}/parse_mirna_targets.py \
        --E_cutoff 2 \
        --num_mismatch_seed 2 \
        --hsp_cutoff 18 \
        --maximum_alignment_length 24 \
    > ${sample_id}_identified_targets.tsv
    """
}

process EGGNOG_ANNOTATION {
    tag "$sample_id"
    publishDir "results/go/eggnog/${sample_id}", mode: 'copy'
    memory '30 GB'
    cpus 10
    conda 'bioconda::eggnog-mapper==2.1.13'

    input:
    tuple val(sample_id), path(targets_faa), path(background_faa)

    output:
    tuple val(sample_id), path("${sample_id}_target_go.txt"), path("${sample_id}_background_go.txt"), emit: go_terms
    path "${sample_id}*", emit: eggnog_outputs
    
    script:
    """
    # 1. Run eggNOG-mapper on target proteins
    emapper.py -i $targets_faa \
        -o ${sample_id}_target \
        --output_dir . \
        -m diamond \
        --cpu ${task.cpus} \
        --itype proteins \
        --go_evidence all \
        --tax_scope 33090 \
        --target_taxa 33090

    
    # 2. Run eggNOG-mapper on background proteins
    emapper.py -i $background_faa \
        -o ${sample_id}_background \
        --output_dir . \
        -m diamond \
        --cpu ${task.cpus} \
        --itype proteins \
        --go_evidence all \
        --tax_scope 33090 \
        --target_taxa 33090


    # 3) Extract GO terms from the GOs column (column 10)
    #    IMPORTANT:
    #      - skip comment/header lines starting with '#'
    #      - keep query ID for EVERY GO term
    #      - split comma-separated list into one GO per line
    awk -F'\\t' 'BEGIN{OFS="\\t"}
      \$1 !~ /^#/ && \$10 != "-" && \$10 != "" {
        n=split(\$10, a, /,/)
        for(i=1;i<=n;i++){
          g=a[i]
          sub(/^ +| +\$/, "", g)
          if(g ~ /^GO:[0-9]+/) print \$1, g
        }
      }' ${sample_id}_target.emapper.annotations \\
      | sort -u > ${sample_id}_target_go.txt

    awk -F'\\t' 'BEGIN{OFS="\\t"}
      \$1 !~ /^#/ && \$10 != "-" && \$10 != "" {
        n=split(\$10, a, /,/)
        for(i=1;i<=n;i++){
          g=a[i]
          sub(/^ +| +\$/, "", g)
          if(g ~ /^GO:[0-9]+/) print \$1, g
        }
      }' ${sample_id}_background.emapper.annotations \\
      | sort -u > ${sample_id}_background_go.txt

    echo "Target proteins with >=1 GO: \$(cut -f1 ${sample_id}_target_go.txt | sort -u | wc -l)"
    echo "Background proteins with >=1 GO: \$(cut -f1 ${sample_id}_background_go.txt | sort -u | wc -l)"

    # quick sanity: show any 'GO-only' lines (should be 0)
    echo "GO-only lines in target_go (should be 0): \$(awk -F'\\t' 'NF==1{c++} END{print c+0}' ${sample_id}_target_go.txt)"
    echo "GO-only lines in background_go (should be 0): \$(awk -F'\\t' 'NF==1{c++} END{print c+0}' ${sample_id}_background_go.txt)"
    """
}

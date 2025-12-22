process ANNOTATE_TARGETS {
    tag "$sample_id"
    container "docker://ncbi/edirect:latest"

    publishDir "results/target_annotations/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(target_file)
    path(mrna_fasta)

    output:
    path("${sample_id}_target_accessions.txt"), emit: target_ids
    path("${sample_id}_functional_annotations.txt"), emit: annotations
    path("${sample_id}_defense_targets.txt"), emit: defense_targets
    path("${sample_id}_summary.txt"), emit: summary
    path("${sample_id}_enrichment_analysis.txt"), emit: enrichment_results

    script:
    """
    set -euo pipefail

    # MiRNATarget output: subject_id like:
    # lcl|XM_062233053.1_cds_XP_062089037.1_1_1
    #
    # We will extract XP_* (protein) accessions preferentially for functional annotation.
    # Fallback: extract XM_* if no XP_* present.

    # 1) Extract unique accessions from subject_id
    # Skip header lines starting with '#'
    awk -F '\\t' 'BEGIN{OFS="\\t"} \$0 !~ /^#/ && NF>1 {print \$2}' ${target_file} \\
      | awk '
          {
            # find XP_ accession in the subject_id
            if (match(\$0, /XP_[0-9]+\\.[0-9]+/)) { print substr(\$0, RSTART, RLENGTH); next }
            # fallback: XM_ accession
            if (match(\$0, /XM_[0-9]+\\.[0-9]+/)) { print substr(\$0, RSTART, RLENGTH); next }
          }
        ' \\
      | sort -u > ${sample_id}_target_accessions.txt

    unique_targets=\$(wc -l < ${sample_id}_target_accessions.txt || echo "0")

    # 2) Fetch functional annotations
    # Output columns: Accession \\t Title/Description \\t Organism
    echo -e "Accession\\tDescription\\tOrganism" > ${sample_id}_functional_annotations.txt

    while read acc; do
      [[ -z "\$acc" ]] && continue

      # Decide which NCBI db to query
      db="protein"
      if echo "\$acc" | grep -q '^XM_'; then
        db="nuccore"
      fi

      # Fetch docsum and extract a compact summary
      # For protein: Title + Organism are usually present.
      # For nuccore: Title often includes gene/product context too.
      if out=\$(esummary -db "\$db" -id "\$acc" 2>/dev/null); then
        # Try common docsum fields; if missing, leave blank rather than crash.
        title=\$(echo "\$out" | xtract -pattern DocumentSummary -element Title 2>/dev/null | head -n 1 || true)
        org=\$(echo "\$out" | xtract -pattern DocumentSummary -element Organism 2>/dev/null | head -n 1 || true)

        # Some protein docsums use different fields; try an additional fallback
        if [[ -z "\$title" ]]; then
          title=\$(echo "\$out" | xtract -pattern DocumentSummary -element Description 2>/dev/null | head -n 1 || true)
        fi

        if [[ -z "\$org" ]]; then
          org=\$(echo "\$out" | xtract -pattern DocumentSummary -element Organism/ScientificName 2>/dev/null | head -n 1 || true)
        fi

        # If still empty, flag it
        if [[ -z "\$title" && -z "\$org" ]]; then
          echo -e "\$acc\\tERROR: Failed to parse docsum\\t" >> ${sample_id}_functional_annotations.txt
        else
          # Replace tabs/newlines just in case
          title=\$(echo "\$title" | tr '\\t\\r\\n' ' ' | sed 's/  */ /g')
          org=\$(echo "\$org" | tr '\\t\\r\\n' ' ' | sed 's/  */ /g')
          echo -e "\$acc\\t\$title\\t\$org" >> ${sample_id}_functional_annotations.txt
        fi
      else
        echo -e "\$acc\\tERROR: esummary failed\\t" >> ${sample_id}_functional_annotations.txt
      fi
    done < ${sample_id}_target_accessions.txt

    # Count successfully annotated (exclude header + ERROR)
    total_annotated_targets=\$(tail -n +2 ${sample_id}_functional_annotations.txt | grep -v '^.*\\tERROR' | wc -l || echo "0")
    failed_annotations=\$(tail -n +2 ${sample_id}_functional_annotations.txt | grep -c '^.*\\tERROR' || echo "0")

    # 3) Defense keyword scan (based on annotation description)
    echo "Defense-Related Targets" > ${sample_id}_defense_targets.txt
    echo "======================" >> ${sample_id}_defense_targets.txt
    echo -e "Accession\\tDescription\\tOrganism" >> ${sample_id}_defense_targets.txt

    if [[ -s ${sample_id}_functional_annotations.txt ]]; then
      # Keywords are broad; adjust for your plant-pathogen context as needed
      grep -Ei "\\t.*(defen|resistan|pathogenesis|PR[0-9]|chitinase|glucanase|NBS|LRR|RLK|RLP|WRKY|jasmon|salicyl|ethylene|phenylpropanoid|lignin|cell wall|pectin|peroxidase).*" \\
        ${sample_id}_functional_annotations.txt \\
        | grep -v '^Accession\\t' >> ${sample_id}_defense_targets.txt || true
    fi

    defense_targets_count=\$(tail -n +4 ${sample_id}_defense_targets.txt | wc -l || echo "0")

    # 4) "Background" defense-like prevalence in mRNA FASTA headers (rough proxy)
    echo "Enrichment Analysis Results" > ${sample_id}_enrichment_analysis.txt
    echo "===========================" >> ${sample_id}_enrichment_analysis.txt
    echo "" >> ${sample_id}_enrichment_analysis.txt

    echo "1. mRNA FASTA header keyword scan (proxy background):" >> ${sample_id}_enrichment_analysis.txt
    echo "---------------------------------------------------" >> ${sample_id}_enrichment_analysis.txt

    total_mrna_transcripts=\$(grep -c "^>" ${mrna_fasta} || echo "0")
    echo "Total transcripts in mRNA FASTA: \$total_mrna_transcripts" >> ${sample_id}_enrichment_analysis.txt

    defense_in_mrna=\$(grep "^>" ${mrna_fasta} \\
      | grep -Eic "(defen|resistan|pathogenesis|PR[0-9]|chitinase|glucanase|NBS|LRR|RLK|RLP|WRKY|jasmon|salicyl|ethylene|phenylpropanoid|lignin|cell wall|pectin|peroxidase)" \\
      || echo "0")
    echo "Defense-keyword headers in mRNA FASTA: \$defense_in_mrna" >> ${sample_id}_enrichment_analysis.txt

    if [[ \$total_mrna_transcripts -gt 0 ]]; then
      mrna_defense_prop=\$(awk -v d="\$defense_in_mrna" -v t="\$total_mrna_transcripts" 'BEGIN {printf "%.4f", d/t}')
    else
      mrna_defense_prop="0.0000"
    fi
    echo "Proportion (background): \$mrna_defense_prop" >> ${sample_id}_enrichment_analysis.txt
    echo "" >> ${sample_id}_enrichment_analysis.txt

    echo "2. Target annotation keyword scan:" >> ${sample_id}_enrichment_analysis.txt
    echo "-------------------------------" >> ${sample_id}_enrichment_analysis.txt
    echo "Total successfully annotated targets: \$total_annotated_targets" >> ${sample_id}_enrichment_analysis.txt
    echo "Defense-related targets (keyword hits): \$defense_targets_count" >> ${sample_id}_enrichment_analysis.txt

    if [[ \$total_annotated_targets -gt 0 ]]; then
      target_defense_prop=\$(awk -v d="\$defense_targets_count" -v t="\$total_annotated_targets" 'BEGIN {printf "%.4f", d/t}')
    else
      target_defense_prop="0.0000"
    fi
    echo "Proportion (targets): \$target_defense_prop" >> ${sample_id}_enrichment_analysis.txt
    echo "" >> ${sample_id}_enrichment_analysis.txt

    echo "3. Enrichment ratio (targets/background):" >> ${sample_id}_enrichment_analysis.txt
    echo "----------------------------------------" >> ${sample_id}_enrichment_analysis.txt

    enrichment_ratio=""
    if [[ \$(awk -v m="\$mrna_defense_prop" 'BEGIN{print (m>0)?1:0}') -eq 1 ]]; then
      enrichment_ratio=\$(awk -v t="\$target_defense_prop" -v m="\$mrna_defense_prop" 'BEGIN {printf "%.4f", t/m}')
      echo "Enrichment ratio: \$enrichment_ratio" >> ${sample_id}_enrichment_analysis.txt
    else
      echo "Background proportion is 0; enrichment ratio not defined." >> ${sample_id}_enrichment_analysis.txt
    fi

    # 5) Summary
    echo "Target Annotation Summary for ${sample_id}" > ${sample_id}_summary.txt
    echo "======================================" >> ${sample_id}_summary.txt
    echo "Input MiRNATarget file: ${target_file}" >> ${sample_id}_summary.txt
    echo "Unique target accessions found: \$unique_targets" >> ${sample_id}_summary.txt
    echo "Successfully annotated: \$total_annotated_targets" >> ${sample_id}_summary.txt
    echo "Failed annotations: \$failed_annotations" >> ${sample_id}_summary.txt
    echo "Defense-related targets: \$defense_targets_count" >> ${sample_id}_summary.txt
    echo "" >> ${sample_id}_summary.txt
    echo "Enrichment (keyword-based proxy):" >> ${sample_id}_summary.txt
    echo "Background defense: \$defense_in_mrna/\$total_mrna_transcripts (\$mrna_defense_prop)" >> ${sample_id}_summary.txt
    echo "Target defense: \$defense_targets_count/\$total_annotated_targets (\$target_defense_prop)" >> ${sample_id}_summary.txt
    if [[ -n "\$enrichment_ratio" ]]; then
      echo "Enrichment ratio: \$enrichment_ratio" >> ${sample_id}_summary.txt
    fi
    """
}

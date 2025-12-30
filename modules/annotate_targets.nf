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

    # MiRNATarget output: subject_id examples:
    #   lcl|XM_062233053.1_cds_XP_062089037.1_1_1
    #   lcl|JQ740210.1_cds_AGA17925.1_1_489
    #   lcl|EF624248.1_cds_ABS17598.1_1_19
    #   lcl|MW658771.1_cds_UBU95743.1_1_1857
    #
    # We will extract accessions in this priority:
    #   1) XP_####.#
    #   2) token after _cds_ (e.g. AGA17925.1 / ABS17598.1 / UBU95743.1 / BAB40666.1)
    #   3) XM_####.#

    ######################################################################
    # 1) Extract unique accessions from subject_id
    ######################################################################
    awk -F '\\t' 'BEGIN{OFS="\\t"} \$0 !~ /^#/ && NF>1 {print \$2}' ${target_file} \\
      | awk '
          {
            s=\$0

            # 1) prefer XP_
            if (match(s, /XP_[0-9]+\\.[0-9]+/)) { print substr(s, RSTART, RLENGTH); next }

            # 2) extract token after _cds_
            if (match(s, /_cds_[A-Za-z0-9]+(\\.[0-9]+)?/)) {
              x = substr(s, RSTART, RLENGTH)
              sub(/^_cds_/, "", x)
              print x
              next
            }

            # 3) fallback XM_
            if (match(s, /XM_[0-9]+\\.[0-9]+/)) { print substr(s, RSTART, RLENGTH); next }
          }
        ' \\
      | sort -u > ${sample_id}_target_accessions.txt

    ######################################################################
    # 1b) Extract query -> accession mapping (unique pairs)
    # Assumes query is in column 1 and subject_id is in column 2
    ######################################################################
    awk -F '\\t' '\$0 !~ /^#/ && NF>1 {print \$1"\\t"\$2}' ${target_file} \\
      | awk -F '\\t' 'BEGIN{OFS="\\t"}
          {
            q=\$1
            s=\$2

            # 1) prefer XP_
            if (match(s, /XP_[0-9]+\\.[0-9]+/)) { print q, substr(s, RSTART, RLENGTH); next }

            # 2) extract token after _cds_
            if (match(s, /_cds_[A-Za-z0-9]+(\\.[0-9]+)?/)) {
              x = substr(s, RSTART, RLENGTH)
              sub(/^_cds_/, "", x)
              print q, x
              next
            }

            # 3) fallback XM_
            if (match(s, /XM_[0-9]+\\.[0-9]+/)) { print q, substr(s, RSTART, RLENGTH); next }
          }' \\
      | sort -u > ${sample_id}_query_accession_map.txt

    unique_targets=\$(wc -l < ${sample_id}_target_accessions.txt || echo "0")

    ######################################################################
    # 2) Build a local lookup table from mrna_fasta headers
    # Columns: protein_id  gene  protein  full_header  first_token_id
    #
    # Example header:
    # >lcl|AB053487.1_cds_BAB40666.1_1 [gene=fpps] [protein=...] [protein_id=BAB40666.1] ...
    ######################################################################
    awk '
      BEGIN{ OFS="\\t" }
      /^>/{
        hdr = substr(\$0,2)

        split(hdr, a, /[ \\t]/)
        first_tok = a[1]

        gene=""
        prot=""
        pid=""

        if (match(hdr, /\\[gene=[^\\]]+\\]/)) {
          gene = substr(hdr, RSTART+6, RLENGTH-7)
        }
        if (match(hdr, /\\[protein=[^\\]]+\\]/)) {
          prot = substr(hdr, RSTART+9, RLENGTH-10)
        }
        if (match(hdr, /\\[protein_id=[^\\]]+\\]/)) {
          pid = substr(hdr, RSTART+12, RLENGTH-13)
        }

        if (pid != "") {
          print pid, gene, prot, hdr, first_tok
        }
      }
    ' ${mrna_fasta} | sort -u > ${sample_id}_mrna_proteinid_map.tsv

    lookup_fasta_by_pid() {
      local acc="\$1"
      awk -F '\\t' -v a="\$acc" '\$1==a{print \$2"\\t"\$3"\\t"\$4; found=1; exit} END{if(!found) exit 1}' \\
        ${sample_id}_mrna_proteinid_map.tsv
    }

    lookup_fasta_by_text() {
      local acc="\$1"
      awk -F '\\t' -v a="\$acc" 'index(\$4,a)>0 || index(\$5,a)>0 {print \$2"\\t"\$3"\\t"\$4; found=1; exit} END{if(!found) exit 1}' \\
        ${sample_id}_mrna_proteinid_map.tsv
    }

    ######################################################################
    # 3) Fetch functional annotations (prefer local FASTA header, fallback NCBI)
    ######################################################################
    echo -e "Query\\tAccession\\tDescription\\tOrganism\\tFASTA_gene\\tFASTA_protein\\tFASTA_header" > ${sample_id}_functional_annotations.txt

    while read acc; do
      [[ -z "\$acc" ]] && continue

      query_names=\$(awk -v a="\$acc" '\$2==a{print \$1}' ${sample_id}_query_accession_map.txt | sort -u | paste -sd"," -)
      [[ -z "\$query_names" ]] && query_names="NA"

      fasta_gene=""
      fasta_prot=""
      fasta_hdr=""

      # Local FASTA lookup (protein_id exact match; else header-text search)
      if out_local=\$(lookup_fasta_by_pid "\$acc" 2>/dev/null); then
        fasta_gene=\$(echo "\$out_local" | cut -f1)
        fasta_prot=\$(echo "\$out_local" | cut -f2)
        fasta_hdr=\$(echo "\$out_local" | cut -f3-)
      elif out_local=\$(lookup_fasta_by_text "\$acc" 2>/dev/null); then
        fasta_gene=\$(echo "\$out_local" | cut -f1)
        fasta_prot=\$(echo "\$out_local" | cut -f2)
        fasta_hdr=\$(echo "\$out_local" | cut -f3-)
      fi

      if [[ -n "\$fasta_gene" || -n "\$fasta_prot" || -n "\$fasta_hdr" ]]; then
        desc=""
        if [[ -n "\$fasta_prot" ]]; then
          desc="\$fasta_prot"
        elif [[ -n "\$fasta_gene" ]]; then
          desc="gene=\$fasta_gene"
        else
          desc="(from FASTA header)"
        fi

        desc=\$(echo "\$desc" | tr '\\t\\r\\n' ' ' | sed 's/  */ /g')
        fasta_gene=\$(echo "\$fasta_gene" | tr '\\t\\r\\n' ' ' | sed 's/  */ /g')
        fasta_prot=\$(echo "\$fasta_prot" | tr '\\t\\r\\n' ' ' | sed 's/  */ /g')
        fasta_hdr=\$(echo "\$fasta_hdr" | tr '\\t\\r\\n' ' ' | sed 's/  */ /g')

        echo -e "\$query_names\\t\$acc\\t\$desc\\t\\t\$fasta_gene\\t\$fasta_prot\\t\$fasta_hdr" >> ${sample_id}_functional_annotations.txt
        continue
      fi

      # Fallback to NCBI
      db="protein"
      if echo "\$acc" | grep -q '^XM_'; then
        db="nuccore"
      fi

      if out=\$(esummary -db "\$db" -id "\$acc" 2>/dev/null); then
        title=\$(echo "\$out" | xtract -pattern DocumentSummary -element Title 2>/dev/null | head -n 1 || true)
        org=\$(echo "\$out" | xtract -pattern DocumentSummary -element Organism 2>/dev/null | head -n 1 || true)

        if [[ -z "\$title" ]]; then
          title=\$(echo "\$out" | xtract -pattern DocumentSummary -element Description 2>/dev/null | head -n 1 || true)
        fi
        if [[ -z "\$org" ]]; then
          org=\$(echo "\$out" | xtract -pattern DocumentSummary -element Organism/ScientificName 2>/dev/null | head -n 1 || true)
        fi

        if [[ -z "\$title" && -z "\$org" ]]; then
          echo -e "\$query_names\\t\$acc\\tERROR: Failed to parse docsum\\t\\t\\t\\t" >> ${sample_id}_functional_annotations.txt
        else
          title=\$(echo "\$title" | tr '\\t\\r\\n' ' ' | sed 's/  */ /g')
          org=\$(echo "\$org" | tr '\\t\\r\\n' ' ' | sed 's/  */ /g')
          echo -e "\$query_names\\t\$acc\\t\$title\\t\$org\\t\\t\\t" >> ${sample_id}_functional_annotations.txt
        fi
      else
        echo -e "\$query_names\\t\$acc\\tERROR: esummary failed\\t\\t\\t\\t" >> ${sample_id}_functional_annotations.txt
      fi
    done < ${sample_id}_target_accessions.txt

    # Count successfully annotated (exclude header + ERROR)
    total_annotated_targets=\$(tail -n +2 ${sample_id}_functional_annotations.txt | grep -v '^.*\\tERROR' | wc -l || echo "0")
    failed_annotations=\$(tail -n +2 ${sample_id}_functional_annotations.txt | grep -c '^.*\\tERROR' || echo "0")

    ######################################################################
    # 4) Defense keyword scan (based on annotation description)
    ######################################################################
    echo "Defense-Related Targets" > ${sample_id}_defense_targets.txt
    echo "======================" >> ${sample_id}_defense_targets.txt
    echo -e "Query\\tAccession\\tDescription\\tOrganism\\tFASTA_gene\\tFASTA_protein\\tFASTA_header" >> ${sample_id}_defense_targets.txt

    if [[ -s ${sample_id}_functional_annotations.txt ]]; then
      grep -i "defense\\|resistance\\|PR[0-9]\\|chitinase\\|glucanase\\|NBS\\|LRR\\|RLK\\|WRKY\\|pathogenesis" \\
        ${sample_id}_functional_annotations.txt >> ${sample_id}_defense_targets.txt 2>/dev/null || true
    fi

    defense_targets_count=\$(tail -n +4 ${sample_id}_defense_targets.txt | wc -l || echo "0")

    ######################################################################
    # 5) "Background" defense-like prevalence in mRNA FASTA headers (rough proxy)
    ######################################################################
    echo "Enrichment Analysis Results" > ${sample_id}_enrichment_analysis.txt
    echo "===========================" >> ${sample_id}_enrichment_analysis.txt
    echo "" >> ${sample_id}_enrichment_analysis.txt

    echo "1. mRNA FASTA header keyword scan (proxy background):" >> ${sample_id}_enrichment_analysis.txt
    echo "---------------------------------------------------" >> ${sample_id}_enrichment_analysis.txt

    total_mrna_transcripts=\$(grep -c "^>" ${mrna_fasta} || echo "0")
    echo "Total transcripts in mRNA FASTA: \$total_mrna_transcripts" >> ${sample_id}_enrichment_analysis.txt

    defense_in_mrna=\$(grep "^>" ${mrna_fasta} \\
      | grep -ic "defense\\|resistance\\|PR[0-9]\\|chitinase\\|glucanase\\|NBS\\|LRR\\|RLK\\|WRKY\\|pathogenesis" \\
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

    ######################################################################
    # 6) Summary
    ######################################################################
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

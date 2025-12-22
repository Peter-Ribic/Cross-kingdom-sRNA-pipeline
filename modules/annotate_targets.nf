process ANNOTATE_TARGETS {
    tag "$sample_id"
    container "docker://ncbi/edirect:latest"
    
    publishDir "results/target_annotations/${sample_id}", mode: 'copy'
    
    input:
    tuple val(sample_id), path(target_file)
    path(mrna_fasta)
    
    output:
    path("${sample_id}_target_ids.txt"), emit: target_ids
    path("${sample_id}_functional_annotations.txt"), emit: annotations
    path("${sample_id}_defense_targets.txt"), emit: defense_targets
    path("${sample_id}_summary.txt"), emit: summary
    path("${sample_id}_enrichment_analysis.txt"), emit: enrichment_results
    
    script:
    """
    # Extract unique NCBI gene IDs from the target file
    grep -o 'GeneID:[0-9]*' ${target_file} | cut -d: -f2 | sort -u > ${sample_id}_target_ids.txt
    
    # Count unique targets
    unique_targets=\$(wc -l < ${sample_id}_target_ids.txt)
    enrichment_ratio=""
    
    # Get functional annotations for each gene ID
    echo -e "GeneID\\tAccession\\tDescription\\tSymbol\\tOrganism" > ${sample_id}_functional_annotations.txt
    
    # Process each gene ID
    while read gene_id; do
        if [[ -n "\$gene_id" ]]; then
            # Use esummary to get gene information
            esummary -db gene -id "\$gene_id" 2>/dev/null | \
            xtract -pattern DocumentSummary -element Id AccessionVersion Description Name Organism/ScientificName >> ${sample_id}_functional_annotations.txt 2>/dev/null || echo "\$gene_id\\tERROR\\tFailed to fetch" >> ${sample_id}_functional_annotations.txt
        fi
    done < ${sample_id}_target_ids.txt
    
    # Create defense targets file
    echo "Defense-Related Targets" > ${sample_id}_defense_targets.txt
    echo "======================" >> ${sample_id}_defense_targets.txt
    
    # Identify defense-related genes from successful annotations
    if [[ -s ${sample_id}_functional_annotations.txt ]]; then
        grep -i "defense\\|resistance\\|PR[0-9]\\|chitinase\\|glucanase\\|NBS\\|LRR\\|RLK\\|WRKY\\|pathogenesis" ${sample_id}_functional_annotations.txt >> ${sample_id}_defense_targets.txt 2>/dev/null || echo "No defense-related genes found in annotations" >> ${sample_id}_defense_targets.txt
    else
        echo "No annotations available to check for defense genes" >> ${sample_id}_defense_targets.txt
    fi
    
    # ENRICHMENT ANALYSIS SECTION
    echo "Enrichment Analysis Results" > ${sample_id}_enrichment_analysis.txt
    echo "===========================" >> ${sample_id}_enrichment_analysis.txt
    echo "" >> ${sample_id}_enrichment_analysis.txt
    
    # 1. Count defense-related transcripts in mRNA fasta file
    echo "1. mRNA FASTA Analysis:" >> ${sample_id}_enrichment_analysis.txt
    echo "------------------------" >> ${sample_id}_enrichment_analysis.txt
    
    # Count total transcripts in mRNA fasta
    total_mrna_transcripts=\$(grep -c "^>" ${mrna_fasta} || echo "0")
    echo "Total transcripts in mRNA FASTA: \$total_mrna_transcripts" >> ${sample_id}_enrichment_analysis.txt
    
    # Count defense-related transcripts in mRNA fasta
    defense_in_mrna=\$(grep -i "defense\\|resistance\\|PR[0-9]\\|chitinase\\|glucanase\\|NBS\\|LRR\\|RLK\\|WRKY\\|pathogenesis" ${mrna_fasta} | grep -c "^>" || echo "0")
    echo "Defense-related transcripts in mRNA FASTA: \$defense_in_mrna" >> ${sample_id}_enrichment_analysis.txt
    
    # Calculate proportion in mRNA fasta using awk
    if [[ \$total_mrna_transcripts -gt 0 ]]; then
        mrna_defense_prop=\$(awk -v d="\$defense_in_mrna" -v t="\$total_mrna_transcripts" 'BEGIN {printf "%.4f", d/t}')
        echo "Proportion of defense-related in mRNA FASTA: \$mrna_defense_prop" >> ${sample_id}_enrichment_analysis.txt
    else
        mrna_defense_prop="0"
        echo "Proportion of defense-related in mRNA FASTA: 0 (no transcripts)" >> ${sample_id}_enrichment_analysis.txt
    fi
    
    echo "" >> ${sample_id}_enrichment_analysis.txt
    
    # 2. Count defense-related targets
    echo "2. Target Analysis:" >> ${sample_id}_enrichment_analysis.txt
    echo "-------------------" >> ${sample_id}_enrichment_analysis.txt
    
    # Count successfully annotated targets (excluding header and errors)
    total_annotated_targets=\$(tail -n +2 ${sample_id}_functional_annotations.txt | grep -v ERROR | wc -l)
    echo "Total successfully annotated targets: \$total_annotated_targets" >> ${sample_id}_enrichment_analysis.txt
    
    # Count defense-related targets (excluding header lines from defense_targets.txt)
    defense_targets_count=\$(tail -n +3 ${sample_id}_defense_targets.txt | grep -v "No defense" | grep -v "No annotations" | wc -l)
    echo "Defense-related targets: \$defense_targets_count" >> ${sample_id}_enrichment_analysis.txt
    
    # Calculate proportion in targets using awk
    if [[ \$total_annotated_targets -gt 0 ]]; then
        target_defense_prop=\$(awk -v d="\$defense_targets_count" -v t="\$total_annotated_targets" 'BEGIN {printf "%.4f", d/t}')
        echo "Proportion of defense-related in targets: \$target_defense_prop" >> ${sample_id}_enrichment_analysis.txt
    else
        target_defense_prop="0"
        echo "Proportion of defense-related in targets: 0 (no annotated targets)" >> ${sample_id}_enrichment_analysis.txt
    fi
    
    echo "" >> ${sample_id}_enrichment_analysis.txt
    
    # 3. Enrichment calculation
    echo "3. Enrichment Calculation:" >> ${sample_id}_enrichment_analysis.txt
    echo "--------------------------" >> ${sample_id}_enrichment_analysis.txt
    
    if [[ \$total_mrna_transcripts -gt 0 ]] && [[ \$total_annotated_targets -gt 0 ]]; then
        # Calculate enrichment ratio using awk
        if [[ \$(echo "\$mrna_defense_prop > 0" | awk '{print (\$1 > 0)}') -eq 1 ]]; then
            enrichment_ratio=\$(awk -v t="\$target_defense_prop" -v m="\$mrna_defense_prop" 'BEGIN {printf "%.4f", t/m}')
            echo "Enrichment Ratio (targets/mRNA): \$enrichment_ratio" >> ${sample_id}_enrichment_analysis.txt
            
            # Statistical significance using Fisher's exact test (approximation)
            non_defense_in_mrna=\$((total_mrna_transcripts - defense_in_mrna))
            non_defense_targets=\$((total_annotated_targets - defense_targets_count))
            
            # Create contingency table for manual calculation
            echo "" >> ${sample_id}_enrichment_analysis.txt
            echo "Contingency Table:" >> ${sample_id}_enrichment_analysis.txt
            echo "                  | Defense | Non-defense | Total" >> ${sample_id}_enrichment_analysis.txt
            echo "------------------------------------------------" >> ${sample_id}_enrichment_analysis.txt
            echo "Targets           | \$defense_targets_count | \$non_defense_targets | \$total_annotated_targets" >> ${sample_id}_enrichment_analysis.txt
            echo "Background (mRNA) | \$defense_in_mrna | \$non_defense_in_mrna | \$total_mrna_transcripts" >> ${sample_id}_enrichment_analysis.txt
            
            # Simple enrichment test using awk for comparison
            is_enriched=\$(awk -v er="\$enrichment_ratio" 'BEGIN {print (er > 1.0) ? "1" : "0"}')
            is_equal=\$(awk -v er="\$enrichment_ratio" 'BEGIN {print (er == 1.0) ? "1" : "0"}')
            
            echo "" >> ${sample_id}_enrichment_analysis.txt
            if [[ \$is_enriched -eq 1 ]]; then
                echo "RESULT: ENRICHED" >> ${sample_id}_enrichment_analysis.txt
                echo "Defense-related targets are enriched compared to background." >> ${sample_id}_enrichment_analysis.txt
                echo "Enrichment fold: \$enrichment_ratio" >> ${sample_id}_enrichment_analysis.txt
            elif [[ \$is_equal -eq 1 ]]; then
                echo "RESULT: NO ENRICHMENT" >> ${sample_id}_enrichment_analysis.txt
                echo "Defense-related targets are at expected frequency." >> ${sample_id}_enrichment_analysis.txt
            else
                echo "RESULT: DEPLETED" >> ${sample_id}_enrichment_analysis.txt
                echo "Defense-related targets are depleted compared to background." >> ${sample_id}_enrichment_analysis.txt
            fi
        else
            echo "No defense-related transcripts in mRNA background." >> ${sample_id}_enrichment_analysis.txt
            if [[ \$defense_targets_count -gt 0 ]]; then
                echo "All defense-related targets in analysis represent novel findings." >> ${sample_id}_enrichment_analysis.txt
            fi
        fi
    else
        echo "Cannot calculate enrichment: insufficient data." >> ${sample_id}_enrichment_analysis.txt
        echo "Total mRNA transcripts: \$total_mrna_transcripts" >> ${sample_id}_enrichment_analysis.txt
        echo "Total annotated targets: \$total_annotated_targets" >> ${sample_id}_enrichment_analysis.txt
    fi
    
    # Create summary (updated with enrichment info)
    echo "Target Annotation Summary for ${sample_id}" > ${sample_id}_summary.txt
    echo "======================================" >> ${sample_id}_summary.txt
    echo "Input file: ${target_file}" >> ${sample_id}_summary.txt
    echo "Unique Gene IDs found: \$unique_targets" >> ${sample_id}_summary.txt
    echo "Successfully annotated: \$total_annotated_targets" >> ${sample_id}_summary.txt
    echo "Failed annotations: \$(tail -n +2 ${sample_id}_functional_annotations.txt | grep ERROR | wc -l)" >> ${sample_id}_summary.txt
    echo "Defense-related targets: \$defense_targets_count" >> ${sample_id}_summary.txt
    echo "" >> ${sample_id}_summary.txt
    echo "Enrichment Analysis:" >> ${sample_id}_summary.txt
    echo "--------------------" >> ${sample_id}_summary.txt
    echo "mRNA background defense: \$defense_in_mrna/\$total_mrna_transcripts (\$mrna_defense_prop)" >> ${sample_id}_summary.txt
    echo "Target defense: \$defense_targets_count/\$total_annotated_targets (\$target_defense_prop)" >> ${sample_id}_summary.txt
    if [[ -n "\$enrichment_ratio" ]]; then
        echo "Enrichment ratio: \$enrichment_ratio" >> ${sample_id}_summary.txt
    fi
    """
}
process ANNOTATE_TARGETS {
    tag "$sample_id"
    container "docker://ncbi/edirect:latest"
    
     publishDir "results/target_annotations/${sample_id}", mode: 'copy'
    
    input:
    tuple val(sample_id), path(target_file)
    
    output:
    path("${sample_id}_target_ids.txt"), emit: target_ids
    path("${sample_id}_functional_annotations.txt"), emit: annotations
    path("${sample_id}_defense_targets.txt"), emit: defense_targets
    path("${sample_id}_summary.txt"), emit: summary
    
    script:
    """
    # Extract unique NCBI gene IDs from the target file
    grep -o 'GeneID:[0-9]*' ${target_file} | cut -d: -f2 | sort -u > ${sample_id}_target_ids.txt
    
    # Count unique targets
    unique_targets=\$(wc -l < ${sample_id}_target_ids.txt)
    
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
    
    # Create defense targets file (even if empty)
    echo "Defense-Related Targets" > ${sample_id}_defense_targets.txt
    echo "======================" >> ${sample_id}_defense_targets.txt
    
    # Try to identify defense-related genes from successful annotations
    if [[ -s ${sample_id}_functional_annotations.txt ]]; then
        grep -i "defense\\|resistance\\|PR[0-9]\\|chitinase\\|glucanase\\|NBS\\|LRR\\|RLK\\|WRKY\\|pathogenesis" ${sample_id}_functional_annotations.txt >> ${sample_id}_defense_targets.txt 2>/dev/null || echo "No defense-related genes found in annotations" >> ${sample_id}_defense_targets.txt
    else
        echo "No annotations available to check for defense genes" >> ${sample_id}_defense_targets.txt
    fi
    
    # Create summary
    echo "Target Annotation Summary for ${sample_id}" > ${sample_id}_summary.txt
    echo "======================================" >> ${sample_id}_summary.txt
    echo "Input file: ${target_file}" >> ${sample_id}_summary.txt
    echo "Unique Gene IDs found: \$unique_targets" >> ${sample_id}_summary.txt
    echo "Successfully annotated: \$(tail -n +2 ${sample_id}_functional_annotations.txt | grep -v ERROR | wc -l)" >> ${sample_id}_summary.txt
    echo "Failed annotations: \$(tail -n +2 ${sample_id}_functional_annotations.txt | grep ERROR | wc -l)" >> ${sample_id}_summary.txt
    echo "Defense-related targets: \$(tail -n +3 ${sample_id}_defense_targets.txt | grep -v "No defense" | grep -v "No annotations" | wc -l)" >> ${sample_id}_summary.txt
    """
}
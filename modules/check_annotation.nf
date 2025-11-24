process CHECK_ANNOTATION {
    tag "$sample_id"
    
    container "quay.io/biocontainers/bedtools:2.31.1--h13024bc_3"
    
    publishDir "results/shortstack_annotated/${sample_id}", mode: 'copy'
    
    input:
    tuple val(sample_id), path(shortstack_dir)
    path genome_gff3
    
    output:
    tuple val(sample_id), path("${sample_id}_annotated_regions.txt"), emit: annotated_regions
    path "*.log", emit: logs
    
    script:
    """
    # Convert ShortStack Results.txt to BED format
    awk 'BEGIN {OFS="\\t"} NR>1 {print \$3, \$4, \$5, \$2, \$7, \$10}' ${shortstack_dir}/Results.txt > ${sample_id}_clusters.bed
    
    # Sort BED file
    sort -k1,1 -k2,2n ${sample_id}_clusters.bed > ${sample_id}_clusters_sorted.bed
    
    # Annotate with GFF3 features
    bedtools intersect -a ${sample_id}_clusters_sorted.bed -b ${genome_gff3} -wa -wb > ${sample_id}_annotated_regions.txt
    
    # Log statistics
    echo "Total clusters: \$(wc -l < ${sample_id}_clusters_sorted.bed)" > annotation_stats.log
    echo "Annotated clusters: \$(wc -l < ${sample_id}_annotated_regions.txt)" >> annotation_stats.log
    """
}
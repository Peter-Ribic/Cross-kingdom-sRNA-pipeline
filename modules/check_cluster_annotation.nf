process CHECK_CLUSTER_ANNOTATION {
    tag "$sample_id"

    conda 'bioconda::bedtools=2.30.0 bioconda::samtools=1.15.1'

    publishDir "results/shortstack/shortstack_annotated/${sample_id}", mode: 'copy'

    input:
    // shortstack_dir still used for MajorRNA.fa (for non-protein-coding fasta output)
    tuple val(sample_id), path(shortstack_dir), path(majorrna_exact_all)
    path genome_gff3

    output:
    tuple val(sample_id), path("${sample_id}_majorRNA_annotated_regions.txt"), emit: annotated_regions
    path "${sample_id}_majorRNA_sites.bed", emit: majorrna_bed
    path "majorRNA_feature_type_counts.tsv", emit: feature_counts
    tuple val(sample_id), path("${sample_id}_MajorRNA_non_protein_coding.fa"), emit: non_protein_coding_fasta
    path "*.log", emit: logs
    path "${task.process}_${sample_id}.tsv", emit: log_info

    script:
    """
    set -euo pipefail

    fasta="${shortstack_dir}/${sample_id}_MajorRNA.fa"
    test -s "\$fasta"
    test -s "${genome_gff3}"
    test -s "${majorrna_exact_all}"


    awk -F'\\t' 'BEGIN{OFS="\\t"}
      NR==1{next}
      {
        chrom=\$3
        start1=\$8
        end1=\$9
        strand=\$10
        mapq=\$11
        cluster=\$2
        maj=\$7
        xw=\$13

        if(start1=="" || end1=="") next
        start0=start1-1
        end0=end1

        name=cluster "|" maj "|XW:" xw
        print chrom, start0, end0, name, mapq, strand
      }' "${majorrna_exact_all}" \
      | sort -k1,1 -k2,2n > ${sample_id}_majorRNA_sites.bed

    total_sites=\$(wc -l < ${sample_id}_majorRNA_sites.bed || true)
    echo "MajorRNA sites from exact table: \$total_sites" > annotation_stats.log

    bedtools intersect \
      -a ${sample_id}_majorRNA_sites.bed \
      -b "${genome_gff3}" \
      -wa -wb > ${sample_id}_majorRNA_annotated_regions.txt

    annotated_rows=\$(wc -l < ${sample_id}_majorRNA_annotated_regions.txt || true)
    echo "Annotated MajorRNA hits (rows in intersect): \$annotated_rows" >> annotation_stats.log

    uniq_loci=\$(awk 'BEGIN{OFS=":"}{print \$1,\$2"-"\$3,\$6}' ${sample_id}_majorRNA_sites.bed | sort -u | wc -l || true)
    echo "Unique MajorRNA loci: \$uniq_loci" >> annotation_stats.log

    total_clusters=\$(cut -f4 ${sample_id}_majorRNA_sites.bed | awk -F'|' '{print \$1}' | sort -u | wc -l || true)
    annotated_clusters=\$(cut -f4 ${sample_id}_majorRNA_annotated_regions.txt | awk -F'|' '{print \$1}' | sort -u | wc -l || true)
    unannotated_clusters=\$(( total_clusters - annotated_clusters ))
    echo "Clusters with NO annotation overlaps: \$unannotated_clusters (of \$total_clusters)" >> annotation_stats.log

    awk 'BEGIN{OFS="\\t"} {cnt[\$9]++}
         END{
           print "feature_type","n_overlaps"
           for(t in cnt) print t,cnt[t]
         }' ${sample_id}_majorRNA_annotated_regions.txt | sort -k2,2nr > majorRNA_feature_type_counts.tsv

    echo "" >> annotation_stats.log
    echo "Top feature types (overlaps):" >> annotation_stats.log
    head -n 20 majorRNA_feature_type_counts.tsv >> annotation_stats.log

    awk 'BEGIN{FS="\\t"} \$0 !~ /^#/ && (\$3=="CDS" || \$3=="exon")' "${genome_gff3}" > protein_coding_features.gff3 || true

    if [ -s protein_coding_features.gff3 ]; then
      bedtools intersect \
        -a ${sample_id}_majorRNA_sites.bed \
        -b protein_coding_features.gff3 \
        -wa > majorrna_overlaps_protein.bed || true
    else
      : > majorrna_overlaps_protein.bed
    fi

    cut -f4 majorrna_overlaps_protein.bed | awk -F'|' '{print \$1}' | sort -u > coding_clusters.txt || true
    coding_n=\$(wc -l < coding_clusters.txt || true)
    echo "Clusters with any site overlapping CDS/exon: \$coding_n" >> annotation_stats.log

    awk '
      NR==FNR {bad[\$1]=1; next}
      /^>/{
        hdr=substr(\$0,2)
        keep = !(hdr in bad)
      }
      keep {print}
    ' coding_clusters.txt "\$fasta" > ${sample_id}_MajorRNA_non_protein_coding.fa

    num_records=\$(awk 'END{print NR/2}' ${sample_id}_MajorRNA_non_protein_coding.fa)
    echo "Non-protein-coding MajorRNA records emitted: \$num_records" >> annotation_stats.log

    echo -e "${task.process}\\t${sample_id}\\t\$num_records" > ${task.process}_${sample_id}.tsv
    """
}

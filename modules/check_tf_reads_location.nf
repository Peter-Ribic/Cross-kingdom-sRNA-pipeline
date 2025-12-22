process CHECK_TF_READS_LOCATION {
  tag "$sample_id"
  container "quay.io/biocontainers/bedtools:2.31.1--h13024bc_3"
  publishDir "results/tf_read_gff3/${sample_id}", mode: 'copy'

  input:
  tuple val(sample_id), path(tf_txt), path(bam)
  path genome_gff3

  output:
  tuple val(sample_id), path("${sample_id}_tf_annotated.txt"), emit: tf_annotated
  path("qids.txt"), emit: qids
  path("tf_reads.bed"), emit: tf_reads_bed
  path("overlaps.tsv"), emit: overlaps

  script:
 """
  # 1) read IDs from TargetFinder (query=<ID>, ...)
  grep -E '^query=' "${tf_txt}" | sed -E 's/^query=([^,]+),.*/\\1/' | sort -u > qids.txt

  # 3) BAM -> BED, then keep only reads in qids.txt
  bedtools bamtobed -i ${bam} > all_reads.bed

  awk 'NR==FNR { keep[\$1]=1; next }
       keep[\$4]==1 { print }' qids.txt all_reads.bed > tf_reads.bed

  # 4) overlap reads with GFF3 (must be a GFF3/GTF/BED, NOT FASTA)
  bedtools intersect -a tf_reads.bed -b "${genome_gff3}" -wa -wb > overlaps.tsv

  # 5) collapse overlaps per read to one string (type;ID or Parent)
  awk 'BEGIN{OFS="\\t"}
       {
         rid=\$4; type=\$9; attr=\$15;
         id=""
         if (match(attr,/ID=[^;]+/)) id=substr(attr,RSTART,RLENGTH)
         else if (match(attr,/Parent=[^;]+/)) id=substr(attr,RSTART,RLENGTH)
         else id=substr(attr,1,60)
         ann=type ";" id
         if (m[rid]=="" ) m[rid]=ann
         else if (index(" | " m[rid] " | "," | " ann " | ")==0) m[rid]=m[rid] " | " ann
       }
       END{ for (r in m) print r, m[r] }' overlaps.tsv | sort -k1,1 > read2ann.tsv

  # 6) append annotation ONLY to query= lines (preserve blocks)
  awk '
    BEGIN{ while((getline < "read2ann.tsv")>0) ann[\$1]=\$2 }
    /^query=/{
      q=\$0; sub(/^query=/,"",q); sub(/,.*/,"",q)
      print \$0 "  gff3=" (q in ann ? ann[q] : "NA")
      next
    }
    {print}
  ' "${tf_txt}" > ${sample_id}_tf_annotated.txt
  """
}

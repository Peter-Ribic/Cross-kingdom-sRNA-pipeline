process CHECK_ANNOTATION {
    tag "$sample_id"

    conda 'bioconda::samtools=1.15.1 bioconda::bedtools=2.30.0'

    publishDir "results/shortstack/shortstack_annotated/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(shortstack_dir)
    path genome_gff3

    output:
    tuple val(sample_id), path("${sample_id}_majorRNA_annotated_regions.txt"), emit: annotated_regions
    path "${sample_id}_majorRNA_sites.bed", emit: majorrna_bed
    path "majorRNA_feature_type_counts.tsv", emit: feature_counts
    tuple val(sample_id), path("${sample_id}_MajorRNA_non_protein_coding.fa"), emit: non_protein_coding_fasta
    path "*.log", emit: logs

    script:
    """
    set -euo pipefail

    bam="${shortstack_dir}/${sample_id}_pathogen_specific_condensed.bam"
    fasta="${shortstack_dir}/${sample_id}_MajorRNA.fa"

    test -s "\$bam"
    test -s "\$fasta"
    test -s "${genome_gff3}"

    # ------------------------------------------------------------------------------
    # 1) Parse FASTA -> majorrna_records.tsv
    #    Columns: normseq \\t original_header(no '>') \\t original_sequence(concatenated)
    #    Then derive majorrna_seqs.txt = unique normseq list
    # ------------------------------------------------------------------------------
    awk '
      function norm(s,   t) {
        t = toupper(s)
        gsub(/U/,"T",t)
        gsub(/[^ACGTN]/,"",t)
        return t
      }

      BEGIN { hdr=""; seq=""; }

      /^>/ {
        if (hdr != "") {
          n = norm(seq)
          if (n != "") print n "\\t" hdr "\\t" seq
        }
        hdr = substr(\$0,2)
        seq = ""
        next
      }

      {
        gsub(/\\r/,"")
        gsub(/[ \\t]/,"")
        seq = seq \$0
      }

      END {
        if (hdr != "") {
          n = norm(seq)
          if (n != "") print n "\\t" hdr "\\t" seq
        }
      }
    ' "\$fasta" > majorrna_records.tsv

    test -s majorrna_records.tsv

    cut -f1 majorrna_records.tsv | sort -u > majorrna_seqs.txt
    echo "Unique MajorRNA sequences in FASTA (normalized): \$(wc -l < majorrna_seqs.txt)" > annotation_stats.log

    # ------------------------------------------------------------------------------
    # 2) Exact match BAM SEQ to normseq list; write BED6
    #    BED name column = normalized sequence (so we can map back to headers)
    # ------------------------------------------------------------------------------
    samtools view -F 4 "\$bam" \\
      | awk 'BEGIN{OFS="\\t"}
             function norm(s,   t) {
               t=toupper(s)
               gsub(/U/,"T",t)
               gsub(/[^ACGTN]/,"",t)
               return t
             }
             NR==FNR { want[\$1]=1; next }
             {
               nseq = norm(\$10)
               if(!(nseq in want)) next

               rname=\$3
               pos0=\$4-1
               mapq=\$5
               flag=\$2+0
               strand = (int(flag/16)%2 ? "-" : "+")

               start=pos0
               end=pos0 + length(nseq)

               print rname, start, end, nseq, mapq, strand
             }' majorrna_seqs.txt - \\
      | sort -k1,1 -k2,2n > ${sample_id}_majorRNA_sites.bed

    total_majorrna_hits=\$(wc -l < ${sample_id}_majorRNA_sites.bed || true)
    echo "MajorRNA hits (exact SEQ match; normalized U->T): \$total_majorrna_hits" >> annotation_stats.log

    # ------------------------------------------------------------------------------
    # 3) Intersect with full annotation (for your existing output)
    # ------------------------------------------------------------------------------
    bedtools intersect \\
      -a ${sample_id}_majorRNA_sites.bed \\
      -b "${genome_gff3}" \\
      -wa -wb > ${sample_id}_majorRNA_annotated_regions.txt

    annotated_majorrna_hits=\$(wc -l < ${sample_id}_majorRNA_annotated_regions.txt || true)
    echo "Annotated MajorRNA hits (rows in intersect): \$annotated_majorrna_hits" >> annotation_stats.log

    uniq_loci=\$(awk 'BEGIN{OFS=":"}{print \$1,\$2"-"\$3,\$6}' ${sample_id}_majorRNA_sites.bed | sort -u | wc -l || true)
    echo "Unique MajorRNA loci: \$uniq_loci" >> annotation_stats.log

    awk 'BEGIN{OFS="\\t"} {cnt[\$9]++}
         END{
           print "feature_type","n_overlaps"
           for(t in cnt) print t,cnt[t]
         }' ${sample_id}_majorRNA_annotated_regions.txt | sort -k2,2nr > majorRNA_feature_type_counts.tsv

    echo "" >> annotation_stats.log
    echo "Top feature types (overlaps):" >> annotation_stats.log
    head -n 20 majorRNA_feature_type_counts.tsv >> annotation_stats.log

    # ------------------------------------------------------------------------------
    # 4) Build FASTA of MajorRNA sequences that do NOT overlap CDS or exon
    #    - Identify normseqs that overlap CDS/exon
    #    - noncoding_normseqs.txt = majorrna_seqs - coding_normseqs
    #    - Output FASTA using ORIGINAL headers/sequences from majorrna_records.tsv
    # ------------------------------------------------------------------------------
    awk 'BEGIN{FS="\\t"} \$0 !~ /^#/ && (\$3=="CDS" || \$3=="exon")' "${genome_gff3}" > protein_coding_features.gff3 || true

    if [ -s protein_coding_features.gff3 ]; then
      bedtools intersect \\
        -a ${sample_id}_majorRNA_sites.bed \\
        -b protein_coding_features.gff3 \\
        -wa > majorrna_overlaps_protein.bed || true
    else
      : > majorrna_overlaps_protein.bed
    fi

    cut -f4 majorrna_overlaps_protein.bed | sort -u > coding_normseqs.txt || true

    # comm requires both sorted
    sort -u coding_normseqs.txt -o coding_normseqs.txt
    sort -u majorrna_seqs.txt -o majorrna_seqs.txt

    comm -23 majorrna_seqs.txt coding_normseqs.txt > noncoding_normseqs.txt || true

    noncoding_n=\$(wc -l < noncoding_normseqs.txt || true)
    coding_n=\$(wc -l < coding_normseqs.txt || true)
    echo "Normalized MajorRNA seqs overlapping CDS/exon: \$coding_n" >> annotation_stats.log
    echo "Normalized MajorRNA seqs NOT overlapping CDS/exon: \$noncoding_n" >> annotation_stats.log

    # Emit FASTA with original headers + original sequences for noncoding normseqs
    awk 'BEGIN{FS="\\t"}
         NR==FNR { keep[\$1]=1; next }
         {
           n=\$1; hdr=\$2; seq=\$3
           if (n in keep) {
             print ">" hdr
             print seq
           }
         }' noncoding_normseqs.txt majorrna_records.tsv > ${sample_id}_MajorRNA_non_protein_coding.fa
    """
}

process SHORTSTACK {

    tag "$sample_id"

    conda "bioconda::shortstack==4.1.2 bioconda::samtools"

    publishDir "results/shortstack/shortstack_results/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(trimmed_reads)
    path genome_fasta

    output:
    tuple val(sample_id), path("ShortStack_out"), emit: shortstack_out
    tuple val(sample_id), path("ShortStack_out/*condensed.bam"), emit: bam
    tuple val(sample_id), path("ShortStack_out/*condensed.fa"), emit: fasta
    tuple val(sample_id), path("ShortStack_out/Results.txt"), emit: results
    tuple val(sample_id), path("ShortStack_out/${sample_id}_MajorRNA.fa"), emit: majorrna_fasta
    tuple val(sample_id), path("ShortStack_out"), path("ShortStack_out/${sample_id}_MajorRNA_exact_locations.all.tsv"), emit: majorrna_exact_all
    tuple val(sample_id), path("ShortStack_out/${sample_id}_MajorRNA_exact_locations.best.tsv"), emit: majorrna_exact_best

    path "${task.process}_${sample_id}.tsv", emit: log_info

    script:
    """
    set -euo pipefail

    ShortStack \
        --readfile ${trimmed_reads} \
        --genomefile ${genome_fasta} \
        --outdir ShortStack_out \
        --threads 10 \
        --mincov 400 \
        --dicermin 20

    awk -F \$'\\t' 'NR>1 && \$11 != "" {
        printf(">%s\\n%s\\n", \$2, \$11)
    }' ShortStack_out/Results.txt > ShortStack_out/${sample_id}_MajorRNA.fa

    num_clusters=\$(awk 'END{print NR}' ShortStack_out/Results.gff3)
    echo -e "${task.process}\\t${sample_id}\\t\$num_clusters" > ${task.process}_${sample_id}.tsv



    bam_file=\$(ls ShortStack_out/*condensed.bam | head -n 1)
    test -s "\$bam_file"
    samtools index "\$bam_file" || true

    out_all="ShortStack_out/${sample_id}_MajorRNA_exact_locations.all.tsv"
    out_best="ShortStack_out/${sample_id}_MajorRNA_exact_locations.best.tsv"

    printf "sample_id\\tcluster\\tchrom\\tcluster_start\\tcluster_end\\tcluster_strand\\tmajorRNA_norm\\taln_start\\taln_end\\taln_strand\\tmapq\\tlen\\tXW\\tCIGAR\\n" > "\$out_all"

    awk -F \$'\\t' 'NR>1 && \$11!="" {print \$2"\\t"\$3"\\t"\$4"\\t"\$5"\\t"\$10"\\t"\$11}' ShortStack_out/Results.txt \\
    | while IFS=\$'\\t' read -r cluster chrom start end strand maj; do

        maj_norm=\$(echo "\$maj" | awk '
            function norm(s, t){t=toupper(s); gsub(/U/,"T",t); gsub(/[^ACGTN]/,"",t); return t}
            {print norm(\$0)}')

        [ -n "\$maj_norm" ] || continue

        samtools view -F 4 "\$bam_file" "\${chrom}:\${start}-\${end}" \\
        | awk -v sample="${sample_id}" -v cl="\$cluster" -v chr="\$chrom" -v cs="\$start" -v ce="\$end" -v cstrand="\$strand" -v maj="\$maj_norm" '
            BEGIN{OFS="\\t"}
            function norm(s, t){t=toupper(s); gsub(/U/,"T",t); gsub(/[^ACGTN]/,"",t); return t}
            function rc(s, i,c,r){
              r=""
              for(i=length(s); i>0; i--){
                c=substr(s,i,1)
                if(c=="A") c="T"; else if(c=="T") c="A"; else if(c=="C") c="G"; else if(c=="G") c="C"; else c="N"
                r=r c
              }
              return r
            }
            function getXW(   i){
              for(i=12;i<=NF;i++){
                if(\$i ~ /^XW:i:/){ sub(/^XW:i:/,"",\$i); return \$i+0 }
              }
              return 1
            }
            {
              # SAM fields: RNAME \$3, POS \$4 (1-based), MAPQ \$5, CIGAR \$6, SEQ \$10, FLAG \$2
              seq = norm(\$10)
              if(seq == "") next
              if(seq != maj && rc(seq) != maj) next

              flag=\$2+0
              astrand = (int(flag/16)%2 ? "-" : "+")
              len = length(seq)

              aln_start = \$4
              aln_end   = \$4 + len - 1

              xw = getXW()

              print sample, cl, chr, cs, ce, cstrand, maj, aln_start, aln_end, astrand, \$5, len, xw, \$6
            }' >> "\$out_all"
    done

    awk -F'\\t' '
      NR==1 {print; next}
      {
        key=\$2
        xw=\$13+0
        mapq=\$11+0
        start=\$8+0
        score = xw*1000000000 + mapq*1000000 - start
        if(!(key in best) || score > best[key]) { best[key]=score; line[key]=\$0 }
      }
      END{
        for(k in line) print line[k]
      }' "\$out_all" \\
    | (head -n 1; tail -n +2 | sort -k2,2) > "\$out_best"
    """
}

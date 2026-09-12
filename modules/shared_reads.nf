process SHARED_READS {
    tag "${sample1} vs ${sample2}"
    container "quay.io/biocontainers/bioawk:1.0--h577a1d6_13"
    publishDir "results/library_similarity/shared_reads", mode: 'copy'

    input:
    tuple val(sample1), path(sample1_fq), val(sample2), path(sample2_fq)

    output:
    path "shared_reads_${sample1}_${sample2}.txt"
    path "shared_reads_${sample1}_${sample2}.table.tsv"

    script:
    """
    set -euo pipefail
    export LC_ALL=C

    total1=\$(zcat "${sample1_fq}" | bioawk -c fastx 'END{print NR}')
    total2=\$(zcat "${sample2_fq}" | bioawk -c fastx 'END{print NR}')

    zcat "${sample1_fq}" | bioawk -c fastx '{print \$seq}' | sort -S 1G | uniq -c > s1.counts.txt
    zcat "${sample2_fq}" | bioawk -c fastx '{print \$seq}' | sort -S 1G | uniq -c > s2.counts.txt

    read shared_reads jaccard_w < <(
      awk '
        NR==FNR {
          c1[\$2]=\$1
          maxsum += \$1  
          next
        }
        {
          seq=\$2
          c=\$1
          if (seq in c1) {
            minsum += (c < c1[seq] ? c : c1[seq])

            if (c > c1[seq]) maxsum += (c - c1[seq])
          } else {
            maxsum += c
          }
        }
        END {
          shared = minsum + 0
          if (maxsum > 0) jac = minsum / maxsum
          else jac = 0
          printf "%d %.6f\\n", shared, jac
        }
      ' s1.counts.txt s2.counts.txt
    )

    pct1_in_2=\$(awk -v s="\$shared_reads" -v t="\$total1" 'BEGIN{if(t>0) printf "%.2f", (s/t)*100; else printf "0.00"}')
    pct2_in_1=\$(awk -v s="\$shared_reads" -v t="\$total2" 'BEGIN{if(t>0) printf "%.2f", (s/t)*100; else printf "0.00"}')

    {
      echo "sample1=${sample1}"
      echo "sample2=${sample2}"
      echo "total_reads_sample1=\$total1"
      echo "total_reads_sample2=\$total2"
      echo "shared_reads=\$shared_reads"
      echo "pct_reads_sample1_in_sample2=\${pct1_in_2}%"
      echo "pct_reads_sample2_in_sample1=\${pct2_in_1}%"
      echo "jaccard_weighted=\$jaccard_w"
    } > "shared_reads_${sample1}_${sample2}.txt"

    {
      printf "sample1\\tsample2\\ttotal_reads_sample1\\ttotal_reads_sample2\\tshared_reads\\tpct_reads_sample1_in_sample2\\tpct_reads_sample2_in_sample1\\tjaccard_weighted\\n"
      printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n" \\
        "${sample1}" "${sample2}" "\$total1" "\$total2" "\$shared_reads" "\$pct1_in_2" "\$pct2_in_1" "\$jaccard_w"
    } > "shared_reads_${sample1}_${sample2}.table.tsv"
    """
}

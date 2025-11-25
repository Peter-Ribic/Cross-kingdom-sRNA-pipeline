process CHECK_ALIGNMENT_DISTRIBUTION {
    tag "$sample_id"
    container "quay.io/biocontainers/samtools:1.22--h96c455f_0"
    publishDir "results/virus_alignment_dist/${sample_id}", mode: 'symlink'

    input:
    tuple val(sample_id), path(alignment_table), path(alignment_raw)
    
    output:
    tuple val(sample_id), path("${sample_id}_viruses_0mm.stats.txt"), path("${sample_id}_idxstats_raw.txt"), emit: results
    path("${sample_id}_viruses_0mm.percent_row.txt"), emit: percent_row 

    script:
    """
    # convert SAM to BAM
    samtools view -b ${sample_id}_viruses_0mm.sam > ${sample_id}_viruses_0mm.bam

    # sort BAM
    samtools sort ${sample_id}_viruses_0mm.bam -o ${sample_id}_viruses_0mm.sorted.bam

    # generate idxstats
    samtools idxstats ${sample_id}_viruses_0mm.sorted.bam > ${sample_id}_idxstats_raw.txt

    # get total mapped reads
    total=\$(awk '{sum+=\$3} END{print sum}' ${sample_id}_idxstats_raw.txt)

    # calculate percent mapped for each reference
    awk -v total="\$total" 'BEGIN{OFS="\\t"; print "ref_name","mapped_reads","percent_mapped"} 
         {if(total>0){percent = (\$3/total)*100} else {percent=0}; 
          print \$1, \$3, percent}' \
         ${sample_id}_idxstats_raw.txt > ${sample_id}_viruses_0mm.stats.txt

    # create single-row file with header: first column sample_id, rest are percent mapped
    header=\$(awk 'BEGIN{OFS="\\t"} {printf "%s\\t", \$1} END{print ""}' ${sample_id}_idxstats_raw.txt)
    percents=\$(awk -v total="\$total" 'BEGIN{OFS="\\t"} {if(total>0){percent = (\$3/total)*100} else {percent=0}; printf "%.2f\\t", percent} END{print ""}' ${sample_id}_idxstats_raw.txt)
    
    # write header and percent row
    (echo -e "sample_id\\t\$header"; echo -e "${sample_id}\\t\$percents") > ${sample_id}_viruses_0mm.percent_row.txt
    """
}

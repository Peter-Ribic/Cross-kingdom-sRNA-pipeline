process COMPARE_SRNA_LIBRARIES_SIMILAR {
    tag "$sample1 vs $sample2"
    container "quay.io/biocontainers/vsearch:2.30.1--hd6d6fdc_0"
    publishDir "results/library_similarity/srna_overlap", mode: 'copy'

    input:
    tuple val(sample1), path(sample1_fq), val(sample2), path(sample2_fq)

    output:
    path "similarity_${sample1}_${sample2}.txt"

    script:
    """
    # Collapse FASTQ reads into unique sequences with counts
    zcat $sample1_fq | awk 'NR%4==2{print}' | sort | uniq -c | awk '{print ">s1_"NR"_c"\$1"\\n"\$2}' > ${sample1}.fa
    zcat $sample2_fq | awk 'NR%4==2{print}' | sort | uniq -c | awk '{print ">s2_"NR"_c"\$1"\\n"\$2}' > ${sample2}.fa
    cat ${sample1}.fa ${sample2}.fa > combined.fa

    # Cluster sequences by 95% identity (~1 nt difference for 21 nt reads)
    vsearch --cluster_fast combined.fa --id 0.95 --uc clusters.uc --minseqlength 15 --threads ${task.cpus}

    # Count clusters with sequences from both samples
    awk '
      BEGIN{s1=0;s2=0;both=0}
      /^H|^S/ {
        if(\$9 ~ /^s1_/) s1c[\$2]=1
        if(\$9 ~ /^s2_/) s2c[\$2]=1
      }
      END{
        for(c in s1c) if(c in s2c) both++
        total=length(s1c)+length(s2c)-both
        print "clusters_sample1="length(s1c)
        print "clusters_sample2="length(s2c)
        print "shared_clusters="both
        print "jaccard_clusters="both/total
      }' clusters.uc > similarity_${sample1}_${sample2}.txt
    """
}

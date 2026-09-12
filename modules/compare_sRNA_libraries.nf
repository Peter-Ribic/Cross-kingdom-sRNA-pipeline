process COMPARE_SRNA_LIBRARIES_SIMILAR {
    tag "$sample1 vs $sample2"
    container "quay.io/biocontainers/vsearch:2.30.1--hd6d6fdc_0"
    publishDir "results/library_similarity/srna_overlap", mode: 'copy'

    input:
    tuple val(sample1), path(sample1_fq), val(sample2), path(sample2_fq)

    output:
    path "similarity_${sample1}_${sample2}.txt"
    path "similarity_${sample1}_${sample2}.table.tsv"

    script:
    """
    zcat $sample1_fq | awk 'NR%4==2{print}' | sort | uniq -c | awk '{print ">s1_"NR"_c"\$1"\\n"\$2}' > ${sample1}.fa
    zcat $sample2_fq | awk 'NR%4==2{print}' | sort | uniq -c | awk '{print ">s2_"NR"_c"\$1"\\n"\$2}' > ${sample2}.fa
    cat ${sample1}.fa ${sample2}.fa > combined.fa

    vsearch --cluster_fast combined.fa --id 0.95 --uc clusters.uc --minseqlength 15 --threads ${task.cpus}

    awk '
      BEGIN{s1=0;s2=0;both=0}
      /^H|^S/ {
        if(\$9 ~ /^s1_/) s1c[\$2]=1
        if(\$9 ~ /^s2_/) s2c[\$2]=1
      }
      END{
        for(c in s1c) if(c in s2c) both++
        n1=length(s1c)
        n2=length(s2c)
        total=n1+n2-both
        jac=(total>0)?both/total:0

        print "clusters_sample1="n1
        print "clusters_sample2="n2
        print "shared_clusters="both
        print "jaccard_clusters="jac
      }' clusters.uc > similarity_${sample1}_${sample2}.txt

    awk -F'=' '
      BEGIN{
        n1=""; n2=""; both=""; jac="";
      }
      \$1=="clusters_sample1"{n1=\$2}
      \$1=="clusters_sample2"{n2=\$2}
      \$1=="shared_clusters"{both=\$2}
      \$1=="jaccard_clusters"{jac=\$2}
      END{
        OFS="\\t"
        print "sample1","sample2","clusters_sample1","clusters_sample2","shared_clusters","jaccard_clusters"
        print "'${sample1}'","'${sample2}'",n1,n2,both,jac
      }' similarity_${sample1}_${sample2}.txt > similarity_${sample1}_${sample2}.table.tsv
    """
}

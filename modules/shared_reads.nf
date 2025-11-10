process SHARED_READS {
    tag "$sample1 vs $sample2"
    container "quay.io/biocontainers/bioawk:1.0--h577a1d6_13"
    publishDir "results/shared_reads", mode: 'copy'

    input:
    tuple val(sample1), path(sample1_fq), val(sample2), path(sample2_fq)

    output:
    path "shared_reads_${sample1}_${sample2}.txt"

    script:
    """
       # Total and unique sequences in sample1
    total1=\$(zcat $sample1_fq | bioawk -c fastx '{print \$seq}' | wc -l)
    unique1=\$(zcat $sample1_fq | bioawk -c fastx '{print \$seq}' | sort | uniq | wc -l)

    # Total and unique sequences in sample2
    total2=\$(zcat $sample2_fq | bioawk -c fastx '{print \$seq}' | wc -l)
    unique2=\$(zcat $sample2_fq | bioawk -c fastx '{print \$seq}' | sort | uniq | wc -l)

    # Count shared unique sequences
    shared=\$(comm -12 <(zcat $sample1_fq | bioawk -c fastx '{print \$seq}' | sort | uniq) \
                        <(zcat $sample2_fq | bioawk -c fastx '{print \$seq}' | sort | uniq) | wc -l)

    # Write statistics to file
    echo "Sample1: $sample1" > shared_reads_${sample1}_${sample2}.txt
    echo "Sample2: $sample2" >> shared_reads_${sample1}_${sample2}.txt
    echo "Total reads Sample1: \$total1" >> shared_reads_${sample1}_${sample2}.txt
    echo "Unique reads Sample1: \$unique1" >> shared_reads_${sample1}_${sample2}.txt
    echo "Total reads Sample2: \$total2" >> shared_reads_${sample1}_${sample2}.txt
    echo "Unique reads Sample2: \$unique2" >> shared_reads_${sample1}_${sample2}.txt
    echo "Shared unique reads: \$shared" >> shared_reads_${sample1}_${sample2}.txt
    echo "Percentage of Sample1 reads in Sample2: \$((100 * shared / unique1))%" >> shared_reads_${sample1}_${sample2}.txt
    """
}

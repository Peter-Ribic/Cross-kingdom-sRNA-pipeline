process FETCH_SRA {
    tag "$sample_id"
    conda 'bioconda::sra-tools=3.0.3'

    input:
    tuple val(sample_id), val(accessions)

    output:
    tuple val(sample_id), path("*.fastq")

    script:
    """
    for acc in ${accessions.join(' ')}; do
        prefetch \$acc
        fasterq-dump \$acc --split-files -O .
    done
    """
}

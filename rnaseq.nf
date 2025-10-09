#!/usr/bin/env nextflow

// Module INCLUDE statements
include { FETCH_SRA } from './modules/fetch_sra.nf'
include { FASTQC } from './modules/fastqc.nf'
include { TRIM_GALORE } from './modules/trim_galore.nf'

include { HISAT2_BUILD } from './modules/hisat2_build.nf'
include { HISAT2_ALIGN } from './modules/hisat2_align.nf'
include { MULTIQC } from './modules/multiqc.nf'


// params.input_csv = "data/single-end.csv"
params.report_id = "all_single-end"
//params.genome_fasta = "data/hops_genome/hops_genome.fa"
params.genome_fasta = "data/verticilium_genome/T2/VnonalfalfaeT2.fasta"
params.reads = "data/sra_data/*/*.fastq"
params.sra_csv = "data/sra_accessions.csv"


workflow {

    read_ch = Channel
        .fromPath(params.sra_csv)
        .splitCsv(header: true)
        .map { row ->
            // Collect all non-empty SRA IDs from the row
            def sra_list = row.values()
                            .findAll { it =~ /^SRR/ }   // keep only SRA IDs
            tuple(row.sample_id, sra_list)
        }
    
    FETCH_SRA(read_ch)
    
    //read_ch = Channel.fromPath(params.reads)
    
    FASTQC(FETCH_SRA.out)
    TRIM_GALORE(FETCH_SRA.out)
    index_ch = HISAT2_BUILD(file(params.genome_fasta))
    HISAT2_ALIGN(TRIM_GALORE.out.trimmed_reads, index_ch)

    qc_ch = FASTQC.out.zip
    .mix(
        FASTQC.out.html,
        TRIM_GALORE.out.trimming_reports,
        TRIM_GALORE.out.fastqc_reports,
        HISAT2_ALIGN.out.log
    )
    .collect()

    MULTIQC(tuple(params.report_id, qc_ch))
}

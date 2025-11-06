#!/usr/bin/env nextflow

// Module INCLUDE statements
include { FETCH_SRA } from './modules/fetch_sra.nf'
include { FETCH_SRA as FETCH_SRA_PATHOGEN } from './modules/fetch_sra.nf'

include { FASTQC } from './modules/fastqc.nf'
include { TRIM_GALORE } from './modules/trim_galore.nf'

include { HISAT2_ALIGN } from './modules/hisat2_align.nf'
include { MULTIQC } from './modules/multiqc.nf'
include { MERGE_READS } from './modules/merge_reads.nf'
include { SHORTSTACK } from './modules/shortstack.nf'
include { BOWTIE_ALIGN_TO_PATHOGEN } from './modules/bowtie_align_to_pathogen.nf'
include { BOWTIE_ALIGN_TO_HOST } from './modules/bowtie_align_to_host.nf'
include { FILTER_PATHOGEN_READS } from './modules/filter_pathogen_reads.nf'
include { BOWTIE_BUILD_PATHOGEN } from './modules/bowtie_build_pathogen.nf'
include { BOWTIE_BUILD_HOST } from './modules/bowtie_build_host.nf'
include { LIST_PATHOGEN_READS } from './modules/list_pathogen_reads.nf'
include { PREDICT_HOP_TARGETS } from './modules/predict_hop_targets.nf'
include { HISAT2_BUILD } from './modules/hisat2_build.nf'
include { FILTER_SRNA_LENGTH } from './modules/filter_srna_length.nf'
include { TARGETFINDER } from './modules/targetfinder.nf'
include { BOWTIE_ALIGN_TO_VIRUSES } from './modules/bowtie_align_to_viruses.nf'
include { BOWTIE_BUILD_VIRUSES } from './modules/bowtie_build_viruses.nf'

// params.input_csv = "data/single-end.csv"
params.report_id = "all_single-end"
params.host_genome_fasta = "data/hops_genome/hops_genome.fa"
params.pathogen_genome_fasta = "data/verticilium_genome/T2/VnonalfalfaeT2.fasta"
params.reads = "data/sra_data/*/*.fastq"
params.sra_csv = "data/sra_accessions.csv"
params.sra_csv_pathogen = "data/sra_accessions_pathogen.csv"
params.host_transcriptome_fasta = "data/hops_transcriptome/combinedGeneModels.fullAssembly.transcripts.fasta"
params.viruses_genome_fasta = "data/viruses_genome/3viroidi_2virusa.fa"


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

    read_ch_pathogen = Channel
        .fromPath(params.sra_csv_pathogen)
        .splitCsv(header: true)
        .map { row ->
            // Collect all non-empty SRA IDs from the row
            def sra_list = row.values()
                            .findAll { it =~ /^SRR/ }   // keep only SRA IDs
            tuple(row.sample_id, sra_list)
        }
    
    FETCH_SRA(read_ch)
    FETCH_SRA_PATHOGEN(read_ch_pathogen)
    
    FASTQC(FETCH_SRA.out)
    TRIM_GALORE(FETCH_SRA.out)
    merged_ch = MERGE_READS(TRIM_GALORE.out.trimmed_reads)
    FILTER_SRNA_LENGTH(merged_ch)

    // CHECK FOR VIRUS INFECTIONS
    viruses_index_ch = BOWTIE_BUILD_VIRUSES(file(params.viruses_genome_fasta))
    BOWTIE_ALIGN_TO_VIRUSES(merged_ch, viruses_index_ch)

    pathogen_index_ch = BOWTIE_BUILD_PATHOGEN(file(params.pathogen_genome_fasta))
    host_index_ch = BOWTIE_BUILD_HOST(file(params.host_genome_fasta))
    

    BOWTIE_ALIGN_TO_PATHOGEN(FILTER_SRNA_LENGTH.out.filtered_reads, pathogen_index_ch)
    BOWTIE_ALIGN_TO_HOST(BOWTIE_ALIGN_TO_PATHOGEN.out.results, host_index_ch)
    LIST_PATHOGEN_READS(BOWTIE_ALIGN_TO_HOST.out.list_input)
    
    filtered_reads = FILTER_PATHOGEN_READS(LIST_PATHOGEN_READS.out.filter_input)

    SHORTSTACK(filtered_reads, file(params.pathogen_genome_fasta))
    //PREDICT_HOP_TARGETS(filtered_reads, file(params.host_transcriptome_fasta))
    //TARGETFINDER(filtered_reads, file(params.host_transcriptome_fasta))
    qc_ch = FASTQC.out.zip
        .mix(
            FASTQC.out.html,
            TRIM_GALORE.out.trimming_reports,
            TRIM_GALORE.out.fastqc_reports,
            BOWTIE_ALIGN_TO_VIRUSES.out.log,
            BOWTIE_ALIGN_TO_PATHOGEN.out.log,
            BOWTIE_ALIGN_TO_HOST.out.log
        )
        .collect()

    MULTIQC(params.report_id, qc_ch)
}

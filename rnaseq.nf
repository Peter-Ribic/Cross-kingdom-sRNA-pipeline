#!/usr/bin/env nextflow

// Module INCLUDE statements
include { FETCH_SRA } from './modules/fetch_sra.nf'
include { FETCH_SRA as FETCH_SRA_PATHOGEN } from './modules/fetch_sra.nf'

include { FASTQC } from './modules/fastqc.nf'
include { FASTQC as FASTQC_PATHOGEN } from './modules/fastqc.nf'
include { TRIM_GALORE } from './modules/trim_galore.nf'
include { TRIM_GALORE as TRIM_GALORE_PATHOGEN } from './modules/trim_galore.nf'
include { MERGE_READS } from './modules/merge_reads.nf'
include { MERGE_READS as MERGE_READS_PATHOGEN} from './modules/merge_reads.nf'
include { FILTER_SRNA_LENGTH } from './modules/filter_srna_length.nf'
include { KEEP_ONLY_PATHOGEN_READS } from './modules/keep_only_pathogen_reads.nf'

include { HISAT2_ALIGN } from './modules/hisat2_align.nf'
include { MULTIQC } from './modules/multiqc.nf'
include { SHORTSTACK } from './modules/shortstack.nf'
include { BOWTIE_ALIGN_TO_PATHOGEN } from './modules/bowtie_align_to_pathogen.nf'
include { BOWTIE_ALIGN_TO_HOST } from './modules/bowtie_align_to_host.nf'
include { BOWTIE_BUILD_PATHOGEN } from './modules/bowtie_build_pathogen.nf'
include { FILTER_PATHOGEN_READS } from './modules/filter_pathogen_reads.nf'
include { BOWTIE_BUILD_HOST } from './modules/bowtie_build_host.nf'
include { LIST_PATHOGEN_READS } from './modules/list_pathogen_reads.nf'
include { HISAT2_BUILD } from './modules/hisat2_build.nf'
include { TARGETFINDER } from './modules/targetfinder.nf'
include { BOWTIE_ALIGN_TO_VIRUSES } from './modules/bowtie_align_to_viruses.nf'
include { BOWTIE_BUILD_VIRUSES } from './modules/bowtie_build_viruses.nf'
include { COMPARE_SRNA_LIBRARIES_SIMILAR } from './modules/compare_sRNA_libraries.nf'
include { SHARED_READS } from './modules/shared_reads.nf'
include { CHECK_ALIGNMENT_DISTRIBUTION } from './modules/check_alignment_distribution.nf'
include { CONCAT_ALIGNMENT_DISTRIBUTION } from './modules/concat_alignment_distribution.nf'
include { KEEP_TREATED_ONLY} from './modules/keep_treated_only.nf'
include { MATCH_WITH_CANDIDATE_SRNAs } from './modules/match_candidate_srnas.nf'
include { BOWTIE_BUILD_CANDIDATE_SRNAS } from './modules/bowtie_build_candidate_srnas.nf'
include { BOWTIE_ALIGN_TO_CANDIDATE_SRNAS } from './modules/bowtie_align_to_candidate_srnas.nf'
include {CHECK_ANNOTATION } from './modules/check_annotation.nf'
include {BOWTIE_BUILD } from './modules/bowtie_build.nf'
include {BOWTIE_ALIGN } from './modules/bowtie_align.nf'
include {BOWTIE_ALIGN as BOWTIE_ALIGN_PATHOGEN_PRELIMINARY } from './modules/bowtie_align.nf'
include {BOWTIE_ALIGN as BOWTIE_ALIGN_HOST_PRELIMINARY } from './modules/bowtie_align.nf'
include { FASTQ_TO_FASTA } from './modules/fastq_to_fasta.nf'
include { CONCAT_TARGETFINDER_RESULTS } from './modules/concat_targetfinder_results.nf'
include { ANNOTATE_TARGETS } from './modules/annotate_targets.nf'
include { PRELIMINARY_MULTIQC } from './modules/preliminary_multiqc.nf'
include { PRELIMINARY_MULTIQC as FASTQC_MULTIQC} from './modules/preliminary_multiqc.nf'
include { PRELIMINARY_MULTIQC as TRIM_GALORE_MULTIQC} from './modules/preliminary_multiqc.nf'
include { CHECK_TF_READS_LOCATION } from './modules/check_tf_reads_location.nf'
include { MIRNA_TARGET } from './modules/mirna_target.nf'
include { PLOT_SHORTSTACK_CLUSTERS } from './modules/plot_clusters.nf'
include { PLOT_SRNA_LENGTH_DISTRIBUTION } from './modules/plot_length_dist.nf'
include {EGGNOG_ANNOTATION} from './modules/eggnog_annotation.nf'
include {GOATOOLS_GOSLIM} from './modules/goatools_goslim.nf'
include {SUMMARIZE_LOGS} from './modules/summarize_logs.nf'
include {FASTP_TRIM} from './modules/fastp.nf'

// params.input_csv = "data/single-end.csv"
params.report_id = "all_single-end"
params.host_genome_fasta = "data/hops_genome/hops_genome.fa"
params.pathogen_genome_fasta = "data/verticilium_genome/T2/VnonalfalfaeT2.fasta"
params.pathogen_genome_gff = "data/verticilium_genome/T2/Verticillium_nonalfalfae_T2.gff3"
params.reads = "data/sra_data/*/*.fastq"
params.sra_csv = "data/sra_accessions.csv"
params.adapter_table = "data/sra_accessions_adapters.csv"
params.sra_csv_pathogen = "data/sra_accessions_pathogen.csv"
params.host_transcriptome_fasta = "data/hops_transcriptome/combinedGeneModels.fullAssembly.transcripts.fasta"
params.annotated_host_mrnas_fasta = "data/annotated_hop_mrnas/sequence.fasta"
params.mirna_target_repo = "data/MiRNATarget"

params.mock_hop_sample_control = "data/mock_data/hop_sample_control.fq"
params.mock_hop_sample_treated = "data/mock_data/hop_sample_treated.fq"
params.mock_verticillium_sample = "data/mock_data/verticillium_sample.fq"
params.mock_candidate_srnas_fasta = "data/mock_data/candidate_srnas.fa"

params.viruses_genome_fasta = "data/viruses_genome/3viroidi_2virusa.fa"
params.candidate_srnas_fasta = "data/candidate_srnas/candidate_srnas.fa"



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

    // read_ch_pathogen = Channel
    //     .fromPath(params.sra_csv_pathogen)
    //     .splitCsv(header: true)
    //     .map { row ->
    //         // Collect all non-empty SRA IDs from the row
    //         def sra_list = row.values()
    //                         .findAll { it =~ /^SRR/ }   // keep only SRA IDs
    //         tuple(row.sample_id, sra_list)
    //     }
 
 // MOCK DATA TESTING
//    trim_input_ch = Channel.from([
//     ['hop_sample_control', file(params.mock_hop_sample_control)],
//     ['hop_sample_treated', file(params.mock_hop_sample_treated)]
// ])
//     TRIM_GALORE(trim_input_ch)

//     trim_input_ch_pathogen = Channel.from([
//     ['verticillium_sample', file(params.mock_verticillium_sample)]
// ])
//     TRIM_GALORE_PATHOGEN(trim_input_ch_pathogen)
// END MOCK DATA TESTING
 
 
    // FETCH_SRA_PATHOGEN(read_ch_pathogen)
    // FASTQC_PATHOGEN(FETCH_SRA_PATHOGEN.out.reads)
    // TRIM_GALORE_PATHOGEN(FETCH_SRA_PATHOGEN.out)
    FETCH_SRA(read_ch)
    FASTQC(FETCH_SRA.out.reads)
    FASTQC_MULTIQC('fastqc_report', FASTQC.out.html.mix(FASTQC.out.zip).collect())

    //TRIM_GALORE(FETCH_SRA.out.reads, file(params.adapter_table))
    //TRIM_GALORE_MULTIQC('trim_galore_report', TRIM_GALORE.out.trimming_reports.mix(TRIM_GALORE.out.fastqc_reports).collect())
    FASTP_TRIM(FETCH_SRA.out.reads)
    merged_ch = MERGE_READS(FASTP_TRIM.out.trimmed_reads)
    
    // PRELIMINARY ALIGNMENT TO HOST AND PATHOGEN GENOMES TO ASSESS CONTAMINATION LEVELS
        host_index_ch = BOWTIE_BUILD_HOST(file(params.host_genome_fasta))
        pathogen_index_ch = BOWTIE_BUILD_PATHOGEN(file(params.pathogen_genome_fasta))
        BOWTIE_ALIGN_HOST_PRELIMINARY(merged_ch, host_index_ch, 'host_index')
        BOWTIE_ALIGN_PATHOGEN_PRELIMINARY(merged_ch, pathogen_index_ch, 'pathogen_index')

        preliminary_qc_ch = BOWTIE_ALIGN_HOST_PRELIMINARY.out.log
            .mix(BOWTIE_ALIGN_PATHOGEN_PRELIMINARY.out.log)
            .collect()
        PRELIMINARY_MULTIQC('preliminary_report', preliminary_qc_ch)
    //

    //KEEP_ONLY_PATHOGEN_READS(flat_ch)
    // CHECK FOR VIRUS INFECTIONS
        viruses_index_ch = BOWTIE_BUILD_VIRUSES(file(params.viruses_genome_fasta))
        BOWTIE_ALIGN_TO_VIRUSES(merged_ch, viruses_index_ch)
        CHECK_ALIGNMENT_DISTRIBUTION(BOWTIE_ALIGN_TO_VIRUSES.out.results)
        CONCAT_ALIGNMENT_DISTRIBUTION(CHECK_ALIGNMENT_DISTRIBUTION.out.percent_row.collect())
    //
    // MATCHING WITH CANDIDATE sRNAs
        candidate_srna_index_ch = BOWTIE_BUILD_CANDIDATE_SRNAS(file(params.candidate_srnas_fasta))
        BOWTIE_ALIGN_TO_CANDIDATE_SRNAS(FASTP_TRIM.out.trimmed_reads, candidate_srna_index_ch)
    //

    FILTER_SRNA_LENGTH(merged_ch)
    PLOT_SRNA_LENGTH_DISTRIBUTION(FILTER_SRNA_LENGTH.out.length_distribution)

    // COMPARE SRNA LIBRARIES FOR SIMILARITY
        pathogen_only_output_ch = FILTER_SRNA_LENGTH.out.filtered_reads
        all_combinations_ch = pathogen_only_output_ch
        .combine(pathogen_only_output_ch)
        .map { pair ->
            // 'pair' is a LinkedList of the two tuples: [ [sample1_id, file1], [sample2_id, file2] ]
            def sample1 = pair[0]
            def sample2 = pair[1]
            def sample3 = pair[2]
            def sample4 = pair[3]
            tuple(sample1, sample2, sample3, sample4)
        }
        .filter { s1_id, s1_file, s2_id, s2_file ->
            s1_id < s2_id  // avoid self-comparison & duplicates
        }
        COMPARE_SRNA_LIBRARIES_SIMILAR(all_combinations_ch)
        SHARED_READS(all_combinations_ch)
    //

    // FILTERING BASED ON TREATED MINUS CONTROL READS
        filtered_ch = FILTER_SRNA_LENGTH.out.filtered_reads
        keyed_ch = filtered_ch
            .map { sample_id, reads ->
                // remove _treated or _control suffix to get base name
                def base = sample_id.replaceAll(/_(treated|control)$/, '')
                tuple(base, sample_id, reads)
            }
            .groupTuple()

        paired_ch = keyed_ch
            .map { base, sample_ids, reads ->
                def treatedIndex = sample_ids.findIndexOf { it.endsWith('_treated') }
                def controlIndex = sample_ids.findIndexOf { it.endsWith('_control') }
                tuple(base, reads[treatedIndex], reads[controlIndex])
            }
        KEEP_TREATED_ONLY(paired_ch)
    //



    //FILTERING BASED ON READS THAT ALIGN PERFEECTLY TO PATHOGEN AND NOT PERFECTLY TO HOST
        BOWTIE_ALIGN_TO_PATHOGEN(KEEP_TREATED_ONLY.out.filtered_reads, pathogen_index_ch)
        BOWTIE_ALIGN_TO_HOST(BOWTIE_ALIGN_TO_PATHOGEN.out.results, host_index_ch)
        LIST_PATHOGEN_READS(BOWTIE_ALIGN_TO_HOST.out.list_input)
        FILTER_PATHOGEN_READS(LIST_PATHOGEN_READS.out.filter_input)
        filtered_reads = FILTER_PATHOGEN_READS.out.pathogen_specific_reads
    //


    SHORTSTACK(filtered_reads, file(params.pathogen_genome_fasta))
    CHECK_ANNOTATION(SHORTSTACK.out.shortstack_out, file(params.pathogen_genome_gff))
    // map filtered reads to host transcriptome
    BOWTIE_BUILD(file(params.annotated_host_mrnas_fasta), 'host_transcriptome_index')
    BOWTIE_ALIGN(filtered_reads, BOWTIE_BUILD.out.index_files, 'host_transcriptome_index')

    // USE SHORTSTACK MAJORRNA OUTPUT TO PREDICT HOP TARGETS
        MIRNA_TARGET(CHECK_ANNOTATION.out.non_protein_coding_fasta, file(params.annotated_host_mrnas_fasta), file(params.mirna_target_repo))
    //

    // ANNOTATE TARGETS
        ANNOTATE_TARGETS(MIRNA_TARGET.out.results, file(params.annotated_host_mrnas_fasta))
        EGGNOG_ANNOTATION(ANNOTATE_TARGETS.out.fastas)
        
        GOATOOLS_GOSLIM(EGGNOG_ANNOTATION.out.go_terms)
    //
    PLOT_SHORTSTACK_CLUSTERS(SHORTSTACK.out.results, file(params.pathogen_genome_fasta))

    // qc_ch = FASTQC.out.zip
    //     .mix(
    //         FASTQC.out.html,
    //       //  FASTQC_PATHOGEN.out.html,
    //         TRIM_GALORE.out.trimming_reports,
    //         TRIM_GALORE.out.fastqc_reports,
    //       //  TRIM_GALORE_PATHOGEN.out.trimming_reports,
    //       //  TRIM_GALORE_PATHOGEN.out.fastqc_reports,
    //         BOWTIE_ALIGN_TO_VIRUSES.out.log,
    //         // KEEP_ONLY_PATHOGEN_READS.out.log,
    //         KEEP_TREATED_ONLY.out.log,
    //         BOWTIE_ALIGN_TO_PATHOGEN.out.log,
    //         BOWTIE_ALIGN_TO_HOST.out.log
    //     )
    //     .collect()

    // MULTIQC(params.report_id, qc_ch)

    SUMMARIZE_LOGS(
        FILTER_SRNA_LENGTH.out.log_info
        .mix(
            FETCH_SRA.out.log_info,
            FASTP_TRIM.out.log_info,
            //TRIM_GALORE.out.log_info,
            BOWTIE_ALIGN_TO_VIRUSES.out.log_info,
            KEEP_TREATED_ONLY.out.log_info,
            BOWTIE_ALIGN_TO_PATHOGEN.out.log_info,
            BOWTIE_ALIGN_TO_HOST.out.log_info,
            FILTER_PATHOGEN_READS.out.log_info,
            SHORTSTACK.out.log_info,
            MIRNA_TARGET.out.log_info
        )
        .collect()
    )
}

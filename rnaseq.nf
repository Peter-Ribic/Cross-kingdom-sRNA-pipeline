#!/usr/bin/env nextflow

include { FETCH_SRA } from './modules/fetch_sra.nf'
include { FASTQC } from './modules/fastqc.nf'
include { FASTQC as FASTQC_AFTER_TRIMMING} from './modules/fastqc.nf'
include { MERGE_READS } from './modules/merge_reads.nf'
include { FILTER_SRNA_LENGTH } from './modules/filter_srna_length.nf'
include { CONCAT_TREATED_ONLY_TABLES } from './modules/concat_treated_only_tables.nf'
include { SHORTSTACK } from './modules/shortstack.nf'
include { BOWTIE_ALIGN_TO_PATHOGEN } from './modules/bowtie_align_to_pathogen.nf'
include { BOWTIE_ALIGN_TO_HOST } from './modules/bowtie_align_to_host.nf'
include { BOWTIE_BUILD_PATHOGEN } from './modules/bowtie_build_pathogen.nf'
include { FILTER_PATHOGEN_READS } from './modules/filter_pathogen_reads.nf'
include { BOWTIE_BUILD_HOST } from './modules/bowtie_build_host.nf'
include { LIST_PATHOGEN_READS } from './modules/list_pathogen_reads.nf'
include { BOWTIE_ALIGN_TO_VIRUSES } from './modules/bowtie_align_to_viruses.nf'
include { BOWTIE_BUILD_VIRUSES } from './modules/bowtie_build_viruses.nf'
include { COMPARE_SRNA_LIBRARIES_SIMILAR } from './modules/compare_sRNA_libraries.nf'
include { SHARED_READS } from './modules/shared_reads.nf'
include { CHECK_ALIGNMENT_DISTRIBUTION } from './modules/check_alignment_distribution.nf'
include { CONCAT_ALIGNMENT_DISTRIBUTION } from './modules/concat_alignment_distribution.nf'
include { KEEP_TREATED_ONLY} from './modules/keep_treated_only.nf'
include { BOWTIE_BUILD_CANDIDATE_SRNAS } from './modules/bowtie_build_candidate_srnas.nf'
include { BOWTIE_ALIGN_TO_CANDIDATE_SRNAS } from './modules/bowtie_align_to_candidate_srnas.nf'
include {CHECK_CLUSTER_ANNOTATION } from './modules/check_cluster_annotation.nf'
include {BOWTIE_BUILD } from './modules/bowtie_build.nf'
include {BOWTIE_ALIGN } from './modules/bowtie_align.nf'
include {BOWTIE_ALIGN as BOWTIE_ALIGN_PATHOGEN_PRELIMINARY } from './modules/bowtie_align.nf'
include {BOWTIE_ALIGN as BOWTIE_ALIGN_HOST_PRELIMINARY } from './modules/bowtie_align.nf'
include { CONCAT_TARGETFINDER_RESULTS } from './modules/concat_targetfinder_results.nf'
include { ANNOTATE_TARGETS } from './modules/annotate_targets.nf'
include { PRELIMINARY_MULTIQC } from './modules/preliminary_multiqc.nf'
include { PRELIMINARY_MULTIQC as FASTQC_MULTIQC} from './modules/preliminary_multiqc.nf'
include { PRELIMINARY_MULTIQC as FASTP_TRIM_MULTIQC} from './modules/preliminary_multiqc.nf'
include { PRELIMINARY_MULTIQC as FASTP_TRIM_QC} from './modules/preliminary_multiqc.nf'
include { PRELIMINARY_MULTIQC as VIRUSES_MULTIQC } from './modules/preliminary_multiqc.nf'
include { PRELIMINARY_MULTIQC as MAIN_FILTERING_MULTIQC } from './modules/preliminary_multiqc.nf'
include { MIRNA_TARGET } from './modules/mirna_target.nf'
include { PLOT_SHORTSTACK_CLUSTERS } from './modules/plot_clusters.nf'
include { PLOT_SRNA_LENGTH_DISTRIBUTION } from './modules/plot_length_dist.nf'
include {EGGNOG_ANNOTATION} from './modules/eggnog_annotation.nf'
include {GOATOOLS_GOSLIM} from './modules/goatools_goslim.nf'
include {SUMMARIZE_LOGS} from './modules/summarize_logs.nf'
include {FASTP_TRIM} from './modules/fastp.nf'

params.host_genome_fasta = "data/hops_genome/hops_genome.fa"
params.pathogen_genome_fasta = "data/verticilium_genome/T2/VnonalfalfaeT2.fasta"
params.pathogen_genome_gff = "data/verticilium_genome/T2/Verticillium_nonalfalfae_T2.gff3"
params.sra_csv = "data/sra_accessions.csv"
params.annotated_host_mrnas_fasta = "data/annotated_hop_mrnas/sequence.fasta"
params.mirna_target_repo = "data/MiRNATarget"
params.viruses_genome_fasta = "data/viruses_genome/3viroidi_2virusa.fa"
params.candidate_srnas_fasta = "data/candidate_srnas/candidate_srnas.fa"

workflow {
    // READ INPUT CSV WITH SRA ACCESSIONS
        read_ch = Channel
            .fromPath(params.sra_csv)
            .splitCsv(header: true)
            .map { row ->
                // Collect all non-empty SRA IDs from the row
                def sra_list = row.values()
                                .findAll { it =~ /^SRR/ }   // keep only SRA IDs
                tuple(row.sample_id, sra_list)
            }
    //
    
    // QC AND PREPROCESSING
        FETCH_SRA(read_ch)
        FASTQC(FETCH_SRA.out.reads)
        FASTQC_MULTIQC('non_processed_reads_report', FASTQC.out.html.mix(FASTQC.out.zip).collect())

        FASTP_TRIM(FETCH_SRA.out.reads)
        FASTQC_AFTER_TRIMMING(FASTP_TRIM.out.trimmed_reads)
        FASTP_TRIM_MULTIQC('fastp_trim_report', FASTQC_AFTER_TRIMMING.out.html.mix(FASTQC_AFTER_TRIMMING.out.zip).collect())
        FASTP_TRIM_QC('fastp_trim_qc_report', FASTP_TRIM.out.qc_html.mix(FASTP_TRIM.out.qc_json).collect())
        merged_ch = MERGE_READS(FASTP_TRIM.out.trimmed_reads)
    //

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

    // CHECK FOR VIRUS INFECTIONS
        viruses_index_ch = BOWTIE_BUILD_VIRUSES(file(params.viruses_genome_fasta))
        BOWTIE_ALIGN_TO_VIRUSES(merged_ch, viruses_index_ch)
        VIRUSES_MULTIQC('viruses_alignment_report', BOWTIE_ALIGN_TO_VIRUSES.out.log.collect())
        CHECK_ALIGNMENT_DISTRIBUTION(BOWTIE_ALIGN_TO_VIRUSES.out.results)
        CONCAT_ALIGNMENT_DISTRIBUTION(CHECK_ALIGNMENT_DISTRIBUTION.out.percent_row.collect())
    //

    // MATCHING WITH CANDIDATE sRNAs
        candidate_srna_index_ch = BOWTIE_BUILD_CANDIDATE_SRNAS(file(params.candidate_srnas_fasta))
        BOWTIE_ALIGN_TO_CANDIDATE_SRNAS(FASTP_TRIM.out.trimmed_reads, candidate_srna_index_ch)
    //

    // LENGTH FILTERING
        FILTER_SRNA_LENGTH(merged_ch)
        PLOT_SRNA_LENGTH_DISTRIBUTION(FILTER_SRNA_LENGTH.out.length_distribution)
    //

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
        CONCAT_TREATED_ONLY_TABLES("keep_treated_only", KEEP_TREATED_ONLY.out.log.collect())
    //



    //FILTERING BASED ON READS THAT ALIGN PERFEECTLY TO PATHOGEN AND NOT PERFECTLY TO HOST
        BOWTIE_ALIGN_TO_PATHOGEN(KEEP_TREATED_ONLY.out.filtered_reads, pathogen_index_ch)
        BOWTIE_ALIGN_TO_HOST(BOWTIE_ALIGN_TO_PATHOGEN.out.results, host_index_ch)
        MAIN_FILTERING_MULTIQC(
            'main_filtering_report',
            BOWTIE_ALIGN_TO_PATHOGEN.out.log
                .mix(BOWTIE_ALIGN_TO_HOST.out.log)
                .collect()
        )

        LIST_PATHOGEN_READS(BOWTIE_ALIGN_TO_HOST.out.list_input)
        FILTER_PATHOGEN_READS(LIST_PATHOGEN_READS.out.filter_input)
        filtered_reads = FILTER_PATHOGEN_READS.out.pathogen_specific_reads
    //

    // CLUSTER IDENTIFICATION AND ANNOTATION
        SHORTSTACK(filtered_reads, file(params.pathogen_genome_fasta))
        PLOT_SHORTSTACK_CLUSTERS(SHORTSTACK.out.results, file(params.pathogen_genome_fasta))
        CHECK_CLUSTER_ANNOTATION(SHORTSTACK.out.shortstack_out, file(params.pathogen_genome_gff))
    //

    // USE SHORTSTACK MAJORRNA OUTPUT TO PREDICT HOP TARGETS
        MIRNA_TARGET(CHECK_CLUSTER_ANNOTATION.out.non_protein_coding_fasta, file(params.annotated_host_mrnas_fasta), file(params.mirna_target_repo))
    //

    // ANNOTATE TARGETS, GO ANALYSIS
        ANNOTATE_TARGETS(MIRNA_TARGET.out.results, file(params.annotated_host_mrnas_fasta))
        EGGNOG_ANNOTATION(ANNOTATE_TARGETS.out.fastas)
        GOATOOLS_GOSLIM(EGGNOG_ANNOTATION.out.go_terms)
    //


    // COLLECT LOGS
        SUMMARIZE_LOGS(
            FILTER_SRNA_LENGTH.out.log_info
            .mix(
                FETCH_SRA.out.log_info,
                FASTP_TRIM.out.log_info,
                BOWTIE_ALIGN_TO_VIRUSES.out.log_info,
                KEEP_TREATED_ONLY.out.log_info,
                BOWTIE_ALIGN_TO_PATHOGEN.out.log_info,
                BOWTIE_ALIGN_TO_HOST.out.log_info,
                FILTER_PATHOGEN_READS.out.log_info,
                SHORTSTACK.out.log_info,
                CHECK_CLUSTER_ANNOTATION.out.log_info,
                MIRNA_TARGET.out.log_info
            )
            .collect()
        )
    //
}

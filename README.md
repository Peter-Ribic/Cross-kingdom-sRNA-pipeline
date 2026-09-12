# Cross Kingdom sRNA Pipeline

## Overview
A Nextflow pipeline for the analysis of small RNA (sRNA) sequencing data from a cross-kingdom RNAi study. The data originate from hop (*Humulus lupulus*) samples — cultivars Celeia and Wye Target — infected with the fungal pathogen *Verticillium nonalfalfae*, along with matching non-infected (control) samples.

The main goals of the pipeline are to:

1. **Filter out fungal sRNA** — identify reads that align to the *V. nonalfalfae* genome but not to the hop genome, thereby isolating putatively pathogen-derived sRNAs.
2. **Identify sRNA clusters** — use ShortStack to group pathogen-specific reads into genomic clusters and characterise them.
3. **Scan the hop transcriptome for potential RNAi targets** — predict hop mRNA targets of the candidate pathogen sRNAs.
4. **Perform enrichment analysis** — annotate predicted targets (eggNOG, GO slim) and assess functional enrichment among them.

Additional quality-control steps are included along the way:

- Raw read retrieval from SRA and initial QC (FastQC / MultiQC).
- Adapter and quality trimming (fastp) with post-trimming QC.
- Preliminary alignment to host and pathogen genomes to assess contamination levels.
- Screening for virus/viroid co-infections in the samples.
- sRNA length filtering with length-distribution plots.
- Pairwise similarity comparison of sRNA libraries (shared reads).
- Filtering of reads present in treated but not control samples, and a final summary of all filtering logs.

The pipeline expects an input CSV file listing SRA accession numbers per sample, with sample IDs suffixed `_treated` or `_control`.

## Prerequisites
- Linux
- Nextflow version 25.04.7
- Git
- DVC

## Setup

Clone the repository and pull the data files:

```bash
git clone https://github.com/Peter-Ribic/Cross-kingdom-sRNA-pipeline.git
cd Cross-kingdom-sRNA-pipeline
dvc pull
```

## How to run
While in Cross-kingdom-sRNA-pipeline directory:
```bash
nextflow run rnaseq.nf
```
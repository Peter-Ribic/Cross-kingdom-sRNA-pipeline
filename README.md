# Cross Kingdom sRNA Pipeline
Pipeline to analyse sRNA, obtained from infected and non-infected hop samples (cultivars Celeia and Wye Target). Main goal of the pipeline was to filter out fungal sRNA, identify sRNA clusters and scan hop transcriptome for potential RNAi targets. Targets are then subjugated to enrichment analysis.

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
If dvc pull does not work, download data/ and results/ folders manually from https://dagshub.com/peterribic0/Cross-kingdom-sRNA-pipeline.
## How to run
While in Cross-kingdom-sRNA-pipeline directory:
```bash
nextflow run rnaseq.nf
```
Running is not required, as all results are already in the results/ folder at https://dagshub.com/peterribic0/Cross-kingdom-sRNA-pipeline.
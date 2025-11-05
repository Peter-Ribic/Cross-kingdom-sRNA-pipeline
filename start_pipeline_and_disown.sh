#!/bin/bash
nohup nextflow run rnaseq.nf -resume > nextflow.log 2>&1 &
echo "Nextflow pipeline started in background."
echo "Logs are being written to nextflow.log"
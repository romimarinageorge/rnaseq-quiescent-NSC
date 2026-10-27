#!/bin/bash
# Module 1: quality control of raw reads (Experiment B)
set -e
cd ~/rnaseq_qNSC
mkdir -p results/qc/fastqc results/qc/multiqc
fastqc --threads 8 --outdir results/qc/fastqc data/raw/*.fastq.gz
multiqc results/qc/fastqc --outdir results/qc/multiqc
echo "=== QC finished ==="

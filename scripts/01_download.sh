#!/bin/bash
# Download Experiment B (GSE116997, paired-end) and convert to compressed FASTQ
set -e
cd ~/rnaseq_qNSC/data/raw
for srr in SRR7514765 SRR7514766 SRR7514767 SRR7514768 SRR7514769; do
    echo "=== Processing $srr ==="
    prefetch $srr
    fasterq-dump --split-files --threads 8 $srr
    pigz -p 8 ${srr}_1.fastq ${srr}_2.fastq
done
echo "=== All downloads finished ==="

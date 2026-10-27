#!/bin/bash
# Module 2a: build a decoy-aware Salmon index for mouse (GENCODE, GRCm39)
set -e
cd ~/rnaseq_qNSC/ref
BASE=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/latest_release
REL=$(wget -qO- $BASE/ | grep -o 'gencode\.vM[0-9]*\.transcripts\.fa\.gz' | head -1 | sed 's/gencode\.\(vM[0-9]*\)\..*/\1/')
if [ -z "$REL" ]; then echo "Could not detect GENCODE release"; exit 1; fi
echo "=== GENCODE release: $REL ==="
wget -c $BASE/gencode.$REL.transcripts.fa.gz
wget -c $BASE/gencode.$REL.annotation.gtf.gz
wget -c $BASE/GRCm39.primary_assembly.genome.fa.gz
echo "=== Making decoy list ==="
zcat GRCm39.primary_assembly.genome.fa.gz | grep '^>' | cut -d ' ' -f 1 | sed 's/>//' > decoys.txt
echo "=== Combining transcriptome + genome ==="
cat gencode.$REL.transcripts.fa.gz GRCm39.primary_assembly.genome.fa.gz > gentrome.fa.gz
echo "=== Building index ==="
salmon index -t gentrome.fa.gz -d decoys.txt -i salmon_index_$REL --gencode -k 31 -p 8
echo "$REL" > gencode_release.txt
echo "=== Index finished ==="

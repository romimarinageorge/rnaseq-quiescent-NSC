#!/bin/bash
# Module 3a: build the transcript-to-gene table that tximport needs
# Reads the GENCODE annotation (GTF) and writes one row per transcript:
#   transcript_id  gene_id  gene_name  gene_type
# Run:  bash ~/rnaseq_qNSC/scripts/06_tx2gene.sh
# Time: about 1 minute
set -euo pipefail
cd ~/rnaseq_qNSC/ref
REL=$(cat gencode_release.txt)

zcat gencode.$REL.primary_assembly.annotation.gtf.gz \
  | awk -F '\t' '
      BEGIN { print "transcript_id\tgene_id\tgene_name\tgene_type" }
      $3 == "transcript" {
        match($9, /transcript_id "[^"]+"/); tx = substr($9, RSTART + 15, RLENGTH - 16)
        match($9, /gene_id "[^"]+"/);       g  = substr($9, RSTART + 9,  RLENGTH - 10)
        match($9, /gene_name "[^"]+"/);     n  = substr($9, RSTART + 11, RLENGTH - 12)
        match($9, /gene_type "[^"]+"/);     t  = substr($9, RSTART + 11, RLENGTH - 12)
        print tx "\t" g "\t" n "\t" t
      }' > tx2gene_$REL.tsv

echo "Rows (transcripts + 1 header line):"
wc -l tx2gene_$REL.tsv
echo "First rows:"
head -5 tx2gene_$REL.tsv | column -t

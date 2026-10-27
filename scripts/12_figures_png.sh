#!/bin/bash
# Module 5a: make PNG copies of the key figures so GitHub can show them inside the README
# (GitHub shows PDFs only as download links, not as pictures.)
# One-time install:  sudo apt install -y poppler-utils
# Run:  bash scripts/12_figures_png.sh      (from ~/rnaseq_qNSC, any environment)
set -euo pipefail
cd ~/rnaseq_qNSC

for pdf in results/deseq2/figures/PCA.pdf \
           results/deseq2/figures/volcano.pdf \
           results/deseq2/figures/marker_genes.pdf \
           results/enrichment/GO_BP_GSEA_dotplot.pdf
do
  png="${pdf%.pdf}"                          # same path without ".pdf"
  pdftoppm -png -r 150 -singlefile "$pdf" "$png"
  echo "made ${png}.png"
done

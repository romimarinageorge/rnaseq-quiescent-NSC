#!/bin/bash
# Module 2b: quantify every Experiment B sample with Salmon
set -euo pipefail
cd ~/rnaseq_qNSC
INDEX=ref/salmon_index_$(cat ref/gencode_release.txt)
tail -n +2 metadata/samples_expB.tsv | while read sample run condition replicate; do
    if [ -f results/salmon/$sample/quant.sf ]; then
        echo "=== $sample already done, skipping ==="
        continue
    fi
    echo "=== Quantifying $sample ($run, $condition) ==="
    salmon quant -i $INDEX -l A \
        -1 data/raw/${run}_1.fastq.gz -2 data/raw/${run}_2.fastq.gz \
        -p 8 --seqBias --gcBias \
        -o results/salmon/$sample < /dev/null
done
echo "=== All samples quantified ==="

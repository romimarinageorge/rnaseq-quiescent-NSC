#!/bin/bash
# Module 2c: collect Salmon QC numbers for every sample into one table
# Run from anywhere:  bash ~/rnaseq_qNSC/scripts/05_salmon_summary.sh
# Time: a few seconds
set -euo pipefail
cd ~/rnaseq_qNSC

OUT=results/salmon/salmon_summary.tsv
echo -e "sample\tcondition\tmapping_rate\tassigned_pairs\tlibrary_type\twrong_strand_fraction" > $OUT

tail -n +2 metadata/samples_expB.tsv | while read sample run condition replicate; do
    dir=results/salmon/$sample
    rate=$(grep -i "mapping rate" $dir/logs/salmon_quant.log | grep -o '[0-9.]*%')
    lib=$(grep -o '"expected_format": "[A-Z]*"' $dir/lib_format_counts.json | grep -o '[A-Z]\{2,\}')
    assigned=$(grep -o '"num_assigned_fragments": [0-9]*' $dir/lib_format_counts.json | grep -o '[0-9]*$')
    bias=$(grep -o '"strand_mapping_bias": [0-9.e-]*' $dir/lib_format_counts.json | grep -o '[0-9.e-]*$')
    echo -e "$sample\t$condition\t$rate\t$assigned\t$lib\t$bias" >> $OUT
done

column -t $OUT

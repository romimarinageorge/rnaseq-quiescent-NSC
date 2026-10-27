# Print the three key results of script 09 again, without rerunning the analysis
suppressPackageStartupMessages(library(DESeq2))

# 1. DESeq2's own summary (padj < 0.05, any size of change)
dds <- readRDS("results/deseq2/dds_fitted.rds")
res <- results(dds, name = "condition_quiescent_vs_proliferating", alpha = 0.05)
cat("=== 1. Summary at padj < 0.05 ===\n")
summary(res)

# 2. Significant genes by both rules (padj < 0.05 and at least 2-fold)
tab <- read.delim("results/deseq2/DE_quiescent_vs_proliferating_all.tsv")
cat("\n=== 2. Significant genes (padj < 0.05 and at least 2-fold) ===\n")
print(table(tab$direction))

# 3. Marker genes
markers <- c("Mki67", "Top2a", "Cdk1", "Ccnb1", "Id1", "Id2", "Id3", "Id4", "Ascl1")
mk <- tab[tab$gene_name %in% markers,
          c("gene_name", "baseMean", "log2FC_shrunken", "padj", "direction")]
mk <- mk[match(markers, mk$gene_name), ]
mk$fold_change <- round(2^mk$log2FC_shrunken, 2)
mk$baseMean <- round(mk$baseMean)
mk$log2FC_shrunken <- round(mk$log2FC_shrunken, 2)
mk$padj <- signif(mk$padj, 3)
cat("\n=== 3. Marker genes ===\n")
print(mk, row.names = FALSE)

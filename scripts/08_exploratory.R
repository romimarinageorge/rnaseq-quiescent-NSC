# Module 3c: exploratory analysis - do the replicates group by condition?
# Run:  Rscript scripts/08_exploratory.R      (from ~/rnaseq_qNSC, env rnaseq-r)
# Time: under 1 minute
# Output: results/deseq2/figures/*.pdf (open them from Windows File Explorer)

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
})

fig <- "results/deseq2/figures"
dir.create(fig, showWarnings = FALSE, recursive = TRUE)

dds <- readRDS("results/deseq2/dds_unfitted.rds")

# ---- 1. Size factors: how much each library is scaled for depth --------------
dds <- estimateSizeFactors(dds)
cat("Normalisation factors combine library size and tximport's gene lengths.\n")
cat("Reads in the kept genes, per sample (millions):\n")
print(round(colSums(counts(dds)) / 1e6, 1))

# ---- 2. Variance-stabilising transformation -----------------------------------
# blind = TRUE: the transformation ignores the condition labels,
# so the plots below are an unbiased look at the data
vsd <- vst(dds, blind = TRUE)
saveRDS(vsd, "results/deseq2/vsd.rds")

# ---- 3. PCA ----------------------------------------------------------------------
pca <- plotPCA(vsd, intgroup = "condition", ntop = 500, returnData = TRUE)
pv  <- round(100 * attr(pca, "percentVar"))
p <- ggplot(pca, aes(PC1, PC2, colour = condition)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = name), show.legend = FALSE) +
  labs(x = paste0("PC1: ", pv[1], "% of variance"),
       y = paste0("PC2: ", pv[2], "% of variance"),
       title = "PCA of the 500 most variable genes (VST)") +
  theme_bw()
ggsave(file.path(fig, "PCA.pdf"), p, width = 6, height = 4.5)
write.table(pca, "results/deseq2/PCA_coordinates.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# ---- 4. Sample-to-sample distances -----------------------------------------------
d <- dist(t(assay(vsd)))
m <- as.matrix(d)
ann <- data.frame(condition = colData(vsd)$condition, row.names = colnames(vsd))
pheatmap(m, clustering_distance_rows = d, clustering_distance_cols = d,
         annotation_col = ann, display_numbers = TRUE, number_format = "%.0f",
         main = "Euclidean distance between samples (VST)",
         filename = file.path(fig, "sample_distances.pdf"), width = 6, height = 5)

cat("Done: PCA.pdf and sample_distances.pdf in", fig, "\n")

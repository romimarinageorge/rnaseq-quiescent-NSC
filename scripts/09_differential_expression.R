# Module 3d: differential expression, quiescent (FGF2+BMP4) vs proliferating (FGF2)
# Run:  Rscript scripts/09_differential_expression.R   (from ~/rnaseq_qNSC, env rnaseq-r)
# Time: about 2-5 minutes (the apeglm shrinkage step is the slowest)
# Output: results/deseq2/*.tsv and results/deseq2/figures/*.pdf

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
})

fig <- "results/deseq2/figures"
dir.create(fig, showWarnings = FALSE, recursive = TRUE)
padj_cut <- 0.05   # false discovery rate threshold
lfc_cut  <- 1      # |log2 fold change| >= 1 means at least 2-fold

dds   <- readRDS("results/deseq2/dds_unfitted.rds")
vsd   <- readRDS("results/deseq2/vsd.rds")
genes <- read.delim("results/deseq2/gene_annotation.tsv")

# ---- 1. Fit the model ----------------------------------------------------------
# DESeq() = size factors + dispersion (gene-wise variability) + negative binomial fit + Wald test
dds <- DESeq(dds)
saveRDS(dds, "results/deseq2/dds_fitted.rds")
print(resultsNames(dds))
coef_name <- "condition_quiescent_vs_proliferating"

pdf(file.path(fig, "dispersion_estimates.pdf"), width = 6, height = 5)
plotDispEsts(dds)
invisible(dev.off())

# ---- 2. Test results -----------------------------------------------------------
res <- results(dds, name = coef_name, alpha = padj_cut)
cat("\nSummary at padj <", padj_cut, ":\n")
summary(res)

# Shrunken log2 fold changes (apeglm): more reliable effect sizes for low-count genes.
# Use these for ranking and plotting; p-values come from 'res'.
res_shr <- lfcShrink(dds, coef = coef_name, type = "apeglm", res = res)

# ---- 3. One annotated results table ---------------------------------------------
tab <- data.frame(
  gene_id         = rownames(res),
  baseMean        = res$baseMean,
  log2FC          = res$log2FoldChange,
  log2FC_shrunken = res_shr$log2FoldChange,
  lfcSE           = res$lfcSE,
  stat            = res$stat,
  pvalue          = res$pvalue,
  padj            = res$padj
)
tab <- merge(genes, tab, by = "gene_id", all.y = TRUE)
tab <- tab[order(tab$padj), ]
tab$significant <- !is.na(tab$padj) & tab$padj < padj_cut & abs(tab$log2FC_shrunken) >= lfc_cut
tab$direction <- ifelse(!tab$significant, "ns",
                        ifelse(tab$log2FC_shrunken > 0, "up_in_quiescent", "down_in_quiescent"))

write.table(tab, "results/deseq2/DE_quiescent_vs_proliferating_all.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tab[tab$significant, ], "results/deseq2/DE_quiescent_vs_proliferating_significant.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nSignificant genes (padj <", padj_cut, "and |shrunken log2FC| >=", lfc_cut, "):\n")
print(table(tab$direction))

# Normalised counts, for browsing and plotting
nc <- counts(dds, normalized = TRUE)
nc_tab <- merge(genes[, c("gene_id", "gene_name")],
                data.frame(gene_id = rownames(nc), round(nc, 1), check.names = FALSE), by = "gene_id")
write.table(nc_tab, "results/deseq2/normalized_counts.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# ---- 4. MA plot -------------------------------------------------------------------
pdf(file.path(fig, "MA_plot_shrunken.pdf"), width = 6, height = 5)
plotMA(res_shr, alpha = padj_cut, ylim = c(-8, 8), main = "Quiescent vs proliferating (shrunken LFC)")
invisible(dev.off())

# ---- 5. Volcano plot ---------------------------------------------------------------
v <- tab[!is.na(tab$padj), ]
v$neglog10padj <- -log10(pmax(v$padj, 1e-300))
lab <- head(v[v$significant, ], 20)
p <- ggplot(v, aes(log2FC_shrunken, neglog10padj, colour = direction)) +
  geom_point(size = 0.8, alpha = 0.6) +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", colour = "grey50") +
  geom_hline(yintercept = -log10(padj_cut), linetype = "dashed", colour = "grey50") +
  geom_text_repel(data = lab, aes(label = gene_name), size = 3, max.overlaps = 30, show.legend = FALSE) +
  scale_colour_manual(values = c(down_in_quiescent = "#2166ac", ns = "grey70", up_in_quiescent = "#b2182b")) +
  labs(x = "log2 fold change (quiescent / proliferating, shrunken)", y = "-log10 adjusted p-value",
       title = "Quiescent vs proliferating NSCs") +
  theme_bw()
ggsave(file.path(fig, "volcano.pdf"), p, width = 7, height = 5.5)

# ---- 6. Heatmap of the 50 most significant genes ---------------------------------
top <- head(tab$gene_id[tab$significant], 50)
if (length(top) >= 2) {
  mat <- assay(vsd)[top, ]
  mat <- mat - rowMeans(mat)                     # centre each gene on its mean
  labels <- genes$gene_name[match(top, genes$gene_id)]
  ann <- data.frame(condition = colData(vsd)$condition, row.names = colnames(vsd))
  pheatmap(mat, annotation_col = ann, labels_row = labels, fontsize_row = 7,
           main = "Top 50 genes (VST, centred per gene)",
           filename = file.path(fig, "heatmap_top50.pdf"), width = 6, height = 9)
}

# ---- 7. Sanity check with known genes ---------------------------------------------
# Expected: proliferation genes DOWN in quiescence; BMP target genes (Id family) UP.
# Blomfield et al. 2019 (eLife) report Ascl1 mRNA unchanged in FGF2+BMP4 quiescent NSCs.
markers <- c("Mki67", "Top2a", "Cdk1", "Ccnb1", "Id1", "Id2", "Id3", "Id4", "Ascl1")
mk <- tab[tab$gene_name %in% markers,
          c("gene_name", "baseMean", "log2FC_shrunken", "padj", "direction")]
cat("\nMarker genes:\n")
print(mk[match(markers, mk$gene_name), ], row.names = FALSE)

long <- do.call(rbind, lapply(intersect(markers, mk$gene_name), function(g) {
  id <- genes$gene_id[match(g, genes$gene_name)]
  if (!id %in% rownames(nc)) return(NULL)
  data.frame(gene = g, sample = colnames(nc), count = nc[id, ],
             condition = colData(dds)$condition)
}))
long$gene <- factor(long$gene, levels = markers)
p2 <- ggplot(long, aes(condition, count + 1, colour = condition)) +
  geom_point(size = 2.5) +
  scale_y_log10() +
  facet_wrap(~ gene, scales = "free_y") +
  labs(y = "normalised count + 1 (log scale)", x = NULL, title = "Marker genes") +
  theme_bw() + theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(fig, "marker_genes.pdf"), p2, width = 8, height = 7)

# ---- 8. Record software versions ----------------------------------------------
writeLines(capture.output(sessionInfo()), "results/deseq2/sessionInfo.txt")
cat("\nDone. Tables in results/deseq2, figures in", fig, "\n")

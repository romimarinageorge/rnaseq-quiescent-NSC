# Module 4: which biological processes change between quiescent and proliferating NSCs?
#
# One-time install first (about 10-30 minutes, mostly download):
#   mamba install -n rnaseq-r -c conda-forge -c bioconda -y \
#     bioconductor-clusterprofiler bioconductor-org.mm.eg.db bioconductor-enrichplot
#
# Run:  Rscript scripts/10_enrichment.R   (from ~/rnaseq_qNSC, env rnaseq-r)
# Time: about 5-15 minutes (GSEA is the slowest part)

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(enrichplot)
  library(ggplot2)
})

out <- "results/enrichment"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
set.seed(42)   # GSEA uses random permutations; fixing the seed makes it reproducible

tab <- read.delim("results/deseq2/DE_quiescent_vs_proliferating_all.tsv")
# org.Mm.eg.db uses Ensembl IDs without the version suffix (".2" etc.)
tab$ensembl <- sub("\\..*$", "", tab$gene_id)
tested <- tab[!is.na(tab$padj), ]

# ---- 1. Over-representation analysis (ORA) --------------------------------------
# Question: are GO terms over-represented among the significant genes,
# compared with all genes that were tested (the "universe")?
universe <- unique(tested$ensembl)
up   <- unique(tested$ensembl[tested$direction == "up_in_quiescent"])
down <- unique(tested$ensembl[tested$direction == "down_in_quiescent"])
cat("Genes: universe", length(universe), "| up", length(up), "| down", length(down), "\n")

run_ora <- function(ids, label) {
  ego <- enrichGO(gene = ids, universe = universe, OrgDb = org.Mm.eg.db,
                  keyType = "ENSEMBL", ont = "BP", pAdjustMethod = "BH",
                  pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    cat("No enriched GO terms for", label, "\n")
    return(invisible(NULL))
  }
  write.table(as.data.frame(ego), file.path(out, paste0("GO_BP_ORA_", label, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  ggsave(file.path(out, paste0("GO_BP_ORA_", label, ".pdf")),
         dotplot(ego, showCategory = 20) + ggtitle(paste("GO BP, genes", label)),
         width = 8, height = 8)
  cat("ORA", label, ":", nrow(as.data.frame(ego)), "enriched terms\n")
}
run_ora(up,   "up_in_quiescent")
run_ora(down, "down_in_quiescent")

# ---- 2. Gene set enrichment analysis (GSEA) --------------------------------------
# Uses ALL tested genes, ranked by the Wald statistic (direction x confidence),
# so no significance cut-off is needed.
ranked <- tested[order(tested$stat, decreasing = TRUE), ]
ranked <- ranked[!duplicated(ranked$ensembl), ]
geneList <- setNames(ranked$stat, ranked$ensembl)

gse <- gseGO(geneList = geneList, OrgDb = org.Mm.eg.db, keyType = "ENSEMBL",
             ont = "BP", minGSSize = 15, maxGSSize = 500,
             pvalueCutoff = 0.05, pAdjustMethod = "BH", eps = 0, seed = TRUE)
gse <- setReadable(gse, OrgDb = org.Mm.eg.db, keyType = "ENSEMBL")
gse_tab <- as.data.frame(gse)
write.table(gse_tab, file.path(out, "GO_BP_GSEA.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("GSEA:", nrow(gse_tab), "significant GO terms",
    "(NES > 0 = higher in quiescent, NES < 0 = higher in proliferating)\n")

if (nrow(gse_tab) > 0) {
  ggsave(file.path(out, "GO_BP_GSEA_dotplot.pdf"),
         dotplot(gse, showCategory = 10, split = ".sign") + facet_grid(. ~ .sign) +
           ggtitle("GSEA, GO biological process"),
         width = 11, height = 8)
  ggsave(file.path(out, "GO_BP_GSEA_top_term.pdf"),
         gseaplot2(gse, geneSetID = 1, title = gse_tab$Description[1]),
         width = 7, height = 5)
}

writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
cat("Done. Results in", out, "\n")

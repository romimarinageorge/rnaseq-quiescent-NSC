# Module 3b: import Salmon results into R and sum transcripts into genes
# Run from the project folder:
#   cd ~/rnaseq_qNSC
#   mamba activate rnaseq-r
#   Rscript scripts/07_import_tximport.R
# Time: about 1-3 minutes (reads 6 tables of 480,180 rows each)

suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
})

dir.create("results/deseq2", showWarnings = FALSE, recursive = TRUE)

# ---- 1. Sample sheet -------------------------------------------------------
samples <- read.delim("metadata/samples_expB.tsv")
# Reference level first: log2 fold changes will be quiescent relative to proliferating
samples$condition <- factor(samples$condition, levels = c("proliferating", "quiescent"))
rownames(samples) <- samples$sample
print(samples)

# ---- 2. Paths to the six quant.sf files -------------------------------------
files <- file.path("results/salmon", samples$sample, "quant.sf")
names(files) <- samples$sample
stopifnot(all(file.exists(files)))

# ---- 3. Transcript-to-gene table ---------------------------------------------
release <- readLines("ref/gencode_release.txt")
tx2gene <- read.delim(file.path("ref", paste0("tx2gene_", release, ".tsv")))
cat("Transcripts in tx2gene:", nrow(tx2gene), "\n")

# ---- 4. Import and summarise to gene level -----------------------------------
# tximport needs only the first two columns: transcript_id, gene_id
txi <- tximport(files, type = "salmon", tx2gene = tx2gene[, c("transcript_id", "gene_id")])
cat("Genes after summarising:", nrow(txi$counts), "\n")

# ---- 5. Build the DESeq2 object ---------------------------------------------
# design = ~ condition : model each gene's counts as depending on condition
dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ condition)

# Remove genes with almost no reads (at least 10 reads in at least 3 samples;
# 3 = the size of the smallest group)
keep <- rowSums(counts(dds) >= 10) >= 3
cat("Genes kept after filtering:", sum(keep), "of", length(keep), "\n")
dds <- dds[keep, ]

# ---- 6. Save --------------------------------------------------------------------
saveRDS(txi, "results/deseq2/txi.rds")
saveRDS(dds, "results/deseq2/dds_unfitted.rds")

# Gene annotation table (one row per gene), used by the later scripts
genes <- unique(tx2gene[, c("gene_id", "gene_name", "gene_type")])
write.table(genes, "results/deseq2/gene_annotation.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Gene-level tables for browsing (not inputs to DESeq2)
counts_tab <- data.frame(gene_id = rownames(txi$counts), round(txi$counts, 2), check.names = FALSE)
tpm_tab    <- data.frame(gene_id = rownames(txi$abundance), round(txi$abundance, 3), check.names = FALSE)
counts_tab <- merge(genes[, c("gene_id", "gene_name")], counts_tab, by = "gene_id")
tpm_tab    <- merge(genes[, c("gene_id", "gene_name")], tpm_tab, by = "gene_id")
write.table(counts_tab, "results/deseq2/gene_counts_estimated.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tpm_tab,    "results/deseq2/gene_TPM.tsv",              sep = "\t", quote = FALSE, row.names = FALSE)

cat("Done: results/deseq2/dds_unfitted.rds and gene tables written\n")

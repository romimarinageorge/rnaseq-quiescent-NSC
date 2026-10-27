library(ggplot2)
qc <- read.delim("results/salmon/salmon_summary.tsv")
qc$mapping_rate <- as.numeric(sub("%", "", qc$mapping_rate))
qc$sample <- factor(qc$sample, levels = qc$sample)
p1 <- ggplot(qc, aes(sample, assigned_pairs / 1e6, fill = condition)) +
  geom_col() +
  labs(x = NULL, y = "Read pairs counted (millions)", title = "Sequencing depth per sample") +
  theme_bw()
p2 <- ggplot(qc, aes(sample, mapping_rate, colour = condition)) +
  geom_point(size = 3) +
  coord_cartesian(ylim = c(80, 100)) +
  labs(x = NULL, y = "Mapping rate (%)", title = "Mapping rate per sample") +
  theme_bw()
pdf("results/salmon/salmon_qc_plots.pdf", width = 7, height = 4.5)
print(p1)
print(p2)
invisible(dev.off())
cat("Saved results/salmon/salmon_qc_plots.pdf\n")

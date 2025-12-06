#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyverse)
  library(glue)
  library(tximport)
  library(DESeq2)
  library(pheatmap)
  library(ggrepel)
})

metadata <- read_csv("analysis/metadata_samples.csv", show_col_types = FALSE)
stopifnot(all(file.exists(glue("salmon/{metadata$sample_id}/quant.sf"))))

# TODO: Replace with a transcript-to-gene mapping matching your reference.
tx2gene <- read_tsv("refs/tx2gene.tsv", col_types = cols())

files <- file.path("salmon", metadata$sample_id, "quant.sf")
names(files) <- metadata$sample_id

# Import Salmon TPM and counts
message("[INFO] Importing Salmon quantifications")
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, countsFromAbundance = "lengthScaledTPM")

# Construct DESeq2 dataset
message("[INFO] Building DESeq2 dataset")
coldata <- metadata %>%
  mutate(condition = factor(condition),
         time_point_dpa = factor(time_point_dpa, levels = sort(unique(time_point_dpa)))) %>%
  column_to_rownames("sample_id")

dds <- DESeqDataSetFromTximport(txi, colData = coldata, design = ~ time_point_dpa + condition + time_point_dpa:condition)

message("[INFO] Running DESeq2")
dds <- DESeq(dds)

# Variance stabilized counts for PCA/heatmap
vsd <- vst(dds, blind = FALSE)
dir.create("results/deseq2", recursive = TRUE, showWarnings = FALSE)

# PCA plot
pca_data <- plotPCA(vsd, intgroup = c("condition", "time_point_dpa"), returnData = TRUE)
p <- ggplot(pca_data, aes(PC1, PC2, color = condition, shape = time_point_dpa, label = name)) +
  geom_point(size = 4) +
  geom_text_repel(max.overlaps = 10) +
  theme_minimal(base_size = 14)
ggsave("results/deseq2/pca.png", p, width = 7, height = 5, dpi = 300)

# Differential expression example: wtap_RNAi vs control at 7 dpa
res_7dpa <- results(dds, contrast = list(c("condition_wtap_RNAi_vs_control", "time_point_dpa7.conditionwtap_RNAi"))) %>%
  as_tibble(rownames = "gene") %>%
  arrange(padj)
write_csv(res_7dpa, "results/deseq2/deseq2_wtap_vs_control_7dpa.csv")

# Volcano plot
volcano <- res_7dpa %>%
  mutate(sig = !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1)
vol <- ggplot(volcano, aes(log2FoldChange, -log10(padj))) +
  geom_point(aes(color = sig), alpha = 0.6) +
  scale_color_manual(values = c("gray70", "firebrick")) +
  theme_minimal(base_size = 14) +
  labs(title = "wtap RNAi vs control at 7 dpa", x = "log2FC", y = "-log10 adj p")
ggsave("results/deseq2/volcano_7dpa.png", vol, width = 6, height = 5, dpi = 300)

# Sample distance heatmap
sample_dists <- dist(t(assay(vsd)))
dist_mat <- as.matrix(sample_dists)
pheatmap(dist_mat, clustering_distance_rows = sample_dists, clustering_distance_cols = sample_dists,
         filename = "results/deseq2/sample_dist_heatmap.png", width = 6, height = 5)

message("[DONE] DESeq2 analysis completed.")

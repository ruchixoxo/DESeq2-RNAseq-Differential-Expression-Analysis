
############################################################
# PROJECT: Differential Gene Expression Analysis of RNA-seq
# USING DESeq2
#
# Dataset: airway (Bioconductor)
# Organism: Homo sapiens
# Comparison: Dexamethasone-treated vs untreated
#
# Author: Ruchi
############################################################


############################################################
# 1. INSTALL AND LOAD REQUIRED PACKAGES
############################################################

# Install BiocManager if required
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install required Bioconductor packages if missing
bioc_packages <- c(
  "airway",
  "DESeq2",
  "clusterProfiler",
  "org.Hs.eg.db",
  "enrichplot",
  "AnnotationDbi"
)

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

# Install CRAN packages if missing
cran_packages <- c("ggplot2", "pheatmap")

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# Load packages
library(airway)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(AnnotationDbi)


############################################################
# 2. CREATE PROJECT DIRECTORIES
############################################################

dir.create("Results", showWarnings = FALSE)
dir.create("Figures", showWarnings = FALSE)


############################################################
# 3. LOAD AND PREPARE DATASET
############################################################

data("airway")

# Examine dataset
airway

# Set untreated samples as reference
colData(airway)$dex <- relevel(
  colData(airway)$dex,
  ref = "untrt"
)

# Construct DESeqDataSet
dds <- DESeqDataSet(
  airway,
  design = ~ cell + dex
)

# Filter low-count genes
dds <- dds[rowSums(counts(dds)) > 1, ]


############################################################
# 4. DIFFERENTIAL GENE EXPRESSION ANALYSIS
############################################################

# Run DESeq2
dds <- DESeq(dds)

# Extract treatment comparison
res <- results(
  dds,
  contrast = c("dex", "trt", "untrt")
)

# Order results by adjusted p-value
res <- res[order(res$padj), ]

# Convert results to data frame
res_df <- as.data.frame(res)

# Export all differential expression results
write.csv(
  res_df,
  "Results/DESeq2_all_results.csv"
)

# Identify significant genes
significant_genes <- subset(
  res_df,
  padj < 0.05 & abs(log2FoldChange) > 1
)

# Separate upregulated and downregulated genes
upregulated <- subset(
  significant_genes,
  log2FoldChange > 0
)

downregulated <- subset(
  significant_genes,
  log2FoldChange < 0
)

# Export results
write.csv(
  significant_genes,
  "Results/DESeq2_significant_genes.csv"
)

write.csv(
  upregulated,
  "Results/Upregulated_Genes.csv"
)

write.csv(
  downregulated,
  "Results/Downregulated_Genes.csv"
)

# Display gene counts
cat("Total significant genes:", nrow(significant_genes), "\n")
cat("Upregulated genes:", nrow(upregulated), "\n")
cat("Downregulated genes:", nrow(downregulated), "\n")


############################################################
# 5. DATA VISUALIZATION
############################################################

# Variance stabilizing transformation
vsd <- vst(dds, blind = FALSE)


# 5A. PCA PLOT

PCA_plot <- plotPCA(
  vsd,
  intgroup = "dex"
) +
  ggtitle("PCA of Airway RNA-seq Samples") +
  theme_minimal()

ggsave(
  "Figures/PCA_Plot.png",
  plot = PCA_plot,
  width = 8,
  height = 6,
  dpi = 300
)


# 5B. VOLCANO PLOT

volcano_data <- as.data.frame(res)

volcano_data$Significance <- "Not Significant"

volcano_data$Significance[
  !is.na(volcano_data$padj) &
  volcano_data$padj < 0.05 &
  abs(volcano_data$log2FoldChange) > 1
] <- "Significant"

Volcano_plot <- ggplot(
  volcano_data,
  aes(
    x = log2FoldChange,
    y = -log10(padj),
    color = Significance
  )
) +
  geom_point(
    alpha = 0.6,
    size = 1.2,
    na.rm = TRUE
  ) +
  scale_color_manual(
    values = c(
      "Significant" = "red",
      "Not Significant" = "grey"
    )
  ) +
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  labs(
    title = "Volcano Plot of Differentially Expressed Genes",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value"
  ) +
  theme_minimal()

ggsave(
  "Figures/Volcano_Plot.png",
  plot = Volcano_plot,
  width = 8,
  height = 6,
  dpi = 300
)


# 5C. HEATMAP OF TOP 30 GENES

top_genes <- rownames(
  head(
    res[!is.na(res$padj), ][order(res[!is.na(res$padj), ]$padj), ],
    30
  )
)

mat <- assay(vsd)[top_genes, ]

# Row-wise scaling
mat_scaled <- t(scale(t(mat)))

# Generate heatmap
pheatmap(
  mat_scaled,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  filename = "Figures/Heatmap.png",
  width = 10,
  height = 8,
  main = "Top 30 Differentially Expressed Genes"
)


############################################################
# 6. GENE ID CONVERSION
############################################################

# Prepare Ensembl gene IDs
clean_ids <- function(x) {
  sub("\\..*$", "", x)
}

sig_ids <- clean_ids(rownames(significant_genes))
up_ids <- clean_ids(rownames(upregulated))
down_ids <- clean_ids(rownames(downregulated))
all_ids <- clean_ids(rownames(dds))

# Convert Ensembl IDs into Entrez IDs and gene symbols
gene_conversion <- bitr(
  unique(all_ids),
  fromType = "ENSEMBL",
  toType = c("ENTREZID", "SYMBOL"),
  OrgDb = org.Hs.eg.db
)

# Prepare gene lists
convert_to_entrez <- function(ids) {
  unique(na.omit(
    gene_conversion$ENTREZID[
      match(ids, gene_conversion$ENSEMBL)
    ]
  ))
}

sig_entrez <- convert_to_entrez(sig_ids)
up_entrez <- convert_to_entrez(up_ids)
down_entrez <- convert_to_entrez(down_ids)
universe_entrez <- convert_to_entrez(all_ids)


############################################################
# 7. GO ENRICHMENT ANALYSIS
############################################################

# All significant genes
ego_all <- enrichGO(
  gene = sig_entrez,
  universe = universe_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

# Upregulated genes
ego_up <- enrichGO(
  gene = up_entrez,
  universe = universe_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

# Downregulated genes
ego_down <- enrichGO(
  gene = down_entrez,
  universe = universe_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

# Export GO results
write.csv(as.data.frame(ego_all),
          "Results/GO_All_DEGs.csv", row.names = FALSE)

write.csv(as.data.frame(ego_up),
          "Results/GO_Upregulated.csv", row.names = FALSE)

write.csv(as.data.frame(ego_down),
          "Results/GO_Downregulated.csv", row.names = FALSE)


############################################################
# 8. KEGG PATHWAY ENRICHMENT ANALYSIS
############################################################

# All significant genes
kegg_all <- enrichKEGG(
  gene = sig_entrez,
  universe = universe_entrez,
  organism = "hsa",
  keyType = "ncbi-geneid",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

# Upregulated genes
kegg_up <- enrichKEGG(
  gene = up_entrez,
  universe = universe_entrez,
  organism = "hsa",
  keyType = "ncbi-geneid",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

# Downregulated genes
kegg_down <- enrichKEGG(
  gene = down_entrez,
  universe = universe_entrez,
  organism = "hsa",
  keyType = "ncbi-geneid",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

# Export KEGG results
write.csv(as.data.frame(kegg_all),
          "Results/KEGG_All_DEGs.csv", row.names = FALSE)

write.csv(as.data.frame(kegg_up),
          "Results/KEGG_Upregulated.csv", row.names = FALSE)

write.csv(as.data.frame(kegg_down),
          "Results/KEGG_Downregulated.csv", row.names = FALSE)


############################################################
# 9. ENRICHMENT VISUALIZATION
############################################################

# GO plots
p1 <- dotplot(
  ego_all,
  showCategory = 15,
  title = "GO Enrichment - All DEGs"
)

p2 <- dotplot(
  ego_up,
  showCategory = 15,
  title = "GO Enrichment - Upregulated Genes"
)

p3 <- dotplot(
  ego_down,
  showCategory = 15,
  title = "GO Enrichment - Downregulated Genes"
)

ggsave("Figures/GO_All_DEGs.png", p1,
       width = 10, height = 7, dpi = 300)

ggsave("Figures/GO_Upregulated.png", p2,
       width = 10, height = 7, dpi = 300)

ggsave("Figures/GO_Downregulated.png", p3,
       width = 10, height = 7, dpi = 300)


# KEGG plots
k1 <- dotplot(
  kegg_all,
  showCategory = 12,
  title = "KEGG Enrichment - All DEGs"
)

k2 <- dotplot(
  kegg_up,
  showCategory = 3,
  title = "KEGG Enrichment - Upregulated Genes"
)

k3 <- dotplot(
  kegg_down,
  showCategory = 15,
  title = "KEGG Enrichment - Downregulated Genes"
)

ggsave("Figures/KEGG_All_DEGs.png", k1,
       width = 10, height = 7, dpi = 300)

ggsave("Figures/KEGG_Upregulated.png", k2,
       width = 10, height = 7, dpi = 300)

ggsave("Figures/KEGG_Downregulated.png", k3,
       width = 10, height = 7, dpi = 300)


############################################################
# 10. TOP 10 ANNOTATED GENES
############################################################

res_df$ensembl_id <- clean_ids(rownames(res_df))

annotation <- gene_conversion

annotated_results <- merge(
  res_df,
  annotation,
  by.x = "ensembl_id",
  by.y = "ENSEMBL",
  all.x = TRUE
)

top10 <- annotated_results[
  !is.na(annotated_results$padj) &
  annotated_results$padj < 0.05 &
  abs(annotated_results$log2FoldChange) > 1,
]

top10 <- top10[order(top10$padj), ]
top10 <- head(top10, 10)

write.csv(
  top10,
  "Results/Top_10_Annotated_Genes.csv",
  row.names = FALSE
)


############################################################
# 11. SAVE SESSION INFORMATION
############################################################

writeLines(
  capture.output(sessionInfo()),
  "Results/SessionInfo.txt"
)


############################################################
# END OF ANALYSIS
############################################################

cat("Analysis completed successfully!\n")
cat("Results saved in the Results folder.\n")
cat("Figures saved in the Figures folder.\n")

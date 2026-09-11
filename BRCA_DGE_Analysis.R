# ============================================================
# TCGA-BRCA Differential Gene Expression Analysis
# Primary Tumor vs Solid Tissue Normal
# ============================================================
#
# Research question:
# Which genes are differentially expressed between breast tumour
# and normal breast tissue in TCGA-BRCA, and what biological
# processes are associated with these changes?
#
# Main tools:
# TCGAbiolinks, DESeq2, EnhancedVolcano, clusterProfiler,
# org.Hs.eg.db, AnnotationDbi, ggplot2, pheatmap, dplyr
#
# ============================================================

# ---------------------------
# 1. Install/load packages
# ---------------------------

cran_packages <- c("ggplot2", "dplyr", "pheatmap")
bioc_packages <- c(
  "TCGAbiolinks", "DESeq2", "EnhancedVolcano",
  "clusterProfiler", "org.Hs.eg.db", "AnnotationDbi"
)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

library(TCGAbiolinks)
library(DESeq2)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ggplot2)
library(dplyr)
library(pheatmap)

# ---------------------------
# 2. Create project folders
# ---------------------------

project_dir <- "TCGA-BRCA-DGE-Analysis"

dir.create(project_dir, showWarnings = FALSE)
dir.create(file.path(project_dir, "figures"), showWarnings = FALSE)
dir.create(file.path(project_dir, "results"), showWarnings = FALSE)

# ---------------------------
# 3. Query TCGA-BRCA RNA-seq
# ---------------------------

query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

# Download the data.
# files.per.chunk = 20 can make a large GDC download more reliable.
GDCdownload(
  query,
  files.per.chunk = 20
)

# Prepare the downloaded data as a SummarizedExperiment object.
data <- GDCprepare(query)

# ---------------------------
# 4. Inspect sample types
# ---------------------------

print(table(data$sample_type))

# Keep only Primary Tumor and Solid Tissue Normal.
# Metastatic samples are excluded.
keep <- data$sample_type %in% c(
  "Primary Tumor",
  "Solid Tissue Normal"
)

data_brca <- data[, keep]

print(table(data_brca$sample_type))

# ---------------------------
# 5. Define experimental condition
# ---------------------------

condition <- factor(
  data_brca$sample_type,
  levels = c(
    "Solid Tissue Normal",
    "Primary Tumor"
  )
)

print(table(condition))

# ---------------------------
# 6. Extract raw unstranded counts
# ---------------------------

counts <- assay(
  data_brca,
  "unstranded"
)

print(dim(counts))

# Create sample metadata.
coldata <- as.data.frame(colData(data_brca))
coldata$condition <- condition

print(table(coldata$condition))

# Verify that sample order is identical.
stopifnot(identical(colnames(counts), rownames(coldata)))

# ---------------------------
# 7. Clean Ensembl gene IDs
# ---------------------------

# Remove version suffixes such as:
# ENSG00000000003.15 -> ENSG00000000003

rownames(counts) <- sub(
  "\\..*$",
  "",
  rownames(counts)
)

# Removing version suffixes can create duplicate Ensembl IDs.
# Combine duplicated rows by summing their counts.

if (anyDuplicated(rownames(counts)) > 0) {
  counts <- rowsum(
    counts,
    group = rownames(counts)
  )
}

print(dim(counts))
print(anyDuplicated(rownames(counts)))

# ---------------------------
# 8. Create DESeq2 object
# ---------------------------

dds <- DESeqDataSetFromMatrix(
  countData = round(counts),
  colData = coldata,
  design = ~ condition
)

# Simple low-count filtering.
dds <- dds[
  rowSums(counts(dds)) >= 10,
]

print(dim(dds))
print(table(dds$condition))

# ---------------------------
# 9. Run differential expression
# ---------------------------

dds <- DESeq(dds)

res <- results(
  dds,
  contrast = c(
    "condition",
    "Primary Tumor",
    "Solid Tissue Normal"
  )
)

res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

# ---------------------------
# 10. Annotate Ensembl IDs
# ---------------------------

gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = res_df$gene,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

res_df$symbol <- gene_symbols[res_df$gene]

# ---------------------------
# 11. Define significant DEGs
# ---------------------------

# Criteria used in the project:
# adjusted p-value < 0.05
# absolute log2 fold change > 1

sig_genes <- res_df[
  !is.na(res_df$padj) &
    res_df$padj < 0.05 &
    abs(res_df$log2FoldChange) > 1,
]

up_genes <- sig_genes[
  sig_genes$log2FoldChange > 1,
]

down_genes <- sig_genes[
  sig_genes$log2FoldChange < -1,
]

# Sort significant results by adjusted p-value.
sig_genes <- sig_genes[
  order(sig_genes$padj),
]

cat("\nTotal significant DEGs:", nrow(sig_genes), "\n")
cat("Upregulated DEGs:", nrow(up_genes), "\n")
cat("Downregulated DEGs:", nrow(down_genes), "\n")

# ---------------------------
# 12. Save differential-expression results
# ---------------------------

write.csv(
  res_df,
  file.path(
    project_dir,
    "results",
    "BRCA_DESeq2_all_results.csv"
  ),
  row.names = FALSE
)

write.csv(
  sig_genes,
  file.path(
    project_dir,
    "results",
    "BRCA_significant_DEGs.csv"
  ),
  row.names = FALSE
)

write.csv(
  up_genes,
  file.path(
    project_dir,
    "results",
    "BRCA_upregulated_genes.csv"
  ),
  row.names = FALSE
)

write.csv(
  down_genes,
  file.path(
    project_dir,
    "results",
    "BRCA_downregulated_genes.csv"
  ),
  row.names = FALSE
)

# ---------------------------
# 13. Volcano plot
# ---------------------------

volcano_plot <- EnhancedVolcano(
  res_df,
  lab = res_df$symbol,
  x = "log2FoldChange",
  y = "padj",
  title = "TCGA-BRCA: Primary Tumor vs Normal",
  subtitle = "Differential Gene Expression Analysis",
  pCutoff = 0.05,
  FCcutoff = 1,
  pointSize = 2.0,
  labSize = 3.0
)

print(volcano_plot)

ggsave(
  file.path(
    project_dir,
    "figures",
    "BRCA_volcano_plot.png"
  ),
  plot = volcano_plot,
  width = 10,
  height = 8,
  dpi = 300
)

# ---------------------------
# 14. Variance-stabilizing transformation
# ---------------------------

vsd <- vst(
  dds,
  blind = FALSE
)

vsd_mat <- assay(vsd)

print(dim(vsd_mat))

# ---------------------------
# 15. PCA plot
# ---------------------------

pca_plot <- plotPCA(
  vsd,
  intgroup = "condition"
)

print(pca_plot)

ggsave(
  file.path(
    project_dir,
    "figures",
    "BRCA_PCA_plot.png"
  ),
  plot = pca_plot,
  width = 9,
  height = 7,
  dpi = 300
)

# ---------------------------
# 16. Sample annotation for heatmap
# ---------------------------

annotation_col <- data.frame(
  Condition = coldata$condition
)

rownames(annotation_col) <- rownames(coldata)

# ---------------------------
# 17. Top 30 DEG heatmap
# ---------------------------

top30 <- head(
  sig_genes[
    !is.na(sig_genes$symbol),
  ],
  30
)

top30_ids <- top30$gene

mat_top30 <- vsd_mat[
  top30_ids,
]

rownames(mat_top30) <- top30$symbol

# Save pheatmap as a high-resolution PNG.
png(
  filename = file.path(
    project_dir,
    "figures",
    "BRCA_heatmap_top30.png"
  ),
  width = 2400,
  height = 1800,
  res = 300
)

pheatmap(
  mat_top30,
  scale = "row",
  show_rownames = TRUE,
  show_colnames = FALSE,
  annotation_col = annotation_col,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  fontsize_row = 8,
  main = "Top 30 Differentially Expressed Genes — TCGA-BRCA"
)

dev.off()

# ---------------------------
# 18. GO enrichment: upregulated genes
# ---------------------------

up_entrez <- bitr(
  up_genes$gene,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

up_entrez <- up_entrez[
  !duplicated(up_entrez$ENTREZID),
]

cat(
  "\nMapped upregulated genes:",
  nrow(up_entrez),
  "\n"
)

ego_up <- enrichGO(
  gene = up_entrez$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

go_up_plot <- dotplot(
  ego_up,
  showCategory = 20,
  title = "GO Biological Process Enrichment — Upregulated Genes"
)

print(go_up_plot)

ggsave(
  file.path(
    project_dir,
    "figures",
    "GO_upregulated.png"
  ),
  plot = go_up_plot,
  width = 10,
  height = 8,
  dpi = 300
)

write.csv(
  as.data.frame(ego_up),
  file.path(
    project_dir,
    "results",
    "GO_upregulated_genes.csv"
  ),
  row.names = FALSE
)

# ---------------------------
# 19. GO enrichment: downregulated genes
# ---------------------------

down_entrez <- bitr(
  down_genes$gene,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

down_entrez <- down_entrez[
  !duplicated(down_entrez$ENTREZID),
]

cat(
  "\nMapped downregulated genes:",
  nrow(down_entrez),
  "\n"
)

ego_down <- enrichGO(
  gene = down_entrez$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

go_down_plot <- dotplot(
  ego_down,
  showCategory = 20,
  title = "GO Biological Process Enrichment — Downregulated Genes"
)

print(go_down_plot)

ggsave(
  file.path(
    project_dir,
    "figures",
    "GO_downregulated.png"
  ),
  plot = go_down_plot,
  width = 10,
  height = 8,
  dpi = 300
)

write.csv(
  as.data.frame(ego_down),
  file.path(
    project_dir,
    "results",
    "GO_downregulated_genes.csv"
  ),
  row.names = FALSE
)

# ---------------------------
# 20. Key breast-cancer genes
# ---------------------------

key_genes <- c(
  "BRCA1",
  "BRCA2",
  "TP53",
  "ESR1",
  "ERBB2"
)

key_results <- res_df[
  res_df$symbol %in% key_genes,
  c(
    "gene",
    "symbol",
    "baseMean",
    "log2FoldChange",
    "pvalue",
    "padj"
  )
]

key_results$significant <- (
  !is.na(key_results$padj) &
    key_results$padj < 0.05 &
    abs(key_results$log2FoldChange) > 1
)

print(key_results)

write.csv(
  key_results,
  file.path(
    project_dir,
    "results",
    "BRCA_key_cancer_genes.csv"
  ),
  row.names = FALSE
)

# ---------------------------
# 21. Key-gene log2FC plot
# ---------------------------

key_gene_plot <- ggplot(
  key_results,
  aes(
    x = symbol,
    y = log2FoldChange
  )
) +
  geom_col() +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  labs(
    title = "Expression Changes of Key Breast Cancer Genes",
    x = "Gene",
    y = "Log2 Fold Change\n(Primary Tumor vs Normal)"
  ) +
  theme_minimal()

print(key_gene_plot)

ggsave(
  file.path(
    project_dir,
    "figures",
    "key_breast_cancer_genes_log2FC.png"
  ),
  plot = key_gene_plot,
  width = 9,
  height = 6,
  dpi = 300
)

# ---------------------------
# 22. Project summary
# ---------------------------

project_summary <- data.frame(
  Category = c(
    "Total samples",
    "Primary Tumor samples",
    "Solid Tissue Normal samples",
    "Genes analyzed after filtering",
    "Significant DEGs",
    "Upregulated DEGs",
    "Downregulated DEGs",
    "PCA PC1 variance",
    "PCA PC2 variance"
  ),
  Result = c(
    ncol(counts),
    sum(data_brca$sample_type == "Primary Tumor"),
    sum(data_brca$sample_type == "Solid Tissue Normal"),
    nrow(dds),
    nrow(sig_genes),
    nrow(up_genes),
    nrow(down_genes),
    "19%",
    "12%"
  )
)

print(project_summary)

write.csv(
  project_summary,
  file.path(
    project_dir,
    "results",
    "BRCA_project_summary.csv"
  ),
  row.names = FALSE
)

# ---------------------------
# 23. Save session information
# ---------------------------

capture.output(
  sessionInfo(),
  file = file.path(
    project_dir,
    "results",
    "sessionInfo.txt"
  )
)

# ---------------------------
# 24. Final message
# ---------------------------

cat("\n============================================================\n")
cat("TCGA-BRCA DGE ANALYSIS COMPLETED\n")
cat("============================================================\n")
cat("Significant DEGs:", nrow(sig_genes), "\n")
cat("Upregulated DEGs:", nrow(up_genes), "\n")
cat("Downregulated DEGs:", nrow(down_genes), "\n")
cat("\nFigures and results saved in:\n")
cat(project_dir, "\n")
cat("============================================================\n")

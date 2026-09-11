# Differential Gene Expression Analysis of Breast Cancer Using TCGA-BRCA

## Overview

This project presents a bioinformatics analysis of RNA-seq gene expression data from The Cancer Genome Atlas Breast Invasive Carcinoma (TCGA-BRCA) project. The analysis compares primary breast tumour samples with solid tissue normal samples to identify differentially expressed genes and investigate the biological processes associated with these changes.

The complete analysis was performed using R and Bioconductor packages, with visualisation and functional enrichment analysis used to interpret the results.

## Research Question

Which genes are differentially expressed between breast tumour and normal breast tissue in TCGA-BRCA, and what biological processes are associated with these changes?

## Dataset

- **Cancer type:** Breast Invasive Carcinoma (TCGA-BRCA)
- **Data type:** RNA-seq Gene Expression Quantification
- **Workflow:** STAR - Counts
- **Primary Tumour samples:** 1,111
- **Solid Tissue Normal samples:** 113
- **Total samples analysed:** 1,224

Raw TCGA data are not included in this repository because of their large size.

## Tools and Packages

- R
- TCGAbiolinks
- DESeq2
- EnhancedVolcano
- clusterProfiler
- org.Hs.eg.db
- ggplot2
- dplyr
- pheatmap

## Analysis Workflow

1. Downloaded TCGA-BRCA RNA-seq data using TCGAbiolinks.
2. Prepared the gene expression count matrix.
3. Selected Primary Tumour and Solid Tissue Normal samples.
4. Removed Ensembl version suffixes and handled duplicated gene identifiers.
5. Applied low-count filtering.
6. Performed differential gene expression analysis using DESeq2.
7. Identified significant differentially expressed genes.
8. Generated a volcano plot.
9. Performed variance-stabilizing transformation and PCA.
10. Generated a heatmap of the top differentially expressed genes.
11. Performed Gene Ontology Biological Process enrichment analysis using clusterProfiler.
12. Examined expression changes in selected breast cancer-related genes.

## Differential Expression Results

Using the criteria:

- Adjusted p-value (padj) < 0.05
- Absolute log2 fold change > 1

the analysis identified:

- **12,304 significant differentially expressed genes**
- **8,186 upregulated genes**
- **4,118 downregulated genes**

Positive log2 fold-change values represent higher expression in Primary Tumour samples relative to Solid Tissue Normal samples.

## Key Findings

Several genes showed substantial expression differences between tumour and normal breast tissue. Examples include:

- MMP11
- COL10A1
- COL11A1
- NEK2
- KIF4A
- PKMYT1
- IQGAP3
- FHL1
- CAV1
- LYVE1

Selected breast cancer-related genes were also examined, including **BRCA1, BRCA2, TP53, ESR1, and ERBB2**.

BRCA1, BRCA2, ESR1, and ERBB2 met the project's differential-expression criteria, while TP53 showed statistically significant expression change but did not meet the absolute log2 fold-change threshold of 1.

These results represent **gene expression differences** between tumour and normal samples and do not represent mutation analysis.

## Principal Component Analysis

PCA was performed using variance-stabilized expression data.

- **PC1 explained 19% of the variance**
- **PC2 explained 12% of the variance**

The PCA showed a clear overall difference in gene expression profiles between Primary Tumour and Solid Tissue Normal samples.

## Functional Enrichment Analysis

Gene Ontology Biological Process enrichment analysis was performed separately for upregulated and downregulated genes.

### Upregulated Genes

Enriched biological processes included processes related to:

- Chromosome and centromere organization
- Nucleosome assembly and organization
- Nuclear division
- Chromosome segregation
- Cell-cycle checkpoint signalling
- Mitotic sister chromatid segregation
- Humoral immune response

### Downregulated Genes

Enriched biological processes included processes related to:

- Blood circulation
- Muscle system processes
- Response to wounding
- Chemotaxis
- Muscle contraction
- Actin filament-based processes
- Regulation of membrane potential
- Vasculature development
- Regulation of angiogenesis
- Blood pressure regulation

Because this is a bulk RNA-seq analysis, some enriched processes may reflect differences in tissue composition and stromal or vascular components rather than tumour-cell-intrinsic changes.

## Visualizations

The `figures/` directory contains:

- PCA plot
- Volcano plot
- Top 30 DEG heatmap
- GO enrichment plot for upregulated genes
- GO enrichment plot for downregulated genes
- Key breast cancer gene expression plot

## Repository Contents

- `BRCA_DGE_Analysis.R` — complete R analysis script
- `BRCA_DGE_Report_Final.pdf` — final project report
- `figures/` — generated visualizations
- `BRCA_DESeq2_all_results.csv` — complete DESeq2 results
- `BRCA_significant_DEGs.csv` — significant differentially expressed genes
- `BRCA_upregulated_genes.csv` — upregulated genes
- `BRCA_downregulated_genes.csv` — downregulated genes
- `BRCA_key_cancer_genes.csv` — selected cancer-related genes
- `BRCA_project_summary.csv` — project summary statistics
- `GO_upregulated_genes.csv` — GO enrichment results for upregulated genes
- `GO_downregulated_genes.csv` — GO enrichment results for downregulated genes

## Reproducibility

The complete R code used for data preparation, differential expression analysis, visualization, gene annotation, and functional enrichment is provided in:

`BRCA_DGE_Analysis.R`

The analysis can be reproduced using the R/Bioconductor packages listed above and the corresponding TCGA-BRCA RNA-seq data.

## Limitations

- This analysis uses bulk RNA-seq data, which contains signals from multiple cell types.
- The tumour and normal groups are unbalanced in sample number.
- The analysis compares tumour and normal samples using an unpaired design.
- Mutation status of BRCA1 or BRCA2 was not analysed.
- The results represent associations in gene expression and do not establish causal relationships.
- Additional analyses such as survival analysis and protein-protein interaction analysis could further extend the project.

## Conclusion

This project demonstrates an end-to-end RNA-seq differential gene expression workflow using TCGA-BRCA data. The analysis identified thousands of genes with significant expression differences between breast tumour and normal tissue and revealed biological processes associated with these changes.

The project demonstrates practical skills in R programming, transcriptomic data analysis, differential expression analysis, data visualization, gene annotation, and functional enrichment analysis.

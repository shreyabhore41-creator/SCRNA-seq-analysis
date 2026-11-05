
# Single-Cell RNA-seq Analysis of Endometriosis

## Project Overview

This mini-project focuses on the computational analysis of endometriosis using publicly available single-cell RNA sequencing (scRNA-seq) data from the NCBI Gene Expression Omnibus (GEO), accession GSE203191.

The study explores cellular heterogeneity, cell-type composition, differential gene expression and biological pathway enrichment to understand cell-specific transcriptional patterns associated with endometriosis.

## Objectives

- Identify major cell populations using single-cell transcriptomic data.
- Perform quality control and doublet detection.
- Integrate sequencing libraries and identify cell clusters.
- Annotate cell populations using marker genes.
- Compare cell-type composition across clinical and preparation groups.
- Investigate ambient RNA contamination using DecontX.
- Identify candidate differentially expressed genes.
- Explore enriched biological pathways using GSEA.

## Tools and Technologies

- R programming
- Seurat
- Harmony
- scDblFinder
- DecontX
- DESeq2
- fgsea
- speckle

## Workflow

1. Dataset retrieval from GEO.
2. Quality control and filtering.
3. Doublet detection.
4. Normalization and PCA.
5. Harmony integration and UMAP visualization.
6. Clustering and cell-type annotation.
7. Cell composition analysis.
8. Ambient RNA correction.
9. Differential expression analysis.
10. Hallmark pathway enrichment.

## Key Findings

- Identified distinct stromal, epithelial, endothelial and immune-associated cell populations.
- Sample preparation influenced perivascular and endothelial cell proportions.
- No cell-type proportion differences between diagnosed cases and controls were significant after FDR correction.
- DecontX correction reduced significant differential-expression genes from 120 to 28.
- Identified candidate genes including PARD6G-AS1, SNHG22, PASK, CXCL5 and KLRC1.
- Identified immune-related and tissue-remodelling pathway signatures.

## Data Source

NCBI GEO: [GSE203191](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE203191\&utm_source=chatgpt.com)

## Limitations

The dataset consists of pooled, non-demultiplexed libraries. Therefore, patient-level differential-expression analysis was not possible. Findings are exploratory and require independent validation.

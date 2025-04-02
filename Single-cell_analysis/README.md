# Single-cell RNA-seq Analysis

This folder contains scripts for processing and analyzing single-cell RNA sequencing data from the retina of chick and zebrafish species.

## Scripts

- **04.0-create_seurat_objects_chick.R**  
  Loads and filters chick E12, E16, and E18 count matrices, then creates Seurat objects and performs basic QC.

- **04.1-scRNAseq_chick_retina_analysis.R**  
  Performs clustering, UMAP dimensionality reduction, marker gene analysis, and cell type annotation for the chick retina.

- **04.2-scRNAseq_zebrafish_analysis.R**  
  Conducts clustering and cell type identification for zebrafish retina cells, focusing on opsin gene expression.

## Dependencies

- [Seurat](https://satijalab.org/seurat/) (v4+)
- dplyr
- Matrix
- ggplot2

All scripts were run in R using Seurat pipelines for normalization, clustering, and visualization.

## Output

Each script generates:
- Annotated Seurat objects
- Cluster UMAP plots
- DotPlots/FeaturePlots for opsin and retinal cell markers

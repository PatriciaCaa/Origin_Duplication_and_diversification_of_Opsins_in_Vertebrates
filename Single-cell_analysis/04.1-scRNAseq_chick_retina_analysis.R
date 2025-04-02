## Chick Retina scRNA-seq Analysis

# Description: Seurat pipeline for E12, E16, and E18 chick retina datasets
# Dataset source: GEO (Gallus gallus - GSE159107)

# Set custom R library path (optional, for HPC environments)
.libPaths("/home/v/vpc5/Single_cell/Rlibrary")

library(dplyr)
library(Seurat)
library(Matrix)

# Load pre-saved Seurat objects
seuratObject12 <- readRDS("/home/v/vpc5/Single_cell/seuratObject12.rds")
seuratObject16 <- readRDS("/home/v/vpc5/Single_cell/seuratObject16.rds")
seuratObject18 <- readRDS("/home/v/vpc5/Single_cell/seuratObject18.rds")

# Merge into one object
combinedSeuratObject <- merge(seuratObject12, y = c(seuratObject16, seuratObject18), add.cell.ids = c("E12", "E16", "E18"))

# Normalize and cluster
combinedSeuratObject <- FindVariableFeatures(combinedSeuratObject, selection.method = "vst", nfeatures = 2000)
combinedSeuratObject <- ScaleData(combinedSeuratObject, features = rownames(combinedSeuratObject))
combinedSeuratObject <- RunPCA(combinedSeuratObject, features = VariableFeatures(object = combinedSeuratObject))
combinedSeuratObject <- FindNeighbors(combinedSeuratObject, dims = 1:10)
combinedSeuratObject <- FindClusters(combinedSeuratObject, resolution = 0.5)
combinedSeuratObject <- RunUMAP(combinedSeuratObject, dims = 1:10)

# Save combined object
saveRDS(combinedSeuratObject, file = "/home/v/vpc5/Single_cell/combinedSeuratObjectchick.rds")

# Initial visualizations
DimPlot(combinedSeuratObject, label = TRUE)
DimPlot(combinedSeuratObject, group.by = "orig.file", cells = sample(colnames(combinedSeuratObject)))
VlnPlot(combinedSeuratObject, "nCount_RNA", pt.size = 0) + RotatedAxis()
VlnPlot(combinedSeuratObject, "nFeature_RNA", pt.size = 0) + RotatedAxis()

# Differential expression
combinedSeuratObject_genes <- FindAllMarkers(combinedSeuratObject)
saveRDS(combinedSeuratObject_genes, file = "/home/v/vpc5/Single_cell/combinedSeuratObjectchick_genes.rds")

# Cluster annotation using marker genes
Idents(combinedSeuratObject) <- "seurat_clusters"
combinedSeuratObject <- DendroOrder(combinedSeuratObject)

RGC_markers <- c("RBPMS", "POU4F1", "THY1")
BC_markers <- c("VSX1", "VSX2")
AC_markers <- c("PAX6", "SLC32A1", "GAD1", "GAD2", "SLC6A9")
MG_markers <- c("SLC1A3", "RLBP1")
HC_markers <- c("ONECUT1", "ONECUT2", "ONECUT3")
Cone_markers <- c("ARR3", "RS1", "GNGT2")
Rod_markers <- c("RHO", "PDE6B", "PDC")
OL_markers <- c("OLIG2")

DotPlot(combinedSeuratObject, features = c(RGC_markers, BC_markers, AC_markers, HC_markers, Cone_markers, Rod_markers, MG_markers, OL_markers)) + RotatedAxis()

# Annotate cell classes manually based on cluster IDs
combinedSeuratObject@meta.data$cell_class <- "N/A"
combinedSeuratObject@meta.data[WhichCells(combinedSeuratObject, idents = c(0, 3, 5, 7, 12, 29, 45)),]$cell_class <- "RGC"
combinedSeuratObject@meta.data[WhichCells(combinedSeuratObject, idents = c(9, 13, 14, 16, 21:27, 31, 33, 35, 37, 40, 41)),]$cell_class <- "BC"
combinedSeuratObject@meta.data[WhichCells(combinedSeuratObject, idents = c(6, 8, 18, 20, 28, 30, 32, 34, 38, 39, 43)),]$cell_class <- "AC"
combinedSeuratObject@meta.data[WhichCells(combinedSeuratObject, idents = c(11, 17)),]$cell_class <- "HC"
combinedSeuratObject@meta.data[WhichCells(combinedSeuratObject, idents = c(4, 10, 19, 36)),]$cell_class <- "PR"
combinedSeuratObject@meta.data[WhichCells(combinedSeuratObject, idents = c(2, 44)),]$cell_class <- "MG"
combinedSeuratObject@meta.data[WhichCells(combinedSeuratObject, idents = c(42)),]$cell_class <- "OL"

# Remove clusters 1 and 15 (low quality)
idents_remove <- c(1, 15)
combinedSeuratObject <- subset(combinedSeuratObject, idents = setdiff(levels(Idents(combinedSeuratObject)), idents_remove))

# Final UMAP with annotations
DimPlot(combinedSeuratObject, group.by = "cell_class")
DimPlot(combinedSeuratObject, label = TRUE)

saveRDS(combinedSeuratObject, "/home/v/vpc5/Single_cell/combinedSeuratObjectchick_annotated.rds")

# Barplot of cell class composition per sample
counts <- table(combinedSeuratObject@meta.data$cell_class, combinedSeuratObject@meta.data$orig.file)
counts <- t(t(counts) / colSums(counts))
barplot(counts, legend = rownames(counts), col= rev(c("red","orange","yellow","white","green","blue", "black")))



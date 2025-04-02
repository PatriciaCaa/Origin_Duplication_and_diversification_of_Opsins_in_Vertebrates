# Zebrafish Retina scRNA-seq Analysis
# Project: Origin, duplication, and diversification of opnsin genes in vertebrates
# Description: Full pipeline for zebrafish retinal scRNA-seq dataset using Seurat

.libPaths("/home/v/vpc5/Single_cell/Rlibrary")

library(dplyr)
library(Seurat)
library(ggplot2)

# Load preprocessed Seurat object
load("/home/v/vpc5/Single_cell/scell-scripts/new_analyses/sub_sam-zebrafishObj.rdata")

# Subsample to 20,000 cells
zebrafishObj <- subset(zebrafishObj, cells = sample(Cells(zebrafishObj), 20000))
save(zebrafishObj, file = "sub_sam-zebrafishObj.rdata")

# Quality Control
zebrafishObj[["percent.mt"]] <- PercentageFeatureSet(zebrafishObj, pattern = "^MT-")
VlnPlot(zebrafishObj, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), ncol = 3)
FeatureScatter(zebrafishObj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# Filter low-quality cells
zebrafishObj <- subset(zebrafishObj, subset = nCount_RNA > 500 & nCount_RNA < 2500 &
                                         nFeature_RNA > 200 & nFeature_RNA < 4000 &
                                         percent.mt < 5)

# Post-filter visualization
VlnPlot(zebrafishObj, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), ncol = 3)

# Normalize and scale
zebrafishObj <- NormalizeData(zebrafishObj)
zebrafishObj <- FindVariableFeatures(zebrafishObj, selection.method = "vst", nfeatures = 2000)
plot1 <- VariableFeaturePlot(zebrafishObj)
print(plot1)

zebrafishObj <- ScaleData(zebrafishObj)
zebrafishObj <- RunPCA(zebrafishObj, features = VariableFeatures(object = zebrafishObj))

# Visualize PCA
DimHeatmap(zebrafishObj, dims = 1:6, cells = 500, balanced = TRUE)
ElbowPlot(zebrafishObj)

# Clustering and UMAP
zebrafishObj <- FindNeighbors(zebrafishObj, dims = 1:15)
zebrafishObj <- FindClusters(zebrafishObj, resolution = 0.5)
zebrafishObj <- RunUMAP(zebrafishObj, dims = 1:10)

umap_plot <- DimPlot(zebrafishObj, reduction = "umap", label = TRUE)
print(umap_plot)

# Save for future use
save(zebrafishObj, file = "zebrafishObj.rdata")

# Marker gene expression visualization
Idents(zebrafishObj) <- "seurat_clusters"
levels(Idents(zebrafishObj)) <- c(0:13)

RGC_markers = c("RBPMS2B", "POU6F2", "THY1")
BC_markers = c("VSX2", "PRKCA", "OTX2", "VSX1")
AC_markers = c("TFAP2A", "GAD2", "SLC6A9")
HC_markers = c("ONECUT1", "CALB1", "TPM3")
Cone_markers = c("ARR3A", "PDE6H")
Rod_markers = c("RHO")
MG_markers = c("GLULA", "APOEB")
MicroG_markers = c("C1QA", "C1QB")

DotPlot(zebrafishObj, features = c(RGC_markers, BC_markers, AC_markers, HC_markers, Cone_markers, Rod_markers, MG_markers, MicroG_markers, "CLDN5A", "CLDN5B", "IGFBP7")) + RotatedAxis()

# Find markers per cluster
cluster_markers <- FindAllMarkers(zebrafishObj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Explore top markers from cluster 1
cluster1_markers <- cluster_markers %>% filter(cluster == 1) %>% top_n(10, avg_log2FC)
cluster1_gene_names <- cluster1_markers$gene
FeaturePlot(zebrafishObj, features = cluster1_gene_names)

# Annotate clusters
X <- c("RGC1", "RGC2", "RGC3", "RGC4", "RGC5", "MG", "RGC6", "RGC7", "BC", "cone/MG", "cone/MG", "RGC8", "RGC9", "CONE/Rod")
names(X) <- c(0:13)
zebrafishObj2 <- RenameIdents(zebrafishObj, X)
zebrafishObj2$CellType <- Idents(zebrafishObj2)
save(zebrafishObj2, file = "zebrafishObj2.rdata")

# Annotated UMAP
umap_plot_annotated <- DimPlot(zebrafishObj2, reduction = "umap", label = TRUE, pt.size = 0.5) +
  ggtitle("UMAP of zebrafish marker genes") + theme(plot.title = element_text(hjust = 0.5))
print(umap_plot_annotated)

# Opsin gene detection
zebrafish_opsin <- read.csv("/home/v/vpc5/Interpro-hit-sequences/ActinopteriCypriniformesDaniorerio.faa.acc-Hitseq.7TMDseqs.names", header = FALSE)
zebrafish_opsin$NAMES <- gsub("^.*_", "", zebrafish_opsin$V1)
zebrafish_opsin$NAMES <- gsub("\\.1", "", zebrafish_opsin$NAMES)

genes_to_check <- c("opn7c", "tmtops3a", "opn1sw1", "rhol", "opn4xb", "tmtops2b",
                    "opn3", "tmtopsb", "tmtops3b", "opn7d", "rrh", "tmtops2a",
                    "valopb", "opn1mw1", "opn1mw2", "opn1mw3", "opn1mw4", "parapinopsinb",
                    "opn1sw2", "opn1lw1", "opn1lw2", "rho", "opn4xa", "opn7b",
                    "parapinopsina", "valopa", "opn7a", "opn4a", "opn4", "opn4b")

genes_to_check <- toupper(genes_to_check)
X <- rownames(zebrafishObj2)
genes_present <- genes_to_check %in% X
result_df <- data.frame(Gene = genes_to_check, IsPresent = genes_present)
print(result_df)

matching_genes <- genes_to_check[genes_present]
if (length(matching_genes) > 0) {
  DotPlot(zebrafishObj2, features = matching_genes) + RotatedAxis()
  FeaturePlot(zebrafishObj2, features = matching_genes)
} else {
  print("No matching genes found for plotting.")
}

DotPlot(zebrafishObj2, features = unique(zebrafish_opsin$NAMES)) + RotatedAxis()
FeaturePlot(zebrafishObj2, features = unique(zebrafish_opsin$NAMES))


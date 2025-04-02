### Chick Retina scRNA-seq: Seurat Object Creation
## Project: Origin, duplication, and diversification of opnsin genes in vertebrates

# Description: Load raw chick retinal data (E12, E16, E18), apply basic filtering, and save Seurat objects as .rds

.libPaths("/home/v/vpc5/Single_cell/Rlibrary")

library(dplyr)
library(Seurat)
library(Matrix)

# Load raw count matrices
scellchicken12 <- read.csv("/home/v/vpc5/Single_cell/sc_seq_papers/Gallus-gallus/GSE159107/Gallus_E12chick_count.matrix.csv", row.names = 1)
scellchicken16 <- read.csv("/home/v/vpc5/Single_cell/sc_seq_papers/Gallus-gallus/GSE159107/Gallus_E16chick_count.matrix.csv", row.names = 1)
scellchicken18 <- read.csv("/home/v/vpc5/Single_cell/sc_seq_papers/Gallus-gallus/GSE159107/Gallus_E18chick_count.matrix.csv", row.names = 1)

# Convert to sparse matrices
scellchicken12 <- as(as.matrix(scellchicken12), "sparseMatrix")
scellchicken16 <- as(as.matrix(scellchicken16), "sparseMatrix")
scellchicken18 <- as(as.matrix(scellchicken18), "sparseMatrix")

# Create Seurat object for E12 and subsample to 10,000 cells
seuratObject12 <- CreateSeuratObject(counts = scellchicken12, min.cells = 3, min.features = 200)
seuratObject12 <- subset(seuratObject12, cells = sample(Cells(seuratObject12), 10000))

# Create Seurat objects for E16 and E18
seuratObject16 <- CreateSeuratObject(counts = scellchicken16, min.cells = 3, min.features = 200)
seuratObject18 <- CreateSeuratObject(counts = scellchicken18, min.cells = 3, min.features = 200)

# Add mitochondrial content (percent.mt) to metadata
seuratObject12[["percent.mt"]] <- PercentageFeatureSet(seuratObject12, pattern = "^MT-")
seuratObject16[["percent.mt"]] <- PercentageFeatureSet(seuratObject16, pattern = "^MT-")
seuratObject18[["percent.mt"]] <- PercentageFeatureSet(seuratObject18, pattern = "^MT-")

# Filter cells with high mitochondrial content (>5%)
seuratObject12 <- subset(seuratObject12, subset = percent.mt < 5)
seuratObject16 <- subset(seuratObject16, subset = percent.mt < 5)
seuratObject18 <- subset(seuratObject18, subset = percent.mt < 5)

# Save objects for later use
saveRDS(seuratObject12, file = "/home/v/vpc5/Single_cell/seuratObject12.rds")
saveRDS(seuratObject16, file = "/home/v/vpc5/Single_cell/seuratObject16.rds")
saveRDS(seuratObject18, file = "/home/v/vpc5/Single_cell/seuratObject18.rds")


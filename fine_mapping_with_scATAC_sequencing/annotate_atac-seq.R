# Script information ------------------------------------------------------

# title: Annotate ATAC-seq data according to OneK1K cell types
# author: Jose Alquicira Hernandez
# date: 2021-03-08
# description: None

# Import libraries --------------------------------------------------------

# Primary
library("tidyverse")
library("dsLib")

# BiocManager::install(c("ggbio", "biovizBase", 
#                        "GenomeInfoDb", "EnsDb.Hsapiens.v75"))
# 
# install.packages("Signac")


# Secondary
library("Seurat")
library("Signac")
library("GenomeInfoDb")
library("EnsDb.Hsapiens.v75")
library("ggplot2")
library("patchwork")


set.seed(1234)

# Set output --------------------------------------------------------------

output <- set_output("2021-03-08", "atac_anno")

# Read data ---------------------------------------------------------------

atac <- readRDS("results/2021-03-13_atac_nextgem/atac.RDS")
pool1 <- readRDS("data/pools/pool_1.RDS")
pool54 <- readRDS("data/pools/pool_54.RDS")


data <- list(pool1, pool54)


# Integrate OneK1K pools --------------------------------------------------

# normalize and identify variable features for each dataset independently
data <- lapply(data, function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(data)
data <- lapply(data, function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

# Find integration anchors
anchors <- FindIntegrationAnchors(data, 
                                  anchor.features = features, 
                                  reduction = "rpca")

# Integrate data
combined <- IntegrateData(anchors)


# specify that we will perform downstream analysis on the corrected data note that the original
# unmodified data still resides in the 'RNA' assay
DefaultAssay(combined) <- "integrated"

# Run the standard workflow for dimensionality reduction
combined <- ScaleData(combined)
combined <- RunPCA(combined, npcs = 30)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:30)

# Transfer labels from OneK1K data to ATAC-seq data -----------------------

transfer.anchors <- FindTransferAnchors(
  reference = combined, 
  query = atac,
  reduction = 'cca', 
  reference.assay = "integrated",
  query.assay = "RNA"
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = combined$cell_type %>% fct_drop(),
  weight.reduction = atac[['lsi']],
  dims = 2:30
)

atac <- AddMetaData(object = atac, metadata = predicted.labels)

DimPlot(atac, group.by = "predicted.id", label = TRUE)

table(atac$predicted.id)

# Save data ---------------------------------------------------------------

saveRDS(atac, here(output, "atac.RDS"))

# Session info ------------------------------------------------------------

print_session(here(output))



#   ____________________________________________________________________________
#   Script information                                                      ####

# title: Subset B cells from OneK1K cohort
# author: Jose Alquicira Hernandez
# date: 2021-03-31
# description: Subsets naive/immature and memory B cells

#   ____________________________________________________________________________
#   HPC details                                                             ####

# screen -S subset
# qrsh -N subset -l mem_requested=150G -q short.q
# conda activate r-4.0.3

#   ____________________________________________________________________________
#   Import libraries                                                        ####

# Primary
library("tidyverse")
library("dsLib")

# Secondary
library("Seurat")

#   ____________________________________________________________________________
#   Set output                                                              ####

output <- set_output("2021-03-31", "b_cells_subset")

#   ____________________________________________________________________________
#   Import data                                                             ####

inicio("Read data")
data <- readRDS(here("results", "2019-08-13_cell_type_anno", "cell_type.RDS"))
fin()

data@assays$SCT@misc$vst.out <- NULL

#   ____________________________________________________________________________
#   Subset data                                                             ####


i <- data$cell_type %in% c("TCL1A- FCER2- B cell", "TCL1A+ FCER2+ B cell")
b_cells <- data[, i]

rm(data)

#   ____________________________________________________________________________
#   Pre-process data                                                        ####

b_cells <- SCTransform(b_cells, 
                       vars.to.regress = c("pool", "percent.mt"), 
                       conserve.memory = TRUE)

b_cells <- RunPCA(b_cells)
b_cells <- RunUMAP(b_cells, dims = 1:30)

pca <- DimPlot(b_cells, reduction = "pca", group.by = "cell_type")
umap <- DimPlot(b_cells, reduction = "umap", group.by = "cell_type")

ggsave(here(output, "pca.png"), pca, width = 7, height = 5)
ggsave(here(output, "umap.png"), umap, width = 7, height = 5)

pca_pool <- DimPlot(b_cells, reduction = "pca", group.by = "pool") + no_legend()
umap_pool <- DimPlot(b_cells, reduction = "umap", group.by = "pool") + no_legend()

ggsave(here(output, "pca_pool.png"), pca_pool, width = 7, height = 5)
ggsave(here(output, "umap_pool.png"), umap_pool, width = 7, height = 5)

#   ____________________________________________________________________________
#   Export data                                                             ####

saveRDS(b_cells, here(output, "b_cells.RDS"))

#   ____________________________________________________________________________
#   Session info                                                            ####

print_session(here(output))



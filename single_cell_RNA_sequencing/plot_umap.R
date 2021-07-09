# Script information ------------------------------------------------------

# title: Generate plots for publication
# author: Jose Alquicira Hernandez
# date: 2020/10/01
# description: Creates UMAP plot from log-normalized data and QC metrics

# screen -S noind
# qrsh -N noind -l mem_requested=200G,tmp_requested=200G -pe smp 3
# conda activate r-4.0.2 

# Import libraries --------------------------------------------------------

# Primary
library("tidyverse")
library("dsLib")

# Secondary
library("Seurat")
library("Nebulosa")
library("colorspace")

# Set output --------------------------------------------------------------

output <- set_output("2020-10-01", "final_plots_noind")

# Helper functions --------------------------------------------------------

plot_feature <- function(data, x, id = "", neb = TRUE, width = 6, height = 5,
                         device = "pdf", ...){
  if(neb){
    p <- plot_density(data, x, size = 0.5, ...) + 
      scale_color_continuous_sequential(palette = "Sunset")
    ggsave(here(output, x %p% "_" %p% id  %p% "_neb." %p%  device), p, 
           width = width, height = height)
  }else{
    p <- FeaturePlot(data, x, ...) + scale_color_viridis_c()
    ggsave(here(output, x %p% "_" %p% id  %p% "." %p%  device), p,
           width = width, height = height)
  }
}


# Read data ---------------------------------------------------------------

# Input
path <- file.path("results", "2019-08-13_cell_type_anno")
filename <- "cell_type.RDS"
input <-  here(path, filename)

# Read file
inicio("Reading gene expression data")
data <- readRDS(input)
fin()

# 1,272,489


# Remove vst output
data@assays$SCT@misc$vst.out <- NULL

# Remove biological outliers ----------------------------------------------
inicio("Detect biological outliers")
cells <- !data$cell_type %in% c("Platelets", "Erythrocytes")
data <- data[, cells]
fin()

# 1,267,768




# Future settings ---------------------------------------------------------

options(future.globals.maxSize = 200 * 1024^3)
future::plan("multisession", workers = 3)

# Normalize data using log transformation ---------------------------------

DefaultAssay(data) <- "RNA"
data@assays$SCT <- NULL

data <- NormalizeData(data)
data <- FindVariableFeatures(data)
data <- ScaleData(data, vars.to.regress = c("percent.mt", "pool"))
data <- RunPCA(data)
p <- ElbowPlot(data, ndims = 25)
ggsave(here(output, "elbow_rna.png"), p, width = 8, height = 3.8)
data <- RunUMAP(data, dims = 1:15)

pal <- c("#882E72", "#B178A6", "#D6C1DE", "#1965B0", "#5289C7", "#7BAFDE", 
         "#4EB265", "#90C987", "#CAE0AB", "#F7EE55", "#F6C141", "#F1932D", 
         "#E8601C", "#DC050C")


relabel_cells <- function(x){
  switch(x, 
         "CD4+ KLRB1- T cell" = "CD4 Naive&CM",
         "CD4+ KLRB1+ T cell" = "CD4 EM&TEMRA",
         "CD4+ SOX4+ T cell" = "CD4 SOX4",
         "CD8+ GNLY+ NKG7+ T cell" = "CD8 EM&TEMRA",
         "CD8+ LTB+ T cell" = "CD8 Naive&CM",
         "CD8+ S100B+ T cell" = "CD8 S100B",
         "XCL1- NK" = "NK",
         "XCL1+ NK" = "NK recruting",
         "IgJ+ B cell" = "Plasma",
         "TCL1A- FCER2- B cell" = "B Mem",
         "TCL1A+ FCER2+ B cell" = "B Imm&Naive",
         "Monocyte CD14+" = "Mono C",
         "Monocyte FCGR3A+" = "Mono NC",
         "Dendritic cell" = "DC",
         "B cell" = "B cell",
         "CD4+ T cell" = "CD4 T cell",
         "CD8+ T cell" = "CD8 T cell",
         "NK cell" = "NK cell",
         "T cell" = "T cell", 
         "pseudo_bulk" = "Pseudo-bulk",
         x
  )
}

cell_type <- sapply(as.character(data$cell_type), relabel_cells)

level_order <- c("CD4 Naive&CM",
                 "CD4 EM&TEMRA",
                 "CD4 SOX4",
                 "CD8 EM&TEMRA",
                 "CD8 Naive&CM",
                 "CD8 S100B",
                 "NK",
                 "NK recruting",
                 "Plasma",
                 "B Mem",
                 "B Imm&Naive",
                 "Mono C",
                 "Mono NC",
                 "DC")

data$cell_type1 <- factor(cell_type, levels = level_order)
table(data$cell_type1)


p <- DimPlot(data, group.by = "cell_type1", cols = pal) + xlab("UMAP 1") + ylab("UMAP 2")
ggsave(here(output, "umap.png"), p, width = 6.7, height = 4.6)


data_10 <- RunUMAP(data, dims = 1:10)
p2 <- DimPlot(data_10, group.by = "cell_type1", cols = pal) + xlab("UMAP 1") + ylab("UMAP 2")
ggsave(here(output, "umap_10pcs.png"), p2, width = 6.7, height = 4.6)


data_10 <- RunUMAP(data, dims = 1:20)
p3 <- DimPlot(data_10, group.by = "cell_type1", cols = pal) + xlab("UMAP 1") + ylab("UMAP 2")
ggsave(here(output, "umap_20pcs.png"), p3, width = 6.7, height = 4.6)

saveRDS(data_10, here(output, "umap.RDS"))


data_10 <- readRDS(here(output, "umap.RDS"))


plot_feature(data_10, "CD14", id = "all", reduction = "umap", neb = TRUE, device = "png")
plot_feature(data_10, "CD4", id = "all", reduction = "umap", neb = TRUE, device = "png")
plot_feature(data_10, "CD8A", id = "all", reduction = "umap", neb = TRUE, device = "png")
plot_feature(data_10, "GZMK", id = "all", reduction = "umap", neb = TRUE, device = "png")
plot_feature(data_10, "CCR7", id = "all", reduction = "umap", neb = TRUE, device = "png")
plot_feature(data_10, "S100B", id = "all", reduction = "umap", neb = TRUE, device = "png")
plot_feature(data_10, "MS4A1", id = "all", reduction = "umap", neb = TRUE, device = "png")
plot_feature(data_10, "KLRB1", id = "all", reduction = "umap", neb = TRUE, device = "png")
plot_feature(data_10, "XCL1", id = "all", reduction = "umap", neb = TRUE, device = "png")


plot_feature(data_10, "XCL2", id = "all", reduction = "umap", neb = TRUE, device = "png")
plot_feature(data_10, "CD3D", id = "all", reduction = "umap", neb = TRUE, device = "png")
plot_feature(data_10, "NKG7", id = "all", reduction = "umap", neb = TRUE, device = "png")
plot_feature(data_10, "GNLY", id = "all", reduction = "umap", neb = TRUE, device = "png")
plot_feature(data_10, "FCGR3A", id = "all", reduction = "umap", neb = TRUE, device = "png")
plot_feature(data_10, "IGJ", id = "all", reduction = "umap", neb = TRUE, device = "png")
plot_feature(data_10, "FCER1A", id = "all", reduction = "umap", neb = TRUE, device = "png")
plot_feature(data_10, "SERPINF1", id = "all", reduction = "umap", neb = TRUE, device = "png")
plot_feature(data_10, "TCL1A", id = "all", reduction = "umap", neb = TRUE, device = "png")
plot_feature(data_10, "CD27", id = "all", reduction = "umap", neb = TRUE, device = "png")
plot_feature(data_10, "SELL", id = "all", reduction = "umap", neb = TRUE, device = "png")
plot_feature(data_10, "LTB", id = "all", reduction = "umap", neb = TRUE, device = "png")
plot_feature(data_10, "GZMB", id = "all", reduction = "umap", neb = TRUE, device = "png")
plot_feature(data_10, "TNFSF13B", id = "all", reduction = "umap", neb = TRUE, device = "png")



get_max <- function(gene){
  p <- plot_density(data_10, gene)
  max(p$data$feature)
}

genes <- c("CD14", "CD4", "CD8A", "GZMK", "CCR7", "S100B", "MS4A1", "KLRB1", "XCL1")
res <- map(genes, get_max)

data.frame(genes, max = flatten_dbl(res) %>% round(2)) %>% 
  knitr::kable()





# Additional marker plots -------------------------------------------------

markers <- c(
"LTB",
"PASK",
"IL7R",
"CD3D",
"XCL2",
"NKG7",
"GNLY",
"IGJ",
"SERPINF1",
"CST3",
"FCER1A",
"TCL1A",
"CD27",
"SELL",
"GZMB",
"GZMA",
"IL4R",
"TNFRSF13B",
"FCER2",
"FCGR3A")


p <- plot_density(data_10, markers, joint = FALSE, combine = FALSE, direction = -1)
q <- patchwork::wrap_plots(p, ncol = 4) 
ggsave(here(output, "markers.png"), q, width = 13, height = 13)


prev <- map(p, ~ . + scale_color_continuous_sequential(palette = "Inferno"))
q <- patchwork::wrap_plots(prev, ncol = 4)
ggsave(here(output, "markers_inferno.png"), q, width = 13, height = 13)



prev <- map(p, ~ . + scale_color_continuous_sequential(palette = "BluYl"))
q <- patchwork::wrap_plots(prev, ncol = 4)
ggsave(here(output, "markers_BluYl.png"), q, width = 13, height = 13)


prev <- map(p, ~ . + scale_color_continuous_sequential(palette = "ag_GrnYl"))
q <- patchwork::wrap_plots(prev, ncol = 4)
ggsave(here(output, "markers_ag_GrnYl.png"), q, width = 13, height = 13)


prev <- map(p, ~ . + scale_color_continuous_sequential(palette = "YlGnBu"))
q <- patchwork::wrap_plots(prev, ncol = 4)
ggsave(here(output, "markers_YlGnBu.png"), q, width = 13, height = 13)

prev <- map(p, ~ . + scale_color_continuous_sequential(palette = "Blues 2"))
q <- patchwork::wrap_plots(prev, ncol = 4)
ggsave(here(output, "markers_blues2.png"), q, width = 13, height = 13)


prev <- map(p, ~ . + scale_color_continuous_sequential(palette = "Purple-Blue"))
q <- patchwork::wrap_plots(prev, ncol = 4)
ggsave(here(output, "markers_Purple-Blue.png"), q, width = 13, height = 13)


prev <- map(p, ~ . + scale_color_continuous_sequential(palette = "ag_Sunset"))
q <- patchwork::wrap_plots(prev, ncol = 4)
ggsave(here(output, "markers_ag_Sunset.png"), q, width = 13, height = 13)


prev <- map(p, ~ . + scale_color_continuous_sequential(palette = "BluGrn"))
q <- patchwork::wrap_plots(prev, ncol = 4)
ggsave(here(output, "markers_BluGrn.png"), q, width = 13, height = 13)



prev <- map(p, ~ . + scale_color_continuous_sequential(palette = "Sunset"))
q <- patchwork::wrap_plots(prev, ncol = 4)
ggsave(here(output, "markers_Sunset.png"), q, width = 13, height = 13)

# Plot pools --------------------------------------------------------------


data_10$pool_i <- data_10$pool %>% 
  as.character() %>% 
  str_remove("pool_") %>% 
  as.numeric()

data_10$pool_i <- factor(as.character(data_10$pool_i), 
                         unique(data_10$pool_i) %>% sort())


p4 <- DimPlot(data_10, split.by = "pool_i", cols = alpha("black", 0.9), ncol = 9, 
              pt.size = 0.001) + 
  xlab("UMAP 1") + ylab("UMAP 2") + theme_pub() + NoLegend() 
ggsave(here(output, "umap_pool.png"), p4, width = 11, height = 13)
ggsave(here(output, "umap_pool.pdf"), p4, width = 11, height = 13)




# Plot PCA ----------------------------------------------------------------

data_10$cell_type <- fct_drop(data_10$cell_type)

level_order <- c("CD4 NC", "CD4 TE", "CD4 SOX4", 
                 "CD8 TE", "CD8 NC","CD8 S100B", 
                 "NK", "NK R", 
                 "Plasma", "B Mem", "B IN",  
                 "Mono C", "Mono NC", 
                 "DC")

cell_type <- recode(data_10$cell_type,
  `CD4+ KLRB1+ T cell`      = "CD4 TE",
  `CD4+ KLRB1- T cell`      = "CD4 NC",
  `CD4+ SOX4+ T cell`       = "CD4 SOX4",
  `CD8+ LTB+ T cell`        = "CD8 NC",
  `CD8+ GNLY+ NKG7+ T cell` = "CD8 TE",
  `CD8+ S100B+ T cell`      = "CD8 S100B",
  `XCL1- NK`                = "NK",
  `XCL1+ NK`                = "NK R",
  `TCL1A- FCER2- B cell`    = "B Mem",
  `TCL1A+ FCER2+ B cell`    = "B IN",
  `IgJ+ B cell`             = "Plasma",
  `Monocyte CD14+`          = "Mono C",
  `Monocyte FCGR3A+`        = "Mono NC",
  `Dendritic cell`          = "DC"
) 

data_10$new_cell_type <- factor(cell_type, level_order)

p5 <- DimPlot(data_10, group.by = "new_cell_type", cols = pal, reduction = "pca") + 
  xlab("PC 1") + ylab("PC 2")
ggsave(here(output, "pca.png"), p5, width = 6.7, height = 5)







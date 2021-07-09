# Script information ------------------------------------------------------

# title: Cell-type annotation
# author: Jose Alquicira Hernandez
# date: 2019/08/13
# description: Aggregates metadata information from classification files
# and labels cells with singularities

# Import libraries --------------------------------------------------------

# Primary
library("tidyverse")
library("dsLib")

# Secondary
library("Seurat")
library("BiocParallel")
library("scExplore")
library("future")
library("reticulate")
use_python("/share/ClusterShare/software/contrib/josalq/python3.7.4/packages/bin/python3", 
           required = TRUE)

# Set output --------------------------------------------------------------

output <- set_output("2019-08-13", "cell_type_anno")

# Read data ---------------------------------------------------------------

# Set up directories
dirs <- c("2019-07-12_myeloid_cell_annotation", "2019-07-12_B_cell_re-annotation", 
  "2019-08-12_cd8_annotation",
  "2019-07-12_myeloid_cell_annotation", "2019-07-12_B_cell_re-annotation", 
  "2019-07-15_Natural_killer_annotation", "2019-07-15_Natural_killer_annotation", 
  "2019-08-12_cd8_annotation", "2019-08-12_cd8_annotation", 
  "2019-08-12_cd8_annotation", "2019-08-12_cd4_annotation", 
  "2019-08-12_cd4_annotation", "2019-08-12_cd4_annotation", 
  "2019-07-12_myeloid_cell_annotation", "2019-07-12_cytotoxic_annotation", 
  "2019-07-16_non_cytotoxic_annotation", "2019-07-12_B_cell_re-annotation", 
  "2019-07-12_cytotoxic_annotation", "2019-07-16_non_cytotoxic_annotation", 
  "2019-07-15_Natural_killer_annotation", "2019-08-12_cd4_annotation", 
  "2019-07-12_B_cell_re-annotation")


files <- c("myeloid.RDS", "Dendritic_cell.RDS", "dc.RDS", "myeloid.RDS", 
  "B_cell.RDS", "nk.RDS", "nk_xcl1.RDS", "GNLY_NKG7_pos.RDS", "LTB_pos.RDS", 
  "S100B_pos.RDS", "KLRB1_pos.RDS", "KLRB1_neg.RDS", "SOX4_pos.RDS", 
  rep("platelets.RDS", 4), rep("erythrocytes.RDS", 5))


# Read metadata slots

readMeta <- function(d, f){
  data <- readRDS(here("results", d, f))
  data@meta.data
}


res <- bpmapply(readMeta, 
                d = dirs, f = files, 
                SIMPLIFY = FALSE, 
                BPPARAM = MulticoreParam(workers = 20)
                )


# File 1
data <- res[[1]][, "class", drop = FALSE]

# File 2 and 3
i <- c(rownames(res[[2]]), rownames(res[[3]]))
dc <- data.frame(row.names = i, class = rep("Dendritic_cell", length(i)))
data <- rbind(data, dc)

# file 4 and 5
subData <- res[4:5] %>% 
  map(~data.frame(id = rownames(.x), class = .x$class)) %>%
  do.call(rbind, .) %>% 
  `rownames<-`(NULL) %>% 
  column_to_rownames("id")

data <- rbind(data, subData)

# file 6
subData <- data.frame(row.names = rownames(res[[6]]), class = rep("XCL1- NK", nrow(res[[6]])))
data <- rbind(data, subData)

# file 7
subData <- data.frame(row.names = rownames(res[[7]]), class = rep("XCL1+ NK", nrow(res[[7]])))
data <- rbind(data, subData)

# file 8 - 13

labels <- c("CD8+ GNLY+ NKG7+ T cell", "CD8+ LTB+ T cell", "CD8+ S100B+ T cell", "CD4+ KLRB1+ T cell", "CD4+ KLRB1- T cell", "CD4+ SOX4+ T cell")

subData <- map2(res[8:13], labels, ~data.frame(id = rownames(.x), class = .y))  %>% 
  do.call(rbind, .) %>% 
  `rownames<-`(NULL) %>% 
  column_to_rownames("id")

data <- rbind(data, subData)


# file 14 - 17

subData <- map2(res[14:17], as.list("Platelets"), ~data.frame(id = rownames(.x), class = .y))  %>% 
  do.call(rbind, .) %>% 
  `rownames<-`(NULL) %>% 
  column_to_rownames("id")

data <- rbind(data, subData)


# file 18 - 22
subData <- map2(res[18:22], as.list("Erythrocytes"), ~data.frame(id = rownames(.x), class = .y))  %>% 
  do.call(rbind, .) %>% 
  `rownames<-`(NULL) %>% 
  column_to_rownames("id")

data <- rbind(data, subData)


# Read gene expression data -----------------------------------------------

# expData <- readRDS(here("results", "2019-08-06_dimRed", "dim_red.RDS"))
expData <- readRDS(here("results", "2019-08-23_norm_sctransform", "norm.RDS"))


# Create function to plot gene expression data ----------------------------

plot_gene_exp <- function(gene1, width = 5.3, height = 4.5, filename = NULL, ...){
  inicio("Plotting " %p% gene1 %p% " gene expression")
  p <- plotExp(singularities, gene1, size = 0.00000001, ...)
  
  if(is.null(filename)){
    filename <- here(output, "sing_" %p% gene1 %p% ".png")
  }else{
    filename <- here(output, "sing_" %p% filename %p% ".png")
  }
  
  ggsave(filename = filename, 
         plot = p, 
         width = width, 
         height = height, 
         dpi = 300)
  fin()
}


# Annotate singularities --------------------------------------------------

i <- (!colnames(expData) %in% rownames(data))
singularities <- expData[,i]


# Classify singularities
plan("multiprocess", workers = 2)
options(future.globals.maxSize = 1024**3 * 25)
singularities <- SCTransform(singularities, vars.to.regress = c("pool", "latent", "percent.mt"), conserve.memory = TRUE)

singularities <- RunPCA(singularities)
singularities <- RunUMAP(singularities, dims = 1:30)

singularities <- FindNeighbors(singularities, dims = 1:30)
singularities <- FindClusters(singularities, resolution = 0.1)


inicio("Plotting UMAP")
p <- plotFeature(singularities, "SCT_snn_res.0.1", size = 0.2)
ggsave(filename = here(output, "umap_sing_0-1.png"), 
       plot = p, 
       width = 6, 
       height = 4.5, 
       dpi = 300)
fin()


# plot_gene_exp(gene1 = "CD3D")
# plot_gene_exp(gene1 = "CD8A")
# plot_gene_exp(gene1 = "FCGR3A")
# plot_gene_exp(gene1 = "HBB")
# plot_gene_exp(gene1 = "PPBP")
# plot_gene_exp(gene1 = "GNLY")
# plot_gene_exp(gene1 = "MS4A1")
# plot_gene_exp(gene1 = "IGJ")
# plot_gene_exp(gene1 = "FCER1A")
# plot_gene_exp(gene1 = "CD14")
# plot_gene_exp(gene1 = "CD3D", gene2 = "CD8A", width = 8.5, height = 2.5, filename = "CD3D-CD8A")


# c1_markers <- FindMarkers(singularities,
#                           ident.1 = 1,
#                           logfc.threshold = 0.5, 
#                           min.pct = 0.8)


# plot_gene_exp(gene1 = "TCL1A")
# plot_gene_exp(gene1 = "CXCR4")
# plot_gene_exp(gene1 = "S100A4")
# plot_gene_exp(gene1 = "KLRB1")
# plot_gene_exp(gene1 = "SOX4")


label_cells <- function(x){
  switch(x,
         "0" = "CD4+ KLRB1- T cell",
         "1" = "IgJ+ B cell",
         "2" = "CD4+ KLRB1+ T cell",
         "3" = "CD4+ KLRB1- T cell",
         "4" = "Dendritic_cell",
         "5" = "CD4+ KLRB1- T cell",
         "6" = "CD4+ KLRB1- T cell",
         "7" = "TCL1A- FCER2- B cell",
         "8" = "CD4+ KLRB1- T cell",
         "9" = "Erythrocytes",
         "10" = "CD4+ KLRB1- T cell",
         "11" = "Platelets")
  
}


singularities$class <-  sapply(as.character(singularities$SCT_snn_res.0.1), label_cells) %>% as.character()
sing_meta <- singularities@meta.data[, "class", drop = FALSE]
data <- rbind(data, sing_meta)
all(colnames(expData) %in% rownames(data))

expData <- AddMetaData(expData, data, col.name = "cell_type")

expData@meta.data$cell_type <- expData@meta.data$cell_type %>% 
  str_replace("_", " ") %>% 
  as.factor()

meta.data <- expData@meta.data

meta.data %>% 
  group_by(pool, cell_type) %>% 
  tally() %>%
  ungroup() %>% 
  group_by(pool) %>% 
  mutate(perc = n / sum(n) * 100) %>% 
  ungroup() %>% 
  mutate(pool = str_remove(pool, "pool_")) %>% 
  mutate(pool = as.numeric(as.character(pool))) %>% 
  mutate(pool = factor(pool, levels = sort(unique(pool)))) -> meta_summary
  
pal <- c("#1f77b4","#ff7f0e","#2ca02c","#d62728", "#25ca7a", "#9467bd","#8c564b","#e377c2","#7f7f7f","#bcbd22","#17becf","#c5b0d5","#c49c94","#f7b6d2","#c7c7c7","#dbdb8d","#9edae5")
                  


p <- ggplot(meta_summary, aes(x = pool, y = perc, fill = cell_type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = pal) +
  theme_classic() + 
  labs(fill = "Cell type") +
  xlab("Pool") +
  ylab("Percentage")

ggsave(filename = here(output, "proportions.png"), 
       plot = p, 
       width = 15, 
       height = 5, 
       dpi = 300)


inicio("Saving data")
saveRDS(expData, file = here(output, "cell_rype.RDS"))
fin()

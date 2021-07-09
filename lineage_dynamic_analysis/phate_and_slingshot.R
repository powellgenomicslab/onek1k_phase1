#   ____________________________________________________________________________
#   Script information                                                      ####

# title: Run PHATE and Slingshot on B cells 
# author: Jose Alquicira Hernandez
# date: 2021-04-11
# description: None

#   ____________________________________________________________________________
#   HPC details                                                             ####

# screen -S phate
# qrsh -N phate -l mem_requested=100G -pe smp 2 -q short.q
# conda activate r-4.0.3

#   ____________________________________________________________________________
#   Import libraries                                                        ####

# Primary
library("tidyverse")
library("dsLib")

# Secondary
library("Seurat")
library("future")
library("phateR")
library("slingshot")
library("dittoSeq")


#   ____________________________________________________________________________
#   Future settings                                                         ####

options(future.globals.maxSize = 1024^3*100, future.seed = TRUE)
plan(multicore(workers = 2))



#   ____________________________________________________________________________
#   Helper functions                                                        ####

rename_cells2 <- function(cell_type){
  
  
  recode(cell_type,
         `CD4+ KLRB1+ T cell`      = "CD4 ET",
         `CD4+ KLRB1- T cell`      = "CD4 NC",
         `CD4+ SOX4+ T cell`       = "CD4 SOX4",
         `CD8+ LTB+ T cell`        = "CD8 NC",
         `CD8+ GNLY+ NKG7+ T cell` = "CD8 ET",
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
  
}


#   ____________________________________________________________________________
#   Set output                                                              ####

output <- set_output("2021-04-11", "b_cells_phate")

#   ____________________________________________________________________________
#   Import data                                                             ####

data <- readRDS(here("results", "2021-03-31_b_cells_subset", "b_cells.RDS"))

#   ____________________________________________________________________________
#   Normalize data                                                          ####

data$cell_type <- rename_cells2(data$cell_type)
data$cell_type <- factor(data$cell_type, levels = c("B Mem", "B IN"))

data <- SCTransform(data, 
                    vars.to.regress = c("pool", "percent.mt"), 
                    do.correct.umi = FALSE, 
                    variable.features.n = 500)
data <- RunPCA(data)
data <- RunUMAP(data, dims = 1:30)


data <- FindNeighbors(data, dims = 1:30)
data <- FindClusters(data, resolution = 0.1)


saveRDS(data, here(output, "sct_500.RDS"))

plot_gene <- function(object, gene, tag = "", reduction.use = "umap", ...){
  p <- dittoDimPlot(object, var = gene, reduction.use = reduction.use, 
                    size = 0.5, ...) + 
    scale_color_viridis_c(name = "Exp") +
    theme(panel.grid = element_blank())
  
  ggsave(here(output, tag %p% "_" %p% gene %p% ".png"), p, 
         height = 5.5, width = 7)
  
}

plot_meta <- function(object, 
                      var, 
                      tag = "", 
                      reduction.use = "umap", 
                      dim.1 = 1, 
                      dim.2 = 2,
                      save = TRUE, 
                      format = "png",
                      height = 5.5, 
                      width = 8.5, 
                      ...){
  
  red_label <- str_to_upper(reduction.use)
  xlab <- paste(red_label, dim.1)
  ylab <- paste(red_label, dim.2)
  
  
  p <- dittoDimPlot(object, 
                    var = var, 
                    reduction.use = reduction.use, 
                    size = 0.5, 
                    main = "", 
                    xlab = xlab,
                    ylab = ylab, 
                    dim.1 = dim.1,
                    dim.2 = dim.2,
                    ...) + 
    theme_pub(panel.grid = element_blank()) 
  
  if(save){
    for(f in format){
      message("Saving ", f, "...")
      ggsave(here(output, tag %p% "_" %p% var %p% "." %p% f), p, 
             height = height, width = width)      
    }
  }else{
    p
  }
}


plot_gene(data, gene = "MS4A1", tag = "all")
plot_gene(data, gene = "CD19", tag = "all")
plot_gene(data, gene = "CD34", tag = "all")
plot_gene(data, gene = "CD3D", tag = "all")
plot_gene(data, gene = "CD3E", tag = "all")
plot_gene(data, gene = "GZMB", tag = "all")
plot_gene(data, gene = "CD8A", tag = "all")

plot_meta(data, var = "seurat_clusters", tag = "all")

Idents(data) <- "seurat_clusters"

markers <- FindAllMarkers(data, only.pos = TRUE)


markers %>% 
  group_by(cluster) %>% 
  arrange(desc(avg_log2FC)) %>% 
  select(-p_val, -p_val_adj) %>% 
  slice(1:8) %>% 
  as.data.frame()




subdata <- data[, !data$seurat_clusters %in% c(4,7)]
subdata$individual[which(subdata$individual == "870_871" & subdata$latent == "b1")] <- "966_967"

subdata <- SCTransform(subdata, 
                       vars.to.regress = c("pool", "percent.mt"), 
                       do.correct.umi = FALSE, 
                       variable.features.n = 500)
subdata <- RunPCA(subdata)
subdata <- RunUMAP(subdata, dims = 1:30)

plot_meta(subdata, var = "cell_type", tag = "sub", 
          dim.1 = 2, dim.2 = 1, format = c("png", "pdf")) 

#   ____________________________________________________________________________
#   Run PHATE                                                               ####

em <- GetAssayData(subdata, slot = "scale.data") %>% t()

p <- ElbowPlot(data, ndims = 50)
ggsave(here(output, "elbow.png"), p, width = 9)

inicio("Run PHATE")
ph <- phate(em, npca = 10)
fin()

saveRDS(ph, file = here(output, "phate_500hvg_10pcs.RDS"))


subdata[["phate"]] <- CreateDimReducObject(ph$embedding, 
                                           key = "PHATE_", 
                                           assay = "SCT")

plot_meta(subdata, var = "cell_type", reduction.use = "phate", 
          tag = "sub_500hvg_10pcs")

plot_meta(subdata, var = "cell_type", tag = "sub_500hvg_10pcs", reduction.use = "phate",
          dim.1 = 2, dim.2 = 1, format = c("png", "pdf")) 

#   ____________________________________________________________________________
#   Pseudotime analysis                                                     ####

subdata <- FindNeighbors(subdata, dims = 1:10)
subdata <- FindClusters(subdata, resolution = 0.1)


plot_meta(subdata, 
          var = "seurat_clusters", 
          tag = "sub_500hvg_10pcs_res0.1", 
          reduction.use = "phate",
          dim.1 = 2, dim.2 = 1, format = c("png", "pdf"))


em <- Embeddings(subdata, reduction = "phate")

inicio("Running pseudotime")
sds <- slingshot(em, clusterLabels = subdata$seurat_clusters, 
                 start.clus = 2, end.clus = 1)
fin()
# ★  Elapsed time: 3.542 hrs

saveRDS(sds, file = here(output, "sds_500hvg_10pcs_0.1_fixed.RDS"))

i <- 1
curve_1 <- slingCurves(sds)[[i]]
pt <- curve_1$lambda %>% as.data.frame() %>% set_names("pt")
subdata <- AddMetaData(subdata, pt)
curve_1 <- curve_1$s[curve_1$ord, 1:2]
colnames(curve_1) <- c("PHATE_1", "PHATE_2")


p <- plot_meta(subdata, 
               var = "pt", 
               reduction.use = "phate", 
               tag = "curve1", 
               dim.1 = 2, 
               dim.2 = 1, 
               save = FALSE) +
  geom_path(aes(PHATE_2, PHATE_1), data = as.data.frame(curve_1),  size = 1) +
  xlab("PHATE 1") +
  ylab("PHATE 2") +
  scale_color_viridis_c(name = "Pseudotime")

ggsave(here(output, "curve1.png"), p, height = 5.5, width = 8.5)



i <- 2
curve_2 <- slingCurves(sds)[[i]]
pt <- curve_2$lambda %>% as.data.frame() %>% set_names("pt2")
subdata <- AddMetaData(subdata, pt)

curve_2 <- curve_2$s[curve_2$ord, 1:2]
colnames(curve_2) <- c("PHATE_1", "PHATE_2")

p <- plot_meta(subdata, 
               var = "pt2", 
               reduction.use = "phate", 
               tag = "curve1", 
               dim.1 = 2, 
               dim.2 = 1, 
               save = FALSE) +
  geom_path(aes(PHATE_2, PHATE_1), data = as.data.frame(curve_2),  size = 1) +
  xlab("PHATE 1") +
  ylab("PHATE 2") +
  scale_color_viridis_c(name = "Pseudotime")

ggsave(here(output, "curve2.png"), p, height = 5.5, width = 8.5)


subdata@misc[["slingshot"]] <- sds


#   ____________________________________________________________________________
#   Define bins                                                             ####

pt <- slingPseudotime(sds, na = FALSE) %>% 
  as.data.frame() %>% 
  rownames_to_column("barcode") %>% 
  select(barcode, curve = curve2) %>% 
  as_tibble()


md <- subdata[[]] %>% 
  rownames_to_column("barcode") %>% 
  select(barcode, individual, pool, latent, cell_type) %>% 
  as_tibble()

pt <- full_join(md, pt, by = "barcode")


pt$Q4 <- Hmisc::cut2(pt$curve, 
                     g = 4) %>% 
  `levels<-`("Q" %p% seq_len(4))


pt$Q5 <- Hmisc::cut2(pt$curve, 
                     g = 5) %>% 
  `levels<-`("Q" %p% seq_len(5))


pt$Q6 <- Hmisc::cut2(pt$curve, 
                     g = 6) %>% 
  `levels<-`("Q" %p% seq_len(6))


bins <- pt %>% 
  as.data.frame() %>% 
  column_to_rownames("barcode") %>% 
  select(Q4, Q5, Q6)


subdata <- AddMetaData(subdata, bins)

plot_meta(subdata, 
          var = "Q4", 
          reduction.use = "phate", 
          tag = "pt", 
          dim.1 = 2, 
          dim.2 = 1,
          format = c("png", "pdf"))

plot_meta(subdata, 
          var = "Q5", 
          reduction.use = "phate", 
          tag = "pt", 
          dim.1 = 2, 
          dim.2 = 1,
          format = c("png", "pdf"))

plot_meta(subdata, 
          var = "Q6", 
          reduction.use = "phate", 
          tag = "pt", 
          dim.1 = 2, 
          dim.2 = 1,
          format = c("png", "pdf"))



get_stats <- function(q){
  
  subdata[[]] %>% 
    group_by(individual, !!sym(q)) %>% 
    tally() %>% 
    ungroup() %>% 
    filter(n >= 3) %>% 
    group_by(!!sym(q)) %>% 
    tally()
  
}


map(c("Q4", "Q5", "Q6"), get_stats) %>% 
  map(knitr::kable)

plot_gene(subdata, 
          gene = "IGJ", 
          tag = "phate", 
          reduction.use = "phate",
          dim.1 = 2, 
          dim.2 = 1, 
          order = "increasing")



norm_data <- GetAssayData(subdata, "scale.data", assay = "SCT")

norm_data <- norm_data %>% t() %>% as.data.frame()

fit_model <- function(gene){
  gene <- data.frame(expression = gene, pseudotime = subdata$pt2)
  lm(expression ~ pseudotime, gene) %>% 
    broom::tidy() %>% 
    filter(term == "pseudotime") %>% 
    select(-term)
}

fit <- map(norm_data, fit_model)

bind_rows(fit, .id = "gene") %>% 
  mutate(p_adj = p.adjust(p.value)) %>% 
  filter(p_adj < 0.05) %>% 
  arrange(p_adj, desc(abs(estimate)), std.error) %>% 
  print(n = 25)


#   ____________________________________________________________________________
#   Export data                                                             ####

saveRDS(subdata, here(output, "b_cells_phate_slingshot.RDS"))
# subdata <- readRDS(here(output, "b_cells_phate_slingshot.RDS"))
md <- subdata[[]] %>% 
  rownames_to_column("barcode") %>% 
  as_tibble() %>% 
  select(barcode, cell_type, individual, pool, latent, SCT_snn_res.0.1, 
         pseudotime = pt2,
         Q4, Q5, Q6)

saveRDS(md, file = here(output, "slingshot_pseudotime.RDS"))

#   ____________________________________________________________________________
#   Session info                                                            ####

print_session(here(output))



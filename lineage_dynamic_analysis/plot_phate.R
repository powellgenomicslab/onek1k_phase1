#   ____________________________________________________________________________
#   Script information                                                      ####

# title: Plot phate and pseudotime results 
# author: Jose Alquicira Hernandez
# date: 2021-04-26
# description: None

#   ____________________________________________________________________________
#   HPC details                                                             ####

# screen -S phate
# qrsh -N phate -l mem_requested=100G
# conda activate r-4.0.3

#   ____________________________________________________________________________
#   Import libraries                                                        ####

# Primary
library("tidyverse")
library("dsLib")

# Secondary
library("Seurat")
library("slingshot")
library("dittoSeq")

#   ____________________________________________________________________________
#   Helper functions                                                        ####

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


#   ____________________________________________________________________________
#   Set output                                                              ####

output <- set_output("2021-04-26", "b_cell_traj_plots")

#   ____________________________________________________________________________
#   Import data                                                             ####

input <- file.path("results", "2021-04-11_b_cells_phate")

inicio("Read object")
subdata <- readRDS(here(input, "b_cells_phate_slingshot.RDS"))
fin()


sds <- subdata@misc$slingshot

#   ____________________________________________________________________________
#   Plot phate (curve)                                                      ####

curve_2 <- slingCurves(sds)[[2]]
pt <- curve_2$lambda %>% as.data.frame() %>% set_names("pt2")
curve_2 <- curve_2$s[curve_2$ord, 1:2]
colnames(curve_2) <- c("PHATE_1", "PHATE_2")

head(subdata[[]])

p <- plot_meta(subdata, 
               var = "pt2", 
               reduction.use = "phate", 
               tag = "curve2", 
               dim.1 = 2, 
               dim.2 = 1, 
               save = FALSE, 
               do.raster = TRUE) +
  geom_path(aes(PHATE_2, PHATE_1), data = as.data.frame(curve_2), size = 1) +
  xlab("PHATE 2") +
  ylab("PHATE 1") +
  scale_color_viridis_c(name = "Pseudotime")



ggsave(here(output, "phate_curve.pdf"), height = 5, width = 7.5)



#   ____________________________________________________________________________
#   Plot quantiles                                                          ####

p <- plot_meta(subdata, 
               var = "Q6", 
               reduction.use = "phate", 
               tag = "quantile", 
               dim.1 = 2, 
               dim.2 = 1, 
               save = FALSE, 
               do.raster = TRUE) +
  geom_path(aes(PHATE_2, PHATE_1), data = as.data.frame(curve_2), size = 1) +
  xlab("PHATE 2") +
  ylab("PHATE 1") +
  guides(color = guide_legend(title = "Quantile", 
                              override.aes = list(size = 5)))


ggsave(here(output, "phate_quantiles.pdf"), height = 5, width = 7.5)



#   ____________________________________________________________________________
#   Plot cell type                                                          ####


p <- plot_meta(subdata, 
               var = "cell_type", 
               reduction.use = "phate", 
               tag = "quantile", 
               dim.1 = 2, 
               dim.2 = 1, 
               save = FALSE, 
               do.raster = TRUE) +
  geom_path(aes(PHATE_2, PHATE_1), data = as.data.frame(curve_2), size = 1) +
  scale_color_manual(values = c("#F7EE55", "#F6C141"), name = "Cell type") +
  xlab("PHATE 2") +
  ylab("PHATE 1")


ggsave(here(output, "phate_cell-type.pdf"), height = 5, width = 7.5)



#   ____________________________________________________________________________
#   Session info                                                            ####

print_session(here(output))



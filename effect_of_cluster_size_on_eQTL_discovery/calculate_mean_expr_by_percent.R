##############################################################################
# Script information                                                      
# Title: Calculate mean expression by percentage
# Author: Seyhan Yazar
# Date: 2021-03-15
# Description: None
##############################################################################

# Import libraries
library(Seurat)
library(dplyr)

# Read the cell type file
df <- readRDS ("/onek1k/cell_type.RDS")

# Correct the id number for a person
df$individual[which(df$individual=="870_871" & df$latent=="b1")] <- "966_967"
df@meta.data$cell <- gsub('\\s+', '', df@meta.data$cell_type)

# Make sure the default assay is set to SCT
DefaultAssay(df) <- "SCT"
Idents(df) <- "cell"

# Subset CD4 cells
df <- subset(df, idents="CD4+KLRB1-Tcell")

# Change the identifier
Idents(df) <- "individual"

# Count number of cells per person (100%)
meta_df <- df@meta.data
meta_df %>% group_by(individual) %>% count()

# Add barcodes as a column
meta_df$barcode=rownames(meta_df)

# Subset meta data 75, 50, 25, 10, 5, 1 percent cells per person
## 75% ##
meta.75 <- meta_df%>%
    group_by(individual) %>% 
    sample_frac(.75)    
meta.75 %>% group_by(individual) %>% count() 
df.75 <- subset(df, cells=meta.75$barcode)
saveRDS(df.75 , "onek1k/effect_of_cluster_size_on_eQTL_discovery/seurat_objs/cd4.75.rds")

## 50% ##
meta.50 <- meta.75 %>%
    group_by(individual) %>% 
    sample_frac((50/75))
meta.50 %>% group_by(individual) %>% count()    
df.50 <- subset(df.75, cells=meta.50$barcode)
saveRDS(df.50 , "/onek1k/effect_of_cluster_size_on_eQTL_discovery/seurat_objs/cd4.50.rds")

## 25% ##
meta.25 <- meta.50 %>%
    group_by(individual) %>% 
    sample_frac(.5)
meta.25 %>% group_by(individual) %>% count()    
df.25 <- subset(df.50, cells=meta.25$barcode)
saveRDS(df.25, "/onek1k/effect_of_cluster_size_on_eQTL_discovery/seurat_objs/cd4.25.rds")

## 10% ##
meta.10 <- meta.25 %>%
    group_by(individual) %>% 
    sample_frac((10/25))
meta.10 %>% group_by(individual) %>% count()    
df.10 <- subset(df.25, cells=meta.10$barcode)
saveRDS(df.10, "/onek1k/effect_of_cluster_size_on_eQTL_discovery/seurat_objs/cd4.10.rds")

## 5% ##
meta.5 <- meta.10 %>%
    group_by(individual) %>% 
    sample_frac(.5)
meta.5 %>% group_by(individual) %>% count() 
df.5 <- subset(df.10, cells=meta.5$barcode)
saveRDS(df.5, "/onek1k/effect_of_cluster_size_on_eQTL_discovery/seurat_objs/cd4.5.rds")

## 1% ##
meta.1 <- meta.5 %>%
    group_by(individual) %>% 
    sample_frac((.2))
meta.1 %>% group_by(individual) %>% count() 
df.1 <- subset(df.5, cells=meta.1$barcode)
saveRDS(df.1, "/onek1k/effect_of_cluster_size_on_eQTL_discovery/seurat_objs/cd4.1.rds")
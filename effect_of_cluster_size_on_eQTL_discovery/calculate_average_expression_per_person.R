##############################################################################
# Script information                                                      
# Title: Calculate average expression per person
# Author: Seyhan Yazar
# Date: 2021-03-15
# Description: None
##############################################################################

# Import libraries
library(Seurat)
library(dplyr)
library(data.table)

# Call the arguments 
args = commandArgs(trailingOnly=TRUE)
percent <- args[1]

# Set the working directory
setwd("/onek1k/effect_of_cluster_size_on_eQTL_discovery")

# Read the cell type file
df <- readRDS (sprintf("seurat_objs/cd4.%s.rds", percent))

# Make sure the default assay is set to SCT
DefaultAssay(df) <- "SCT"
Idents(object=df) <- 'individual'

person.averages <- AverageExpression(object= df, assays = "SCT", features = NULL,
  return.seurat = FALSE, add.ident = NULL, slot = "counts",
  use.scale = FALSE, use.counts = FALSE, verbose = TRUE)

fwrite(person.averages$SCT, sprintf("count_averages/count_averages_per_person.%s.tsv", percent), 
    quote=F, row.names=TRUE, sep="\t")
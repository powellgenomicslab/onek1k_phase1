##############################################################################
# Script information                                                      
# Title: Concatenate coloc raw associations
# Author: Seyhan Yazar
# Date: 2021-03-15
# Description: None
##############################################################################

# Import libraries
library(data.table)
library(dplyr)
library(coloc)

args = commandArgs(trailingOnly=TRUE)
condition <- args[1]
# condition <- "as"
celltype <- args[2]
# celltype <- "BimmNaive"

# Main directory
main_dir <- "/onek1k"
cur_dir <- sprintf("%s/colocalisation_using_coloc/associations", main_dir)

setwd(cur_dir)
files <- list.files(path=".", pattern="*.tsv", full.name=T)

assoc_df_list <- lapply(files, function (x) {tryCatch(fread(x), error=function(e) NULL)})
assoc_df <- bind_rows(assoc_df_list)

assoc_df_85 <- assoc_df %>% filter(SNP.PP.H4 > 0.85)

summary_table <- table(assoc_df_85$condition, assoc_df_85$celltype)

write.table(summary_table, sprintf("%s/colocalisation_using_coloc/associations_summary_table.tsv", main_dir), sep="\t", quote=F)
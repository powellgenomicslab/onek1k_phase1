##############################################################################
# Script information                                                      
# Title: Identify unique eSNP-eGene pairs in the Matrix eQTL results
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
# condition <- "t1dm"
celltype <- args[2]
# celltype <- "Plasma"

# Main directory
main_dir <- "/onek1k"
cur_dir <- sprintf("%s/colocalisation_using_coloc/results/%s/%s", main_dir, condition, celltype)

setwd(cur_dir)
files <- list.files(path=".", pattern="*.rds", full.name=T)

extract_snp <- function(file) {
    df <- readRDS(file)
    result_df <- df$results %>% filter(SNP.PP.H4 > 0.5)
    if (nrow(result_df)==0) {
        NULL
    } else {
    fname <- sub(".rds", "", file)
    gname <- sub(".*_", "", fname)
    result_df$condition <- condition
    result_df$celltype <- celltype
    result_df$gene <- gname
    result_df
    }  
}

summaries <- lapply(files, extract_snp)
summaries_df <- bind_rows(summaries)

fwrite(summaries_df, sprintf("%s/colocalisation_using_coloc/associations/%s_%s_SNP_associations.tsv", main_dir, condition, celltype), sep="\t")
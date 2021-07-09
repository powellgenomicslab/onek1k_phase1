##############################################################################
# Script information                                                      
# Title: Concatenate association summaries
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
cur_dir <- sprintf("%s/colocalisation_using_coloc/results/%s/%s", main_dir, condition, celltype)

setwd(cur_dir)
files <- list.files(path=".", pattern="*.rds", full.name=T)

extract_h4 <- function(file) {
    df <- readRDS(file)
    result <- t(df$summary)
    result_df <- data.frame(result)
    fname <- sub(".rds", "", file)
    gname <- sub(".*_", "", fname)
    result_df$condition <- condition
    result_df$celltype <- celltype
    result_df$gene <- gname
    result_df
}

summaries <- lapply(files, extract_h4)
summaries_df <- bind_rows(summaries)

fwrite(summaries_df, sprintf("%s/colocalisation_using_coloc/summaries/%s_%s_summary.tsv", main_dir, condition, celltype), sep="\t")
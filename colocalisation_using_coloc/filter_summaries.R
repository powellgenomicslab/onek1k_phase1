##############################################################################
# Script information                                                      
# Title: Filter summaries
# Author: Seyhan Yazar
# Date: 2021-03-15
# Description: None
##############################################################################

# Import libraries
library(data.table)
library(dplyr)
library(coloc)

# Main directory
main_dir <- "/onek1k"
cur_dir <- sprintf("%s/colocalisation_using_coloc/summaries", main_dir)

setwd(cur_dir)
files <- list.files(path=".", pattern="*.tsv", full.name=T)

extract_h4 <- function(file) {
    df <- fread(file)
    df
}

summaries <- lapply(files, extract_h4)
summaries_df <- bind_rows(summaries)
summaries_df <- summaries_df %>% filter(PP.H4.abf > 0.75)

fwrite(summaries_df, sprintf("%s/colocalisation_using_coloc/summaries/%s_%s_summary.tsv", main_dir, condition, celltype), sep="\t")
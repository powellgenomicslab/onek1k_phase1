##############################################################################
# Script information                                                      
# Title: Combine gene lists for all conditions
# Author: Seyhan Yazar
# Date: 2021-03-15
# Description: None
##############################################################################

# Import libraries
library(data.table)
library(dplyr)
library(stringr)

# List of condition
condition_list <- c("as","crohns","ibd","ms","ra","sle","t1dm")

# Make an output directory
dir.create("/onek1k/colocalisation_using_SMR/smr_results/plot_gene_list_condition_level")

# Set the working directory
setwd("/onek1k/colocalisation_using_SMR/smr_results/plot_gene_list_condition_cell_level")

# Function to combine results
condition_plot_list <- function (x) {
    condition <- x
    files <- list.files(path=".", pattern=sprintf("^%s_*", condition))
    dataset <- lapply(files, function(x) { tryCatch(fread(x) , error=function(e) NULL)})
    filenames <- sub("*_gene_list", "", files)
    filenames <- sub(".*\\_", "", filenames)
    names(dataset) <- filenames
    dataset2<- bind_rows(dataset, .id="cell_type")
    dataset2 <- dataset2 %>% select(cell_type, ProbeChr, probeID)
    fwrite(dataset2,sprintf("/onek1k/colocalisation_using_SMR/smr_results/plot_gene_list_condition_level/%s_smr.lst", condition), 
        sep=" ", quote=FALSE, col.names=F)
}

# Apply the function over the condition_list vector
lapply (condition_list, condition_plot_list)

# Make an output directory
dir.create("/onek1k/colocalisation_using_SMR/smr_results/plot_summary_data_condition_level")

# Set the working directory
setwd("/onek1k/colocalisation_using_SMR/smr_results/plot_summary_data_condition_cell_level")

# Function to combine summary results
condition_summary <- function(y) {
    condition <- y
    files <- list.files(path=".", pattern=sprintf("^%s_*", condition))
    dataset <- lapply(files, function(x) { tryCatch(fread(x) , error=function(e) NULL)})
    filenames <- sub("*_summary.tsv", "", files)
    filenames <- sub(".*\\_", "", filenames)
    names(dataset) <- filenames
    dataset2<- bind_rows(dataset, .id="cell_type")
    fwrite(dataset2,sprintf("/onek1k/colocalisation_using_SMR/smr_results/plot_summary_data_condition_level/%s_smr_summary.tsv", 
        condition), sep="\t", quote=FALSE)
}

# Apply the function over the condition_list vector
lapply (condition_list, condition_summary)
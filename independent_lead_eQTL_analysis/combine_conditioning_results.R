##############################################################################
# Script information                                                      
# Title: Combine conditioning results
# Author: Seyhan Yazar
# Date: 2021-03-15
# Description: None
##############################################################################

# Call variables
args = commandArgs(trailingOnly=TRUE)

celltype1 <- args[1]
celltype2 <- args[2]

# celltype1 <- "CD4all"
# celltype2 <- "CD4effCM"

# celltype1 <- "Plasma"
# celltype2 <- "MonoNC"

# import libraries
library(data.table)
library(dplyr)
library(stringr)

workdir = '/onek1k/independent_lead_eQTL_analysis/conditional_spearman_results'
setwd(workdir)

files <- list.files(path=".", pattern="*.tsv")
files <- files[str_detect(files, paste(celltype1,"_vs_", celltype2, sep=""))==TRUE]
dataset <- lapply(files, function (x) { tryCatch(fread(x), error=function(e) NULL)})
filenames <- sub("*.tsv","",files)
names(dataset) <- filenames
dataset2 <- bind_rows(dataset, .id="list_name")

fwrite(dataset2,sprintf("/onek1k/independent_lead_eQTL_analysis/cs_all_results/%s_vs_%s_ALL_correlation_summary.tsv", celltype1, celltype2), sep="\t")

print("JOB IS DONE!")

quit()

##############################################################################
# Script information                                                      
# Title: Collate results and prepare plot data
# Author: Seyhan Yazar
# Date: 2021-03-15
# Description: None
##############################################################################

# Import libraries
library(data.table)
library(dplyr)
library(stringr)
library(qvalue)
library(tibble)

args = commandArgs(trailingOnly=TRUE)
condition <- args[1]
celltype <- args[2]

# condition <- "as"
# celltype <- "Bmem"

setwd(sprintf("/onek1k/colocalisation_using_SMR/smr_output/%s/%s", condition, celltype))

files <- list.files(path=".", pattern="*.smr")
dataset <- lapply(files, function(x) { tryCatch(fread(x) , error=function(e) NULL)})
# Get file names and add to the list
dataset2<- bind_rows(dataset)
nrow(dataset2)
# Calculate the qvalues for pvalues
pvalues <- dataset2$p_SMR
qobj <- qvalue(p = pvalues, pi0=1)
dataset2 <- dataset2 %>% add_column(qvalue=qobj$qvalues,localFDR=qobj$lfdr)
# Write out the raw results
fwrite(dataset2, sprintf("/onek1k/colocalisation_using_SMR/smr_results/smr_raw_results_condition_cell_level/%s_%s_raw_results.tsv", condition, celltype), sep="\t", quote=F)

# Identify the significantly assiociated pairs (summary data for plots)
smr <- dataset2
smr <- smr %>% filter (localFDR<0.05)
smr <- smr[which(-log10(smr$p_GWAS)>8),]
smr <- smr[which(-log10(smr$p_eQTL)>4),]
fwrite(smr, sprintf("/onek1k/colocalisation_using_SMR/smr_results/plot_summary_data_condition_cell_level/%s_%s_summary.tsv", condition, celltype), sep="\t", quote=F)

# Generate gene list for plots
smr <- smr %>% select(probeID, ProbeChr)
fwrite(smr, sprintf("/onek1k/colocalisation_using_SMR/smr_results/plot_gene_list_condition_cell_level/%s_%s_gene_list", condition, celltype), sep="\t", quote=F)

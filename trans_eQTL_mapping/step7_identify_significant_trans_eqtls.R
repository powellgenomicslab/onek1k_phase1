##############################################################################
# Script information                                                      
# Title: Identify significant trans eQTLs
# Author: Seyhan Yazar
# Date: 2021-03-05
# Description: None
##############################################################################

# Call the cell type
args = commandArgs(trailingOnly=TRUE)
celltype <- args[1]
print(celltype)

# celltype <- "CD4TGFbStim"

# Upload libraries
library(qvalue)
library(matrixStats)
library(data.table)
library(tidyverse)
library(broom)

cell_type <- c("CD4all","CD4effCM","CD4TGFbStim","CD8eff","CD8all","CD8unknown",
                "NKmat","NKact","Plasma","Bmem","BimmNaive","MonoC","MonoNC","DC")
cell_labels <- c("CD4 Naïve&CM","CD4 EM&TEMRA","CD4 SOX4","CD8 EM&TEMRA","CD8 Naïve&CM","CD8 S100B","NK", "NK Recruiting","Plasma",
                 "B Mem","B Imm&Naïve","Mono C","Mono NC", "DC")
new_cell_type <- c("CD4 NC", "CD4 ET", "CD4 SOX4", "CD8 ET","CD8 NC", "CD8 S100B", "NK", "NK R", "Plasma","B Mem", "B IN","Mono C", "Mono NC","DC")

cells <- data.frame(cell_type,cell_labels,new_cell_type)

new_celltype <- cells$new_cell_type[cells$cell_type==celltype]

# Set the directory path
wd <- '/onek1k/trans_eQTL_mapping'
setwd(wd)

# Get all the result files (significant eQTls)
files <- list.files(path=sprintf("./output/%s", celltype), pattern="*correlation_results.tsv", full.name=T)
dataset <- lapply(files, function(x) { fread(x)})
dataset_df <- bind_rows(dataset)
dataset_df <- dataset_df %>% filter(localFDR < 0.05)
dim(dataset_df)
head(dataset_df)

# cis-eqtl results
cis_df <- fread("../single_cell_cis_eQTL_mapping/supplementary/OneK1K_eQTLs_Results_03022021.tsv")
cis_df <- cis_df %>% filter(cell_type==new_celltype & eSNP_rank=='eSNP1')
cis_df$snp <- paste0(cis_df$Chromosome,":",cis_df$Position)
cis_df <- cis_df %>% select(cell_type, GeneID, rsID, Chromosome, Position, snp)

# join cis and trans tables
combined_df <- left_join(dataset_df, cis_df, by='snp')
head(combined_df)
results_to_save <- combined_df %>% select(cell_type, rsID, Chromosome, Position, GeneID, gene, corr.coeff, s_stat,
    pval, qvalue, localFDR)
colnames(results_to_save) <- c("celltype", "rsID", "SNP_chr", "SNP_pos", 
    "eGene_cis", "eGene_trans", "rho_corr_coeff", "s_statistics", "pvalue", "qvalue", "FDR")
results_to_save %>% head()

fwrite(results_to_save, file=sprintf("./trans_eqtls_celltype_specific_output/%s_trans_eqtls.tsv", gsub("\\s+","",new_celltype)))
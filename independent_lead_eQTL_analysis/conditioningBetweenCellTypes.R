##############################################################################
# Script information                                                      
# Title: Identify unique eSNP-eGene pairs in the Matrix eQTL results
# Author: Seyhan Yazar
# Date: 2021-03-15
# Description: None
##############################################################################


# Call variables
args = commandArgs(trailingOnly=TRUE)
celltype1 <- args[1]
celltype2 <- args[2]
chrNumber <- args[3]

# celltype1 <- "Plasma"
# celltype2 <- "NKmat"
# chrNumber <- "2"

# Import libraries
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(broom)
library(future)

setwd("/onek1k/independent_lead_eQTL_analysis")

# Input filenames
celltype1_filename <- sprintf("/onek1k/single_cell_cis_eQTL_mapping/%s/round1/%s_chr%s_lead_eSNP.tsv", celltype1, celltype1, chrNumber)
celltype2_filename <- sprintf("/onek1k/single_cell_cis_eQTL_mapping/%s/round1/%s_chr%s_lead_eSNP.tsv", celltype2, celltype2, chrNumber)
log_residuals_filename <- sprintf("/onek1k/single_cell_cis_eQTL_mapping/%s/round1/%s_chr%s_log_residuals.tsv", celltype1, celltype1, chrNumber)
genotype_filename <- sprintf("/onek1k/data/Genotype_Files/genotype_chr%s.tsv", chrNumber)

# Read in files
celltype1_df <- fread(celltype1_filename)
celltype1_df %>% head()
celltype2_df <- fread(celltype2_filename)
celltype2_df %>% head()
log_residuals_df <- fread(log_residuals_filename)

# Add prefix to colnames
colnames(celltype1_df) <- paste("celltype1", colnames(celltype1_df), sep = "_")
colnames(celltype2_df) <- paste("celltype2", colnames(celltype2_df), sep = "_")
colnames(celltype1_df)[1] <- "geneid"
colnames(celltype2_df)[1] <- "geneid"

# Merge top snps to see the overlap
top_snps_df <- full_join(celltype1_df, celltype2_df, by="geneid")
top_snps_df %>% head()

# Check if celltype1_snpid and celltype2_snpid are the same
top_snps_df <- top_snps_df %>% mutate(agreement = celltype1_snpid==celltype2_snpid)

# Number of eGenes with same lead SNP in two cell types 
print(nrow(top_snps_df %>% filter(agreement==TRUE))) 

# Subset those agree=0 to calculate new betas
different_top_snps_df <- top_snps_df %>% filter(agreement==FALSE)

No_eGenes <- print(nrow(top_snps_df %>% filter(agreement==FALSE)))

print(sprintf("There is %s overlapping SNPs!", No_eGenes))

stopifnot(No_eGenes!= 0)

# Read the genotype file
genotype_df <- fread(genotype_filename)

# Rename the files for conditioning function
# c1_exp=expression1_df
# df_snps=different_top_snps_df
# c2_geno=genotype2_df
# c1_covs=covariates1_df

calculate_pval_snp2 <- function(geneid) {
    genename <- geneid
    c1 <- celltype1
    c2 <- celltype2
    c1_exp <- log_residuals_df
    df_snps <- different_top_snps_df 
    df_geno <- genotype_df
      # print(gene)
    # get expression of the gene in cell type1
    expression_celltype1 <- c1_exp %>% select("sampleid",all_of(genename))
    colnames(expression_celltype1)[2]<- "expression"
    # identify the snps for testing
    celltype1_lead_snp <- df_snps$celltype1_snpid[df_snps$geneid==genename]
    celltype2_lead_snp <- df_snps$celltype2_snpid[df_snps$geneid==genename]
    # get the genotype of celltype2_snp in celltype1
    genotype_celltype2_snp <- df_geno %>% select("sampleid",all_of(celltype2_lead_snp))
    colnames(genotype_celltype2_snp)[2] <- "genotype_c2snp"
    df <- merge(expression_celltype1,genotype_celltype2_snp, by="sampleid")
    head(df)
    model <- lm(expression ~ genotype_c2snp, data=df)
    beta_c2snp <- as.numeric(coefficients(model)[2])
    t_stat_c2snp <- as.numeric(summary(model)$coefficients[,3][2])
    pvalue_c2snp <- as.numeric(summary(model)$coefficients[,4][2])
    summary <- data.frame(genename,celltype1,celltype2,celltype1_lead_snp,celltype2_lead_snp,beta_c2snp,t_stat_c2snp,pvalue_c2snp)
    summary
}

# Apply regression function to all eGenes
options(future.globals.maxSize = 2 * 1024^3)
plan("multicore", workers = 4)
exp_c1_snp_c2 <- future(lapply(different_top_snps_df$geneid,calculate_pval_snp2))
exp_c1_snp_c2 <- value(exp_c1_snp_c2)
results_exp_c1_snp_c2<- do.call(rbind.data.frame,exp_c1_snp_c2)
fwrite(results_exp_c1_snp_c2, sprintf("expC1_snpC2_regression_results/%s_vs_%s_chr%s_linear_model_summary.tsv", celltype1, celltype2, chrNumber), sep="\t", row.names=F, col.names=T, quote=F)

# Calculate residuals for each individual
calculate_residuals <- function(geneid) {
 genename <- geneid
    c1 <- celltype1
    c2 <- celltype2
    c1_exp <- log_residuals_df
    df_snps <- different_top_snps_df 
    df_geno <- genotype_df
      # print(gene)
    # get expression of the gene in cell type1
    expression_celltype1 <- c1_exp %>% select("sampleid",all_of(genename))
    colnames(expression_celltype1)[2]<- "expression"
    # identify the snps for testing
    celltype1_lead_snp <- df_snps$celltype1_snpid[df_snps$geneid==genename]
    celltype2_lead_snp <- df_snps$celltype2_snpid[df_snps$geneid==genename]
    # get the genotype of celltype2_snp in celltype1
    genotype_celltype2_snp <- df_geno %>% select("sampleid",all_of(celltype2_lead_snp))
    colnames(genotype_celltype2_snp)[2] <- "genotype_c2snp"
    df <- merge(expression_celltype1,genotype_celltype2_snp, by="sampleid")
    model <- lm(expression ~ genotype_c2snp, data=df)
    residuals_all <- data.frame(id=df$sampleid, fitted=residuals(model))
    residuals_all 
}

options(future.globals.maxSize = 2 * 1024^3)
plan("multicore", workers = 4)
res_expression <- future(lapply(different_top_snps_df$geneid, calculate_residuals))
res_expression <- value(res_expression)
names(res_expression) <- different_top_snps_df$geneid
r.expression <- bind_rows(res_expression, .id = "geneid")
colnames(r.expression)[2] <- "sampleid"
fwrite(r.expression, sprintf("residuals/%s_vs_%s_chr%s_linear_model_residuals.tsv", celltype1, celltype2, chrNumber), sep="\t", row.names=F, col.names=T, quote=F)

# Spearman's rank correlation 
# x is the data frame with chr and pos, y is snpid and geneid
spearman_correlation <- function (x,y) {
  gene <- y$geneid
  snp <- y$celltype1_snpid

  # Select values to test
  res_val <- r.expression %>% filter(geneid==gene)
  genotype_val <- genotype_df %>% select("sampleid",all_of(snp))
  
  # Create a test matrix
  test_df <- left_join(res_val,genotype_val,by="sampleid")
  colnames(test_df) <- c("geneid", "sampleid", "residual", "SNP")
  
  # Generate model
  model <- cor.test(x=test_df$SNP, y=test_df$residual, method = 'spearman', exact=T)
  model_table <- tidy(model) 
  model_table
}

options(future.globals.maxSize = 2 * 1024^3)
plan("multicore", workers = 4)
adjusted_spearman_df <- future(different_top_snps_df %>% group_by(celltype1_snpid,geneid) %>% 
    group_modify(spearman_correlation))
adjusted_spearman_df <- value(adjusted_spearman_df)
colnames(adjusted_spearman_df)[3] <- "c1_new_rho"
colnames(adjusted_spearman_df)[4] <- "c1_new_S_statistics"
colnames(adjusted_spearman_df)[5] <- "c1_new_p.value"
adjusted_spearman_df <- adjusted_spearman_df %>% ungroup() %>% select("geneid", "c1_new_rho",
    "c1_new_S_statistics", "c1_new_p.value")

different_top_snps_df <- full_join(different_top_snps_df,adjusted_spearman_df,by="geneid")
write.table(different_top_snps_df, sprintf("conditional_spearman_results/%s_vs_%s_chr%s_correlation_summary.tsv", celltype1, celltype2, chrNumber), sep="\t", row.names=F, col.names=T, quote=F)

print('JOB IS DONE!')

quit()


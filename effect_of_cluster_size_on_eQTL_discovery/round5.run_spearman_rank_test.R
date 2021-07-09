##############################################################################
# Script information                                                      
# Title: Conditional cis-eQTL mapping for CD4 NC cells - Round 5
# Author: Seyhan Yazar
# Date: 2020-12-23
# Description: This R script was written to run an array job per chromosome 
# for five subsets using "round5.run_spearman_rank_test.sh" script
##############################################################################

# Call variables
args = commandArgs(trailingOnly=TRUE)
percent <- args[1]
print(percent)
chrNumber <- args[2]
print(chrNumber)

# cellLabel <- "NKact"
# chrNumber <- "22"

# Import libraries
library(tidyverse)
library(broom)
library(ggplot2)
library(qvalue)
library(matrixStats)
library(reshape2)
library(data.table)
library(dplyr)
library(magrittr)
library(future)

# Set the directory path
wd <- '/onek1k/effect_of_cluster_size_on_eQTL_discovery/'
analysis_dir <- paste0(wd, "percent", percent, collapse='')
setwd(analysis_dir)

data.dir <- '/onek1k/data/'

# Input filenames
r4_residual_filename <- sprintf("round4/percent%s_chr%s_eSNP3_adjusted_residuals.tsv", percent, chrNumber)
r4_significant_SNPs <- sprintf("round4/percent%s_chr%s_round4_significant_correlation_results.tsv", percent, chrNumber)
genotype_filename <- sprintf("%sGenotype_Files/genotype_chr%s.tsv", data.dir, chrNumber)
geneLoc_filename <- sprintf("%sGene_Location_Files/geneloc_chr%s.tsv", data.dir, chrNumber)
snpLoc_filename <- sprintf("%sSNP_Location_Files/snpsloc_chr%s.tsv", data.dir, chrNumber)

# Read in files
## Count matrix
r4_residual_df <- fread(r4_residual_filename)
print("r4_residaul_df info ...")
dim(r4_residual_df)
print(r4_residual_df %>% head())

## Significant SNPs file
r4_significant_SNPs_df <- fread(r4_significant_SNPs)
print("r4_significant_SNPs_df info...")
dim(r4_significant_SNPs_df)

## Genotype file
genotype_df <- fread(genotype_filename)
print("genotype_df info...")
dim(genotype_df) 

## Gene location file
geneLoc_df <- fread(geneLoc_filename)
print("geneLoc_df info...")
dim(geneLoc_df)

## SNP location file
snpLoc_df <- fread(snpLoc_filename)
print("snpLoc_df info...")
dim(snpLoc_df)

# Identify the top eSNP for each eGene 
eSNP4 <- r4_significant_SNPs_df %>%
  group_by(geneid) %>%
  arrange(qvalue) %>%
  filter(row_number()==1)

print("eSNP4 info...")
dim(eSNP4)

fwrite(eSNP4,sprintf("round4/percent%s_chr%s_eSNP4.tsv",percent,chrNumber),sep="\t", quote=F)

# Remaning eSNPs to test
eSNPs_to_test <- r4_significant_SNPs_df %>% 
    group_by(geneid) %>%
    arrange(qvalue) %>%
    filter(row_number()!=1)  

## Subset residuals for the genes to be tested
dim(r4_residual_df)
sample_ids <- r4_residual_df$sampleid
gene_ids <- eSNP4$geneid
r4_residual_df <- r4_residual_df %>% select (all_of(gene_ids))
r4_residual_df$sampleid <- sample_ids

# Subset genotype file for the significant SNPs
dim(genotype_df)
genotype_df <- genotype_df[(genotype_df$sampleid %in% sample_ids),]
dim(genotype_df)
genotype_df <- genotype_df %>% select(sampleid, r4_significant_SNPs_df$snpid)
genotype_df <- genotype_df %>% ungroup()
snp_ids <- colnames(genotype_df[-1])

# Find residuals after adjustment of lead SNP
calculate_adjusted_residuals <- function (x) {
  gene <- x

  # select gene to regress
  exprs_val <- r4_residual_df %>% select("sampleid", all_of(gene))

  # select SNP to add
  snp = as.character(eSNP4$snpid[(eSNP4$geneid==gene)])
  snp_genotype = genotype_df %>% select(sampleid, all_of(snp))

  # Create a test df by adding covariates
  test_df <- left_join(exprs_val,snp_genotype,by="sampleid")
  colnames(test_df)[2] <- "expression"
  colnames(test_df)[3] <- "genotype"

  # Generate model
  model <- lm(expression ~ genotype , data=test_df)
  residuals=resid(model)
  residuals  
}

options(future.globals.maxSize = 1 * 1024^3)
plan("multicore", workers = 4)
adjusted_residual_mat <- future(sapply(gene_ids,calculate_adjusted_residuals))
adjusted_residual_mat <- value(adjusted_residual_mat)
rownames(adjusted_residual_mat) <- sample_ids
adjusted_residual_df <- data.frame(adjusted_residual_mat, check.names=FALSE)
adjusted_residual_df$sampleid <- sample_ids
dim(adjusted_residual_df)
fwrite(adjusted_residual_df,sprintf("round5/percent%s_chr%s_eSNP4_adjusted_residuals.tsv",percent,chrNumber),sep="\t", quote=F)

# Spearman's rank correlation 
# x is the data frame with chr and pos, y is snpid and geneid
spearman_correlation <- function (x,y) {
  gene <- y$geneid
  snp <- y$snpid

  # Select values to test
  res_val <- adjusted_residual_df %>% select("sampleid",all_of(gene))
  genotype_val <- genotype_df %>% select("sampleid", all_of(snp))
  
  # Create a test matrix
  test_df <- left_join(res_val,genotype_val,by="sampleid")
  colnames(test_df) <- c("sampleid","residual", "SNP")
  
  # Generate model
  model <- cor.test(x=test_df$SNP, y=test_df$residual, method = 'spearman', exact=T)
  model_table <- tidy(model)
  model_table
}

gene_snp_test_df <- eSNPs_to_test %>% select(snpid,geneid)
options(future.globals.maxSize = 1 * 1024^3)
plan("multicore", workers = 4)
adjusted_spearman_df <- future(gene_snp_test_df %>% group_by(snpid,geneid) %>% group_modify(spearman_correlation))
adjusted_spearman_df <- value(adjusted_spearman_df)

# Calculate the qvalues for pvalues
pvalues <- adjusted_spearman_df$p.value
qobj <- qvalue(p = pvalues, pi0=1)

# Save regardless of significance
adjusted_spearman_df <- adjusted_spearman_df %>% add_column(qvalue=qobj$qvalues,localFDR=qobj$lfdr)
nrow(adjusted_spearman_df)
fwrite(adjusted_spearman_df, sprintf("round5/percent%s_chr%s_round5_correlation_results.tsv", percent, chrNumber), sep="\t", quote=F)

# Save only significant
adjusted_spearman_df_significant <- adjusted_spearman_df[(adjusted_spearman_df$localFDR < 0.05),]
nrow(adjusted_spearman_df_significant)
fwrite(adjusted_spearman_df_significant, sprintf("round5/percent%s_chr%s_round5_significant_correlation_results.tsv", percent, chrNumber), sep="\t",  quote=F)

print('JOB IS DONE!')

quit()
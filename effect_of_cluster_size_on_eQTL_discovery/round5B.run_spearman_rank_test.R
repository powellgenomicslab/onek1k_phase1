##############################################################################
# Script information                                                      
# Title: Conditional cis-eQTL mapping for CD4 NC cells - Round 5B
# Author: Seyhan Yazar
# Date: 2020-12-23
# Description: This R script was written to run an array job per chromosome 
# for five subsets using "round5B.run_spearman_rank_test.sh" script
##############################################################################

# Call variables
args = commandArgs(trailingOnly=TRUE)
percent <- args[1]
print(percent)
chrNumber <- args[2]
print(chrNumber)

# percent <- "5"
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
analysis_dir <- paste0(wd,"percent", percent, collapse='')
setwd(analysis_dir)

data.dir <- '/onek1k/data/'

# Input filenames
r5_residual_filename <- sprintf("round5/percent%s_chr%s_eSNP4_adjusted_residuals.tsv", percent, chrNumber)
r5_significant_SNPs <- sprintf("round5/percent%s_chr%s_round5_significant_correlation_results.tsv", percent, chrNumber)
genotype_filename <- sprintf("%sGenotype_Files/genotype_chr%s.tsv", data.dir, chrNumber)
geneLoc_filename <- sprintf("%sGene_Location_Files/geneloc_chr%s.tsv", data.dir, chrNumber)
snpLoc_filename <- sprintf("%sSNP_Location_Files/snpsloc_chr%s.tsv", data.dir, chrNumber)

# Read in files
## Count matrix
r5_residual_df <- fread(r5_residual_filename)
print("r5_residaul_df info ...")
dim(r5_residual_df)
print(r5_residual_df %>% head())

## Significant SNPs file
r5_significant_SNPs_df <- fread(r5_significant_SNPs)
print("r5_significant_SNPs_df info...")
dim(r5_significant_SNPs_df)

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
eSNP5 <- r5_significant_SNPs_df %>%
  group_by(geneid) %>%
  arrange(qvalue) %>%
  filter(row_number()==1)

print("eSNP5 info...")
dim(eSNP5)

fwrite(eSNP5,sprintf("round5/percent%s_chr%s_eSNP5.tsv",percent,chrNumber),sep="\t", quote=F)

print('JOB IS DONE!')

quit()
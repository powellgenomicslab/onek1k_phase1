##############################################################################
# Script information                                                      
# Title: Conditional cis-eQTL mapping - round 5B (to save eSNP5)
# Author: Seyhan Yazar
# Date: 2020-12-23
# Description: This R script was written to run an array job per chromosome 
# for 14 cell types using "round5B.run_spearman_rank_test.sh" script
##############################################################################

# Call the cell type
args = commandArgs(trailingOnly=TRUE)
cellLabel <- args[1]
print(cellLabel)
chrNumber <- args[2]
print(chrNumber)

# Example
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
wd <- '/onek1k/single_cell_cis_eQTL_mapping/'
analysis_dir <- paste0(wd, cellLabel, collapse='')
setwd(analysis_dir)

data.dir <- '/onek1k/data/'

# Input filenames
r5_residual_filename <- sprintf("round5/%s_chr%s_eSNP4_adjusted_residuals.tsv", cellLabel, chrNumber)
r5_significant_SNPs <- sprintf("round5/%s_chr%s_round5_significant_correlation_results.tsv", cellLabel, chrNumber)
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

fwrite(eSNP5,sprintf("round5/%s_chr%s_eSNP5.tsv",cellLabel,chrNumber),sep="\t", quote=F)

print('JOB IS DONE!')

quit()
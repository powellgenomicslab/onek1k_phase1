##############################################################################
# Script information                                                      
# Title: Convert .besd files to Matrix eQTL format
# Author: Seyhan Yazar
# Date: 2021-03-15
# Description: None
##############################################################################

# Import libraries
library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(gridExtra)
library(car)
library(ggthemes)

# Call variables
args = commandArgs(trailingOnly=TRUE)
cellLabel <- args[1]
print(cellLabel)
chrNumber <- args[2]
print(chrNumber)

# cellLabel <- "BimmNaive"
# chrNumber <- "1"

# Set the directory path
wd <- '/onek1k/colocalisation_using_SMR/'
analysis_dir <- paste0(wd, cellLabel, collapse='')
setwd(analysis_dir)

df_correlation <- fread(sprintf("linear_regression/%s_chr%s_regression_results.tsv",cellLabel,chrNumber))
df_correlation$SNP <- gsub("*_.", "", df_correlation$snpid)
# df_correlation$chr <- gsub(":.*", "", df_correlation$SNP)
# df_correlation$pos <- gsub(".*:", "", df_correlation$SNP)

df_snp_info <- fread(sprintf("/onek1k/imputed_data/snps_with_maf_greaterthan0.05/chr%s.SNPs.txt", chrNumber))
colnames(df_snp_info)
df_snp_info <- df_snp_info %>% select(SNP, `REF(0)`, `ALT(1)`, ALT_Frq)
colnames(df_snp_info)

df_correlation2 <- left_join(df_correlation, df_snp_info, by="SNP")
colnames(df_correlation2)

df_correlation2 <- df_correlation2 %>% select(SNP, geneid, estimate, statistic, p.value, localFDR)
colnames(df_correlation2) <- c("SNP", "gene", "beta", "t-stat", "p-value", "FDR")
# df_correlation2$n_miss <- "0"
#colnames(df_correlation2) <- c("chr", "rs", "ps", "allel1", "allel0", "af", "beta", "se", "l_remle", "p_wald", "n_miss")
head(df_correlation2)
# df_correlation2 <- df_correlation2 %>% select(c("chr", "rs", "ps", "n_miss", "allel1", "allel0", "af", "beta", "se", "l_remle", "p_wald"))

fwrite(df_correlation2, sprintf("/onek1k/colocalisation_using_SMR/onek1kdata/%s/%s_chr%s_matrixeqtl_format.tsv",cellLabel,cellLabel,chrNumber), sep="\t", quote=F)

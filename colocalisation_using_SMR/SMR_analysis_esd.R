##############################################################################
# Script information                                                      
# Title: Generate .esd files for SMR analysis
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
# chrNumber <- "21"

# Set the directory path
wd <- '/onek1k/colocalisation_using_SMR  /matrix_eQTL/'
analysis_dir <- paste0(wd, cellLabel, collapse='')
setwd(analysis_dir)

df_correlation <- fread(sprintf("%s_chr%s_regression_results.tsv",cellLabel,chrNumber))
df_correlation$SNP <- gsub("*_.", "", df_correlation$snpid)
df_correlation$Chr <- gsub(":.*", "", df_correlation$SNP)
df_correlation$Bp <- gsub(".*:", "", df_correlation$SNP)

df_snp_info <- fread(sprintf("/onek1k/imputed_data//snps_with_maf_greaterthan0.05/chr%s.SNPs.txt", chrNumber))
colnames(df_snp_info)
df_snp_info <- df_snp_info %>% select(SNP, `REF(0)`, `ALT(1)`, ALT_Frq)
colnames(df_snp_info)

df_correlation2 <- left_join(df_correlation, df_snp_info, by="SNP")
colnames(df_correlation2)

df_correlation3 <- df_correlation2 %>% select(Chr, SNP, Bp, "ALT(1)", "REF(0)", ALT_Frq, estimate, std.error, p.value)
colnames(df_correlation3) <- c("Chr", "SNP", "Bp", "A1", "A2", "Freq", "Beta", "se", "p")
head(df_correlation3)

outfile_names <- sprintf("/onek1k/colocalisation_using_SMR/smr_data/%s/%s_chr%s_%d.esd",cellLabel,cellLabel,chrNumber, seq(1:length(df_correlation$geneid)))

for (i in 1:length(outfile_names)) {
fwrite(df_correlation3, outfile_names[i] ,sep="\t", quote=F)
}

df_gene_info <- fread(sprintf("/onek1k/Gene_Location_Files/geneloc_chr%s.tsv", chrNumber))
df_gene_info <- df_gene_info %>% filter(geneid %in% df_correlation$geneid)
colnames(df_gene_info)
df_gene_info$Chr <- df_gene_info$chr
df_gene_info$ProbeID <- df_gene_info$geneid
df_gene_info$Gene <- df_gene_info$geneid
df_gene_info$ProbeBp <- df_gene_info$start
df_gene_info$Orientation <- "NA"
df_gene_info$PathOfEsd <- sprintf("/onek1k/colocalisation_using_SMR/smr_data/%s/%s_chr%s_%d.esd",cellLabel,cellLabel,chrNumber, seq(1:length(df_gene_info$ProbeID)))
df_gene_info$GeneticDistance <- "0"
colnames(df_gene_info)
df <- df_gene_info %>% select(Chr, ProbeID, GeneticDistance, ProbeBp, Gene, Orientation, PathOfEsd)

fwrite(df, sprintf("/onek1k/colocalisation_using_SMR/smr_data/%s/%s_chr%s.flist",cellLabel,cellLabel,chrNumber), sep="\t", quote=F)

##############################################################################
# Script information                                                      
# Title: Make .esi file for SMR analysis
# Author: Seyhan Yazar
# Date: 2021-03-15
# Description: None
##############################################################################

# Import libraries
library(data.table)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)
chrNumber <- args[1]

# chrNumber <- "16"

df_bim <- fread(sprintf("/onek1k/imputed_data/filter_vcf_r08_maf005/plink_chr%s.bim", chrNumber))
df_bim$snpid <- paste0(df_bim$V2,"_",df_bim$V5)
df_bim %>% filter(V2=="16:28510044")

df_freq <- fread(sprintf("/onek1k/imputed_data/snps_with_maf_greaterthan0.05/chr%s.SNPs.txt", chrNumber))
df_freq %>% filter(SNP=="16:28510044")

df_freq$snpid <- paste0(df_freq$SNP,"_",df_freq$'ALT(1)')
df_freq$effect_allele <- df_freq$'REF(0)'
df_freq$other_allele <- df_freq$'ALT(1)'

df_joined <- left_join(df_bim, df_freq, by="snpid")
# df_joined$realsnpid <- paste0(df_joined$SNP,"_",df_joined$'ALT(1)')
# df_joined <- df_joined %>% select(V1, realsnpid, V3, V4, effect_allele,  other_allele, ALT_Frq)
# colnames(df_joined)[2] <- "snpid"
df_joined$Freq <- 1-df_joined$ALT_Frq
df_joined <- df_joined %>% select(V1, snpid, V3, V4, effect_allele, other_allele, Freq)

fwrite(df_joined, sprintf("/colocalisation_using_SMR/esi_files/chr_%s.esi", chrNumber), sep="\t", col.names=F, quote=F)

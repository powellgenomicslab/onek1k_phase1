##############################################################################
# Script information                                                      
# Title: Identify unique eSNP-eGene pairs in the Matrix eQTL results
# Author: Seyhan Yazar
# Date: 2021-03-15
# Description: None
##############################################################################

# Import libraries
library("dplyr")
library("valr")
library("refGenome")
library("BEDMatrix")
library("data.table")

# Call variables
args = commandArgs(trailingOnly=TRUE)
cellType <- args[1]
chrNumber <- args[2]

# cellType <- "CD4all"
# chrNumber <- "10"

# Read cis eqtls
cis <- fread(sprintf("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/science_revision_Dec20/smr_analysis/matrix_eQTL/%s/output_files/%s_chr%s_cis_eqtls_210211.tsv", cellType, cellType, chrNumber))
nrow(cis)

# Identify cis-eQTLs with unique SNP-gene pairs
cis_clean <- cis  %>% distinct(SNP, gene, .keep_all = T)
nrow(cis_clean)

fwrite(cis_clean, sprintf("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/science_revision_Dec20/smr_analysis/matrix_eQTL/%s/output_files/%s_chr%s_cis_eqtls_210211_uniq.tsv", cellType, cellType, chrNumber), sep="\t", quote=F)
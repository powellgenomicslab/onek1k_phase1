##############################################################################
# Script information                                                      
# Title: Generate genotype text files from bed files
# Author: Seyhan Yazar
# Date: 2021-03-02
# Description: None
##############################################################################

# Import libraries
library("dplyr")
library("valr")
library("refGenome")
library("BEDMatrix")
library("data.table")

# Set the directory path
wd <- '/onek1k/trans_eqtl_mapping'
setwd(wd)

# Prepare genotype file
bed_df <- BEDMatrix("snps_20210302.bed")
chrMatrix <- as.matrix (bed_df)
chrMatrix <- chrMatrix[, !duplicated(colnames(chrMatrix))]
sampleid <- gsub("^0_*","", rownames(chrMatrix))
rownames (chrMatrix) <- sampleid
chrMatrix <- data.frame(chrMatrix, check.names=F)
chrMatrix$sampleid <- rownames(chrMatrix)
genotype <- chrMatrix %>% select (sampleid, everything())
dim(genotype)
print(genotype[1:5,1:5])

write.table(genotype, file="genotype_20210302.tsv", quote=F,row.names=F, sep="\t")
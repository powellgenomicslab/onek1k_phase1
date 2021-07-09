##############################################################################
# Script information                                                      
# Title: Run spearman's rank correlation tests for trans eqtl mapping
# Author: Seyhan Yazar
# Date: 2021-03-02
# Description: None
##############################################################################

# Call variables
args = commandArgs(trailingOnly=TRUE)
celltype <- args[1]
print(celltype)

# celltype <- "CD4TGFbStim"

# Import libraries
library(qvalue)
library(matrixStats)
library(data.table)
library(tidyverse)
library(broom)

# Set the directory path
wd <- '/onek1k'
setwd(wd)

# Get all the result files (significant eQTls)
files <- list.files(path="./single_cell_cis_eQTL_mapping/All_Results", pattern="*_eSNP1.tsv", full.name=T)
dataset <- lapply(files, function(x) { fread(x)})
# Get file names and add to the list
filenames <- sub("./single_cell_cis_eQTL_mapping/All_Results/*", "", files)
filenames <- sub("*_eSNP1.tsv", "", filenames)
names(dataset) <- filenames
df_eqtl <- bind_rows(dataset, .id="cell_type")

cis_eqtl_cell_df <- df_eqtl %>% filter(cell_type==celltype)
cis_eqtl_cell_pairs <- cis_eqtl_cell_df %>% select(snpid, geneid)
cis_eqtl_cell_pairs$snpid <- cis_eqtl_cell_pairs$snpid %>% gsub("\\_.*", "", .)
head(cis_eqtl_cell_pairs)

# Correct snpids to match previous outputs 
# bim_df <- fread("trans_eQTL_mapping/snps_20210302.bim")
# bim_df$ids <- paste0(bim_df$V2, "_", bim_df$V5)
# bim_df$snpid <- paste0(bim_df$V2, "_", bim_df$V6)

## Genotype file
genotype_df <- fread("trans_eQTL_mapping/genotype_20210302.tsv")
dim(genotype_df)
print(genotype_df[1:5,1:5])
# Note bedmatrix get A1 from the bim file. snpid has the assessed allele(usually A2)
colnames(genotype_df) <- colnames(genotype_df) %>% gsub("\\_.*", "", .)
print(genotype_df[1:5,1:5])

log_res <- fread(sprintf("./single_cell_cis_eQTL_mapping/%s/round1/%s_chr1_log_residuals.tsv",celltype, celltype))
print(log_res[1:5,1:5])
gene_list <- colnames(log_res)[!colnames(log_res)=="sampleid"]

results <- list() 

func_corr <- function (x) {
    snp <- x
    # snp <- "10:71963608"
    genotype_val <- genotype_df[,c("sampleid", snp), with = FALSE]
    egene <- cis_eqtl_cell_pairs$geneid[cis_eqtl_cell_pairs$snpid==snp]

    for (i in gene_list[!gene_list %in% egene]) {
        gene <- i
        # gene="PPA1"
        res_val <- log_res[,c("sampleid",gene), with = FALSE]

        # Create a test matrix
        test_df <- left_join(res_val,genotype_val,by=c("sampleid"))
        colnames(test_df) <- c("sampleid","residual", "SNP")
        test_df$SNP <- as.numeric(test_df$SNP)
  
        # Generate model
        model <- cor.test(x=test_df$SNP, y=test_df$residual, method = 'spearman', exact=F)
        corr.coeff <- model$estimate
        pval <- model$p.value
        s_stat <- model$s
        results[[i]] <- data.frame(gene, snp, corr.coeff, s_stat, pval)
    }
    results_df <- bind_rows(results)
    qobj <- qvalue(p=results_df$pval, pi0=1)
    results_df <- results_df %>% mutate(qvalue=qobj$qvalues, localFDR=qobj$lfdr)
    
    fwrite(results_df, file=sprintf("%s/trans_eQTL_mapping/output/%s/%s_%s_correlation_results.tsv", wd, celltype, celltype, snp), sep="\t", quote=F)
}

lapply(cis_eqtl_cell_pairs$snpid, func_corr)

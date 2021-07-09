##############################################################################
# Script information                                                      
# Title: Run Matrix eQTL
# Author: Seyhan Yazar
# Date: 2021-03-15
# Description: None
##############################################################################

# Import libraries
library(data.table)
library(dplyr)

# Call variables
args = commandArgs(trailingOnly=TRUE)
cellType <- args[1]
chrNumber <- args[2]

# cellType <- "Plasma"
# chrNumber <- "6"

setwd(sprintf('/onek1k/colocalisation_using_SMR/matrix_eQTL/%s/input_files',cellType))

library("MatrixEQTL")

## Settings

# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

# Genotype file name
SNP_file_name = paste0('SNPs_chr',chrNumber,'_210211.tsv', collapse='');
snps_location_file_name = paste0('snpsloc_chr',chrNumber,'_210211.tsv', collapse='');

# Gene expression file name
expression_file_name = paste0('expression_chr',chrNumber,'.tsv', collapse='');
gene_location_file_name = paste0('geneloc_chr',chrNumber,'.tsv', collapse='');

# Covariates file name
# Set to character() for no covariates
covariates_file_name = character();

# Output file name
output_file_name_cis = tempfile();
output_file_name_tra = tempfile();

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 1;
pvOutputThreshold_tra = 0;

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");

# Distance for local gene-SNP pairs
cisDist = 1e6;

## Load genotype data

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

## Load covariates

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
cvrt$LoadFile(covariates_file_name);
}

## Run the analysis
snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

me = Matrix_eQTL_main(
snps = snps,
gene = gene,
cvrt = cvrt,
output_file_name     = output_file_name_tra,
pvOutputThreshold     = pvOutputThreshold_tra,
useModel = useModel,
errorCovariance = errorCovariance,
verbose = TRUE,
output_file_name.cis = output_file_name_cis,
pvOutputThreshold.cis = pvOutputThreshold_cis,
snpspos = snpspos,
genepos = genepos,
cisDist = cisDist,
pvalue.hist = "qqplot",
min.pv.by.genesnp = FALSE,
noFDRsaveMemory = FALSE);

unlink(output_file_name_tra);
unlink(output_file_name_cis);

## Results:

cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat('Detected local eQTLs:', '\n');
# show(me$cis$eqtls)
print(head(me$cis$eqtls))
print(nrow(me$cis$eqtls))

cisEQTLs <- me$cis$eqtls %>% select(snps, gene, beta, statistic, pvalue, FDR)
colnames(cisEQTLs) <- c("SNP", "gene", "beta", "t-stat", "p-value", "FDR")

fwrite(cisEQTLs, sprintf("/onek1k/colocalisation_using_SMRs/matrix_eQTL/%s/output_files/%s_chr%s_cis_eqtls_210211.tsv", cellType, cellType, chrNumber), quote=F,row.names=F,sep="\t")

print("JOB IS DONE!")

quit()
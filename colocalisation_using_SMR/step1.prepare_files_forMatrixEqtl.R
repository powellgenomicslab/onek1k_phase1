##############################################################################
# Script information                                                      
# Title: Identify effect of age on cell counts per cell type
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

# rs ids
df_hrc <- readRDS ("/onek1k/hrc_ids_all.rds")

# chrNumber <- "22"
# cellType="Plasma"

# Prepare genotype file
# library(BEDMatrix)

# Set the directory path
wd <- '/onek1k/single_cell_cis_eQTL_mapping/'
setwd(wd)

# Prepare expression file 
expression <- fread(sprintf("%s/round1/%s_chr%s_log_residuals.tsv", cellType, cellType, chrNumber))
sampleids <- expression$sampleid
t_expression <- data.frame(t(expression[,-c("sampleid")]), check.names=F)
colnames(t_expression) <- sampleids

# Prepare gene location file
setwd("/onek1k/reference_data/refdata-cellranger-hg19-3.0.0/genes")
gtf <- ensemblGenome()
read.gtf(gtf, "genes.gtf")
my_genes <- gtf@ev$genes
my_genes <- my_genes[(my_genes$seqid %in% c(1:22)),]
my_genes <- my_genes[order(my_genes$gene_source),]
my_genes <- my_genes[!duplicated(my_genes$gene_name),]
my_genes <- my_genes[,c("seqid","start","end","gene_name", "strand")]
table(duplicated(my_genes$gene_name))
setwd("/onek1k/single_cell_cis_eQTL_mapping")
my_genes$chr <- paste0('chr', my_genes$seqid)
my_genes <- my_genes[(my_genes$gene_name %in% rownames(t_expression)),]
head(my_genes)
my_genes <- my_genes[(my_genes$chr==paste0('chr',chrNumber)),]
my_genes <- my_genes[,c(4,6,2,3)]
colnames(my_genes) <- c("geneid","chr","start","end")
fwrite(my_genes, sprintf("/onek1k/colocalisation_using_SMR/matrix_eQTL/%s/input_files/geneloc_chr%s.tsv",  cellType, chrNumber), quote=F,row.names=F,sep="\t")

# subset expression dataset
t_expression  <- t_expression[c(row.names(t_expression) %in% my_genes$geneid),]
t_expression$id <- rownames(t_expression)
t_expression <- t_expression %>% select(id, everything())

fwrite(t_expression, sprintf('onek1k/colocalisation_using_SMR/matrix_eQTL/%s/input_files/expression_chr%s.tsv', cellType, chrNumber), quote=F,row.names=F, sep="\t")

# Prepare Genotype file
chr <- BEDMatrix(paste0('/onek1k/imputed_data/filter_vcf_r08_maf005/plink_chr',chrNumber,'.bed',collapse=''))
chrMatrix <- t(as.matrix (chr))
n <- gsub("^0_*","", colnames(chrMatrix))
colnames (chrMatrix) <- n
chrMatrix <- chrMatrix [,c(sampleids)]
chrMatrix <- data.frame(chrMatrix, check.names=F)
chrMatrix$snpid <- rownames(chrMatrix)
chrMatrix2 <- left_join(chrMatrix, df_hrc, by="snpid")
chrMatrix2 <- chrMatrix2 %>% select(-snpid)
genotype <- chrMatrix2 %>% select (ID, everything())
colnames(genotype)[1] <- "id"

fwrite(genotype, sprintf("/onek1k/colocalisation_using_SMR/matrix_eQTL/%s/input_files/SNPs_chr%s_210211.tsv", cellType, chrNumber), quote=F,row.names=F, sep="\t")

# Prepare snp location file
# snp_df <- fread(sprintf("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/cell_specific_eQTL_analysis_April20/%s/linear_regression/input_files/snpsloc_chr%s.tsv", cellType, chrNumber))
snps<- data.frame(chrMatrix$snpid)
snps <- cbind(snps, read.table(text = as.character(snps$chrMatrix.snpid), sep = ":"))
snps <- cbind(snps, read.table(text = as.character(snps$V2), sep = "_"))
snps <- snps [,c(1,2,4)]
snps$V1 <- paste0('chr',chrNumber)
colnames(snps) <- c("snpid","chr","pos")
snp_df <- left_join(snps, df_hrc, by="snpid")
snp_df <- snp_df[c("ID","chr","pos")]
colnames(snp_df)[1] <- "id"

fwrite(snp_df, sprintf('/onek1k/colocalisation_using_SMR/matrix_eQTL/%s/input_files/snpsloc_chr%s_210211.tsv', cellType, chrNumber), quote=F,row.names=F,sep="\t")



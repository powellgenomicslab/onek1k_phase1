##############################################################################
# Script information                                                      
# Title: Supplementary tables
# Author: Seyhan Yazar
# Date: 2020-12-25
# Description: This R script was written to make supplementary tables
##############################################################################

# Import libraries
library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(gridExtra)
library(ggthemes)
library(qvalue)

# Table for genes tested for each cell types

setwd("/onek1k/single_cell_cis_eQTL_mapping/All_Results")

cell_type <- c("BimmNaive", "Bmem", "CD4all", "CD4effCM", "CD4TGFbStim", "CD8all", "CD8eff", "CD8unknown", 
    "DC", "MonoC", "MonoNC", "NKact", "NKmat", "Plasma")

extract_genename <- function(CT) {
    # CT="BimmNaive"
    directory <- sprintf("/onek1k/single_cell_cis_eQTL_mapping/%s/round1", CT)
    setwd(directory)
    files <- list.files(path=".", pattern="*gene_SNP_pairs.tsv", full.name=T)
    dataset <- lapply(files, function(x) {fread(x)})
    dataset2<- bind_rows(dataset)
    unique(dataset2$geneid)
}

dir.create("/onek1k/single_cell_cis_eQTL_mapping/supplementary")
dir.create("/onek1k/single_cell_cis_eQTL_mapping/supplementary/gene_tested")

output <- sapply(cell_type, extract_genename)
setwd("/onek1k/single_cell_cis_eQTL_mapping/supplementary/gene_tested")
for(i in 1:length(output)){
  write.table(output[[i]], paste0(names(output)[i], ".txt"), col.names = F, row.names=F, quote=F)
}

# All the eQTLs

# Cell types
cell_type <- c("CD4all","CD4effCM","CD4TGFbStim","CD8eff","CD8all","CD8unknown",
                "NKmat","NKact","Plasma","Bmem","BimmNaive","MonoC","MonoNC","DC")
cell_labels <- c("CD4 Naïve&CM","CD4 EM&TEMRA","CD4 SOX4","CD8 EM&TEMRA","CD8 Naïve&CM","CD8 S100B","NK", "NK Recruiting","Plasma",
                 "B Mem","B Imm&Naïve","Mono C","Mono NC", "DC")
level_order <- c("CD4 NC", "CD4 ET", "CD4 SOX4", "CD8 ET","CD8 NC", "CD8 S100B", "NK", "NK R", "Plasma","B Mem", "B IN","Mono C", "Mono NC","DC")

cells <- data.frame(cell_type,cell_labels,level_order)

wd <- '/onek1k/single_cell_cis_eQTL_mapping/All_Results'
setwd(wd)
# Get all the result files (significant eQTls)
files <- list.files(path=".", pattern="*.tsv")
dataset <- lapply(files, function(x) { fread(x)})
# Get file names and add to the list
filenames <- sub("*.tsv", "", files)
names(dataset) <- filenames
df_eqtl <- bind_rows(dataset, .id="list.name")
df_eqtl$cell_type <- gsub("_.*","", df_eqtl$list.name)
df_eqtl$eSNP_rank <- gsub(".*_","",df_eqtl$list.name)
df_eqtl <- left_join(df_eqtl , cells, by="cell_type")

head(df_eqtl)
nrow(df_eqtl)

# Read HRC ids and add to the dataframes
hrc_ids <- fread("/onek1k/reference/HRC.r1-1.GRCh37.wgs.mac5.sites.tab")
head(hrc_ids)
hrc_ids$snpid <- sprintf("%s:%s_%s", hrc_ids$'#CHROM', hrc_ids$POS, hrc_ids$ALT)
hrc_ids <- hrc_ids %>% select(snpid, '#CHROM', POS, ID)
colnames(hrc_ids) [2] <- "CHR"

# Add the HRC ids
df_eqtl <- left_join(df_eqtl, hrc_ids, by="snpid")
nrow(df_eqtl)

# Read MAF files to extract allele info 
maf_files <- list.files(path="/onek1k/imputed_data/snps_with_maf_greaterthan0.05", 
    pattern=".txt", full.name=TRUE)
mafs_lst <- lapply(maf_files, function(x) fread(x))
mafs_df <- bind_rows(mafs_lst)
mafs_df$snpid <- sprintf("%s_%s", mafs_df$SNP, mafs_df$'ALT(1)')
mafs_df <- mafs_df %>% select(snpid, 'REF(0)','ALT(1)')
colnames(mafs_df) <- c("snpid", "other allele", "assessed allele")

# Add the allele info 
df_eqtl <- left_join(df_eqtl, mafs_df, by="snpid")

# Gene names 
library("refGenome")
setwd("/onek1k/reference_data/refdata-cellranger-hg19-3.0.0/genes")
gtf <- ensemblGenome()
read.gtf(gtf, "genes.gtf")
my_genes <- gtf@ev$genes
my_genes <- my_genes[(my_genes$seqid %in% c(1:22)),]
my_genes <- my_genes[order(my_genes$gene_source),]
my_genes <- my_genes[!duplicated(my_genes$gene_name),]
my_genes <- my_genes[,c("gene_id","gene_name")]
table(duplicated(my_genes$gene_name))
head(my_genes)
colnames(my_genes) <- c("ensemblid", "geneid")

df_eqtl <- left_join(df_eqtl, my_genes, by="geneid")
str(df_eqtl) 

# df_eqtl_final  <- df_eqtl %>% select(geneid, ensemblid, ID, `assessed allele`, new_names, estimate, statistic, p.value, qvalue, localFDR)
df_eqtl_final  <- df_eqtl %>% select(level_order, geneid, ensemblid, ID, CHR, POS, `assessed allele`, eSNP_rank, estimate, statistic, p.value, qvalue, localFDR)
str(df_eqtl_final)

colnames (df_eqtl_final) <- c("GeneID", "Gene_EnsemblID", "SNP_rsID", "SNP_assessed_allele", "cell_type", "rho_correlation_coefficient", "s-statistics", "pvalue", "qvalue","FDR")
colnames (df_eqtl_final) <- c("cell_type", "GeneID", "Gene_EnsemblID", "rsID", "Chromosome", "Position", "SNP_assessed_allele",  "eSNP_rank", "rho_correlation_coefficient", "S-statistics","pvalue", "qvalue","FDR")
nrow(df_eqtl_final)
specify_decimal <- function(x) trimws(format(round(x, 3), nsmall=3))
df_eqtl_final[,c("rho_correlation_coefficient","S-statistics")] <- apply(df_eqtl_final[,c("rho_correlation_coefficient","S-statistics")], 2, specify_decimal)
scientific_notation <- function (x) formatC(x, format = "e", digits = 3)
df_eqtl_final[,c("pvalue","qvalue","FDR")] <- apply(df_eqtl_final[,c("pvalue","qvalue","FDR")], 2,scientific_notation) 
setwd(wd)

fwrite(df_eqtl_final, "/onek1k/single_cell_cis_eQTL_mapping/supplementary/OneK1K_eQTLs_Results_03022021.tsv", sep="\t", quote=FALSE)


##############################################################################
# Script information                                                      
# Title: Concatenate significant trans eQTLs
# Author: Seyhan Yazar
# Date: 2021-03-15
# Description: None
##############################################################################

# Import libraries
library(data.table)
library(dplyr)

# Set the directory path
wd <- '/onek1k/trans_eQTL_mapping'
setwd(wd)

# Get all the result files (significant eQTls)
files <- list.files(path="./trans_eqtls_celltype_specific_output", pattern="*trans_eqtls.tsv", full.name=T)
dataset <- lapply(files, function(x) { fread(x)})
dataset_df <- bind_rows(dataset)
dim(dataset_df)
head(dataset_df)
dataset_df$SNP_chr <- as.character(dataset_df$SNP_chr)
#dataset_df <- dataset_df %>% select(eGene_cis, eGene_trans, SNP_chr)

### Filter out genes on the same chromosome
same_chr_func <- function(x,y) {
    # cisgene <- "FCRL3"
    # transgene <- "FCRL2"
    cisgene <- y$eGene_cis
    transgene <- y$eGene_trans
    cis_chr <- y$SNP_chr
    print(cis_chr)
    gene_loc_filename <- paste0("/onek1k/single_cell_cis_eQTL_mapping/Gene_Location_Files/geneloc_chr",cis_chr,".tsv")
    gene_loc_df <- fread(gene_loc_filename)
    cis_in <- ifelse(cisgene %in% gene_loc_df$geneid,1,0)
    trans_in <- ifelse(transgene %in% gene_loc_df$geneid,1,0)
    data.frame(cis_in, trans_in)
}

same_chr_results <- dataset_df %>% group_by(eGene_trans, eGene_cis, SNP_chr) %>% group_modify(same_chr_func)
dataset_df2 <- left_join(dataset_df, same_chr_results, by=c("eGene_cis", "eGene_trans", "SNP_chr"))
dataset_diff_chr_df <- dataset_df2 %>% filter(cis_in!=trans_in)
fwrite(dataset_diff_chr_df, file="./trans_eqtls_all_celltypes.tsv", sep="\t")

dataset_diff_chr_df01 <- dataset_diff_chr_df %>% filter(FDR < 0.01)
fwrite(dataset_diff_chr_df01, file="./trans_eqtls_all_celltypes.01.tsv", sep="\t")
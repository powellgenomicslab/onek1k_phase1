#!/usr/bin/env Rscript

# Library path
.libPaths("/directflow/SCCGGroupShare/projects/SeyhanYazar/R_libs_3.6.1")

# Upload libraries
library(data.table)
library(dplyr)

# Set the directory path
wd <- '/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/science_revision_Dec20/trans_eqtl_analysis'
setwd(wd)

trans_df <- fread("./trans_eqtls_all_celltypes.01.tsv")

# Get all the result files (significant eQTls)
setwd("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/cell_specific_eQTL_analysis_April20/Gene_Location_Files")
files <- list.files(path=".", pattern="*.tsv", full.name=T)
dataset <- lapply(files, function(x) { fread(x)})
dataset_df <- bind_rows(dataset)
dim(dataset_df)
head(dataset_df)
colnames(dataset_df) <- c("eGene_trans", "eGene_trans_chr", "eGene_trans_start", "eGene_trans_end", "eGene")
dataset_df <- dataset_df[,-5] 
trans_df_with_coord <- left_join(trans_df, dataset_df, by="eGene_trans")

# Prepare gene location file
setwd("/directflow/SCCGGroupShare/projects/data/reference_data/refdata-cellranger-hg19-3.0.0/genes")
library("refGenome")
gtf <- ensemblGenome()
read.gtf(gtf, "genes.gtf")
my_genes <- gtf@ev$genes
# my_genes <- my_genes[(my_genes$seqid %in% c(1:22,"X")),]
my_genes <- my_genes[order(my_genes$gene_source),]
my_genes <- my_genes[!duplicated(my_genes$gene_name),]
my_genes <- my_genes[,c("seqid","start","end","gene_name","strand")]
table(duplicated(my_genes$gene_name))
my_genes <- my_genes[,c(4,1,2,3)]
colnames(my_genes) <- c("eGene_trans", "eGene_trans_chr", "eGene_trans_start", "eGene_trans_end")
dim(my_genes)
trans_df_with_coord <- left_join(trans_df, my_genes, by="eGene_trans")


### Correct gene names to get the coordinates correct 
trans_to_correct_gene <- trans_df_with_coord[!complete.cases(trans_df_with_coord)]
trans_to_correct_gene$eGene_trans <- gsub("*\\.[0-9]", "", trans_to_correct_gene$eGene_trans)
trans_to_correct_gene <- trans_to_correct_gene [,-c("eGene_trans_chr", "eGene_trans_start", "eGene_trans_end")]
trans_to_correct_gene <- left_join(trans_to_correct_gene,my_genes, by="eGene_trans")

### Bind the two sets to get the complete number of associations
trans_df_with_coord <- trans_df_with_coord[complete.cases(trans_df_with_coord)]
trans_df_with_coord <- bind_rows(trans_df_with_coord, trans_to_correct_gene)

### Check to see if the genes are on the same chromosome ##
trans_df_with_coord <- trans_df_with_coord %>% filter(eGene_trans_chr!=SNP_chr)
trans_df_with_coord <- trans_df_with_coord [,-c("cis_in", "trans_in")]

setwd(wd)
fwrite(trans_df_with_coord, file="./trans_eqtls_all_celltypes.01_with_coordinates.tsv", sep="\t")
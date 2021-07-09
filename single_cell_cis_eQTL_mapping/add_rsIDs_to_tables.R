##############################################################################
# Script information                                                      
# Title: Addition of rsIDs
# Author: Seyhan Yazar
# Date: 2020-12-25
# Description: This R script was written to add the rsIDs to the tables
##############################################################################

# Import libraries
library(data.table)
library(dplyr)
library(tidyverse)

# Order of independent eQTLs
rank_order <- c("eSNP1","eSNP2","eSNP3", "eSNP4","eSNP5")

# Order of cells
cell_order <- c("CD4all","CD4effCM","CD4TGFbStim","CD8eff","CD8all","CD8unknown",
                "NKmat","NKact","Plasma","Bmem","BimmNaive","MonoC","MonoNC","DC")

cell_labels <- c("CD4 all","CD4 Eff&CM","CD4 TGFbstim","CD8 Eff","CD8 all","CD8 other","NK mat", "NK rct","Plasma","B Mem","B imm&Naive","Mono c","Mono nc", "DC")

cells <- data.frame(cell_order,cell_labels)
colnames(cells)[1] <- "celltype"
cells$celltype <- as.character(cells$celltype)

setwd("/onek1k/single_cell_cis_eQTL_mapping/All_Results")
# Get all the result files (significant eQTls)
files <- list.files(path=".", pattern="*.tsv")
dataset <- lapply(files, function(x) { tryCatch(read_tsv(x) , error=function(e) NULL)})
# Get file names and add to the list
filenames <- sub("*.tsv", "", files)
names(dataset) <- filenames

# Combine lists
dataset2<- bind_rows(dataset, .id="list_name")
dataset2$celltype <- gsub("_.*","", dataset2$list_name)
dataset2$eSNP_rank <- gsub(".*_","",dataset2$list_name)
dataset2$celltype <- factor(dataset2$celltype, levels=cell_order)
dataset2 <- left_join(dataset2, cells, by="celltype")

# hrc_filename <- "/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/cell_specific_eQTL_analysis_February20/hrc_ids.tsv"
# df_hrc <- fread(hrc_filename)
# df_hrc %>% colnames()
# df_hrc$CHRPOS2 <- paste0(df_hrc$CHROM,":",df_hrc$POS)
# df_hrc$snpid <- paste0(df_hrc$CHRPOS2,"_",df_hrc$ALT) 
# df_hrc_rs <- df_hrc %>% select(snpid,ID)
# fwrite(df_hrc_rs, "/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/cell_specific_eQTL_analysis_April20/hrc_ids_all.tsv", quote=F, sep="\t")

df_hrc_rs <- fread("/onek1k/single_cell_cis_eQTL_mapping/hrc_ids_all.tsv")
dataset2 <- left_join(dataset2, df_hrc_rs, by="snpid")
dataset2 <- dataset2 %>% select(snpid,ID) %>% distinct()
fwrite(dataset2, "/onek1k/single_cell_cis_eQTL_mapping/hrc_ids_eqtl.tsv", quote=F, sep="\t")

df_leads <- dataset2 %>% filter(eSNP_rank=="eSNP1")
fwrite(df_leads, "/onek1k/single_cell_cis_eQTL_mapping/OneK1K_eqtl_analysis_lead_SNPs.tsv", quote=F, sep="\t")

##############################################################################
# Script information                                                      
# Title: Identify overlapping cis-eGenes and colocalised genes
# Author: Seyhan Yazar
# Date: 2021-03-15
# Description: None
##############################################################################

# Import libraries
library(data.table)
library(dplyr)
library(ggplot2)

celltype <- c("CD4all","CD4effCM","CD4TGFbStim","CD8eff","CD8all","CD8unknown",
                "NKmat","NKact","Plasma","Bmem","BimmNaive","MonoC","MonoNC","DC")
cell_labels <- c("CD4 NC","CD4 ET","CD4 SOX4","CD8 ET","CD8 NC","CD8 S100B","NK", "NK R","Plasma","B Mem","B IN","Mono C","Mono NC", "DC")
cells <- data.frame(celltype,cell_labels)

setwd("/onek1k/science_revision_Dec20/colocalisation_using_coloc/associations")
files <- list.files(path=".", pattern="*.tsv", full.name=T)
assoc_df_list <- lapply(files, function (x) {tryCatch(fread(x), error=function(e) NULL)})
assoc_df <- bind_rows(assoc_df_list)

df_eqtl <- fread("/onek1k/single_cell_cis_eQTL_mapping/supplementary/OneK1K_eQTLs_Results_03022021.tsv")
df_eqtl_esnp1 <- df_eqtl %>% filter(eSNP_rank=="eSNP1") %>% select(GeneID, rsID)
colnames(df_eqtl_esnp1) <- c("gene","snp") 

assoc_ciseqtl <- assoc_df[paste0(assoc_df$snp, assoc_df$gene) %in% paste0(df_eqtl_esnp1$snp, df_eqtl_esnp1$gene), ]

## Filter PP > 80
assoc_ciseqtl <- assoc_ciseqtl %>% filter(SNP.PP.H4 > .8)
dim(assoc_ciseqtl)
assoc_ciseqtl <- left_join(assoc_ciseqtl, cells, by="celltype")
assoc_ciseqtl <- assoc_ciseqtl %>% select("condition", "cell_labels", "gene", everything())
assoc_ciseqtl <- assoc_ciseqtl %>% select(-c("celltype"))
colnames(assoc_ciseqtl) <- c("condition", "cell type", "gene name", "rsID", "V.GWAS", "z.GWAS", "r.GWAS", "lABF.GWAS",
 "V.eQTL", "z.eQTL", "r.eQTL", "lABF.eQTL", "internal.sum.lABF", "SNP.PP.H4")
assoc_ciseqtl <- assoc_ciseqtl %>% mutate(condition = toupper(condition))
fwrite(assoc_ciseqtl, "/onek1k/colocalisation_using_coloc/colocalization_with_eSNP1.tsv", sep="\t")

assoc_ciseqtl %>% group_by(condition) %>% summarise(no_egenes = n_distinct(`gene name`)) %>% data.frame()
table_ciscoloc <- assoc_ciseqtl %>% group_by(condition, `cell type`) %>% summarise(no_egenes = n_distinct(`gene name`)) %>% data.frame()
fwrite(table_ciscoloc, "/onek1k/colocalisation_using_coloc/coloc_summary_eGenes.tsv", sep="\t")

summary_tbl <- table_ciscoloc %>% tidyr::pivot_wider(names_from=condition, values_from=no_egenes)
summary_tbl <- summary_tbl %>% arrange(factor(`cell.type`, levels=cell_labels)) 

# assoc_df <- assoc_df %>% select(condition, celltype, gene)
# assoc_counts <- assoc_df %>% group_by(condition, celltype) %>% summarise(no_egenes = n_distinct(gene)) %>% data.frame()
# write.table(assoc_counts,"/onek1k/colocalisation_using_coloc/association_egenes.tsv", quote=F, sep="\t" )
# assoc_counts$id <- seq(1,98)

# Check against smr results
smr_dir <- "/onek1k/colocalisation_using_SMR/smr_results/plot_summary_data_condition_level"
setwd(smr_dir)
files <- list.files(path=".", pattern="*.tsv")
smr_df_list <- lapply(files, function (x) {tryCatch(fread(x), error=function(e) NULL)})
names(smr_df_list) <- gsub("*_smr_summary.tsv", "", files)
smr_df <- bind_rows(smr_df_list, .id="condition")
smr_df %>% group_by(condition) %>% summarise(no_egenes = n_distinct(Gene)) %>% data.frame()

#!/usr/bin/env Rscript

.libPaths("/directflow/SCCGGroupShare/projects/SeyhanYazar/R_libs_3.6.1")
library("dplyr")
library("data.table")
library("ggplot2")

cellType <- args[1]
cellLabel <- args[2]

cellType <- "Plasma"
cellLabel <- "Plasma"

setwd(sprintf('/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/trans_eQTL_analysis_November20/%s', cellType))
files <- list.files(path=".", pattern="*.tsv")
dataset <- lapply(files, function(x) { tryCatch(fread(x) , error=function(e) NULL)})
# Get file names and add to the list
filenames <- sub("*_correlation_results.tsv", "", files)
names(dataset) <- filenames

# Calculate bonferroni correction
dataset <- lapply(dataset, function (x) { cbind(x, bonferroni= p.adjust(x$pval, "bonferroni"))})

df_eqtl <- bind_rows(dataset)
no_genes_tested <- length(unique(df_eqtl$gene))
df_eqtl2 <- df_eqtl %>% filter(localFDR < 0.05)
length(unique(df_eqtl2$gene))
df_eqtl_egene_esnp_pairs <- df_eqtl2 %>% select(gene,snp) %>% unique()
length(unique(df_eqtl_egene_esnp_pairs$gene))

ms_eqtl <- fread("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/trans_eQTL_analysis_November20/OneK1K_eQTLs_For_MS_200803.tsv")
ms_eqtl$snpid <- paste0(ms_eqtl$Chromosome,":",ms_eqtl$Position,"_",ms_eqtl$SNP_assessed_allele)
ms_eqtl$ciseqtl <- ifelse(ms_eqtl$cell_type==cellLabel, "Yes", "No")
ms_eqtl <- ms_eqtl %>% filter(eSNP_rank=="eSNP1")
length(unique(ms_eqtl$snpid))
ms_eqtl_egene_esnp_pairs <- ms_eqtl %>% select (GeneID, snpid)
colnames(ms_eqtl_egene_esnp_pairs) <- c("gene","snp")
ms_eqtl_egene_esnp_pairs$cis <- "yes"

df_combine_res <- left_join(df_eqtl_egene_esnp_pairs, ms_eqtl_egene_esnp_pairs, by=c("gene","snp"))
length(unique(df_combine_res$gene))
df_combine_res <- df_combine_res %>% filter(is.na(cis))
length(unique(df_combine_res$gene))
trans.list <- split(df_combine_res, f=df_combine_res$gene)
trans.list <- lapply(trans.list, function(x) x <- x %>% select(snp))
no_trans_egenes_FDR <- length(trans.list)

eqtl_cellType <- ms_eqtl %>% filter(ciseqtl=="Yes")
trans.list.cellType <- trans.list[names(trans.list) %in% eqtl_cellType$GeneID] 
no_cis_trans_egenes <- length(trans.list.cellType)

no_eqtl_per_gene  <- data.frame(sapply(trans.list, nrow))
no_eqtl_per_gene$genes <-rownames(no_eqtl_per_gene)
colnames(no_eqtl_per_gene)[1] <- "no of trans eqtl"
trans.list.NOTcellType <- trans.list[!names(trans.list) %in% eqtl_cellType$GeneID] 
length(trans.list.NOTcellType)

length(unique(eqtl_cellType$GeneID))

df_eqtl_b <- df_eqtl %>% filter(bonferroni < 0.05)
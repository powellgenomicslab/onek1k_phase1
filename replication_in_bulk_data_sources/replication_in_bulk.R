##############################################################################
# Script information                                                      
# Title: Replication in Bulk
# Author: Seyhan Yazar
# Date: 2020-12-26
# Description: None
##############################################################################

# Set output ------------------------------------------
output <- set_output("2020-12-26", "replication-in-bulk")

# Import libraries
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(broom)

onek1k <- fread("/onek1k/single_cell_cis_eQTL_mapping/supplementary/OneK1K_eQTLs_Results_03022021.tsv")

cell_labels <- c("CD4 NC", "CD4 ET", "CD4 SOX4", "CD8 ET", "CD8 NC", "CD8 S100B", "NK", "NK R", "Plasma", "B Mem", "B IN", "Mono C", "Mono NC", "DC")

# Overlap with GTEX v8
gtex <- fread("/onek1k/replication_in_bulk_data_sources/eqtls/Whole_Blood.v8.EUR.signif_pairs.txt.gz")
key_table <- fread("/onek1k/replication_in_bulk_data_sources/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt")
gtex <- left_join(gtex, key_table, by="variant_id")
gtex$Gene_EnsemblID <- gsub("\\..*","", gtex$phenotype_id)
head(gtex)
gtex$Gene_EnsemblID <- gsub("\\..*", "", gtex$Gene_EnsemblID)
gtex$rsID <- gtex$rs_id_dbSNP151_GRCh38p7

overlap <- onek1k %>% 
                group_by(cell_type) %>%
                left_join(., gtex, by=c("rsID", "Gene_EnsemblID")) 

# Overlap with eQTLGen
eqtlgen <- fread("/onek1k/replication_in_bulk_data_sources/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt")
colnames(eqtlgen)[2] <- "rsID"
colnames(eqtlgen)[8] <- "Gene_EnsemblID"
overlap_eqtlgen <- onek1k %>% 
                group_by(cell_type) %>%
                left_join(., eqtlgen, by=c("rsID", "Gene_EnsemblID")) 

# Check overlap of eGenes
overlap$gtex_gene <- ifelse(overlap$slope=="NA", 0, 1)
overlap$gtex_gene[is.na(overlap$gtex_gene)] <- 0
overlap$gtex_eGene <- ifelse(overlap$gtex_gene==1 & overlap$pval_beta < 10^-8, 1, 0)
gtex_1k1k <- overlap %>% select("cell_type","rsID", "Gene_EnsemblID", "gtex_gene", "gtex_eGene")

gtex_1k1k$gtexsum <- gtex_1k1k$gtex_gene + gtex_1k1k$gtex_eGene
gtex_1k1k_2 <- gtex_1k1k %>% distinct(cell_type, Gene_EnsemblID, gtexsum, .keep_all=TRUE)

overlap_eqtlgen$eQTLGen_gene <- ifelse(overlap_eqtlgen$Zscore=="NA", 0, 1)
overlap_eqtlgen$eQTLGen_gene[is.na(overlap_eqtlgen$eQTLGen_gene)] <- 0
overlap_eqtlgen$eQTLGen_eGene <- ifelse(overlap_eqtlgen$eQTLGen_gene==1 & overlap_eqtlgen$FDR.y < 0.05, 1, 0)
eQTLGen_1k1k <- overlap_eqtlgen %>% select("cell_type","rsID", "Gene_EnsemblID","eQTLGen_gene", "eQTLGen_eGene")

gtex_eqtlgen <- left_join(eQTLGen_1k1k, gtex_1k1k, by=c("cell_type","rsID", "Gene_EnsemblID"))

gtex_eqtlgen$group <- ifelse (gtex_eqtlgen$gtex_eGene==1 & gtex_eqtlgen$eQTLGen_eGene==1, "GTEx eGene and eQTLGen eGene",
                        ifelse (gtex_eqtlgen$gtex_eGene==1 & gtex_eqtlgen$eQTLGen_gene==1, "GTEx eGene and eQTLGen gene",
                            ifelse (gtex_eqtlgen$gtex_eGene==1 & gtex_eqtlgen$eQTLGen_gene==0, "GTEx eGene but not in eQTLGen",
                                ifelse (gtex_eqtlgen$eQTLGen_eGene==1 & gtex_eqtlgen$gtex_gene==1, "eQTLGen eGene and GTEx gene",
                                    ifelse(gtex_eqtlgen$eQTLGen_eGene==1 & gtex_eqtlgen$gtex_gene==0, "eQTLGen eGene and not in GTEx gene",
                                        ifelse (gtex_eqtlgen$gtex_gene==1 & gtex_eqtlgen$eQTLGen_gene==0,"GTEx gene but not in eQTLGen", 
                                            ifelse (gtex_eqtlgen$gtex_gene==0 & gtex_eqtlgen$eQTLGen_gene==1, "eQTLGen gene but not in GTEx", "Neither in GTEx nor in eQTLGen")))))))

gtex_eqtlgen2 <- gtex_eqtlgen %>% distinct(cell_type, Gene_EnsemblID, group, .keep_all=TRUE)                            

library(viridis)
library(hrbrthemes)
library(ggthemes)

stacked_data <- gtex_eqtlgen2 %>% group_by (cell_type, group) %>% count()
stacked_data$cell_type <- factor(stacked_data$cell_type, levels=cell_labels)

stacked_plot <- 
    ggplot(stacked_data, aes(fill=group, y=n, x=cell_type)) + 
    geom_bar(position="stack", stat="identity") +
    scale_fill_viridis(discrete = T) +
    ggtitle("") +
    theme_clean() +
    xlab("") + ylab("eGenes(n)") +
    theme (legend.position="bottom", axis.text.x = element_text(angle = 45, hjust=1), legend.title = element_blank()) +
    guides(fill=guide_legend(nrow=4,byrow=TRUE))

ggsave ("/onek1k/replication_in_bulk_data_sources/overlap_egenes.png", stacked_plot)

gtex_eqtlgen2$gensum <- as.numeric((gtex_eqtlgen2$eQTLGen_gene + gtex_eqtlgen2$eQTLGen_eGene))
gtex_eqtlgen2$gensum <- factor(gtex_eqtlgen2$gensum, levels=c(2,1,0), labels=rev(c("not tested", "not significant", "significant")))
plot_data_gen <- gtex_eqtlgen2 %>% group_by (cell_type, gensum) %>% count()
colnames(plot_data_gen) <- c("cell type", "group", "count")
plot_data_gen$consortium <- "eQTLGen"

gtex_eqtlgen2$gtexsum <- as.numeric((gtex_eqtlgen2$gtex_gene + gtex_eqtlgen2$gtex_eGene))
gtex_eqtlgen2$gtexsum <- factor(gtex_eqtlgen2$gtexsum, levels=c(2,1,0), labels=rev(c("not tested", "not significant", "significant"))) 
gtex_eqtlgen2 %>% group_by (cell_type, gtexsum) %>% count()
plot_data_gtex <- gtex_eqtlgen2 %>% group_by (cell_type, gtexsum) %>% count()
colnames(plot_data_gtex) <- c("cell type", "group", "count")
plot_data_gtex$consortium <- "GTEx"
plot_data <- bind_rows(plot_data_gen, plot_data_gtex)
plot_data$`cell type` <- factor(plot_data$`cell type`, levels=cell_labels)

gen_plot <-
    ggplot(plot_data, aes(fill=group, y=count, x=`cell type`)) + 
    geom_bar(position="stack", stat="identity") +
    scale_fill_viridis(discrete = T) +
    ggtitle("") +
    theme_clean() +
    xlab("") + ylab("eGenes(n)") +
    theme (legend.position="bottom", axis.text.x = element_text(angle = 45, hjust=1), legend.title = element_blank()) +
    guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
    facet_wrap(~consortium)

ggsave ("/onek1k/replication_in_bulk_data_sources/overlap_egenes_210419.png", gen_plot)

# Helper function for string wrapping. 
# Default 20 character target width.
swr = function(string, nwrap=4) {
  paste(strwrap(string, width=nwrap), collapse="\n")
}
swr = Vectorize(swr)
gen_plot2$`cell type` = swr(gen_plot2$`cell type`)

gen_plot2 <- plot_data %>% 
    ggplot(aes(x=consortium, y=count, fill=group)) + geom_col() + facet_grid(.~`cell type`) +
    scale_fill_viridis(discrete = T) +
    ggtitle("") +
    theme_clean(18) +
    xlab("") + ylab("Genes(n)") +
    theme (legend.position="bottom", axis.text.x = element_text(angle = 45, hjust=1), legend.title = element_blank(), panel.border = element_blank()) +
    guides(fill=guide_legend(nrow=1,byrow=TRUE))

ggsave ("/onek1k/replication_in_bulk_data_sources/overlap_egenes_210419_2.png", gen_plot2, width=17, height=7, units="in")

colnames(overlap)

dataset1 <- overlap %>% select("cell_type", "GeneID", "Gene_EnsemblID", "rsID", "Chromosome", "Position", "SNP_assessed_allele",
     "pvalue", "FDR", "pval_nominal", "pval_nominal_threshold")
colnames(dataset1)[8:11] <- c("OneK1K_pvalue", "OneK1K_FDR", "GTEx_nominal_pvalue", "GTEx_nominal_pvalue_threshold")

dataset2 <- overlap_eqtlgen %>% select("cell_type", "GeneID", "rsID", "Pvalue", "FDR.y")    
colnames(dataset2)[4:5] <- c("eQTLGen_pvalue", "eQTLGen_FDR") 

dataset3 <- full_join(dataset1, dataset2, by=c("cell_type","GeneID", "rsID"))

dim(dataset3)

dataset3 <- dataset3 %>% select("cell_type", "GeneID", "Gene_EnsemblID", "rsID", "Chromosome", "Position", "SNP_assessed_allele", "OneK1K_pvalue", "OneK1K_FDR",
    "eQTLGen_pvalue", "eQTLGen_FDR", "GTEx_nominal_pvalue", "GTEx_nominal_pvalue_threshold")

head(dataset3) %>% data.frame

fwrite(dataset3, "/onek1k/replication_in_bulk_data_sources/eqtlgen_gtex_overlap_with_conditional_eQTLs.tsv", sep="\t" )

dataset1 %>% filter(GTEx_nominal_pvalue < GTEx_nominal_pvalue_threshold) %>% nrow()

dataset2 %>% filter(eQTLGen_FDR < 0.05) %>% nrow()

# Session info ------------------------------------------
print_session(here(output))
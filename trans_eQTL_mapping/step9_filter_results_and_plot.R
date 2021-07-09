##############################################################################
# Script information                                                      
# Title: Filter associations in the MHC region and plot the results (Figure 6)
# Author: Seyhan Yazar
# Date: 2021-07-07
# Description: None
##############################################################################

# Import libraries
Packages <- c("plyr", "dplyr", "ggplot2", "reshape2", "magrittr",
              "purrr", "tidyr", "UpSetR", "RColorBrewer", "cowplot",
              "gridExtra", "ggpubr", "ComplexHeatmap", "tibble",
              "tidyverse", "vctrs", "data.table", "ggthemes", "likert", 
              "epiR", "zoo", "ggpmisc")
invisible(lapply(Packages, library, character.only = TRUE))

# Colour palettee
tol14rainbow=c("#882E72", "#B178A6", "#D6C1DE", "#1965B0", "#5289C7", "#7BAFDE", "#4EB265", "#90C987", 
               "#CAE0AB", "#F7EE55", "#F6C141", "#F1932D", "#E8601C", "#DC050C")

# Set order of plotting
corder <- c("CD4 NC", "CD4 ET", "CD4 SOX4", "CD8 ET","CD8 NC", "CD8 S100B", "NK", "NK R", "Plasma","B Mem", "B IN","Mono C", "Mono NC","DC")

# Folder to save figures nd tables
FigTabDir <- "onek1k/Figures_and_Tables"
# dir.create(FigTabDir)

trans_df <- fread("/onek1k/trans_eQTL_mapping/trans_eqtls_all_celltypes.01_with_coordinates.tsv")
trans_df <- trans_df %>% filter(eGene_trans_chr %in% c(1:22, "MT"))

trans_df2 <- trans_df [!(trans_df$SNP_chr==6 & trans_df$SNP_pos < 5e7),]
trans_df2$rsID[trans_df2$SNP_chr==17 & trans_df2$SNP_pos=="44355634"] <- "rs2696531" 
write.table(trans_df2, sprintf("%s/trans_eqtls_non_hla_and_linkage_20210707.tsv", FigTabDir), 
            quote = F, row.names = F, sep="\t")

# Libraries
library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(viridis)
library(ggrepel)

# The dataset is provided in the gapminder library
data <- trans_df2 %>% 
  group_by(celltype, eGene_cis) %>% 
  count(eGene_cis)

data2 <- data %>% 
  group_by(celltype) %>%
  summarise(no_trans=sum(n), no_cis=n())

### Add cis-eqtl size ###
cis_df <- fread("/onek1k/single_cell_cis_eQTL_mapping/OneK1K_eQTLs_Results_03022021.tsv")
cis_size_df <- cis_df %>% select(cell_type, GeneID) %>%
  group_by(cell_type) %>%
  summarise(cis_size=n())
colnames(cis_size_df)[1] <- "celltype"

data2 <- left_join(data2, cis_size_df, by="celltype")
data2$celltype <- factor(data2$celltype,levels=corder)

buble2 <- data2 %>% arrange(desc(cis_size)) %>%
  ggplot(aes(x=no_cis, y=no_trans,size=cis_size, fill=celltype)) +
  geom_point(shape=21, colour="white") +
  scale_size(range = c(5, 24), name="Population (M)") +
  scale_fill_manual(values=tol14rainbow) +
  theme_ipsum() +
  theme(legend.position="bottom") +
  ylab("trans eGenes (n)") +
  xlab("cis eGenes (n)") +
  theme(legend.position = "none") +
  geom_text_repel(aes(label=celltype), size=4) 

buble2

ggsave(sprintf("%s/Figure6A_cis_trans_bubble_plot_based_on_cis_eqtl_discovery_size_20210707_updated.eps", FigTabDir),
       width = 12, height = 8, dpi = 150, units = "in", device='eps')

data2$population <- c(436528,61786,4065,205077,133482,34528,159820,9677,3625,48023,82068,38233,15166,8690)

cor.test(data2$no_trans, data2$population, method="spearman")
cor.test(data2$no_trans, data2$cis_size, method="spearman")

## Number of eGenes shared among cell type
sharegene <-  trans_df2 %>% group_by(celltype) %>% count(eGene_trans) 
sharegene2 <- sharegene %>% group_by(eGene_trans) %>% count()
sharegene2$sharedcells <- factor(sharegene2$n)
ggplot(sharegene2) +
  geom_bar(stat="count", aes(sharedcells), fill="white", color="black")+   
  theme_minimal(15) +
  theme_ipsum() +
  xlab("Number of shared cell types") +
  ylab("Number of trans eGenes")
ggsave(sprintf("%s/eGenes_in_multiple_celltypes_20210707_updated.eps", FigTabDir))

chr21 <- trans_df2 %>% filter(SNP_chr=="21")




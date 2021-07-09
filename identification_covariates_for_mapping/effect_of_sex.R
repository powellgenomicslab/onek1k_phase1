##############################################################################
# Script information                                                      
# Title: Identify effect of sex on cell counts per cell type
# Author: Seyhan Yazar
# Date: 2020-12-22
# Description: None
##############################################################################

# Set output ------------------------------------------
output <- set_output("2020-12-22", "sex_effect")

# Set the library path --------------------------------
.libPaths("/directflow/SCCGGroupShare/projects/SeyhanYazar/R_libs_3.6.1")

# Import libraries
library(tidyverse)
library(broom)
library(ggplot2)
library(data.table)
library(dplyr)
library(forcats)
library(hrbrthemes)
library(ggpubr)

set.seed(1234)

# Original cell labels --------------------------------
cell <- c("TCL1A+FCER2+Bcell", "TCL1A-FCER2-Bcell", "CD4+KLRB1-Tcell", "CD4+KLRB1+Tcell", 
    "CD4+SOX4+Tcell", "CD8+LTB+Tcell", "CD8+GNLY+NKG7+Tcell", "CD8+S100B+Tcell", 
    "Dendriticcell", "MonocyteCD14+", "MonocyteFCGR3A+", "XCL1+NK", "XCL1-NK", "IgJ+Bcell",
    "Erythrocytes","Platelets")

# Updated cell labels ---------------------------------
new_names <- c("B IN", "B Mem", "CD4 NC", "CD4 ET", "CD4 SOX4", "CD8 NC", "CD8 ET",
    "CD8 S100B", "DC", "Mono C", "Mono NC", "NK R", "NK", "Plasma",
    "Erythrocytes","Platelets")

cells <- data.frame(cell, new_names)

# Color palette to use in plots -----------------------
tol14rainbow=c("#882E72", "#B178A6", "#D6C1DE", "#1965B0", "#5289C7", "#7BAFDE", "#4EB265", "#90C987", 
               "#CAE0AB", "#F7EE55", "#F6C141", "#F1932D", "#E8601C", "#DC050C")

######### READ DATA ##########
covar_file <- "/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/cell_specific_eQTL_analysis_April20/CD4all/round1/CD4all_covariates.tsv"
df_covar <- fread(covar_file)
colnames(df_covar)[1] <-"individual"

CellsMeta <- readRDS("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/azimuth_cell_specific_eQTL/metadata.RDS")
# Remove spaces in cell label 
CellsMeta$cell <- gsub('\\s+', '', CellsMeta$cell_type)

# Add the updated cell labels to metadata
CellsMeta <- left_join(CellsMeta, cells, by="cell")
head(CellsMeta)

# Count cells per individual
cell_counts <- CellsMeta %>%
    # remove erythorcyes and platelets
    filter(new_names!=c("Erythrocytes","Platelets")) %>%
    droplevels %>%
    group_by(individual) %>%
    count(new_names) 

# Reorder the cells to match the color palette
order_names <- c("CD4 NC", "CD4 ET", "CD4 SOX4", "CD8 ET","CD8 NC", "CD8 S100B", "NK", "NK R", "Plasma","B Mem", "B IN","Mono C", "Mono NC","DC")

# Merge cell counts and sex dataframes
cell_counts <- left_join(cell_counts, df_covar, by="individual")
# Add labels to sex info
cell_counts$SEX <- factor(cell_counts$SEX, labels=c("Male","Female"))

# Plot cell counts between males and females 
sex_vs_counts_ps <- cell_counts %>%
    mutate(new_names=fct_relevel(new_names, order_names)) %>%
    ggboxplot(x="SEX", y="n", fill="new_names") +
    theme_ipsum(plot_title_face = "bold", plot_title_size = 18) +
    theme(legend.position="none", strip.text = element_text(face="bold")) +
    xlab("") + ylab("cell counts") + 
    stat_compare_means(aes(label=paste0("p=",..p.format..)), label.x=1.3, label.y.npc=0.9) +
    scale_fill_manual(values=tol14rainbow, labels=new_names) +
    facet_wrap(~new_names, scales="free", ncol=3)
ggsave(sex_vs_counts_ps, file='/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/science_revision_Dec20/effect_of_age/sex_vs_counts_20210330_ps.png',
                            width = 9, height = 12, dpi = 300, units = "in")

# Session info ------------------------------------------
print_session(here(output))
##############################################################################
# Script information                                                      
# Title: Identify effect of age on cell counts per cell type
# Author: Seyhan Yazar
# Date: 2020-12-22
# Description: None
##############################################################################

# Set output ------------------------------------------
output <- set_output("2020-12-22", "age_effect")

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
age_df <- fread("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/cell_specific_eQTL_analysis_April20/age_covariate.tsv")
colnames(age_df)[1] <- "individual"

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

# Merge cell counts and age dataframes
cell_counts <- left_join(cell_counts, age_df, by="individual")

# Plot correlation of cell counts with age
age_vs_counts <- cell_counts %>%
    mutate(new_names=fct_relevel(new_names, order_names)) %>%
    ggplot(aes(x=age, y=n, color=new_names), stat="identity") +
    geom_point() +
    geom_smooth(method='lm', formula= y~x, color='black', se = FALSE, show.legend = FALSE) +
    labs(x="age", y="cell counts(n)", fill = "Cell Type") +
    theme_ipsum(plot_title_face = "bold", plot_title_size = 18) +
    theme(legend.position="none", strip.text = element_text(face="bold")) +
    scale_color_manual(values=tol14rainbow, labels=new_names) +
    facet_wrap(~new_names, scales="free", ncol=3) 
ggsave(age_vs_counts, file='/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/science_revision_Dec20/age_vs_counts_20210112.png',
                            width = 9, height = 12, dpi = 300, units = "in")

# Get the correlation coefficients and p-values
corr_coeffs <- cell_counts %>%
                 group_by(new_names) %>%
                group_map(~ broom::tidy(lm(n ~ age, data = .x)))

# Save data --------------------------------------------
saveRDS(corr_coeffs, here(output, "age_effect.RDS"))

# Session info ------------------------------------------
print_session(here(output))

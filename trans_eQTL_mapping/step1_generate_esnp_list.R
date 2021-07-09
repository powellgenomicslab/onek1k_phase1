##############################################################################
# Script information                                                      
# Title: Make a eSNP list to test for trans effects
# Author: Seyhan Yazar
# Date: 2021-03-02
# Description: None
##############################################################################

# Import libraries
library("dplyr")
library("data.table")

setwd("/onek1k")

cis_eqtl_df <- fread("single_cell_cis_eQTL_mapping/supplementary/OneK1K_eQTLs_Results_03022021.tsv")

cis_eqtl_df$newId <- paste0(cis_eqtl_df$Chromosome,":",cis_eqtl_df$Position)

cis_eqtl_list <- cis_eqtl_df %>% select(newId)

fwrite(cis_eqtl_list, "trans_eqtl_mapping/OneK1K_cis_eSNP_list.txt", quote=F, col.names=F)
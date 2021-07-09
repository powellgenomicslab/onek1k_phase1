##############################################################################
# Script information                                                      
# Title: Make epi files 
# Author: Seyhan Yazar
# Date: 2021-03-15
# Description: None
##############################################################################

# Import libraries
library(data.table)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)
cellType <- args[1]
chrNumber <- args[2]

# cellType <- "BimmNaive"
# chrNumber <- "22"

input <- sprintf("/onek1kcolocalisation_using_SMR/matrix_eQTL/%s/input_files/geneloc_chr%s.tsv", cellType, chrNumber)
df_input <- fread(input)
df_input$chromosome <- gsub("chr*","", df_input$chr)
df_input$probeID <- df_input$geneid
df_input$distance <- 0
df_input$orientation <- "+"

df_input <- df_input %>% select(chromosome, probeID, distance, start, geneid, orientation)
fwrite(df_input, sprintf("/onek1k/colocalisation_using_SMR/smr_data/%s/epi_files/%s_chr%s.epi", cellType, cellType, chrNumber), sep="\t", quote=F, col.names=F)


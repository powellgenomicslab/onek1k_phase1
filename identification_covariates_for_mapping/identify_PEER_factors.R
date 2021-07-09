##############################################################################
# Script information                                                      
# Title: Identify PEER factors
# Author: Seyhan Yazar
# Date: 2020-12-22
# Description: This R script was written to run an array job for 14 cell types
# using "identify_peer_factors.sh" script
##############################################################################

# Set output ------------------------------------------
output <- set_output("2020-12-22", "peer_factors")

# Set the library path for PEER package ---------------
.libPaths("/directflow/SCCGGroupShare/projects/SeyhanYazar/.conda/envs/peerEnv/lib/R/library")

# Import libraries
library(peer)

set.seed(1234)

# Set the arguments to call from the bash script -----
args = commandArgs(trailingOnly=TRUE)
cellLabel <- args[1]

# Get expression file --------------------------------
expr = read.csv(sprintf("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/cell_specific_eQTL_analysis_April20/%s/round1/%s_expression.tsv", cellLabel, cellLabel), sep="\t", header=T)
dim(expr)
rnames = as.character(expr$sampleid)
expr2 = expr[,-1]
expr2[1:5,1:5]
cnames = colnames(expr2)

# Get covariate files ---------------------------------
covs = read.csv('/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/cell_specific_eQTL_analysis_April20/covariates.tsv', sep='\t')
age = read.csv('/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/cell_specific_eQTL_analysis_April20/age_covariate.tsv', sep='\t')
colnames(age)[1] <- "sampleid"

# Merge covarite and age file
covs = merge(covs, age, by="sampleid")
dim(covs)

# Check all the individuals exist in both files
table(expr$sampleid %in% covs$sampleid)

# ***IMPORTANT*** Match the order of individuals 
covs = covs[match(expr$sampleid, covs$sampleid),]
covs2 = covs[,-1]

# Set PEER paramaters based on the instructions from PEER package website

model = PEER()

PEER_setPhenoMean(model, as.matrix(expr2))

PEER_setCovariates(model, as.matrix(covs2))

dim(PEER_getPhenoMean(model))

# PEER_setAdd_mean(model, TRUE)

# PEER_setNmax_iterations(model, 100)

PEER_setNk(model,10) # Set to generate 10 PEER factors

PEER_getNk(model)

PEER_update(model)

# Set a directory for each cell type to save outputs
dir.create(sprintf("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/science_revision_Dec20/peer_factors/%s", cellLabel))

# Calculate and save the PEER factors
factors = PEER_getX(model)
dim(factors)
factors_df = data.frame(factors)
factors_df$sampleid = rnames
factors_df = factors_df[c(19,1:18)]
colnames(factors_df) = c("sampleid", "sex", "pc1", "pc2", "pc3", "pc4", "pc5", "pc6", "age",
    "pf1", "pf2","pf3", "pf4", "pf5",
    "pf6", "pf7","pf8", "pf9", "pf10")
write.table(factors_df, file=sprintf("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/science_revision_Dec20/peer_factors/%s/%s_peer_factors.tsv", cellLabel, cellLabel),
    sep="\t", col.names=T, row.names=F, quote=F)    

# Calculate and save the weights for each factor
weights = PEER_getW(model)
dim(weights)
weights_df = data.frame(weights)
weights_df$geneid = cnames
weights_df = weights_df[c(19,1:18)]
colnames(weights_df) = c("geneid", "sex", "pc1", "pc2", "pc3", "pc4", "pc5", "pc6", "age",
    "pf1", "pf2","pf3", "pf4", "pf5",
    "pf6", "pf7","pf8", "pf9", "pf10")
write.table(weights_df, file=sprintf("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/science_revision_Dec20/peer_factors/%s/%s_peer_factor_weights.tsv", cellLabel, cellLabel),
    sep="\t", col.names=T, quote=F, row.names=F)   

# Calculate and save the precision values
precision = PEER_getAlpha(model)
dim(precision)
precision_df = data.frame(precision)
precision_df$covariate = c("sex", "pc1", "pc2", "pc3", "pc4", "pc5", "pc6", "age",
    "pf1", "pf2","pf3", "pf4", "pf5",
    "pf6", "pf7","pf8", "pf9", "pf10")
write.table(precision_df, file=sprintf("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/science_revision_Dec20/peer_factors/%s/%s_peer_factor_precision.tsv", cellLabel, cellLabel),
    sep="\t", col.names=T, quote=F, row.names=F)   

# Calculate and save the residuals
residuals = PEER_getResiduals(model)
dim(residuals)
residuals_df = data.frame(residuals)
colnames(residuals_df) = cnames
residuals_df = cbind(sampleid=rnames, residuals_df)
write.table(residuals_df, file=sprintf("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/science_revision_Dec20/peer_factors/%s/%s_peer_factor_residuals.tsv", cellLabel, cellLabel),
    sep="\t", col.names=T, quote=F, row.names=F)

# Session info ------------------------------------------
print_session(here(output))
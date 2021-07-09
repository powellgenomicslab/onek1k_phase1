##############################################################################
# Script information                                                      
# Title: Pi analysis for effect_of_cluster_size_on_eQTL_discovery
# Author: Seyhan Yazar
# Date: 2021-03-15
# Description: None
##############################################################################

# Import libraries
library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(magrittr)
library(qvalue)

wd <- '/onek1k/science_revision_Dec20/effect_of_cluster_size_on_eQTL_discovery/percent75/round1'
setwd(wd)
# Get all the result files (significant eQTls)
files <- list.files(path=".", pattern="*_correlation_results.tsv")
files <- files[seq(1,length(files),by=2)]
dataset <- lapply(files, function(x) { fread(x)})
df_eqtl <- bind_rows(dataset)
df_eqtl_2 <- df_eqtl %>% filter(localFDR <= 0.05)
df_eqtl_2 <- df_eqtl_2 %>% select("geneid", "snpid")
df_eqtl_2$sign <- 1

wd <- '/onek1k/single_cell_cis_eQTL_mapping/CD4all/round1'
setwd(wd)
# Get all the result files (significant eQTls)
files <- list.files(path=".", pattern="*_correlation_results.tsv")
files <- files[seq(1,length(files),by=2)]
dataset <- lapply(files, function(x) { fread(x)})
df_eqtl_100 <- bind_rows(dataset)

data <- left_join(df_eqtl_100, df_eqtl_2, by=c("geneid", "snpid") )


pvalues=data$p.value
pobj=pi0est(pvalues,lambda = seq(0.05,0.95,0.05), pi0.method = c("smoother"), smooth.log.pi0="TRUE")
png("/onek1k/effect_of_cluster_size_on_eQTL_discovery/pi_test_smooth_2303.png")
plot(pobj$lambda,pobj$pi0.smooth)
dev.off()

png("/onek1k/effect_of_cluster_size_on_eQTL_discovery/pi_test_hist_75_2303.png")
hist(df_eqtl$p.value)
dev.off()
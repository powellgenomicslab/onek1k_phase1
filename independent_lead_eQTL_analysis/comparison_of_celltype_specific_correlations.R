##############################################################################
# Script information                                                      
# Title: Compare correlations between cell pairs
# Author: Seyhan Yazar
# Date: 2021-03-15
# Description: None
##############################################################################

# Import libraries
library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(gridExtra)
# library(car)
library(ggthemes)

# Call variables
args = commandArgs(trailingOnly=TRUE)
celltypeA <- args[1]
celltypeB <- args[2]

# celltypeA <- "Plasma"
# celltypeB <- "Bmem"

# celltypeB <- "Plasma"
# celltypeA <- "Bmem"

# Get the eQTL results
wd <- '/onek1k/single_cell_RNA_sequencing/All_Results'
setwd(wd)
# Get all the result files (significant eQTls)
files <- list.files(path=".", pattern="*_eSNP1.tsv")
dataset <- lapply(files, function(x) { tryCatch(fread(x) , error=function(e) NULL)})
# Get file names and add to the list
filenames <- sub("*.tsv", "", files)
names(dataset) <- filenames
df_eqtl <- bind_rows(dataset, .id="list.name")
df_eqtl$celltype <- sub("*_eSNP1", "", df_eqtl$list.name)

df_eqtl <- df_eqtl %>% select("celltype","geneid", "snpid", "estimate", "localFDR")
df_eqtl <- df_eqtl[!duplicated(df_eqtl),]
df_eqtl_celltypeA <- df_eqtl %>% filter(celltype==celltypeA)

df_eqtl_celltypeB <- df_eqtl %>% filter(celltype==celltypeB)  
df_eqtl_celltypeB <- df_eqtl_celltypeB %>% select("geneid", "snpid", "estimate", "localFDR")
  
df_together <- left_join(df_eqtl_celltypeA, df_eqtl_celltypeB, by=c("geneid", "snpid"))
df_together$significance_status <- ifelse(df_together$localFDR.y < 0.05, 1, 0)
df_together[is.na(df_together)] <- 0
df_together$significance_status <- factor(df_together$significance_status, levels=c(1,0), labels = c("significant", "non-significant"))
  
outfile <- sprintf("/onek1k/independent_lead_eQTL_analysis/%s_vs_%s_comparison.pdf",
              celltypeA, celltypeB) 
ggplot(df_together, aes(x=estimate.x, y=estimate.y, group=significance_status)) + 
    geom_point(aes(shape=significance_status), size=2) +
    scale_shape_manual(values=c(16, 3)) +
    geom_hline(yintercept=0) +
    geom_vline(xintercept=0) + 
    theme_minimal(15) +
    xlab(sprintf("Correlation coefficients of lead eQTLs in %s cells", celltypeA)) +
    ylab(sprintf("Correlation coefficients of same eSNPs in %s cells", celltypeB)) +
    theme(legend.position="top", legend.title=element_blank()) +
    geom_rug(col=rgb(.5,0,0,alpha=.2)) 
ggsave(filename = outfile)





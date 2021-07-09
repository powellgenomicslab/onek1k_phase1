#!/usr/bin/env Rscript

# Library path
.libPaths("/directflow/SCCGGroupShare/projects/SeyhanYazar/R_libs_3.6.1")

# Upload libraries
library(data.table)
library(dplyr)
library(stringr)

setwd("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/science_revision_Dec20/main_analysis/All_Results")
files <- list.files(path=".", pattern=".tsv")
dataset <- lapply(files, function(x) { fread(x) })
filenames <- sub("*.tsv", "", files)
names(dataset) <- filenames
dataset2<- bind_rows(dataset, .id="list_name")
dataset2$celltype <- gsub("_.*","",dataset2$list_name)
dataset2$eSNPrank <- gsub(".*_\\s*","",dataset2$list_name)
dataset2 <- dataset2 %>% filter(eSNPrank=="eSNP1")

# get the chr, bp, and allele seperated
eQTL <- dataset2
chr <- array(0, nrow(eQTL))
bp <- array(0, nrow(eQTL))

for(i in 1:nrow(eQTL)) {
	chr[i] <- strsplit(eQTL$snpid[i], split=":")[[1]][1]
	bp[i] <- strsplit(strsplit(eQTL$snpid[i], split=":")[[1]][2], split="_")[[1]][1]
}
​
info <- as.data.frame(cbind(chr, bp))
eQTL <- cbind(eQTL, info)
​
ms <- fread("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/cell_specific_eQTL_analysis_April20/SMR_analysis/ms/discovery_metav3.0.meta")
# MS
gwas_ms <- ms 
index <- which(gwas_ms$P!=0)
gwas_ms <- gwas_ms[index,]
index <- which(-log10(gwas_ms$P)<25)
gwas_ms <- gwas_ms[index,]
​
#GWAS / eQTL overlap
gwas_ms_sig <- gwas_ms[which(-log10(gwas_ms$P)>8),]
overlap <- match(gwas_ms_sig$BP, eQTL$bp)
​
gwas_ms_sig_overlap <- gwas_ms_sig[!is.na(overlap),]
eQTL_ms_overlap <- eQTL[overlap[!is.na(overlap)],]

MS_results <- cbind(gwas_ms_sig_overlap, eQTL_ms_overlap)
colnames(MS_results)
MS_results$newnames <- recode(MS_results$celltype,
                      "CD4effCM" = "CD4 ET",
                      "CD4all" = "CD4 NC",
                      "CD4TGFbStim" = "CD4 SOX4",
                      "CD8all" = "CD8 NC",
                      "CD8eff" = "CD8 ET",
                      "CD8unknown" = "CD8 S100B",
                      "NKmat" = "NK",
                      "NKact" = "NK R",
                      "Bmem" = "B Mem",
                      "BimmNaive" = "B IN",
                      "Plasma" = "Plasma",
                      "MonoC" = "Mono C",
                      "MonoNC" = "Mono NC",
                      "DC" = "DC"
) 
MS_results <- MS_results %>% select("newnames", "geneid", "SNP", "chr", "bp", "A1", "A2", "OR", "P", "estimate", "p.value", "localFDR")
colnames(MS_results) <- c("CELL TYPE", "GENE ID", "SNP rsID", "CHR", "BP","A1","A2", "GWAS OR", "GWAS P", "eQTL CORR COEFF", "eQTL P", "eQTL FDR")
fwrite(MS_results, "/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/science_revision_Dec20/additional_info/MS_GWAS_eQTL_overlap_210416.tsv", sep="\t", quote=F)

# Number of risk genes
length(unique(MS_results$`GENE ID`))
# risk genes in a single cell type
MS_results %>% group_by(`GENE ID`) %>% count()
MS_results %>% group_by(`GENE ID`) %>% count() %>% filter(n==1) %>% nrow()
# risk genes in two cell types
MS_results %>% group_by(`GENE ID`) %>% count() %>% filter(n==2) %>% nrow()
# risk genes in three cell types
MS_results %>% group_by(`GENE ID`) %>% count() %>% filter(n==3) %>% nrow()
# risk genes in four cell types
MS_results %>% group_by(`GENE ID`) %>% count() %>% filter(n==4) %>% nrow()

MS_results %>% filter(`GENE ID`=="RMI2")
MS_results %>% filter(`GENE ID`=="METTL21B")


celltypes=c("Bmem", "CD4all", "CD4effCM", "CD4TGFbStim", "CD8all", "CD8eff", "CD8unknown", "Plasma",
    "MonoC", "MonoNC", "NKact", "NKmat", "DC")

for (i in 1:13) {
celltype <- celltypes[i]
print(celltype)
setwd(sprintf("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/cell_specific_eQTL_analysis_April20/%s/round1", celltype))
files <- list.files(path=".", pattern="*correlation_results.tsv")
files <- files[seq(1, length(files), 2)]
dataset <- lapply(files, function(x) { fread(x) })

filenames <- sub("*.tsv", "", files)
names(dataset) <- filenames
dataset2<- bind_rows(dataset, .id="list_name")

fwrite(dataset2, sprintf("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/cell_specific_eQTL_analysis_April20/Data_For_Enrichment_Analysis/correlation_results_%s.tsv", celltype), sep="\t", quote=F)
}



fwrite(dataset2, "/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/cell_specific_eQTL_analysis_April20/Data_For_Enrichment_Analysis/OneK1K_eQTL_results.tsv", sep="\t", quote=F)

# Read HRC ids and add to the dataframes
hrc_ids <- fread("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/reference/HRC.r1-1.GRCh37.wgs.mac5.sites.tab")
head(hrc_ids)
hrc_ids$snpid <- sprintf("%s:%s_%s", hrc_ids$'#CHROM', hrc_ids$POS, hrc_ids$ALT)
hrc_ids <- hrc_ids %>% select(snpid, '#CHROM', POS, ID)
colnames(hrc_ids) [2] <- "CHR"

# Add the HRC ids
df_eqtl <- left_join(dataset2, hrc_ids, by="snpid")
nrow(df_eqtl)
df_eqtl$SNP <- df_eqtl$ID
df_eqtl <- df_eqtl %>% select(SNP, geneid) 
df_eqtl <- df_eqtl[!duplicated(df_eqtl$SNP),]



# IBD
cd <- fread("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/cell_specific_eQTL_analysis_April20/SMR_analysis/gwas_results/EUR.CD.gwas_info03_filtered.assoc")
gwas_ms <- cd
index <- which(gwas_ms$P!=0)
gwas_ms <- gwas_ms[index,]
index <- which(-log10(gwas_ms$P)<25)
gwas_ms <- gwas_ms[index,]
​
#GWAS / eQTL overlap
gwas_ms_sig <- gwas_ms[which(-log10(gwas_ms$P)>9,]
overlap <- match(gwas_ms_sig$BP, eQTL$bp)
​
gwas_ms_sig_overlap <- gwas_ms_sig[!is.na(overlap),]
eQTL_ms_overlap <- eQTL[overlap[!is.na(overlap)],]

MS_results <- cbind(gwas_ms_sig_overlap, eQTL_ms_overlap)
colnames(MS_results)
MS_results$newnames <- recode(MS_results$celltype,
                      "CD4effCM" = "CD4 TE",
                      "CD4all" = "CD4 NC",
                      "CD4TGFbStim" = "CD4 SOX4",
                      "CD8all" = "CD8 NC",
                      "CD8eff" = "CD8 ET",
                      "CD8unknown" = "CD8 S100B",
                      "NKmat" = "NK",
                      "NKact" = "NK R",
                      "Bmem" = "B Mem",
                      "BimmNaive" = "B IN",
                      "Plasma" = "Plasma",
                      "MonoC" = "MonoC",
                      "MonoNC" = "MonoNC",
                      "DC" = "DC"
) 
MS_results <- MS_results %>% select("newnames", "geneid", "SNP", "chr", "bp", "A1", "A2", "OR", "P", "estimate", "p.value", "localFDR")

fwrite(MS_results, "/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/cell_specific_eQTL_analysis_April20/Supplementary_Material/CD_results.tsv", sep="\t", quote=F)

rm(list=ls())

# ra
ra <- fread("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/cell_specific_eQTL_analysis_April20/SMR_analysis/ra/RA_GWASmeta_European_v2.txt.gz")
colnames(ra)[2] <- "chr"
colnames(ra)[3] <- "bp"
colnames(ra)[9] <- "P"

gwas_ms <- ra
index <- which(gwas_ms$P!=0)
gwas_ms <- gwas_ms[index,]
index <- which(-log10(gwas_ms$P)<25)
gwas_ms <- gwas_ms[index,]
​
#GWAS / eQTL overlap
gwas_ms_sig <- gwas_ms[which(-log10(gwas_ms$P)>9),]
overlap <- match(gwas_ms_sig$BP, eQTL$bp)
​
gwas_ms_sig_overlap <- gwas_ms_sig[!is.na(overlap),]
eQTL_ms_overlap <- eQTL[overlap[!is.na(overlap)],]

MS_results <- cbind(gwas_ms_sig_overlap, eQTL_ms_overlap)
fwrite(MS_results, "/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/cell_specific_eQTL_analysis_April20/Supplementary_Material/IBD_results.tsv", sep="\t", quote=F)

##############################################################################
# Script information                                                      
# Title: Run colocalisation analysis with coloc
# Author: Seyhan Yazar
# Date: 2021-03-15
# Description: None
##############################################################################

# Import libraries
library(data.table)
library(dplyr)
library(coloc)

# Call variables
args = commandArgs(trailingOnly=TRUE)
condition <- args[1]
# condition <- "crohns"
celltype <- args[2]
# celltype <- "MonoNC"

# Main directory
main_dir <- "/onek1k"

# Get gwas file
gwas_file <- sprintf("%s/colocalisation_using_SMR/gwas_data/%s.ma", main_dir, condition)
gwas_df <- fread(gwas_file)
# gwas_df$pseudoid <- paste0(gwas_df$SNP,"_",gwas_df$A2)
head(gwas_df)
dim(gwas_df)
gwas_df <- gwas_df[!duplicated(gwas_df$SNP)]


# Calculate case/control proportions in gwas
sas <- 146/(146+361048)
sms <- 47429/(47429+68374)
scd <- 5956/(5956+14927)
sibd <- 12882/(12882+21770)
sra <- 14261/(14261+43923)
st1dm <- 844/(844+360350)
ssle <- 4036/(4036+6959)

cond <- c("as", "crohns", "ibd", "ms", "ra", "sle", "t1dm")
prop <- c(sas, scd, sibd, sms, sra, ssle, st1dm)
size <- c("361194", "20883", "34652", "115803", "58184", "10995", "361194")
prop_df <- data.frame(cond, prop, size)

# Get matrix eqtl outputs
eqtl_files <- list.files(path=sprintf("%s/colocalisation_using_SMR/matrix_eQTL/%s/output_files", main_dir, celltype), 
    pattern="*_uniq.tsv", full.names=T)
eqtl_df_list <- lapply(eqtl_files, function (x) {tryCatch(fread(x), error=function(e) NULL)})
eqtl_df <- bind_rows(eqtl_df_list)
eqtl_df$beta_se <- eqtl_df$beta/eqtl_df$'t-stat'
eqtl_df$pseudoid <- eqtl_df$SNP

# Get onek1k mafs 
maf_files <- list.files(path=sprintf("%s/colocalisation_using_SMR/smr_data/%s", main_dir, celltype),
    pattern="*11.esi", full.names=T)
maf_df_list <- lapply(maf_files, function (x) {tryCatch(fread(x), error=function(e) NULL)})
maf_df <- bind_rows(maf_df_list)
head(maf_df)
colnames(maf_df) <- c("chr","SNP", "ignore","position","a1","a2", "maf")

# Number of individuals included in the eqtl analysis per cell type
ctype <- c("BimmNaive", "Bmem", "CD4all", "CD4effCM", "CD4TGFbStim",
    "CD8all", "CD8eff", "CD8unknown", "DC", "MonoC", "MonoNC",
    "NKact", "NKmat", "Plasma")
cno <- c("982", "982", "982", "982", "858",
    "982", "982", "981", "968", "969", "934",
    "969", "982", "795")

cno_df <- data.frame(ctype,cno=as.numeric(cno))

# Add mafs to eqtl dataframe
eqtl_df <- left_join(eqtl_df, maf_df, by="SNP")

gene_ids <- unique(eqtl_df$gene)
print(length(gene_ids))

egene_coloc <- function(region) {
    # region <- "CTSS"
    print(region)
    print(gene_ids %>% { which(. == region) })
    # Subset matrix output per gene
    eqtl_df_gene <- eqtl_df %>% filter(gene==region)
    print(dim(eqtl_df_gene))
    
    # Add gwas data and remove rows with NAs
    df <- left_join(eqtl_df_gene, gwas_df, by="SNP")
    df <- df %>% select(-n)
    dim(df)
    df <- na.omit(df)
    if (nrow(df)==0) {
        NULL 
        } else { 
         # Calculate variance
        df$se2 <- (df$se)^2
        df$beta_se2 <- (df$beta_se)^2
        df <- df %>% filter(se2!=0)
        df <- df %>% arrange(position)
        # print(head(df))  
        # Run coloc analysis
        my.res <- coloc.abf(dataset1=list(snp=df$SNP, beta=df$b, varbeta=df$se2, N=prop_df$size[prop_df$cond==condition], s=prop_df$prop[prop_df$cond==condition], type="cc"),
        dataset2=list(snp=df$SNP, beta=df$beta, varbeta=df$beta_se2, N=cno_df$cno[cno_df$ctype==celltype],type="quant"),
                        MAF=df$maf)
        # print(my.res)
        # Save coloc output as a rds file
        saveRDS(my.res, file=sprintf("%s/coloc_analysis/results/%s/%s/coloc_result_%s_%s_%s.rds", main_dir, condition, celltype, condition, celltype, region))

        # Print summary results
        gene_summary <- my.res$summary
        gene_summary$gene <- region
        gene_summary  
    }
} 

results <- lapply(gene_ids,egene_coloc)
results_df <- bind_rows(results)

# Save summary hypthesis results
fwrite(results_df, file=sprintf("%s/coloc_analysis/summary_results/summary_hypothesis_results_%s_%s.tsv", main_dir, condition, celltype),
    quote=F, sep="\t")
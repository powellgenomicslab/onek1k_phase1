#   ____________________________________________________________________________
#   Script information                                                      ####

# title: Run dynamic eQTL analysis
# author: Jose Alquicira Hernandez
# date: 2021-09-28
# description: Perform dynamic eQTL mapping in across naive-memory B cells

#   ____________________________________________________________________________
#   HPC details                                                             ####

# screen -S dyn
# qrsh -N dyn -l mem_requested=50G
# conda activate r-4.1.0

#   ____________________________________________________________________________
#   Import libraries                                                        ####

library("dsLib")
library("data.table")
library("lme4")
library("qvalue")


#   ____________________________________________________________________________
#   Set output                                                              ####

output <- set_output("2021-10-06", "dyn-eqtl")

#   ____________________________________________________________________________
#   Import data                                                             ####

path <- "/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/science_revision_Dec20/pseudotime-analysis_2021-09"

eqtl <- fread(here(path, "gene_snps_to_test_20211005.csv"))
data <- fread(here(path, "LMM_data.tsv"))

# Extract covariates
md <- data[, .(quantile, sampleid, pf1, pf2, pc1, pc2, pc3, pc4, pc5, pc6, sex, age)]

cols <- colnames(md)

# Extract genes and SNP columns
gene_snp <- data[, !..cols]

# Split eqtl data by SNP-gene pair
eqtl <- split(eqtl, seq_len(nrow(eqtl)))


# gene <- "CD37"
# snp <- "19:49878115_A"
# x <- eqtl[[1]]


fit_lmm <- function(x){
  
  gene <- x$GeneID
  snp <- x$snpid
  
  # Combine covariates with gene and SNP data
  test <- cbind(data[, ..gene], data[,..snp])
  setnames(test, new = c("gene", "snp"))
  test <- cbind(test, md)

  # Apply log transformation to gene expression
  test[, gene := log(gene + 1)]
  
  # Fit null model
  fit0 <- lmer(gene ~ snp + age + sex + pc1 + pc2 + pc3 + 
                 pc4 + pc5 + pc6 + pf1 + pf2 + 
                 (1 | sampleid) + quantile, data = test, REML = FALSE)
  
  # Fit augmented model
  fit <- lmer(gene ~ snp + age + sex + pc1 + pc2 + pc3 + 
                pc4 + pc5 + pc6 + pf1 + pf2  + 
                (1 | sampleid) + snp*quantile, data = test, REML = FALSE)
  
  # Perform LRT
  list(anova = anova(fit, fit0), fit = fit)
  
}

# Fit model for all SNP-gene pairs
inicio("Run LMM")
res <- lapply(eqtl, fit_lmm)
fin()


# Aggregate pairs
dyn_eqtl <- rbindlist(eqtl)

# Store and tidy results
#dyn_eqtl[, fit := lapply(res, function(x) x$fit)]
dyn_eqtl[, anova := lapply(res, function(x) x$anova)]

# Extract p-values and apply multiple testing correction
dyn_eqtl[, p.value := sapply(anova, function(x) x[2, "Pr(>Chisq)"])]
dyn_eqtl[, fdr := qvalue(p.value)$lfdr]


#   ____________________________________________________________________________
#   Run quadratic model                                                     ####

fit_lmm_sq <- function(x){
  
  gene <- x$GeneID
  snp <- x$snpid
  # Combine covariates with gene and SNP data
  test <- cbind(data[, ..gene], data[,..snp])
  setnames(test, new = c("gene", "snp"))
  test <- cbind(test, md)
  # Apply log transformation to gene expression
  test[, gene := log(gene + 1)]
  # Set up quadratic term
  test[, quantile_2 := quantile**2]
  
  # Fit null model
  fit0 <- lmer(gene ~ snp + age + sex + pc1 + pc2 + pc3 + 
                 pc4 + pc5 + pc6 + pf1 + pf2 + 
                 (1 | sampleid) + quantile + quantile_2, data = test, REML = FALSE)
  
  # Fit augmented model
  fit <- lmer(gene ~ snp + age + sex + pc1 + pc2 + pc3 + 
                pc4 + pc5 + pc6 + pf1 + pf2  + 
                (1 | sampleid) + snp*quantile + snp*quantile_2, data = test, REML = FALSE)
  
  # Perform LRT
  list(anova = anova(fit, fit0), fit = fit)
  
}

# Fit model for all SNP-gene pairs
inicio("Run squared LMM")
res_sq <- lapply(eqtl, fit_lmm_sq)
fin()

# Store and tidy results
#dyn_eqtl[, fit_sq := lapply(res, function(x) x$fit)]
dyn_eqtl[, anova_sq := lapply(res_sq, function(x) x$anova)]

# Extract p-values and apply multiple testing correction
dyn_eqtl[, p.value_sq := sapply(anova_sq, function(x) x[2, "Pr(>Chisq)"])]
dyn_eqtl[, fdr_sq := qvalue(p.value_sq)$lfdr]

# Create snp-gene pairs
dyn_eqtl[, id := paste(snpid, GeneID, sep = "-")]

fit <- lapply(res, function(x) x$fit)
fit_sq <- lapply(res_sq, function(x) x$fit)

dyn_eqtl[,singular := sapply(fit, isSingular)]
dyn_eqtl[,singular_sq := sapply(fit_sq, isSingular)]


#   ____________________________________________________________________________
#   Export data                                                             ####

saveRDS(dyn_eqtl, here(output, "dyn_eqtl_linear-sq_data.RDS"))

#   ____________________________________________________________________________
#   Session info                                                            ####

print_session(here(output))

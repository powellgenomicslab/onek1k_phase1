#   ____________________________________________________________________________
#   Script information                                                      ####

# title: Plots OneK1K dynamic eQTL results from B cells
# author: Jose Alquicira Hernandez
# date: 2021-04-16
# description: None

#   ____________________________________________________________________________
#   HPC details                                                             ####

# screen -S dyn-eqtl
# qrsh -N dyn-eqtl -l mem_requested=20G -pe smp 1 -q short.q
# conda activate r-4.1.0

#   ____________________________________________________________________________
#   Import libraries                                                        ####

# Primary
library("dsLib")
library("data.table")

# Secondary
library("stringr")
library("broom")
library("ComplexHeatmap")
library("LDlinkR")

#   ____________________________________________________________________________
#   Set output                                                              ####

output <- set_output("2021-10-06", "dyn_eqtls_summary")

#   ____________________________________________________________________________
#   Import data                                                             ####

path <- file.path("/",
                  "directflow", 
                  "SCCGGroupShare", 
                  "projects", 
                  "SeyhanYazar",
                  "onek1k",
                  "science_revision_Dec20",
                  "pseudotime-analysis_2021-09", 
                  "matrix_output")
quantile_dirs <- list.files(path, full.names = TRUE)

quantile_names <- sapply(str_split(quantile_dirs, "/"), 
                         function(x) str_split(x[length(x)], "_", simplify = TRUE)[1,1])

data <- lapply(quantile_dirs, fread)
names(data) <- quantile_names


# Combine data across quantiles
data <- rbindlist(data, idcol = "quantile")

# Index data by snp-gene pair
data[,id :=  paste(SNP, gene, sep = "-")]
setkey(data, quantile, gene)

#   ____________________________________________________________________________
#   Read original analysis                                                  ####

# Read SNP information with rsids
ref <- fread(here("data", "Revised_Supplementary_Tables.xlsx - Table.S10.tsv"), skip = 1)

# Load dynamic eQTL mapping
dyn <- readRDS(here("results", "2021-10-06_dyn-eqtl", "dyn_eqtl_linear-sq_data.RDS"))
setnames(dyn, "snpid", "SNP")

# Get assessed allele
dyn[, assessed_allele := lapply(str_split(SNP, "_"), "[", 2)]

# Retrieve significant dynamic eQTLs
dyn_sig <- dyn[fdr < 0.05 | fdr_sq < 0.05]

dyn_sig[fdr < 0.05 & !singular| fdr_sq < 0.05 & !singular_sq]
# 299/333

# Retrieve MatrixEQTL results for each dynamic eQTL
data_sig <- data[id %chin% dyn_sig$id, ]

stopifnot(all(dyn_sig$id %in% data_sig$id))

# Add dummy variable to specify pseudotime rank
data_sig[, rank := as.integer(str_remove(quantile, "Q"))]

# Add rs ids
setnames(ref, "SNP", "rsid")
ref[, SNP := paste0(Chromosome,":", Position, "_", `SNP assessed allele`)]
ref <- unique(ref[, .(SNP, rsid)])
data_sig <- merge(data_sig, ref, by = "SNP")
dyn_sig <- merge(dyn_sig, ref, by = "SNP")

# Nest data by SNP-gene pairs
data_sig <- data_sig[,.(data = list(.SD)), by = .(rsid, gene, id)]

table(data_sig[, sapply(data, nrow)])
# 3   4   5   6 
# 18  35  15 265 

# Merge dynamic mapping with MatrixEQTL results
dyn_sig <- merge(dyn_sig, data_sig[, .(id, data)], by = c("id"))

# Get number of quantiles where eGene is expressed
dyn_sig[, n_quantiles := sapply(data, nrow)]

# Define SNP-gene pairs
dyn_sig[, pair := paste0(rsid, "_", assessed_allele, "-", GeneID)]

# Extract betas and transform to matrix format
betas <- as.data.table(tidyr::unnest(dyn_sig[, .(pair, data)], "data"))
betas <- dcast(betas, pair ~ quantile, value.var = "beta")
betas <- betas[dyn_sig$pair]

stopifnot(all(dyn_sig$pair == betas$pair))


# Define linear and quadratic eQTLs
dyn_sig[, linear_significant := fifelse(fdr < 0.05, TRUE, FALSE)]
dyn_sig[, quadratic_significant := fifelse(fdr_sq < 0.05, TRUE, FALSE)]

# Scale betas 
all_scaled_betas <- betas[,.SD,.SDcols = patterns("Q")]
rownames(all_scaled_betas) <- betas$pair


#   ____________________________________________________________________________
#   Import GWAS data                                                        ####

# Define function to create batches
chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))

# Split eSNPs into 20 batches
snps <- chunk2(dyn_sig[, rsid], 20)

# Retrieve SNPs in LD to each eSNP
gwas <- lapply(snps, function(x){
  tryCatch(LDtrait(x,
                   pop = "GBR",
                   r2d = "r2",
                   r2d_threshold = 0.8,
                   win_size = 500000,
                   token = Sys.getenv("LDLINK_TOKEN"),
                   file = FALSE), error = function(e) NA)
})



gwas <- lapply(gwas, as.data.table)
gwas[[18]] <- NULL # Remove empty results
gwas <- rbindlist(gwas)
saveRDS(gwas, here(output, "gwas_r2_0.8.RDS"))
# gwas <- readRDS(here(output, "gwas_r2_0.8.RDS"))

autoimmune <-c("Type 1 diabetes", "Rheumatoid arthritis", "Inflammatory bowel disease", "Multiple sclerosis", "Crohn's disease", "Systemic lupus erythematosus", "Ankylosing spondylitis")

gwas_autoimmune <- gwas[GWAS_Trait %chin% autoimmune]
gwas_autoimmune <- unique(gwas_autoimmune[, .(snp = Query, GWAS_Trait)])

# Validate with GWAS catalog
gwas_catalog <- fread("data/gwas_catalog_11-10-21", sep = "\t", quote = "", fill = TRUE)
snps <- unlist(snps, use.names = FALSE)


gwas_catalog[, snp_list := sapply(SNPS, function(x) str_split(x, ","))]


gwas_catalog_sub <- gwas_catalog[sapply(snp_list, function(x) any(x %chin% snps)), ]
gwas_catalog_sub <- gwas_catalog_sub[, .(snp = SNPS, GWAS_Trait = `DISEASE/TRAIT`)]


gwas_catalog_sub <- unique(gwas_catalog_sub[GWAS_Trait %in% autoimmune, ])

gwas_autoimmune <- unique(rbind(gwas_autoimmune, gwas_catalog_sub))

#   ____________________________________________________________________________
#   Combine GWAS hits with dynamic eQTLs                                    ####

dyn_sig[, gwas := ifelse(rsid %in% gwas_autoimmune$snp, "True", "False")]


dyn_sig[, method := fifelse(linear_significant & !quadratic_significant , "Linear", 
                            fifelse(!linear_significant & quadratic_significant, "Quadratic", 
                                    "Linear & Quadratic"))]


dyn_sig[, method := factor(method, levels = c("Linear", "Quadratic", "Linear & Quadratic"))]
dyn_sig[, method := factor(method, levels = c("Linear", "Linear & Quadratic", "Quadratic"))]
dyn_sig[, gwas := factor(gwas, levels = c("True", "False"))]



ha <-  HeatmapAnnotation(Model = dyn_sig$method,
                         col = list(Model = c(Linear = "#38618c", 
                                              "Linear & Quadratic" = "#F6D65C",
                                              Quadratic = "#ff5964"
                                              )
                         )
)



hb <-  HeatmapAnnotation(`Autoimmune GWAS` = dyn_sig$gwas,
                         col = list(
                           `Autoimmune GWAS` = c(False = "#296A6C", 
                                                 True = "#F89423")
                         )
)


# Convert effect size data.table to matrix
all_scaled_betas_mat <- apply(all_scaled_betas, 1, scale)
colnames(all_scaled_betas_mat) <- dyn_sig$pair
rownames(all_scaled_betas_mat) <- names(all_scaled_betas)


# Define eucledian distance function to deal with NA values
robust_euclidean <- function(x, y){
  x[is.na(x)] <- 0
  y[is.na(y)] <- 0
  sqrt(sum((x - y) ^ 2))
}

postscript(here(output, "heatmap_betas.eps"), height = 3, width = 17)
Heatmap(all_scaled_betas_mat,
        col = rev(RColorBrewer::brewer.pal(11, "RdYlBu")),
        top_annotation = ha,
        bottom_annotation = hb,
        na_col = "gray",
        row_title = "Pseudotime",
        column_title = "eQTL",
        column_title_side = "bottom",
        cluster_rows  = FALSE, 
        clustering_distance_columns = robust_euclidean,
        column_names_gp = grid::gpar(fontsize = 2),
        column_split = dyn_sig$method, 
        column_gap = unit(0, "mm"), 
        row_names_side = "left",
        heatmap_legend_param = list(title = expression("Scaled "*beta*" values")))
dev.off()

# Add GWAS hits
gwas_autoimmune_summary <- gwas_autoimmune[, .(rsid = snp, GWAS_Trait)][, .(GWAS_Trait = paste0(.SD$GWAS_Trait, collapse = ",")), rsid]
dyn_sig <- merge(dyn_sig, gwas_autoimmune_summary, by = "rsid", all.x = TRUE)

dyn_sig_summary <- dyn_sig[, .(SNP = rsid, 
                               `Gene ID`= GeneID, 
                               ID = id, 
                               `Linear p-value` = p.value, 
                               `Linear FDR` = fdr, 
                               `Quadratic p-value` = p.value_sq, 
                               `Quadratic FDR` = fdr_sq, 
                               `Type` = method, 
                               `GWAS Trait` = GWAS_Trait,
                               pair)]

dyn_sig_summary <- merge(dyn_sig_summary, betas, by = "pair")
dyn_sig_summary[, pair := NULL]
dyn_sig_summary <- dyn_sig_summary[order(`Linear FDR`, `Quadratic FDR` )]

#   ____________________________________________________________________________
#   Export results                                                          ####

fwrite(dyn_sig_summary, here(output, "dyn-eqtl_summary.tsv"), sep = "\t")

#   ____________________________________________________________________________
#   Session info                                                            ####

print_session(here(output))


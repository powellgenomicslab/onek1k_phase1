# Script information ------------------------------------------------------

# title: Quality control of OneK1K data
# author: Jose Alquicira Hernandez
# date: 2020-10-01
# description: Removes outlier cells based on the expression of mitochondrial
# gene expression, number of UMIs, and genes detected per cell. Data is read 
# from multiple directories and integrated to create a single gene expression 
# matrix. A Seurat object is created and quality control is applied to remove 
# outlier cells based on number of genes, counts, and mitochondrial expression. 
# For the percentage of mitochondrial gene expression, Cells that deviate 2 negative 
# and 3 positive SDs are considered as outliers and removed. For the total
# number of UMIs and features (genes) per cell, only the low threshold is applied.
# Doublets are filtered out based on demuxlet results prior QC.


# Import libraries --------------------------------------------------------

# Primary
library("tidyverse")
library("dsLib")
library("data.table")
library("Seurat")
library("bestNormalize")

# Set output --------------------------------------------------------------

output <- set_output("2020-10-01", "QC_plots")

# Read data ---------------------------------------------------------------

# Input directory
input <- file.path("",
                   "share", 
                   "ScratchGeneral", 
                   "annsen", 
                   "data", 
                   "experimental_data", 
                   "CLEAN", 
                   "OneK1K_scRNA")

# Get file names for all pools
samples <- list.dirs(input, full.names = FALSE, recursive = FALSE) %>% 
  str_subset("OneK1K_scRNA_Sample")

# rename pools
pools <- samples %>% 
  str_split("_") %>% 
  map(~ .[3]) %>% 
  unlist() %>% 
  str_remove("Sample")

# Create helper function to read data
readData <- function(dirSample, pool){
  # Input
  dirSample <- file.path(input, dirSample, "outs", "filtered_gene_bc_matrices", "hg19")
  # Read file
  x <- Read10X(data.dir = dirSample, gene.column = 2)
  colnames(x) <- colnames(x) %>% str_split("-") %>% map_chr(1)
  # Assign pool id to barcode
  colnames(x) <- paste0(colnames(x),"-", pool)
  x
}


# Read gene expression matrices
inicio("Reading gene expression data")
data <- map2(samples, pools, readData)
fin()

# Get gene and cell counts
data %>% 
  map(dim) %>% 
  reduce(rbind) %>% 
  as.data.frame() %>% 
  set_names(c("genes", "cells")) %>% 
  mutate(pool = paste0("pool_",pools)) -> geneCellInfo

# Save data summary per pool
write_delim(geneCellInfo, path = here(output, "pool_genes_cells.txt"), delim = "\t")

# Get gene intersection from all pools
data %>% 
  map(row.names) %>% 
  reduce(intersect) -> geneIntersection

# Validate that all genes are shared among all pools and in the same order
cat("\nAre all genes shared among all pools?\n")
data %>% 
  map(row.names) %>% 
  map_lgl(~all(geneIntersection == .)) %>% 
  all() %>% 
  cat("\n")

# The intersection of all gene names is equal to the number of genes
# in each pool. Therefore, all genes are shared among all btches. 
# All gene names are in the same order in all matrices.

# Create pool info
inicio("Creating metadata")
pmap(list(geneCellInfo$pool, geneCellInfo$cells), rep) %>% 
  unlist() -> poolInfo

data %>% 
  map(colnames) %>% 
  unlist() %>% 
  data.frame(row.names = ., pool = poolInfo) -> poolInfo
fin()

# Merge cells into a single gene expression matrix
inicio("Integrating pools")
data <- do.call(cbind, data)
fin()

# Create Seurat object
inicio("Creating Seurat object")
data <- CreateSeuratObject(data, 
                           project = "onek1k", 
                           meta.data = poolInfo)
fin()

# Get percentage expression of mitochondrial and ribosomal genes
data <- PercentageFeatureSet(data, pattern = "^MT-", col.name = "percent.mt")


inicio("Saving pre-QC metadata")
saveRDS(data@meta.data, file = here(output, "preQC_metadata.RDS"))
raw_meta.data <- readRDS(file = here(output, "preQC_metadata.RDS"))
fin()


raw_meta.data %>% 
  as_tibble() %>% 
  select(-orig.ident) %>% 
  pivot_longer(-pool, names_to = "metric", values_to = "value") %>% 
  mutate(pool = str_remove(pool, "pool_")) %>% 
  mutate(pool = factor(pool, pool %>% unique() %>% as.numeric %>% sort())) %>%
  mutate(metric = case_when(metric == "nCount_RNA" ~ "Number of counts", 
                            metric == "nFeature_RNA" ~ "Number of features",
                            metric == "percent.mt" ~ "Percentage of mitochondrial gene expression")) %>% 
  ggplot() +
  aes(pool, value) +
  geom_boxplot(color = "gray20", outlier.size = 0.5) +
  ggtitle("Pre-QC") +
  xlab("Pool") +
  ylab("Value") +
  facet_wrap(~metric, ncol = 1, scales = "free_y") +
  rotate_x() +
  theme_pub() +
  theme(axis.text.y = element_text(size = 14), 
        plot.title = element_text(size = 20, hjust = 0.5),
        strip.background = element_rect(fill = "#F1F3BD"),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        strip.text.x = element_text(size = 17)) -> p

ggsave(here(output, "pre-qc.png"), p, width = 12, height = 15)
ggsave(here(output, "pre-qc.pdf"), p, width = 12, height = 15)



# Keep singlets -----------------------------------------------------------

# Keep singlets only as determined by demuxlet
inicio("Extracting singlets")
path <- file.path("data", "singlets", "barcodes_assigned_to_ppl.txt")
singlets <- fread(file = here(path), data.table = FALSE)
singlets$id <- paste(str_remove(singlets$BARCODE, "-.*$"), str_remove(singlets$batch, "sample"), sep = "-")
i <- colnames(data) %in% singlets$id
table(i)
data <- data[,i]
fin()


# Add individul information to metadata
inicio("Adding individual metadata")

i <- singlets$id %in% colnames(data)
singlets <- singlets[i,]

singlets %>% 
  `rownames<-`(NULL) %>% 
  column_to_rownames("id") %>% 
  select(-BARCODE, -batch) -> singlets

data <- AddMetaData(data, metadata = singlets, col.name = "individual")
levs <- fct_drop(data@meta.data$pool)

i <- levs %>% 
  levels() %>% 
  str_split("_") %>% 
  map(2) %>% 
  flatten_chr() %>% 
  as.integer() %>% 
  order()

data$pool <- factor(levs, levels = levels(levs)[i])

fin()


# QC ----------------------------------------------------------------------


inicio("Saving metadata after doublet removal")
saveRDS(data@meta.data, file = here(output, "singlets_metadata.RDS"))
fin()

# Extract metadata
meta.data <- data@meta.data

# Define QC metrics to evaluate to remove outliers
metricsLabels <- list("% Mitochondrial gene expression", "Number of UMIs", "Number of features")
metrics <- list("percent.mt", "nCount_RNA", "nFeature_RNA")
names(metrics) <- metrics


# object <- meta.data 
# group <- "pool"
# metric <- "percent.mt"
# 
# lowSD = 3
# highSD = 2

# Create function to get the cutoffs based on normalized data
getCutoffs <- function(object, group, metric, lowSD = 3, highSD = 2){
  
  # Normalize QC metric across pools
  res <- object %>% 
    split(.[[group]]) %>%                               # Split data by pool
    map(~pull(., {{metric}})) %>%                    # Extract pool information only
    map(quietly(orderNorm)) %>%                         # Transform to a normal distribution
    map("result")                                       # Extract results from `orderNorm`
  
  # Extract higher threshold based on `lowSD` standard deviations from the mean
  res %>% map(~mean(.$x.t) + sd(.$x.t) * highSD) -> higher
  # Extract lower threshold based on `highSD` standard deviations from the mean
  res %>% map(~mean(.$x.t) - sd(.$x.t) * lowSD) -> lower
  
  # Transform Z-scores to original distributions 
  higherTransform <- map2(res, higher, predict, inverse = TRUE)
  lowerTransform <- map2(res, lower, predict, inverse = TRUE)
  
  # Bind results and merge into a tibble
  higherTransform <- bind_cols(higherTransform) %>% 
    pivot_longer(everything(), names_to = "group", values_to = "higher")
  lowerTransform <- bind_cols(lowerTransform) %>% 
    pivot_longer(everything(), names_to = "group", values_to = "lower")
  
  # Combine lower and higher thresholds into a single tibble
  cutoffs <- full_join(higherTransform, lowerTransform, by = "group") %>% 
    pivot_longer(-group, names_to = "type", values_to = "metric")
  
  # Combine lower and higher Z-score thresholds into a single tibble
  higher <- bind_cols(higher) %>%  pivot_longer(everything(), names_to = "group", values_to = "higher")
  lower <- bind_cols(lower) %>%  pivot_longer(everything(), names_to = "group", values_to = "lower")
  cutoffsTrans <- full_join(higher, lower, by = "group") %>% 
    pivot_longer(-group, names_to = "type", values_to = "metric")
  
  # Extract original and Z-score distributions
  res %>% 
    map(~.[c("x", "x.t")]) %>%                                # Extract original and transformed distributions
    map(as.data.frame) %>%                                    # Bind into a data.frame
    map(`names<-`, c(metric, paste0(metric, "_norm"))) %>%    # Rename columns
    bind_rows(.id = group) -> res                             # Combine distributions across pools
  
  # Returns:
  # - cutoffs: A data.frame containing the group (pool), type of threhold (higher or lower), and the value
  # - normCutoffs: Normalized cutoffs
  # - trans: A data.frame including the original and the normalized distribution by group (pool)
  list(cutoffs = as.data.frame(cutoffs), normCutoffs = as.data.frame(cutoffsTrans), trans = res)
}

# Get cutoffs
inicio("Obtaining cutoffs across all pools")
cutoffs <- map(metrics, ~getCutoffs(meta.data, "pool", .))
names(cutoffs) <- names(metrics)
fin()

# Filter data -------------------------------------------------------------

# Get a list of cutoffs by pool and metric
# Sets Inf higher cutoffs for "nCount_RNA" and "nFeature_RNA"

setHighCutoff <- function(x, inf = TRUE){
  
  x %>% 
    map("cutoffs") %>% 
    map(spread, key = "type", value = "metric", 2) -> x
  
  if(inf){
    x <- x %>% map(~mutate(.x, higher = Inf))
  }
  
  x %>% 
    map(~split(., .["group"])) 
}


cutoffsBypool <- setHighCutoff(cutoffs[c("nCount_RNA", "nFeature_RNA")])
cutoffsBypool <- append(cutoffsBypool, setHighCutoff(cutoffs[c("percent.mt")], inf = FALSE))


findOutliers <- function(cutoff, metric, group){
  # Iterate for each value in `group` variable and filter data based on
  # corresponding filters. Data for each value in group is stored in `res` list
  res <- map(names(cutoff), function(val){
    
    meta.data %>% 
      rownames_to_column("barcode") %>% 
      filter(UQ(as.name(group)) == !!quo(val)) -> subMeta
    
    lower <- subMeta %>% filter(!!sym(metric) < cutoff[[val]]$lower)
    higher <-  subMeta %>%  filter(!!sym(metric) > cutoff[[val]]$higher)
    
    list(lower = lower, higher = higher) %>% 
      bind_rows(.id = "type")
  })  
  
  # Combine 
  res %>% 
    map(bind_rows)
  
}

# Get outliers
inicio("Extracting outliers")
res <- map2(cutoffsBypool, names(cutoffsBypool), findOutliers, group = "pool")
fin()

outliers <- res %>% 
  map(bind_rows) %>% 
  bind_rows(.id = "origin") %>% 
  distinct(barcode, .keep_all = TRUE)

# Get number of outliers by QC metric and type of cutoff
outliers %>% 
  group_by(type, origin) %>% 
  summarise(n = n())

# Get number of outliers by QC metric
outliers %>% 
  group_by(origin) %>% 
  summarise(n = n())
 

i <- !colnames(data) %in% outliers$barcode
table(i)

# Filter data
inicio("Removing outliers")
data <- data[,i]
fin()

inicio("Saving metadata after QC")
saveRDS(data[[]], file = here(output, "post-qc_metadata.RDS"))
fin()

data[[]] %>% 
  as_tibble() %>% 
  select(-orig.ident, -individual) %>% 
  pivot_longer(-pool, names_to = "metric", values_to = "value") %>% 
  mutate(pool = str_remove(pool, "pool_")) %>% 
  mutate(pool = factor(pool, pool %>% unique() %>% as.numeric %>% sort())) %>%
  mutate(metric = case_when(metric == "nCount_RNA" ~ "Number of counts", 
                            metric == "nFeature_RNA" ~ "Number of features",
                            metric == "percent.mt" ~ "Percentage of mitochondrial gene expression")) %>% 
  ggplot() +
  aes(pool, value) +
  geom_boxplot(outlier.size = 0.5) +
  ggtitle("Post-QC") +
  xlab("Pool") +
  ylab("Value") +
  facet_wrap(~metric, ncol = 1, scales = "free_y") +
  rotate_x() +
  theme_pub() +
  theme(axis.text.y = element_text(size = 14), 
        plot.title = element_text(size = 20, hjust = 0.5),
        strip.background = element_rect(fill = "#F1F3BD"),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        strip.text.x = element_text(size = 17)) -> p2

ggsave(here(output, "post-qc.png"), p2, width = 12, height = 15)
ggsave(here(output, "post-qc.pdf"), p2, width = 12, height = 15)


library("gridExtra")

ggsave(
  filename = here(output, "qc.pdf"),
  plot = marrangeGrob(list(p, p2), nrow = 1, ncol = 1), 
  height = 15, width = 12
)




print_session(here(output))


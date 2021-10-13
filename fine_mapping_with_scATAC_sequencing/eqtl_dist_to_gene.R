# Script information ------------------------------------------------------

# title: Calculates distance from eQTLs to gene body
# author: Jose Alquicira Hernandez
# date: 2021-09-23
# description: None

# Import libraries --------------------------------------------------------

# Primary
library("data.table")
library("dsLib")

# Secondary
library("ggplot2")
library("stringr")
library("rtracklayer")

# Set output --------------------------------------------------------------

output <- set_output("2021-09-28", "eqtl_distance_to_gene")

#   ____________________________________________________________________________
#   Helper functions                                                        ####


rename_cells <- function(cell_type){
  # Map labels
  dplyr::recode(cell_type,
                `CD4effCM`      = "CD4 ET",
                `CD4all`      = "CD4 NC",
                `CD4TGFbStim`       = "CD4 SOX4",
                `CD8all`        = "CD8 NC",
                `CD8eff` = "CD8 ET",
                `CD8unknown`      = "CD8 S100B",
                `NKmat`                = "NK",
                `NKact`                = "NK R",
                `Bmem`    = "B Mem",
                `BimmNaive`    = "B IN",
                `Plasma`             = "Plasma",
                `MonoC`          = "Mono C",
                `MonoNC`        = "Mono NC",
                `DC`          = "DC"
  )
}

# Map labels
rename_cells2 <- function(cell_type){
  
  
  dplyr::recode(cell_type,
                `CD4+ KLRB1+ T cell`      = "CD4 ET",
                `CD4+ KLRB1- T cell`      = "CD4 NC",
                `CD4+ SOX4+ T cell`       = "CD4 SOX4",
                `CD8+ LTB+ T cell`        = "CD8 NC",
                `CD8+ GNLY+ NKG7+ T cell` = "CD8 ET",
                `CD8+ S100B+ T cell`      = "CD8 S100B",
                `XCL1- NK`                = "NK",
                `XCL1+ NK`                = "NK R",
                `TCL1A- FCER2- B cell`    = "B Mem",
                `TCL1A+ FCER2+ B cell`    = "B IN",
                `IgJ+ B cell`             = "Plasma",
                `Monocyte CD14+`          = "Mono C",
                `Monocyte FCGR3A+`        = "Mono NC",
                `Dendritic cell`          = "DC"
  )
  
}

pal <- c("#882E72", "#B178A6", "#D6C1DE", "#1965B0", "#5289C7", "#7BAFDE", 
         "#4EB265", "#90C987", "#CAE0AB", "#F7EE55", "#F6C141", "#F1932D", 
         "#E8601C", "#DC050C")

level_order <- c("CD4 NC", "CD4 ET", "CD4 SOX4", 
                 "CD8 ET", "CD8 NC","CD8 S100B", 
                 "NK", "NK R", 
                 "Plasma", "B Mem", "B IN",  
                 "Mono C", "Mono NC", 
                 "DC")


#   ____________________________________________________________________________
#   Import gene annotation                                                  ####

anno <- import("data/refdata-cellranger-hg19-3.0.0_genes.gtf")
anno <- as.data.table(as.data.frame(anno))

# Extract autosomal genes
anno <- anno[type == "gene" & seqnames %in% 1:22, 
             .(gene_id, gene_name, seqnames, start, end, strand, gene_biotype)]
# Refactor chromosomes
anno <- anno[, seqnames := factor(seqnames, levels = sort(as.integer(unique(as.character(seqnames)))))]
setkey(anno, seqnames, gene_id)

#   ____________________________________________________________________________
#   Import eQTL results                                                     ####


# Read data
eqtl_raw <- fread(input = here("data", "Revised_Supplementary_Tables.xlsx - Table.S10.tsv"), skip = 2)
eqtl <- eqtl_raw[, .(cell_type = `Cell type`, gene = `Gene Ensembl ID`, rsID = SNP, Chromosome, Position)]
eqtl[, pair := paste(cell_type, rsID, gene, sep = "_")]

eqtl <- split(eqtl, seq_len(nrow(eqtl)))

# Get distances between eSNP and gene body
distances <- sapply(eqtl, function(x){
  gene_pos <- anno[seqnames == x$Chromosome & gene_id == x$gene, .(start, end)]
  
  start_dist <- gene_pos$start - x$Position
  end_dist <- gene_pos$end - x$Position
  
  min(c(abs(start_dist), abs(end_dist)))
})

eqtl_raw$distances <- distances

eqtl_raw[, `Cell type` := factor(`Cell type`, levels = level_order)]

# Extract number of eQTLS per cell type
egenes <- eqtl_raw[, .N, .(`Gene ID`, `Cell type`)]
setnames(egenes, "Gene ID", "geneid")
setkey(egenes, geneid, `Cell type`)

#   ____________________________________________________________________________
#   Create controls                                                         ####


##  ............................................................................
##  Import data                                                             ####

files <- list.files("data/eqtl_inputs", recursive = TRUE)
files <- str_subset(files, ".md", negate = TRUE)
cell_type_list <- sapply(str_split(files, "/"), "[", 1)
cell_type_list <- rename_cells(cell_type_list)

# For each chromosome-cell type pair, sample N number of SNPs for each eGene
# (i.e. same number of eQTLs matched by gene)

read_snps <- function(file, cell_type){
  cat(file, "\n", sep = "")
  
  # Read cis-SNP data for each cell-type/chromosome
  data <- fread(file = here("data", "eqtl_inputs", file),
                colClasses = c("integer", "integer", "character", "character"))   
  
  setkey(data, geneid)
  
  egenes_sub <- egenes[`Cell type` == cell_type, !"Cell type"]
  
  data <- data[geneid %chin% egenes_sub$geneid]
  
  if(!nrow(data)){
    message("No eGenes for ", cell_type, " in this chromosome file")
    return(NA)
  }
  
  data <- merge(data, egenes_sub, by = "geneid")
  
  data <- data[, {
    n <- unique(.SD$N)
    set.seed(66)
    if(nrow(.SD) >= n){
      i <- sample(nrow(.SD), n)
    }else{
      i <- sample(nrow(.SD), n, replace = TRUE)
    }
    .SD[i, ]
  }, geneid]
  
  
  # Add cell type info to results
  data$cell_type <- cell_type
  
  data
}

inicio("Read cis-SNPs")
snps <- mapply(read_snps, files, cell_type_list, SIMPLIFY = FALSE)
fin()

# Remove instances where no eSNPs are contained in a given chromosome
snps <- snps[!is.na(snps)]

# Aggregate results
snps <- rbindlist(snps)
setkey(snps, snpid, chr, pos, cell_type)
saveRDS(snps, here(output, "snps.RDS"))
#snps <- readRDS(here(output, "snps.RDS"))

# Verify all genes from cis mapping are in the annotation data
all(snps$geneid %in% anno$gene_name)

# Merge random cic SNPs with annotation data
snps <- merge(snps, 
              anno[, .(seqnames, start, end, geneid = gene_name, ensembl_id = gene_id)], 
              by = "geneid", allow.cartesian = TRUE)

# Calculate distances to start end end of egene
snps[, `:=`(dist_start = abs(start - pos), dist_end = abs(end - pos))]

# Get minimum distance
snps$distances <- snps[, .(distances = min(dist_start, dist_end)), by = seq_len(nrow(snps))][, distances]

# Combine eQTLs and random SNPs
all_snps <- rbind(eqtl_raw[, .(cell_type = `Cell type`, distances, set = "eqtl")],
                  snps[, .(cell_type, distances, set = "control")])

# Calculate logarithmic distances
all_snps[, distances_log := log10(distances + 1)]
all_snps[, set := factor(set, levels = c("eqtl", "control"), 
                         labels = c("eQTL", "Random SNPs"))]


# Nest data by cell type
all_snps_nest <- all_snps[, .(data = list(.SD)), cell_type]

# Test for difference between eQTLS and random SNPs
all_snps_nest[, p.value := sapply(data, function(x){
  wilcox.test(x[set == "eQTL", distances], x[set != "eQTL", distances], alternative = "less")$p.value
})]

# Perform multiple testing correction
all_snps_nest[, p.adj := p.adjust(p.value, method = "fdr")]

#   ____________________________________________________________________________
#   Plot data                                                               ####

exponents <- c(0, 2, 4, 6)

p <- ggplot(all_snps, aes(cell_type, distances_log, fill = set)) +
  geom_boxplot(size = 0.5, outlier.size = 0.3, notch = TRUE) +
  xlab("Cell type") +
  ylab("Distance to gene body (bp)") +
  scale_fill_manual(values = c("#FFC600", "#A1D7F9"), 
                    limits = c("eQTL", "Random SNPs"), 
                    name = "Loci") +
  scale_y_continuous(labels = format(10 ** exponents, 
                                     scientific = FALSE, 
                                     big.mark = ","), 
                     breaks = exponents) +
  annotation_logticks(sides = "l") +
  rotate_x() + 
  theme_pub()

#   ____________________________________________________________________________
#   Save results                                                            ####

ggsave(here(output, "distances.png"), p, height = 4.5, width = 9.5)
ggsave(here(output, "distances.pdf"), p, height = 4.5, width = 9.5)

# Session info ------------------------------------------------------------

print_session(here(output))

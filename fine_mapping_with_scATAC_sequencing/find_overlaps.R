# Script information ------------------------------------------------------

# title: Calculates distance from eQTLs and random cis-SNPs to ATAC-seq
# author: Jose Alquicira Hernandez
# date: 2021-03-08
# description: None

# Import libraries --------------------------------------------------------

# Primary
library("tidyverse")
library("dsLib")

# BiocManager::install(c("ggbio", "biovizBase", 
#                        "GenomeInfoDb", "EnsDb.Hsapiens.v75"))
# 
# install.packages("Signac")


# Secondary
library("Seurat")
library("Signac")
library("GenomeInfoDb")
library("EnsDb.Hsapiens.v75")
library("ggplot2")
library("patchwork")
library("data.table")

set.seed(1234)

# Set output --------------------------------------------------------------

output <- set_output("2021-03-08", "atac_anno")

#   ____________________________________________________________________________
#   Helper functions                                                        ####


rename_cells <- function(cell_type){
  # Map labels
  recode(cell_type,
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
  
  
  recode(cell_type,
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

get_distance <- function(data, start_offset = 0, end_offset = 0){
  
  map2(data$eqtl_range, data$seurat, function(x, y, start_offset, end_offset){
    
    # Get cell type label
    cell_type <- unique(y$predicted.id) %>% as.character()
    
    # Format cell type label to match with peak results
    cell_type <- str_replace(cell_type, " ", "_")
    
    # Extract peaks from all cell types (this information is the same in
    # all seurat objects for each cell type)
    gr <- y@assays$peaks@misc$cell_type_peaks
    
    # Extract peaks present in cell type
    i <- grepl(cell_type, gr$peak_called_in)
    gr <- gr[i, ]
    
    # Incorporate offsets if necessary
    start(x) <- start(x) - start_offset
    end(x) <- end(x) + end_offset
    
    # Calculate nearest distance to peak for each eQTL
    distance <- distanceToNearest(x, gr)
    
    # Extract distances
    distances <- mcols(distance)$distance
    
    # Return snp, distance to peak, and analyzed cell type
    tibble(snp = names(x@ranges), distances, cell_type) %>% 
      arrange(distances) 
    
  }, start_offset, end_offset)
  
}


find_overlaps <- function(data, start_offset = 0, end_offset = 0){
  
  map2(data$eqtl_range, data$seurat, function(x, y, start_offset, end_offset){
    
    
    cell_type_subject <- unique(y$predicted.id)
    cell_type_subject <- str_replace(cell_type_subject, " ", "_")
    
    gr <- y@assays$peaks@misc$cell_type_peaks
    
    i <- grepl(cell_type, gr$peak_called_in)
    gr <- gr[i, ]
    
    i <- names(x) %>% duplicated() %>% `!`()
    x <- x[i, ]
    
    start(x) <- start(x) - start_offset
    end(x) <- end(x) + end_offset
    
    x$cell_type_query <- x$cell_type
    x$cell_type <- NULL
    
    hits <- findOverlapPairs(x, gr)
    
    res <- S4Vectors::first(hits) %>% as.data.frame() %>% as_tibble() %>% set_names("query_" %p% names(.))
    res_subject <- S4Vectors::second(hits) %>% as.data.frame() %>% as_tibble() %>% set_names("subject_" %p% names(.))
    
    res_subject$subject_cell_type <- str_replace_all(cell_type_subject, "_", " ")
    
    
    bind_cols(res, res_subject)
    
  }, start_offset, end_offset)
  
}



randomize_snps <- function(seed){
  
  # Get number of unique eQTLs/SNPs per cell type
  n <- eqtl$data %>% 
    map(tally) %>%
    map_int(~ pull(., n))
  
  # Create distribution of seeds for all cell types
  # using input seed
  set.seed(seed)
  seeds <- sample(1:1e4, length(n))
  
  # Subsets n cis-SNPs from each cell type
  randomize <- function(size, seed_i, cell_type_i){
    
    # Subsets SNPs for target cell type 
    sub_snps <- snps[cell_type == cell_type_i]
    
    # Subset SNPs
    set.seed(seed_i)
    i <- sample(seq_len(nrow(sub_snps)), size)
    sub_snps <- sub_snps[i, ]
    
    # Converts output
    sub_snps %>% 
      as.data.frame() %>% 
      as_tibble()
  }
  
  # Get random samples per cell type with different sample size and seed values
  # by cell type
  sub_snps <- pmap(list(n, seeds, eqtl$cell_type), randomize)
  
  # Create template to store results based on eqtl tibble
  control <- eqtl
  
  # Create GenomicRanges for each cell type
  sub_snps <- map(sub_snps, ~ {
    
    GRanges(IRanges(start = .$pos, names = .$snpid), 
            seqnames = "chr" %p% .$chr)
  }
  )
  
  # Store ranges in `eqtl_range` column
  control$eqtl_range <- sub_snps
  
  control
}


pal <- c("#882E72",
         "#1965B0", "#5289C7", "#7BAFDE", 
         "#4EB265", "#90C987", 
         "#F7EE55", "#F6C141", 
         "#F1932D", "#E8601C", 
         "#DC050C")


level_order <- c("CD4 NC",
                 "CD8 ET", "CD8 NC","CD8 S100B", 
                 "NK", "NK R", 
                 "B Mem", "B IN",  
                 "Mono C", "Mono NC", 
                 "DC")
#   ____________________________________________________________________________
#   Read data                                                               ####

atac <- readRDS(here(output, "atac.RDS"))

atac <- atac[, atac$predicted.id != "CD4+ KLRB1+ T cell"]

atac$predicted.id <- rename_cells2(atac$predicted.id) 
atac$predicted.id <- factor(atac$predicted.id, levels = level_order) 


DimPlot(atac, group.by = "predicted.id", label = TRUE, repel = TRUE, pt.size = 0.8) +
  scale_color_manual(values = alpha(pal, 0.7)) +
  guides(color = guide_legend(override.aes = list(shape = 15, size = 6, alpha = 1))) +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  ggtitle("") +
  theme_pub()

ggsave(here(output, "umap.png"), p, height = 5.5, width = 7)


DefaultAssay(atac) <- "peaks"

#   ____________________________________________________________________________
#   Call peaks by cell type                                                 ####

# Rename cells
atac$predicted.id <- rename_cells2(atac$predicted.id)

# Call peaks for each cell type

peaks <- CallPeaks(
  object = atac,
  group.by = "predicted.id",
  macs2.path = "/Users/josealquicirahernandez/opt/miniconda3/bin/macs2"
)

saveRDS(peaks, here(output, "peaks_cell_type.RDS"))
# peaks <- readRDS(here(output, "peaks_cell_type.RDS"))

Idents(atac) <- "predicted.id"

# Add cell-type-specific peaks to "peaks" assay
atac@assays$peaks@misc$cell_type_peaks <- peaks

# Split data by cell type
atac <- SplitObject(atac, split.by = "predicted.id")

# Get number of peaks called in each cell type
peaks_tab <- peaks %>% as.data.frame() %>% as_tibble()
peaks_tab$peak_called_in %>% 
  str_split(",") %>% 
  flatten_chr() %>% 
  table() %>% 
  as.data.frame() %>% 
  arrange(desc(Freq)) %>% 
  set_names(c("cell_type", "N")) %>% 
  dplyr::filter(cell_type != "CD4_ET") %>% 
  pull(N) %>% 
  mean()
  


# # Assign 
# atac <- map(atac, function(x){
#   
#   # Get average number of UMIs per peak across all cells
#   data <- GetAssayData(x, slot = "data", assay = "peaks")
#   means <- Matrix::rowMeans(data)
#   
#   # Extract ranges data
#   range_data <- GetAssayData(x, slot = "ranges", assay = "peaks")
#   
#   
#   chr <- range_data %>% seqnames() %>% as.character()
#   range_names <- paste(chr, start(range_data), end(range_data), sep = "-")
#   
#   if(!all(names(means) == range_names)) stop("Range mistmatch")
#   
#   range_data$means <- means
#   
#   x@assays$peaks@ranges <- range_data
#   x
#   
# })



#   ____________________________________________________________________________
#   Import eQTL results                                                     ####


# Read data
eqtl_raw <- read_tsv(here("data", "eqtl.tsv"))

# Select SNP information and remove duplicates (due to SNP-gene pairs)
eqtl <- eqtl_raw %>% 
  dplyr::select(cell_type, rsID, Chromosome, Position) %>% 
  distinct()

# Nest eQTL info by cell type
eqtl <- eqtl %>% group_by(cell_type) %>% nest()

# Create ranges for each eQTL 
eqtl <- eqtl %>% 
  mutate(eqtl_range = map(data, ~{
    GRanges(IRanges(start = .$Position, names = .$rsID), 
            seqnames = "chr" %p% .$Chromosome)
  }))

# Gather scATAC-seq
atac <- tibble(cell_type = names(atac), seurat = atac)
eqtl <- right_join(eqtl, atac, by = "cell_type")


eqtl_res <- get_distance(eqtl) %>% bind_rows()


#   ____________________________________________________________________________
#   Create controls                                                         ####


##  ............................................................................
##  Import data                                                             ####

files <- list.files("data/eqtl_inputs", recursive = TRUE)
files <- files %>% str_subset(".md", negate = TRUE)
cell_type_list <- str_split(files, "/") %>% map_chr(1)
cell_type_list <- rename_cells(cell_type_list)


read_snps <- function(file, cell_type){
  # Read cis-SNP data for each cell-type/chromosome
  data <- fread(file = here("data", "eqtl_inputs", file),
                colClasses = c("integer", "integer", "character", "character"))   
  
  # Extract SNP info only
  data <- data[,c("chr", "pos", "snpid")]  
  
  # Remove duplicates (due to snp-gene testing)
  data <- unique(data)  
  
  # Add cell type info to results
  data$cell_type <- cell_type
  data
}

inicio("Read cis-SNPs")
snps <- map2(files, cell_type_list, read_snps)
fin()

snps <- rbindlist(snps)
setkey(snps, snpid)

saveRDS(snps, here(output, "snps.RDS"))
#snps <- readRDS(here(output, "snps.RDS"))


##  ............................................................................
##  Randomize tested SNPs                                                   ####

set.seed(66)
seeds <- sample(1e6, 10)

control_ranges <- map(seeds, randomize_snps)

controls_res <- map(control_ranges, get_distance) %>% 
  bind_rows()


#   ____________________________________________________________________________
#   Find overlaps with eQTLs                                                ####


# Combine results (eQTL and null distances)
controls_res$set <- "control"
eqtl_res$set <- "eqtl"
res <- bind_rows(controls_res, eqtl_res)

# Format cell type labels
res$cell_type <- str_replace(res$cell_type, "_", " ")

# Set order of cell types for plotting
res  %>% 
  mutate(cell_type = factor(cell_type, level_order)) -> res

res %>% 
  group_by(cell_type, set) %>% 
  summarise(mean_dist = mean(distances)) %>% 
  group_by(set) %>% 
  summarise(max(mean_dist), min(mean_dist))


n <- eqtl$data %>% 
  map(~ dplyr::filter(., !duplicated(rsID))) %>% 
  map(tally) %>%
  map_int(~ pull(., n))



summ <- tibble(cell_type = eqtl$cell_type, n_eqtl = n, n = eqtl$seurat %>% map_int(ncol))

res <- left_join(res, summ, by = "cell_type")

res <- res %>% 
  mutate(label = paste0(cell_type, " | eQTLs = ",  n_eqtl, " | Cells = ", n))

labels <- res$label %>% unique() 

labels_name <- labels %>% str_split(" \\| ") %>% map_chr(1)

i <- map(level_order, `==`, labels_name) %>% map(which) %>% unlist()

res$label <- factor(res$label, labels[i])


#   ____________________________________________________________________________
#   Plot results                                                            ####


res %>% 
  group_by(cell_type, set) %>% 
  summarise(mean_dist = mean(distances)) %>% 
  arrange(set) %>%
  as.data.frame() %>% 
  dplyr::filter(set == "control") %>% 
  ungroup() %>% 
  summarise_at(vars(mean_dist), list(mean = mean, max = max, min = min)) %>% 
  t() %>% 
  apply(2, round)

res_nest <- res %>% 
  group_by(cell_type) %>% 
  nest() 

test_diff <- function(x){
  dist_cont <- x %>% 
    dplyr::filter(set == "control") %>% pull(distances) 
  
  dist_eqtl <- x %>% 
    dplyr::filter(set == "eqtl") %>% pull(distances) 
  
  #wilcox.test(dist_eqtl, dist_cont)
  t.test(dist_eqtl, dist_cont)$p.value
  
}


res_nest <- res_nest %>% 
  mutate(t_test = map_dbl(data, test_diff)) %>% 
  dplyr::select(cell_type, p_value = t_test)


res_nest


exponents <- c(0, 2, 4, 6)

res %>% 
  mutate(set = factor(set, levels = c("eqtl", "control"), 
                      labels = c("eQTL", "Random SNPs"))) %>% 
  dplyr::filter(cell_type != "CD4 ET") %>% 
  mutate(cell_type = factor(cell_type, level_order)) %>% 
  mutate(dist_log = log10(distances + 1)) %>% 
  ggplot() +
  aes(cell_type, dist_log, fill = set) +
  geom_boxplot(size = 0.5, outlier.size = 0.3, notch = TRUE) +
  xlab("Cell type") +
  ylab("eQTL distance to peak (bp)") +
  scale_y_continuous(labels = format(10 ** exponents, 
                                     scientific = FALSE, 
                                     big.mark = ","), 
                     breaks = exponents) +
  theme_pub() +
  annotation_logticks(sides = "l") +
  rotate_x() + 
  scale_fill_manual(values = c("#FFC600", "#A1D7F9"), 
                    limits = c("eQTL", "Random SNPs"), 
                    name = "Loci") -> p


ggsave(here(output, "distances.png"), p, height = 4.5, width = 9.5)
ggsave(here(output, "distances.pdf"), p, height = 4.5, width = 9.5)



# Session info ------------------------------------------------------------

print_session(here(output))



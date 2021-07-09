# Script information ------------------------------------------------------

# title: Process ATAC-seq dataset from 10x genomics website
# author: Jose Alquicira Hernandez
# date: 2021-03-08

# Import libraries --------------------------------------------------------

# Primary
library("tidyverse")
library("dsLib")

# BiocManager::install(c("ggbio", "biovizBase", 
#                        "GenomeInfoDb", "EnsDb.Hsapiens.v75"))

# install.packages("Signac")


# Secondary
library("Seurat")
library("Signac")
library("GenomeInfoDb")
library("EnsDb.Hsapiens.v75")
library("ggplot2")
library("patchwork")
set.seed(1234)

# Set output --------------------------------------------------------------

output <- set_output("2021-03-13", "atac_nextgem")

# Read data ---------------------------------------------------------------

counts <- Read10X_h5(filename = "data/atac-seq_gem1.1/atac_pbmc_10k_nextgem_filtered_peak_bc_matrix.h5")
metadata <- read.csv(
  file = "data/atac-seq_gem1.1/atac_pbmc_10k_nextgem_singlecell.csv",
  header = TRUE,
  row.names = 1
)

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'hg19',
  fragments = 'data/atac-seq_gem1.1/atac_pbmc_10k_nextgem_fragments.tsv.gz',
  min.cells = 10,
  min.features = 200
)

pbmc <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)



# Add genomic annotations -------------------------------------------------

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)

# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg19"

# add the gene information to the object
Annotation(pbmc) <- annotations


# QC ----------------------------------------------------------------------

# compute nucleosome signal score per cell
pbmc <- NucleosomeSignal(pbmc)

# compute TSS enrichment score per cell
pbmc <- TSSEnrichment(pbmc, fast = FALSE)

# add blacklist ratio and fraction of reads in peaks
pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100
pbmc$blacklist_ratio <- pbmc$blacklist_region_fragments / pbmc$peak_region_fragments

pbmc$high.tss <- ifelse(pbmc$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(pbmc, group.by = 'high.tss') + NoLegend()

pbmc$nucleosome_group <- ifelse(pbmc$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = pbmc, group.by = 'nucleosome_group')


pbmc <- subset(
  x = pbmc,
  subset = peak_region_fragments > 3000 &
    peak_region_fragments < 20000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)



# Dimensionality reduction ------------------------------------------------

pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc <- RunSVD(pbmc)

pbmc <- RunUMAP(pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindNeighbors(pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindClusters(pbmc, verbose = FALSE, algorithm = 3)
DimPlot(pbmc, label = TRUE) + NoLegend()



# Gene activity matrix ----------------------------------------------------

gene.activities <- GeneActivity(pbmc)

# add the gene activity matrix to the Seurat object as a new assay and normalize it
pbmc[['RNA']] <- CreateAssayObject(counts = gene.activities)
pbmc <- NormalizeData(
  object = pbmc,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(pbmc$nCount_RNA)
)

DefaultAssay(pbmc) <- 'RNA'

FeaturePlot(
  object = pbmc,
  features = c('MS4A1', 'CD3D', 'LEF1', 'NKG7', 'TREM1', 'LYZ'),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)


# Save data ---------------------------------------------------------------


saveRDS(pbmc, here(output, "atac.RDS"))


# Session info ------------------------------------------------------------

print_session(here(output))


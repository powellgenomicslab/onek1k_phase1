# Script information ------------------------------------------------------

# title: Normalize OneK1k data
# author: Jose Alquicira Hernandez
# date: 2019/07/10
# description: Data normalization is performed using the new SCTransform method
# from Seurat.
# See:
# - https://satijalab.org/seurat/v3.0/sctransform_vignette.html

# Import libraries --------------------------------------------------------

# Primary
library("tidyverse")
library("here")
library("dsLib")

# Secondary
library("Seurat")


# Set output --------------------------------------------------------------

output <- set_output("2019-08-23", "norm_sctransform")

# Read data ---------------------------------------------------------------

inicio("Reading gene expression data")
input <- file.path("results", "2019-06-16_QC")
data <- readRDS(here(input, "QC.RDS"))
fin()

# Normalize data ----------------------------------------------------------

inicio("Normalizing data using SCTransform")
options(future.globals.maxSize = 1024**3 * 2500)
data <- SCTransform(object = data, vars.to.regress = c("percent.mt", "pool"), conserve.memory = TRUE)
fin()

# Save output -------------------------------------------------------------

inicio("Saving normalized data")
saveRDS(data, file = here(output, "norm.RDS"))
fin()

# Session info ------------------------------------------------------------
print_session(output)

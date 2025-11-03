suppressPackageStartupMessages({
    library(methylTFR)
    library(tools)
    library(methylTFRAnnotationHg38)
    library(dplyr)
    library(GenomicRanges)
    library(RnBeads)
})
set.seed(42)
source("/icbb/projects/nitschre/methylTFR/scripts/memoryTcells/plots.R")
source("/icbb/projects/nitschre/methylTFR/R/plot_helpers.R")
source("/icbb/projects/nitschre/methylTFR/R/expected_deviations.R")

plot_dir <- "/icbb/projects/nitschre/methylTFR/figures/figure4/EchoFootprints"
## Pseudobulks from ECHO methylation
# Define the paths to the methylome data for Th and Tc samples
paths <- list.files("/icbb/projects/igunduz/DARPA_analysis/ECHO_Cell_PB_261025/", full.names=TRUE)

cell_type_colors <- c(
  "Observed divided Expected B-cell" = "#AE017E",
  "Observed divided Expected Monocyte" = "#CC4C02",
  "Observed divided Expected NK-cell" = "#A65628",
  "Observed divided Expected Th-Mem" = "#41B6C4",
  "Observed divided Expected Tc-Mem" = "#4292C6",
  "Observed divided Expected Tc-Naive" = "#888FB5",
  "Observed divided Expected Th-Naive" = "#C7E9B4"
)

msites <- GRangesList()
for(path in paths){
    msite <-read_methylome(path, "bismarkcytosine")
    name <- file_path_sans_ext(basename(path))
    msites[[name]] <- msite
}

# Load annotation files
gcfreqs <- getGCfreq(motifSet = "jaspar2020_distal")
tf_bindsites <- getTFbindsites(motifSet = "jaspar2020")
gc_dist <- getGenomeGC()
enhancer <- readRDS("/icbb/projects/share/annotations/methylTFRAnnotationHg38/inst/extdata/distal_regions.RDS")

samples <- msites
cell_type <- names(msites)[1]

# Footprint plots
source("/icbb/projects/nitschre/methylTFR/scripts/blueprint/plot_and_save_difference_covid.r")
plot_and_save_difference_covid(msites, plot_dir, cell_type_colors)

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
th_mem_sample <- "/icbb/projects/igunduz/new_scmeth_pseudobulks_070224/Th-Mem.bedGraph"
th_naive_sample <- "/icbb/projects/igunduz/new_scmeth_pseudobulks_070224/Th-Naive.bedGraph"
tc_naive_sample <- "/icbb/projects/igunduz/new_scmeth_pseudobulks_070224/Tc-Naive.bedGraph"
tc_mem_sample <- "/icbb/projects/igunduz/new_scmeth_pseudobulks_070224/Tc-Mem.bedGraph"

obs_colors <- c(
  "Observed divided Expected Th-Mem" = "#41B6C4",
  "Observed divided Expected Th-Naive" = "#C7E9B4",
  "Observed divided Expected Tc-Mem" = "#4292C6",
  "Observed divided Expected Tc-Naive" = "#888FB5")

paths <- c(th_mem_sample, th_naive_sample, tc_mem_sample, tc_naive_sample)
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

# Get TFS from differentials
diff_cm <- readRDS("/icbb/projects/nitschre/methylTFR/r_objects/diff_cm.RDS")
diff_em <- readRDS("/icbb/projects/nitschre/methylTFR/r_objects/diff_em.RDS")

diff_cm <- rownames(diff_cm[diff_cm$p_value_adjusted < 0.05,])
diff_em <- rownames(diff_em[diff_em$p_value_adjusted < 0.05,])

logical_index <- names(tf_bindsites) %in% c(diff_cm, diff_em)
tf_bindsites_diff <- tf_bindsites[logical_index] 

# Footprint plots
source("/icbb/projects/nitschre/methylTFR/scripts/blueprint/plot_and_save_difference_covid.r")
plot_and_save_difference_covid(msites, plot_dir, obs_colors)

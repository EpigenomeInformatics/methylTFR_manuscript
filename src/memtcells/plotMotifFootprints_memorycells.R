# Adaption of plotMotifFootprints for Memory Tcells
set.seed(13)

suppressPackageStartupMessages({
    library(methylTFR)
    library(RnBeads)
    library(methylTFRAnnotationHg38)
    library(dplyr)
})
source("/icbb/projects/nitschre/methylTFR/scripts/memoryTcells/plots.R")
source("/icbb/projects/nitschre/methylTFR/R/plot_helpers.R")
source("/icbb/projects/nitschre/methylTFR/R/expected_deviations.R")

plot_dir <- "/icbb/projects/nitschre/methylTFR/figures"

# Import Rnbeads object
rnbset <- load.rnb.set("/icbb/projects/nitschre/methylTFR/results/memoryTcells/reports/data_import_data/rnb.set_preprocessed")

# Merging cell type replicates
rnbset_merged <- mergeSamples(rnbset, "cellType")

msites <- rnb.RnBSet.to.GRangesList(rnbset_merged)
save(msites, file ="//icbb/projects/nitschre/methylTFR/r_objects/msites.rda")

# Load annotation files
gcfreqs <- getGCfreq(motifSet = "jaspar2020_distal")
tfset <- "jaspar2020"
tf_bindsites <- getTFbindsites(motifSet = tfset)
gc_dist <- readRDS("/icbb/projects/igunduz/annotation/genome_wide_GC.RDS")
#enhancer <- readRDS("/icbb/projects/share/annotations/methylTFRAnnotationHg38/inst/extdata/distal_regions.RDS")
#filtered_names <- names(tf_bindsites)

# Plot
path <- file.path(plot_dir, "motifFootprint_TN_FOXF2.pdf")
combinded_data <- lapply(names(msites), function(celltype){
    p <- plotMotifFootprint(
    motif = "FOXF2",
    tf_bindsites = tf_bindsites,
    msites = msites[[celltype]],
    sample_name = celltype,
    gc_dist = gc_dist,
    gcfreqs = gcfreqs,
    method = "division"
    )
    return(p)
})



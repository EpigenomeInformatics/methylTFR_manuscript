#!/usr/bin/env Rscript

#####################################################################
# 02_run_mTFR_blueprint.R
# created on 27-10-25 by Irem B Gunduz
# Updated by IBG on 23-08-2026
# Run methylTFR analysis on the preprocessed Blueprint RnBeads object
#####################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(GenomicRanges)
  library(logger)
  library(RnBeads)
  library(methylTFR)
  library(methylTFRAnnotationHg38)
})
set.seed(42)

# Motif sets to run, names are lower case in the new annotation package
motifSetList <- c("jaspar2020", "jaspar2020_distal")
num.cores <- 32
chunk.size <- 10

# Directory where the annotation resources are stored locally
annotation.dir <- "/icbb/projects/share/annotations/methylTFRAnnotationHg38/inst/extdata"
if (!dir.exists(annotation.dir)) {
  stop("Annotation directory does not exist: ", annotation.dir)
}
options(methylTFRAnnotationHg38.datadir = annotation.dir)

# Directory where the RnBeads output was written to by 01
analysis.dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/blueprint/"
rnb.dir <- file.path(analysis.dir, "RnBeads_230826")
rnb.set.path <- file.path(rnb.dir, "reports", "data_import_data", "rnb.set_preprocessed")

# Directory where the deviations should be written to
out.dir <- file.path(analysis.dir, "mTFR_devs_230826")
if (!dir.exists(out.dir)) dir.create(out.dir, recursive = TRUE)

# Distal regions used to restrict the jaspar2020_distal run.
# Ensembl Regulatory Build v104, hg38, filtered to distal. This is the same
# region set that jaspar2020_distal_motif_gcfreq.rds was built against, so it
# must not be regenerated from the LOLA bed files: a different region set would
# make the observed footprints and the GC background disagree.
distal.file <- "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFRAnnotationHg38_old/inst/extdata/distal_regions.RDS"
if (!file.exists(distal.file)) {
  stop("Distal regions file does not exist: ", distal.file)
}
distal <- readRDS(distal.file)
log_info("Loaded ", length(distal), " distal regions from ", distal.file)

# Load the preprocessed RnBeads object
log_info("Loading RnBeads object from ", rnb.set.path)
rnb.set <- RnBeads::load.rnb.set(rnb.set.path)

# Genome wide GC distribution is shared across the motif sets
gc_dist <- getGenomeGC("hg38")

for (motifSet in motifSetList) {
  out.file <- file.path(out.dir, paste0(motifSet, "_deviations.RDS"))
  if (file.exists(out.file)) {
    log_info("Skipping ", motifSet, ", deviations already exist")
    next
  }

  log_info("Running methylTFR for ", motifSet)
  log_info("Loading the TF binding sites and GC freqs")

  # The distal run reuses the JASPAR2020 binding sites, restricted to distal regions
  enhancer <- if (motifSet == "jaspar2020_distal") distal else NULL
  tfset <- if (motifSet == "jaspar2020_distal") "jaspar2020" else motifSet

  gcfreqs <- getGCfreq(motifSet = motifSet)
  tf_bindsites <- getTFbindsites(motifSet = tfset)

  log_info("Number of motifs in gcfreqs: ", length(gcfreqs))
  log_info("Out file: ", out.file)

  # methylTFR reads the calls at single cytosine resolution from the RnBSet
  deviations <- run_methylTFR_RnBeads(
    rnb_set = rnb.set,
    tf_bindsites = tf_bindsites,
    gcfreqs = gcfreqs,
    gc_dist = gc_dist,
    chunkSize = chunk.size,
    threads = num.cores,
    enhancer = enhancer,
    ignoreStrand = TRUE,
    cov_threshold = 1
  )

  saveRDS(deviations, out.file)
  log_success("Finished ", motifSet)
}

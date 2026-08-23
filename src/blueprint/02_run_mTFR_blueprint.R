#!/usr/bin/env Rscript

#####################################################################
# 02_run_mTFR_blueprint.R
# created on 27-10-25 by IBG
# Run methylTFR analysis on BLUEPRINT data
#####################################################################

set.seed(42)
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(methylTFR)
  library(methylTFRAnnotationHg38)
  library(logger)
  library(GenomicRanges)
  library(muLogR)
  library(stringr)
  library(RnBeads)
})

# For no bias correction
# source("/icbb/projects/nitschre/methylTFR/methylTFR_manuscript/src/run_mTFR_RnBeads_noCorrection.R", chdir = TRUE)

motifSetList <- c("JASPAR2020_distal", "JASPAR2020")
main.dir <- "/scratch/icbb/regina/data"
sample_dir <- paste0(main.dir, "RnBeads")
out.dir <- paste0(main.dir, "mTFR_devs/")

if (!dir.exists(out.dir)) {
  dir.create(out.dir)
}
if (!file.exists("/icbb/projects/igunduz/annotation/methylTFRAnnotationHg38/inst/extdata/distal_regions.RDS")) {
  distal <- fread("/icbb/projects/share/annotations/lolaDB/hg38/EnsemblRegBuildBP/regions/regionSet_1.bed", header = FALSE)
  distal$V6 <- str_replace(distal$V6, ".", "*")
  distal <- GRanges(
    seqnames = distal$V1,
    ranges = IRanges(start = distal$V2, end = distal$V3),
    strand = distal$V6
  )
  saveRDS(distal, "/icbb/projects/igunduz/annotation/methylTFRAnnotationHg38/inst/extdata/distal_regions.RDS")
} else {
  distal <- readRDS("/icbb/projects/igunduz/annotation/methylTFRAnnotationHg38/inst/extdata/distal_regions.RDS")
}
logger.info("Loading RnBeads object...")
rnb_set <- RnBeads::load.rnb.set(paste0(sample_dir, "reports/data_import_data/rnb.set_preprocessed"))

for (motifSet in motifSetList) {
  logger.info(paste0("Running methylTFR for ", motifSet))
  logger.info("Loading the TF binding sites, GC freqs and GC dist")
  distal <- if (motifSet != "JASPAR2020_distal") {
    NULL
  } else {
    distal
  }
  tfset <- if (motifSet == "JASPAR2020_distal") {
    "JASPAR2020"
  } else {
    motifSet
  }
  gcfreqs <- getGCfreq(motifSet)
  gc_dist <- getGenomeGC()  
  tf_bindsites <- getTFbindsites(motifSet = tfset)

  logger::log_info("Number of motifs in gcfreqs: ", length(gcfreqs))
  logger::log_info("Out dir: ", out.dir)
  deviations <- run_methylTFR_RnBeads(
    rnb_set = rnb_set,
    threads = 32,
    chunkSize = 15,
    tf_bindsites = tf_bindsites,
    gcfreqs = gcfreqs,
    gc_dist = gc_dist,
    enhancer = distal,
    ignoreStrand = TRUE
  )
  saveRDS(deviations, paste0(out.dir, motifSet, "_deviations.RDS"))
}

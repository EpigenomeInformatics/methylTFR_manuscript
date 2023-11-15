#!/usr/bin/env Rscript

#####################################################################
# 03_run_methylTFR_blueprint.R
# created on 2023-11-12 by Irem Gunduz
# Run methylTFR JASPAR2020 and ALTIUS motifs analysis on BLUEPRINT data
#####################################################################

set.seed(42)
logger::log_info("Loading libraries...")
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(methylTFR)
  library(methylTFRAnnotationHg38)
  library(logger)
  library(muLogR)
})

motifSetList <- c("jaspar2020", "altius")
sample_dir <- "/icbb/projects/igunduz/methylTFR_manuscript/data/BLUEPRINT"


for (motifSet in motifSetList) {
  logger.info(paste0("Running methylTFR for ", motifSet))
  logger::log_info("Loading the TF binding sites, GC freqs and GC dist")
  gcfreqs <- getGCfreq(motifSet = motifSet)
  gc_dist <- getGenomeGC()
  tf_bindsites <- getTFbindsites(motifSet = motifSet)

  logger::log_info("Number of motifs in gcfreqs: ", length(gcfreqs))
  if (length(gcfreqs) != length(tf_bindsites)) {
    logger::log_info("Number of motifs in gcfreqs and tf_bindsites are not equal")
    logger::log_info("Number of motifs in tf_bindsites before filtering: ", length(tf_bindsites))
    tf_bindsites <- tf_bindsites[names(gcfreqs)]
    logger::log_info("Number of motifs in tf_bindsites after filtering: ", length(tf_bindsites))
  }
  logger::log_info("Number of NAs in gc_dist: ", sum(is.na(gc_dist)))

  out.dir <- paste0("/icbb/projects/igunduz/methylTFR_manuscript/results/BLUEPRINT/mtfr_", motifSet, "_121123/")
  if (!dir.exists(out.dir)) {
    dir.create(out.dir)
  }
  deviations <- run_methyltfr(
    sample_ann = "samples.tsv",
    sample_dir = sample_dir,
    full_path = FALSE,
    threads = 30,
    chunkSize = 10,
    tf_bindsites = tf_bindsites,
    gcfreqs = gcfreqs,
    gc_dist = gc_dist,
    filetype = "EPP",
    enhancer=NULL
  )
  saveRDS(deviations, paste0(out.dir, motifSet, "_deviations.RDS"))
}

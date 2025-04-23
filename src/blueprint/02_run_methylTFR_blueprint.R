#!/usr/bin/env Rscript

#####################################################################
# 03_run_methylTFR_blueprint.R
# created on 2023-11-12 by Irem Gunduz
# Run methylTFR analysis on BLUEPRINT data
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
  library(RnBeads)
  library(stringr)
})

motifSetList <- c("jaspar2020_distal","altius")[1]
#sample_dir <- "/icbb/projects/igunduz/methylTFR_manuscript/results/BLUEPRINT/"
sample_dir <- "/scratch/icbb/mtfr_manuscript/BLUEPRINT/"
mtfr_dir <- "/icbb/projects/igunduz/mTFR_bias_fix_v3/"
out.dir <- paste0(mtfr_dir,"blueprint_020724/")

if (!dir.exists(out.dir)) {dir.create(out.dir)}
if (!file.exists("/icbb/projects/share/annotations/methylTFRAnnotationHg38/inst/extdata/distal_regions.RDS")) {
    distal <- fread("/icbb/projects/share/annotations/lolaDB/hg38/EnsemblRegBuildBP/regions/regionSet_1.bed", header = FALSE) 
    distal$V6 <- str_replace(distal$V6, ".", "*")
    distal <- GRanges(seqnames = distal$V1,
                  ranges = IRanges(start = distal$V2, end = distal$V3), 
                  strand = distal$V6)
  saveRDS(distal, "/icbb/projects/igunduz/annotation/methylTFRAnnotationHg38/inst/extdata/distal_regions.RDS")
}else{
  distal <- readRDS("/icbb/projects/share/annotations/methylTFRAnnotationHg38/inst/extdata/distal_regions.RDS")
}
logger.info("Loading RnBeads object...")
rnb_set <- RnBeads::load.rnb.set(paste0(sample_dir,"reports/data_import_data/rnb.set_preprocessed"))
gc_dist <- getGenomeGC()

for (motifSet in motifSetList) {
  logger.info(paste0("Running methylTFR for ", motifSet))
  logger.info("Loading the TF binding sites, GC freqs and GC dist")
  distal <- if(motifSet != "jaspar2020_distal"){NULL}else{distal}
  tfset <- if(motifSet == "jaspar2020_distal"){ "jaspar2020"}else{ motifSet}
  gcfreqs <- getGCfreq(motifSet = motifSet)
  tf_bindsites <- getTFbindsites(motifSet = tfset)
  logger::log_info("Number of motifs in gcfreqs: ", length(gcfreqs))
  logger::log_info("Out dir: ", out.dir)
  deviations <- run_methylTFR_RnBeads(
    rnb_set= rnb_set,
    threads = 32,
    chunkSize = 15,
    tf_bindsites = tf_bindsites,
    gcfreqs = gcfreqs,
    gc_dist = gc_dist,
    enhancer=distal,
    ignoreStrand = TRUE
  )
  saveRDS(deviations, paste0(out.dir, motifSet, "_deviations.RDS"))
}

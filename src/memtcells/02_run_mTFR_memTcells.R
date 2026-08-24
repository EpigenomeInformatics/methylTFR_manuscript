#!/usr/bin/env Rscript

#####################################################################
# 02_run_mTFR_memTcells.R
# created on 29-10-2025 by Irem B Gunduz
# Updated by IBG on 23-08-2026
# Run methylTFR on the preprocessed memory T cell RnBeads object
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
chunk.size <- 15

# TEMRA is unreplicated, it is dropped here so every downstream script
# works on the same samples
drop.samples <- "51_Hf03_BlTR_Ct_WGBS_S_1.MCSv3.20170714.GRCh38.cpg.filtered.CG.bed"

# Distal regions used to restrict the jaspar2020_distal run.
# Ensembl Regulatory Build v104, hg38, filtered to distal. This is the same
# region set that jaspar2020_distal_motif_gcfreq.rds was built against, so it
# must not be regenerated from the LOLA bed files.
distal.file <- "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFRAnnotationHg38_old/inst/extdata/distal_regions.RDS"
if (!file.exists(distal.file)) {
  stop("Distal regions file does not exist: ", distal.file)
}
distal <- readRDS(distal.file)
log_info("Loaded ", length(distal), " distal regions from ", distal.file)

# Directory where the RnBeads output was written to by 01
analysis.dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/memoryTcells"
rnb.set.path <- file.path(analysis.dir, "reports", "data_import_data", "rnb.set_preprocessed")

# Directory where the deviations should be written to
out.dir <- file.path(analysis.dir, "mTFR_devs_230826")
if (!dir.exists(out.dir)) dir.create(out.dir, recursive = TRUE)

# Load the preprocessed RnBeads object once, not once per motif set
log_info("Loading RnBeads object from ", rnb.set.path)
rnb.set <- RnBeads::load.rnb.set(rnb.set.path)

drop.samples <- intersect(drop.samples, samples(rnb.set))
if (length(drop.samples) > 0) {
  log_info("Removing ", length(drop.samples), " sample(s): ", paste(drop.samples, collapse = ", "))
  rnb.set <- remove.samples(rnb.set, drop.samples)
}
log_info(length(samples(rnb.set)), " samples enter methylTFR")

# Genome wide GC distribution is shared across the motif sets
gc_dist <- getGenomeGC()

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

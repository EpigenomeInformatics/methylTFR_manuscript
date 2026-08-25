#!/usr/bin/env Rscript

#####################################################################
# debug_sc_mtfr.R
# Diagnose why computeDeviation fails on a single ECHO allc file.
#
# run_methyltfr hides the real error: mclapply returns a try-error
# object per failed motif, and write_block_to_sink then calls $dev on
# it, which is the "$ operator is invalid for atomic vectors" message.
# Everything below runs serially so the actual error is printed.
#
# Usage:
#   Rscript debug_sc_mtfr.R [/path/to/one.allc.tsv.gz]
#####################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(methylTFR)
  library(methylTFRAnnotationHg38)
})

motifSet <- "jaspar2020_distal"
tfSet <- "jaspar2020"
n.motifs <- 5
cov.threshold <- 1
distal.file <- "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFRAnnotationHg38_old/inst/extdata/distal_regions.RDS"

args <- commandArgs(trailingOnly = TRUE)
allc.file <- if (length(args) > 0) {
  args[1]
} else {
  "/icbb/projects/share/datasets/ECHO/allcFiles/HIV/CGN_allC_StrndMrgd_HIV_1-HIV_3-A1-AD001.CGN-Merge.allc.tsv.gz"
}

hdr <- function(x) cat("\n===", x, "===\n")

#####################################################################
hdr("1. Raw allc file")
#####################################################################

cat("file:", allc.file, "\n")
if (!file.exists(allc.file)) stop("file does not exist")
raw <- fread(allc.file, header = FALSE, nrows = 5, showProgress = FALSE)
cat("columns:", ncol(raw), "\n")
print(raw)
cat("\nparse_allc reads V1=chr V2=pos V3=strand V5=mc V6=cov\n")

#####################################################################
hdr("2. read_methylome")
#####################################################################

msites <- read_methylome(allc.file, type = "allc", cov_threshold = cov.threshold)
cat("sites after cov >=", cov.threshold, ":", length(msites), "\n")
cat("seqlevels:", paste(head(seqlevels(msites), 30), collapse = ", "), "\n")
cat("strand levels:", paste(unique(as.character(strand(msites))), collapse = ", "), "\n")
cat("score range:", paste(range(msites$score, na.rm = TRUE), collapse = " .. "), "\n")
cat("non finite scores:", sum(!is.finite(msites$score)), "\n")
cat("coverage summary:\n")
print(summary(msites$coverage))

#####################################################################
hdr("3. Annotation")
#####################################################################

gcfreqs <- getGCfreq(motifSet = motifSet)
tf_bindsites <- getTFbindsites(motifSet = tfSet)
gc_dist <- getGenomeGC()
enhancer <- readRDS(distal.file)

cat("motifs in gcfreqs:", length(gcfreqs), "\n")
cat("motifs in tf_bindsites:", length(tf_bindsites), "\n")
cat("shared motif names:", length(intersect(names(gcfreqs), names(tf_bindsites))), "\n")
cat("gc_dist ranges:", length(gc_dist), "\n")
cat("distal ranges:", length(enhancer), "\n")

cat("\nseqlevel overlap with msites:\n")
cat("  gc_dist     :", length(intersect(seqlevels(msites), seqlevels(gc_dist))), "shared\n")
cat("  tf_bindsites:", length(intersect(seqlevels(msites), seqlevels(tf_bindsites[[1]]))), "shared\n")
cat("  distal      :", length(intersect(seqlevels(msites), seqlevels(enhancer))), "shared\n")

cat("\nGC_bin values in gc_dist:\n")
print(table(gc_dist$GC_bin, useNA = "ifany"))

#####################################################################
hdr("4. addGCBintoMethylome, the usual failure point")
#####################################################################

# computeExpectations does t(gcfreq) %*% binMsites[, 2], so the number of
# GC bins this cell covers MUST equal nrow(gcfreq). A cell that misses a
# bin, or that picks up an NA bin, gives non conformable arguments for
# EVERY motif at once, which is what "all scheduled cores encountered
# errors" looks like.
bin_meth <- addGCBintoMethylome(msites, gc_dist, TRUE)
cat("bin_meth dim:", paste(dim(bin_meth), collapse = " x "), "\n")
print(bin_meth)

nrow_gcfreq <- nrow(gcfreqs[[intersect(names(gcfreqs), names(tf_bindsites))[1]]])
cat("\nnrow(gcfreq) =", nrow_gcfreq, "  nrow(bin_meth) =", nrow(bin_meth), "\n")
if (nrow(bin_meth) != nrow_gcfreq) {
  cat("*** MISMATCH: this alone makes every motif fail ***\n")
  cat("    bins present:", paste(bin_meth[, 1], collapse = ", "), "\n")
} else {
  cat("    conformable, so the failure is elsewhere\n")
}
if (any(is.na(bin_meth))) cat("*** bin_meth contains NA ***\n")

#####################################################################
hdr("5. Per motif overlaps")
#####################################################################

motifs <- intersect(names(gcfreqs), names(tf_bindsites))[seq_len(n.motifs)]
for (m in motifs) {
  tfbs <- tf_bindsites[[m]]
  tfbs <- resize(tfbs, width(tfbs)[1] + 130, fix = "center")
  n_all <- length(tfbs)
  tfbs_d <- suppressWarnings(subsetByOverlaps(tfbs, enhancer, ignore.strand = TRUE))
  hits <- suppressWarnings(
    findOverlaps(msites, tfbs_d, type = "within", ignore.strand = TRUE)
  )
  cat(sprintf(
    "%-16s tfbs %7d -> distal %7d   covered CpGs %6d   gcfreq %d x %d\n",
    m, n_all, length(tfbs_d), length(hits),
    nrow(gcfreqs[[m]]), ncol(gcfreqs[[m]])
  ))
}

#####################################################################
hdr("6. computeDeviation, serial, real error message")
#####################################################################

for (m in motifs) {
  res <- tryCatch(
    computeDeviation(m, msites, tf_bindsites, gcfreqs,
      enhancer = enhancer, ignoreStrand = TRUE, binMsites = bin_meth
    ),
    error = function(e) e
  )
  if (inherits(res, "error")) {
    cat(sprintf("%-16s ERROR: %s\n", m, conditionMessage(res)))
  } else {
    cat(sprintf("%-16s ok  dev=%s exp_dev=%s\n", m,
      format(res$dev, digits = 4), format(res$exp_dev, digits = 4)))
  }
}

#####################################################################
hdr("7. GC bins UNDER the distal restriction, which is what the run does")
#####################################################################

# Section 4 above used the full gc_dist, but methyltfr_core restricts it to
# the enhancer regions BEFORE binning:
#
#   if (!is.null(enhancer)) gc_dist <- subsetByOverlaps(gc_dist, enhancer, ...)
#   bin_meth <- addGCBintoMethylome(msites, gc_dist, ignoreStrand)
#
# So the bins that matter are the ones this cell covers inside distal
# regions only. computeExpectations still needs exactly nrow(gcfreq) of
# them, so a missing bin breaks every motif at once.
gc_d <- suppressWarnings(subsetByOverlaps(gc_dist, enhancer, ignore.strand = TRUE))
cat("gc_dist ranges after distal restriction:", length(gc_d),
  sprintf("(%.2f%% of genome wide)\n", 100 * length(gc_d) / length(gc_dist)))
cat("GC_bin values available in the distal restricted table:\n")
print(table(gc_d$GC_bin, useNA = "ifany"))

bin_meth_d <- tryCatch(addGCBintoMethylome(msites, gc_d, TRUE), error = function(e) e)
if (inherits(bin_meth_d, "error")) {
  cat("*** addGCBintoMethylome FAILED:", conditionMessage(bin_meth_d), "***\n")
} else {
  cat("\nbin_meth dim under distal restriction:",
    paste(dim(bin_meth_d), collapse = " x "), "\n")
  print(bin_meth_d)
  cat("\nnrow(gcfreq) =", nrow_gcfreq, "  nrow(bin_meth) =", nrow(bin_meth_d), "\n")
  if (nrow(bin_meth_d) != nrow_gcfreq) {
    cat("*** THIS IS THE BUG: t(gcfreq) %*% binMsites[, 2] is non conformable,\n")
    cat("    so every motif fails and mclapply returns try-errors. ***\n")
    cat("    bins covered:", paste(bin_meth_d[, 1], collapse = ", "), "\n")
  } else {
    cat("    conformable here, so this cell would have run\n")
  }
}

cat("\nSites overlapping the distal restricted GC table:\n")
hits_d <- suppressWarnings(findOverlaps(msites, gc_d, ignore.strand = TRUE))
cat("  ", length(unique(hits_d@from)), "of", length(msites), "sites\n")

cat("\nDone.\n")

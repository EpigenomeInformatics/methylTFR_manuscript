#!/usr/bin/env Rscript

#####################################################################
# 01_run_rneads_blueprint.R
# created on 27-10-25 by Irem B Gunduz
# Updated by IBG on 23-08-2026
# Filter the existing Blueprint RnBeads run down to the healthy samples
#
# The full RnBeads import and preprocessing takes too long to repeat, so
# this reuses the preprocessed object of the 291025 run, removes the
# disease samples and writes the result under a new date. The source run
# is never modified. The vanilla import that produced it is in git
# history, before this revision.
#####################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(logger)
  library(RnBeads)
})
set.seed(12)

# Directories, the source run is read only from here on
analysis.dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/blueprint/"
source.tag <- "RnBeads_291025"
target.tag <- "RnBeads_230826"

source.dir <- file.path(analysis.dir, source.tag, "reports", "data_import_data")
source.set <- file.path(source.dir, "rnb.set_preprocessed")

target.dir <- file.path(analysis.dir, target.tag, "reports", "data_import_data")
target.set <- file.path(target.dir, "rnb.set_preprocessed")
target.annot <- file.path(target.dir, "annotation.csv")

# Number of samples the filtered set must contain, set to NULL to skip
expected.samples <- 147

num.cores <- 32

if (!file.exists(source.set)) {
  stop("Source RnBeads object does not exist: ", source.set)
}
if (file.exists(target.set)) {
  stop(
    "Target RnBeads object already exists: ", target.set,
    "\nRemove it or change target.tag, this script never overwrites."
  )
}
if (!dir.exists(target.dir)) dir.create(target.dir, recursive = TRUE)

# Multiprocess
parallel.setup(num.cores)

#####################################################################
# Load the existing preprocessed object
#####################################################################

log_info("Loading RnBeads object from ", source.set)
rnb.set <- load.rnb.set(source.set)
sannot <- pheno(rnb.set)
log_info(source.tag, ": ", nrow(sannot), " samples")

#####################################################################
# Remove the disease samples
#####################################################################

# The Blueprint cancer samples are flagged through the DISEASE column
# (7 Acute Lymphocytic Leukemia Bcell_pre and 3 Multiple Myeloma plasma).
# They are annotated as cellTypeGroup "Bcell" and "plasma", not "cancer",
# so filtering on cellTypeGroup alone would keep every one of them.
disease <- as.character(sannot$DISEASE)
drop_idx <- which(is.na(disease) | disease != "None")

# Kept as a safeguard, no sample currently carries this label
if ("cellTypeGroup" %in% colnames(sannot)) {
  drop_idx <- unique(c(
    drop_idx,
    which(as.character(sannot$cellTypeGroup) == "cancer")
  ))
}

if (length(drop_idx) > 0) {
  log_info("Removing ", length(drop_idx), " samples: ",
    paste(unique(disease[drop_idx]), collapse = ", "))
  rnb.set <- remove.samples(rnb.set, drop_idx)
}

sannot <- pheno(rnb.set)
log_info(target.tag, ": ", nrow(sannot), " samples remaining")

if (!is.null(expected.samples) && nrow(sannot) != expected.samples) {
  stop(
    "Expected ", expected.samples, " samples after filtering but got ",
    nrow(sannot), ". Check the DISEASE column of the source run before ",
    "running anything downstream."
  )
}

#####################################################################
# Save the filtered object and its annotation
#####################################################################

# The annotation is written from the object itself, so it can never
# disagree with the samples that 02 to 06 read.
write.csv(sannot, target.annot, row.names = FALSE)
log_info("Wrote ", target.annot)

save.rnb.set(rnb.set, target.set, archive = FALSE)
log_success("Wrote ", target.set)

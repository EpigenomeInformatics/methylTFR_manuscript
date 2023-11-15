#!/usr/bin/env Rscript

#####################################################################
# 02_bedToEPP_blueprint.R
# created on 2023-11-10 by Irem Gunduz
# Export EPPs from RnBeads object, for BLUEPRINT data
#####################################################################

suppressPackageStartupMessages({
  #library(dplyr)
  library(RnBeads)
  #library(data.table)
  library(muLogR)
  library(methylTFR)
})
set.seed(42)

logger.info("Setting up directories...")
rnbeads.dir <- "/icbb/projects/igunduz/methylTFR_manuscript/results/BLUEPRINT/"
out.dir <- "/icbb/projects/igunduz/methylTFR_manuscript/results/BLUEPRINT/bedToEPP/"
if (!dir.exists(out.dir)) {dir.create(out.dir)}


report.dir <- file.path(rnbeads.dir, "reports")
rnbset <- RnBeads::load.rnb.set(paste0(report.dir, "/data_import_data/rnb.set_preprocessed"))
result <- RnBeads::rnb.RnBSet.to.GRangesList(rnbset, "sites")
rm(rnbset)


# convert sites to EPP object
epps <- methylTFR::convertToEPP(result,
  filePath = out.dir,
  save = TRUE,
  threads = 30,
  verbose = TRUE
)

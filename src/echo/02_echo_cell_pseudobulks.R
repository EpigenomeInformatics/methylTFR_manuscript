#!/usr/bin/env Rscript

#####################################################################
# 02_echo_cell_pseudobulks.R
# created on 25-10-24 by IBG
# Create pseudobulks for the methylation data
#####################################################################

suppressPackageStartupMessages({
  library(muLogR)
  library(dplyr)
  library(data.table)
})
source("/icbb/projects/igunduz/exposure_atlas_manuscript/src/utils/createPseudoBulks.R")
outputDir <- "/icbb/projects/igunduz/DARPA_analysis/ECHO_Cell_PB_261025/"
if (!dir.exists(outputDir)) {
  dir.create(outputDir)
}

# read the sample annotation
sampleAnnot <- data.table::fread("/icbb/projects/igunduz/DARPA_analysis/artemis_031023/sample_annot.tsv") %>%
  dplyr::filter(!cell_type == "Other-cell")

# add full file paths
sampleAnnot$allC_FilePathfull <- gsub(sampleAnnot$allC_FilePath_fixed,
  pattern = "/icbb/projects/igunduz/DARPA/allcFiles/",
  replacement = "/icbb/projects/share/datasets/ECHO/allcFiles/"
)

# Create single base-pair pseudobulks
createPseudoBulks(sampleAnnot,
  filePathCol = "allC_FilePathfull",
  sampleIdCol = "CommonMinID",
  numThreads = 30, mcType = "CGN", groupName = "cell_type", singleBP = TRUE,
  fileType = "allc", singleCovOff = 99999, groupName2 = NULL,
  indexed = FALSE, excludeChr = c("chrX", "chrY", "chrM", "chrL"),
  outputDir = outputDir
)

#####################################################################

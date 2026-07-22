#!/usr/bin/env Rscript

#####################################################################
# 01_run_rneads_memTcells.R
# created on 2023-11-15 by Irem Gunduz
# Run RnBeads vanilla analysis for memorTcell Methylation Data
#####################################################################
suppressPackageStartupMessages({
  library(dplyr)
  library(RnBeads)
  library(grid)
})
set.seed(12)

# Directory where your data is located
data.dir <- "/icbb/projects/igunduz/methylTFR_manuscript/data"
bed.dir <- file.path(data.dir, "memoryTcells")
sample.annotation <- file.path(bed.dir, "samples.tsv")
num.cores <- 30

# Directory where the output should be written to
analysis.dir <- "/icbb/projects/igunduz/methylTFR_manuscript/results"
if (!dir.exists(analysis.dir)) dir.create(analysis.dir)
analysis.dir <- file.path(analysis.dir, "memoryTcells")
if (!dir.exists(analysis.dir)) dir.create(analysis.dir)

# Directory where the report files should be written to
report.dir <- file.path(analysis.dir, "reports")
rnb.initialize.reports(report.dir)
rnb.options(
  analysis.name = "CD4 T memory cell subtypes  Methylation Data",
  assembly = "hg38",
  import.table.separator = "\t",
  disk.dump.big.matrices = TRUE,
  strand.specific = FALSE,
  region.types = c("cpgislands", "promoters", "genes", "tiling"),
  import.default.data.type = "data.dir",
  import.bed.style = "BisSNP",
  filtering.sex.chromosomes.removal = TRUE,
  identifiers.column = "bedFile",
  differential.comparison.columns = "cellType"
)

# Multiprocess
parallel.setup(num.cores)

data.source <- c(bed.dir, sample.annotation, 1)
result <- rnb.run.import(data.source = data.source, data.type = "bs.bed.dir", dir.reports = report.dir)
rnbset <- result$rnb.set

## Quality Control
rnb.run.qc(rnbset, report.dir)

## Preprocessing
rnbset <- rnb.run.preprocessing(rnbset, dir.reports = report.dir)$rnb.set

## save the object
save.rnb.set(rnbset, paste0(report.dir, "/data_import_data/rnb.set_preprocessed"), archive = FALSE)

## Differential methylation
rnb.run.differential(rnbset, report.dir)

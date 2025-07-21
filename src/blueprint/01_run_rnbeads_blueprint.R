#!/usr/bin/env Rscript

#####################################################################
# 01_run_rneads_blueprint.R
# created on 2023-11-09 by Irem Gunduz
# Run RnBeads vanilla analysis for Blueprint Bcell Methylation Data
#####################################################################
suppressPackageStartupMessages({
  library(dplyr)
  library(RnBeads)
  library(grid)
})
set.seed(12)

# Directory where your data is located
data.dir <- "/icbb/projects/share/datasets/"
bed.dir <- file.path(data.dir, "blueprint")
sample.annotation <- file.path(bed.dir, "samples.tsv")
num.cores <- 30

# Directory where the output should be written to
analysis.dir <- "/icbb/projects/igunduz/methylTFR_manuscript/results"
if (!dir.exists(analysis.dir)) dir.create(analysis.dir)
analysis.dir <- file.path(analysis.dir, "BLUEPRINT_080725")
if (!dir.exists(analysis.dir)) dir.create(analysis.dir)

# Directory where the report files should be written to
report.dir <- file.path(analysis.dir, "reports")
rnb.initialize.reports(report.dir)

# Multiprocess
parallel.setup(num.cores)

# Create tiling regions for 1kb
tiling1kb <- muRtools::getTilingRegions("hg38", width=1000L, onlyMainChrs=TRUE)%>%
  data.table::as.data.table() %>%
  dplyr::select(seqnames, start, end) %>%
  as.data.frame()
colnames(tiling1kb) <- c("Chromosome", "Start", "End")
rnb.set.annotation(type = "tiling1kb", regions = tiling1kb, assembly = "hg38")


# Set up the analysis
rnb.options(
  analysis.name = "Blueprint",
  assembly = "hg38",
  import.table.separator = "\t",
  region.aggregation = "sum",
  region.types = c("tiling1kb"),
  import.default.data.type = "data.dir",
  import.bed.style = "EPP",
  disk.dump.big.matrices = TRUE,
  strand.specific = FALSE,
  filtering.sex.chromosomes.removal = TRUE,
  differential.enrichment.lola = FALSE,
  identifiers.column = "bedFile",
  differential.comparison.columns = "cellTypeGroup" # not exclusive cell-types
)
if (!file.exists(paste0(report.dir, "/data_import_data/rnb.set_preprocessed"))) {
  data.source <- c(bed.dir, sample.annotation, 1)
  result <- rnb.run.import(data.source = data.source, data.type = "bs.bed.dir", dir.reports = report.dir)
  rnbset <- result$rnb.set

  ## save the object
  save.rnb.set(rnbset, paste0(report.dir, "/data_import_data/rnb.set_preprocessed"), archive = FALSE)

} else {
  rnbset <- RnBeads::load.rnb.set(paste0(report.dir, "/data_import_data/rnb.set_preprocessed"))
}

## Preprocessing
rnbset <- rnb.run.preprocessing(rnbset, dir.reports = report.dir)$rnb.set

## save the object
save.rnb.set(rnbset, paste0(report.dir, "/data_import_data/rnb.set_filtered_160725"), archive = FALSE)

## Quality Control
##rnb.run.qc(rnbset, report.dir)

## Preprocessing
#rnbset <- rnb.run.preprocessing(rnbset, dir.reports = report.dir)$rnb.set

## Differential methylation
#rnb.run.differential(rnbset, report.dir)

## Exploratory analysis
#rnb.run.exploratory(rnbset, report.dir)


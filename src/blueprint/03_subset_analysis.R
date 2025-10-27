#!/usr/bin/env Rscript

#####################################################################
# 03_subset_analysis.R
# created on 2023-11-09 by Irem Gunduz
# Run RnBeads vanilla analysis for Blueprint Methylation Data
# Subset to only include B cells and T cells
#####################################################################

set.seed(42)
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(logger)
  library(muLogR)
  library(RnBeads)
})

# Directory where your data is located
data.dir <- "/icbb/projects/share/datasets/"
bed.dir <- file.path(data.dir, "blueprint")
sample.annotation <- file.path(bed.dir, "samples_TB.tsv")
if (!file.exists(sample.annotation)) {
  sample.annotation <- file.path(bed.dir, "samples.tsv")
  sannot <- data.table::fread(sample.annotation)
  # Filter the cancer cells
  sannot <- sannot %>% filter(cellTypeGroup %in% c("Bcell", "Tcell"))
  sample.annotation <- file.path(bed.dir, "samples_TB.tsv")
  data.table::fwrite(sannot, sample.annotation, sep = "\t")
}
num.cores <- 30

# Directory where the output should be written to
analysis.dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/blueprint/"
if (!dir.exists(analysis.dir)) dir.create(analysis.dir)
analysis.dir <- file.path(analysis.dir, "TB_RnBeads_271025")
if (!dir.exists(analysis.dir)) dir.create(analysis.dir)
# Multiprocess
parallel.setup(30)

# Directory where the report files should be written to
report.dir <- file.path(analysis.dir, "reports")
rnb.initialize.reports(report.dir)

# Create tiling regions for 1kb
tiling1kb <- muRtools::getTilingRegions("hg38", width = 1000L, onlyMainChrs = TRUE) %>%
  data.table::as.data.table() %>%
  dplyr::select(seqnames, start, end) %>%
  as.data.frame()
colnames(tiling1kb) <- c("Chromosome", "Start", "End")
rnb.set.annotation(type = "tiling1kb", regions = tiling1kb, assembly = "hg38")

# Create distal regions
distal <- readRDS("/icbb/projects/share/annotations/methylTFRAnnotationHg38/inst/extdata/distal_regions.RDS")
distal <- as.data.frame(distal) %>%
  dplyr::select(seqnames, start, end)
colnames(distal) <- c("Chromosome", "Start", "End")
rnb.set.annotation(type = "distal", regions = distal, assembly = "hg38")

# Set up the analysis
rnb.options(
  analysis.name = "Blueprint Healthy Cells",
  assembly = "hg38",
  import.table.separator = "\t",
  region.aggregation = "sum",
  analyze.sites = FALSE,
  region.types = c("tiling1kb", "distal", "promoters"),
  import.default.data.type = "data.dir",
  import.bed.style = "EPP",
  disk.dump.big.matrices = TRUE,
  strand.specific = FALSE,
  filtering.sex.chromosomes.removal = TRUE,
  differential.enrichment.lola = TRUE,
  differential.enrichment.lola.dbs = "/icbb/projects/share/annotations/lolaDB/",
  identifiers.column = "bedFile",
  differential.comparison.columns = "cellTypeGroup" # not exclusive cell-types
)

# Data Import
data.source <- c(bed.dir, sample.annotation, 1)
result <- rnb.run.import(data.source = data.source, data.type = "bs.bed.dir", dir.reports = report.dir)
rnbset <- result$rnb.set

# Quality Control
rnb.run.qc(rnb.set, report.dir)

## Preprocessing
rnbset <- rnb.run.preprocessing(rnbset, dir.reports = report.dir)$rnb.set

## save the object
save.rnb.set(rnb.set, paste0(report.dir, "/data_import_data/rnb.set_preprocessed"), archive = FALSE)

# Differential methylation
rnb.run.differential(rnbset, report.dir)

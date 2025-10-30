#!/usr/bin/env Rscript

#####################################################################
# 01_run_rneads_blueprint.R
# created on 27-10-25 by Irem Gunduz
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
sample.annotation <- file.path(bed.dir, "samples_without_cancer.tsv")
if (!file.exists(sample.annotation)) {
  sample.annotation <- file.path(bed.dir, "samples.tsv")
  sannot <- data.table::fread(sample.annotation)
  # Filter the cancer cells
  sannot <- sannot %>% filter(cellTypeGroup != "cancer")
  sample.annotation <- file.path(bed.dir, "samples_without_cancer.tsv")
  data.table::fwrite(sannot, sample.annotation, sep = "\t")
}
num.cores <- 30

# Directory where the output should be written to
analysis.dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/blueprint/"
if (!dir.exists(analysis.dir)) dir.create(analysis.dir)
analysis.dir <- file.path(analysis.dir, "RnBeads_291025")
if (!dir.exists(analysis.dir)) dir.create(analysis.dir)

# Directory where the report files should be written to
report.dir <- file.path(analysis.dir, "reports")
rnb.initialize.reports(report.dir)

# Multiprocess
parallel.setup(num.cores)

# Create tiling regions for 1kb
tiling1kb <- muRtools::getTilingRegions("hg38", width = 1000L, onlyMainChrs = TRUE) %>%
  data.table::as.data.table() %>%
  dplyr::select(seqnames, start, end) %>%
  as.data.frame()
colnames(tiling1kb) <- c("Chromosome", "Start", "End")
rnb.set.annotation(type = "tiling1kb", regions = tiling1kb, assembly = "hg38")

# Create distal regions
# Distal regions
distal <- readRDS("/icbb/projects/share/annotations/methylTFRAnnotationHg38/inst/extdata/distal_regions.RDS")
distal <- as.data.frame(distal) %>%
  dplyr::select(seqnames, start, end)
colnames(distal) <- c("Chromosome", "Start", "End")
rnb.set.annotation(type = "distal", regions = distal, assembly = "hg38")

# Set up the analysis
rnb.options(
  analysis.name = "Blueprint",
  assembly = "hg38",
  import.table.separator = "\t",
  region.aggregation = "sum",
  region.types = c("tiling1kb", "distal"),
  import.default.data.type = "data.dir",
  import.bed.style = "EPP",
  disk.dump.big.matrices = TRUE,
  strand.specific = FALSE,
  filtering.sex.chromosomes.removal = TRUE,
  differential.enrichment.lola = FALSE,
  identifiers.column = "bedFile",
  differential.comparison.columns = "cellTypeGroup" # not exclusive cell-types
)

# Data Import
data.source <- c(bed.dir, sample.annotation, 1)
result <- rnb.run.import(data.source = data.source, data.type = "bs.bed.dir", dir.reports = report.dir)
rnb.set <- result$rnb.set

# Quality Control
# rnb.run.qc(rnb.set, report.dir)

## Preprocessing
rnb.set <- rnb.run.preprocessing(rnb.set, dir.reports = report.dir)$rnb.set

## save the object
save.rnb.set(rnb.set, paste0(report.dir, "/data_import_data/rnb.set_preprocessed"), archive = FALSE)

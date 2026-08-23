#!/usr/bin/env Rscript

#####################################################################
# 01_run_rneads_blueprint.R
# created on 27-10-25 by Irem B Gunduz
# Updated by IBG on 23-08-2026
# Run RnBeads vanilla analysis for Blueprint Methylation Data
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
  sannot <- sannot %>% filter(cellTypeGroup != "cancer" & DISEASE == "None")
  sample.annotation <- file.path(bed.dir, "samples_without_cancer.tsv")
  data.table::fwrite(sannot, sample.annotation, sep = "\t")
}
num.cores <- 30

# Directory where the output should be written to
analysis.dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/blueprint/"
if (!dir.exists(analysis.dir)) dir.create(analysis.dir)
analysis.dir <- file.path(analysis.dir, "RnBeads_230826")
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
# Ensembl Regulatory Build v104, hg38, filtered to distal. This is the same
# region set that jaspar2020_distal_motif_gcfreq.rds was built against, so
# all downstream scripts must read it from here.
distal.file <- "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFRAnnotationHg38_old/inst/extdata/distal_regions.RDS"
if (!file.exists(distal.file)) {
  stop("Distal regions file does not exist: ", distal.file)
}
distal <- readRDS(distal.file)
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

rnb.set <- rnb.execute.sex.removal(rnb.set.unfiltered)$dataset

# Remove sites that have an exceptionally high coverage
rnb.set <- rnb.execute.highCoverage.removal(rnb.set)$dataset

# Remove sites containing NA for beta values
rnb.set <- rnb.execute.na.removal(rnb.set)$dataset

# Remove sites for which the beta values have low standard deviation
rnb.set <- rnb.execute.variability.removal(rnb.set, 0.005)$dataset

## Preprocessing
rnb.set <- rnb.run.preprocessing(rnb.set, dir.reports = report.dir)$rnb.set

## save the object
save.rnb.set(rnb.set, paste0(report.dir, "/data_import_data/rnb.set_preprocessed"), archive = FALSE)

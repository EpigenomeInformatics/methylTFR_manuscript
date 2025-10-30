#!/usr/bin/env Rscript

#####################################################################
# 01_run_rneads_memTcells.R
# created on 29-10-2025 by Irem Gunduz
# Run RnBeads vanilla analysis for memoryTcell Methylation Data
#####################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(RnBeads)
  library(grid)
})
set.seed(12)

# Directory where your data is located
data.dir <- "/icbb/projects/share/datasets/"
bed.dir <- file.path(data.dir, "memoryTcells")
sample.annotation <- file.path(bed.dir, "samples.tsv")
num.cores <- 30

# Directory where the output should be written to
analysis.dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/"
analysis.dir <- file.path(analysis.dir, "memoryTcells")
if (!dir.exists(analysis.dir)) dir.create(analysis.dir)

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


rnb.options(
  analysis.name = "CD4 T memory cell subtypes  Methylation Data",
  assembly = "hg38",
  import.table.separator = "\t",
  disk.dump.big.matrices = TRUE,
  strand.specific = FALSE,
  region.types = c("tiling1kb", "distal"),
  import.default.data.type = "data.dir",
  import.bed.style = "BisSNP",
  filtering.sex.chromosomes.removal = TRUE,
  identifiers.column = "bedFile",
  differential.comparison.columns.all.pairwise = "cellType" #
)

# Multiprocess
parallel.setup(num.cores)

if (!file.exists(paste0(report.dir, "/data_import_data/rnb.set_preprocessed"))) {
  # Data Import
  data.source <- c(bed.dir, sample.annotation, 1)
  result <- rnb.run.import(data.source = data.source, data.type = "bs.bed.dir", dir.reports = report.dir)
  rnb.set <- result$rnb.set

  # Quality Control
  rnb.run.qc(rnb.set, report.dir)

  # Preprocessing
  rnb.set <- rnb.run.preprocessing(rnb.set, dir.reports = report.dir)$rnb.set

  # save the object
  save.rnb.set(rnb.set, paste0(report.dir, "/data_import_data/rnb.set_preprocessed"), archive = FALSE)
} else {
  # Load the preprocessed object
  rnb.set <- load.rnb.set(paste0(report.dir, "/data_import_data/rnb.set_preprocessed"))
}
if (!file.exists(paste0(report.dir, "/differential_methylation_data/differential_rnbDiffMeth"))) {
  # Remove TEMRA sample
  rnb.set <- remove.samples(rnb.set, "51_Hf03_BlTR_Ct_WGBS_S_1.MCSv3.20170714.GRCh38.cpg.filtered.CG.bed")
  # Differential methylation
  rnb.run.differential(rnb.set, report.dir)
} else {
  # Load differential methylation results
  diffMeth <- load.rnb.diffmeth(paste0(analysis.dir, "/reports/differential_methylation_data/differential_rnbDiffMeth/"))
}

# Run LOLA for differential methylation data
logger.start("Running LOLA")
lolaDb_path <- "/icbb/projects/share/annotations/lolaDB/lolaSubCluster/hg38/"

# Run LOLA
res <- performLolaEnrichment.diffMeth(rnb.set, diffMeth, lolaDb_path)
logger.info("Saving results")
saveRDS(res, paste0(analysis.dir, "/reports/differential_methylation_data/differential_rnbDiffMeth/lola_results.rds"))
logger.completed()

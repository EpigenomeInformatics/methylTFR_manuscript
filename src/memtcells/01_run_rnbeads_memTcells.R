#!/usr/bin/env Rscript

#####################################################################
# 01_run_rnbeads_memTcells.R
# created on 29-10-2025 by Irem B Gunduz
# Updated by IBG on 23-08-2026
# Run RnBeads vanilla analysis for the CD4 memory T cell data,
# then differential methylation and LOLA motif enrichment
#####################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(logger)
  library(LOLA)
  library(RnBeads)
})
set.seed(12)

# Directory where the data is located
data.dir <- "/icbb/projects/share/datasets/"
bed.dir <- file.path(data.dir, "memoryTcells")
sample.annotation <- file.path(bed.dir, "samples.tsv")
num.cores <- 30

# TEMRA is a single unreplicated sample, it cannot enter a group
# comparison and is excluded from every downstream script
drop.samples <- "51_Hf03_BlTR_Ct_WGBS_S_1.MCSv3.20170714.GRCh38.cpg.filtered.CG.bed"

# Directory where the output should be written to
analysis.dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/memoryTcells"
if (!dir.exists(analysis.dir)) dir.create(analysis.dir, recursive = TRUE)

# Directory where the report files should be written to
report.dir <- file.path(analysis.dir, "reports")
rnb.initialize.reports(report.dir)

rnb.set.path <- file.path(report.dir, "data_import_data", "rnb.set_preprocessed")
diffmeth.path <- file.path(report.dir, "differential_methylation_data", "differential_rnbDiffMeth")
lola.file <- file.path(diffmeth.path, "TF_motifs_lola.rds")

# Distal regions, Ensembl Regulatory Build v104 filtered to hg38 distal.
# Same region set the distal GC frequency tables were built against.
distal.file <- "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFRAnnotationHg38_old/inst/extdata/distal_regions.RDS"
if (!file.exists(distal.file)) {
  stop("Distal regions file does not exist: ", distal.file)
}

lolaDb.path <- "/icbb/projects/share/annotations/lolaDB/lolaTFmotifs/hg38/"

# Create tiling regions for 1kb
tiling1kb <- muRtools::getTilingRegions("hg38", width = 1000L, onlyMainChrs = TRUE) %>%
  data.table::as.data.table() %>%
  dplyr::select(seqnames, start, end) %>%
  as.data.frame()
colnames(tiling1kb) <- c("Chromosome", "Start", "End")
rnb.set.annotation(type = "tiling1kb", regions = tiling1kb, assembly = "hg38")

# Create distal regions
distal <- readRDS(distal.file)
distal <- as.data.frame(distal) %>%
  dplyr::select(seqnames, start, end)
colnames(distal) <- c("Chromosome", "Start", "End")
rnb.set.annotation(type = "distal", regions = distal, assembly = "hg38")

# Set up the analysis
rnb.options(
  analysis.name = "CD4 T memory cell subtypes Methylation Data",
  assembly = "hg38",
  import.table.separator = "\t",
  region.types = c("tiling1kb", "distal"),
  import.default.data.type = "data.dir",
  import.bed.style = "BisSNP",
  disk.dump.big.matrices = TRUE,
  strand.specific = FALSE,
  filtering.sex.chromosomes.removal = TRUE,
  identifiers.column = "bedFile",
  differential.comparison.columns.all.pairwise = "cellType"
)

# Multiprocess
parallel.setup(num.cores)

#####################################################################
# Import and preprocessing
#####################################################################

if (!file.exists(rnb.set.path)) {
  log_info("Importing the bed files from ", bed.dir)
  data.source <- c(bed.dir, sample.annotation, 1)
  result <- rnb.run.import(
    data.source = data.source, data.type = "bs.bed.dir",
    dir.reports = report.dir
  )
  rnb.set <- result$rnb.set

  # Quality Control
  rnb.run.qc(rnb.set, report.dir)

  # Preprocessing
  rnb.set <- rnb.run.preprocessing(rnb.set, dir.reports = report.dir)$rnb.set

  # Save the object
  save.rnb.set(rnb.set, rnb.set.path, archive = FALSE)
  log_success("Wrote ", rnb.set.path)
} else {
  log_info("Loading the preprocessed object from ", rnb.set.path)
  rnb.set <- load.rnb.set(rnb.set.path)
}

# Drop the unreplicated samples before any group comparison
drop.samples <- intersect(drop.samples, samples(rnb.set))
if (length(drop.samples) > 0) {
  log_info("Removing ", length(drop.samples), " sample(s): ", paste(drop.samples, collapse = ", "))
  rnb.set <- remove.samples(rnb.set, drop.samples)
}
log_info(length(samples(rnb.set)), " samples enter the differential analysis")

#####################################################################
# Differential methylation
#####################################################################

if (!file.exists(diffmeth.path)) {
  log_info("Running differential methylation")
  rnb.run.differential(rnb.set, report.dir)
}
# Loaded in both branches, the original only assigned diffMeth when the
# results already existed, so a first run reached LOLA with no object
diffMeth <- load.rnb.diffmeth(diffmeth.path)

#####################################################################
# LOLA motif enrichment
#####################################################################

if (!file.exists(lola.file)) {
  log_info("Running LOLA against ", lolaDb.path)
  res <- performLolaEnrichment.diffMeth(rnb.set, diffMeth, lolaDb.path)
  saveRDS(res, lola.file)
  log_success("Wrote ", lola.file)
} else {
  log_info("LOLA results already exist: ", lola.file)
}

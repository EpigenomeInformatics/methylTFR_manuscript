#!/usr/bin/env Rscript

#####################################################################
# 01_run_rneads_blueprint.R
# created on 27-10-25 by Irem Gunduz
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
  # Filter out the cancer cells
  sannot <- sannot %>% filter(cellTypeGroup != "cancer")
  sample.annotation <- file.path(bed.dir, "samples_without_cancer.tsv")
  data.table::fwrite(sannot, sample.annotation, sep = "\t")
}
num.cores <- 30

# Directory where the output should be written to
analysis.dir <- "/scratch/icbb/regina/data/blueprint/"
if (!dir.exists(analysis.dir)) dir.create(analysis.dir)
analysis.dir <- file.path(analysis.dir, "RnBeads")
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

if (!file.exists(paste0(report.dir, "/data_import_data/rnb.set_preprocessed"))) {
  data.source <- c(bed.dir, sample.annotation, 1)
  result <- rnb.run.import(data.source = data.source, data.type = "bs.bed.dir", dir.reports = report.dir)
  rnbset <- result$rnb.set

  ## Quality Control
  rnb.run.qc(rnbset, report.dir)

  ## Preprocessing
  rnbset <- rnb.run.preprocessing(rnbset, dir.reports = report.dir)$rnb.set

  ## save the object
  save.rnb.set(rnbset, paste0(report.dir, "/data_import_data/rnb.set_preprocessed"), archive = FALSE)
} else {
  rnbset <- RnBeads::load.rnb.set(paste0(report.dir, "/data_import_data/rnb.set_preprocessed"))
}

## Differential methylation
rnb.run.differential(rnbset, report.dir)

## Exploratory analysis
rnb.run.exploratory(rnbset, report.dir)



#####################################################################
# run_rnbeads for only T an Bcell samples
#####################################################################

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

# Directory where the output should be written to
analysis.dir <- "/scratch/icbb/regina/data/blueprint/"
if (!dir.exists(analysis.dir)) dir.create(analysis.dir)
analysis.dir <- file.path(analysis.dir, "TB_RnBeads")
if (!dir.exists(analysis.dir)) dir.create(analysis.dir)

# Directory where the report files should be written to
report.dir <- file.path(analysis.dir, "reports")
rnb.initialize.reports(report.dir)

# Set up the analysis
rnb.options(
  analysis.name = "Blueprint Healthy Cells",
  assembly = "hg38",
  import.table.separator = "\t",
  region.aggregation = "sum",
  analyze.sites = FALSE,
  region.types = c("tiling1kb", "distal"),
  import.default.data.type = "data.dir",
  import.bed.style = "EPP",
  disk.dump.big.matrices = TRUE,
  strand.specific = FALSE,
  filtering.sex.chromosomes.removal = TRUE,
  differential.enrichment.lola = FALSE, # Running ERROR
  # differential.enrichment.lola.dbs = "/icbb/projects/share/annotations/lolaDB/",
  identifiers.column = "bedFile",
  differential.comparison.columns = "cellTypeGroup" # not exclusive cell-types
)

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

## Differential methylation
diffMeth <- load.rnb.diffmeth(paste0(analysis.dir, "/reports/differential_methylation_data/differential_rnbDiffMeth/"))

#Subset to only include B cells
idx <- rnb.set@pheno[rnb.set@pheno$cellTypeGroup != "Bcell", ]$bedFile
idxt <- rnb.set@pheno[rnb.set@pheno$cellTypeGroup == "Bcell", ]$bedFile
rnbset_bcells <- remove.samples(rnb.set,idx)
rnbset_tcells <- remove.samples(rnb.set,idxt)

# Set report directories
breport.dir <- paste0(report.dir,"/bcell")
treport.dir <- paste0(report.dir,"/tcell")
if(!dir.exists(breport.dir)){dir.create(breport.dir)}
if(!dir.exists(treport.dir)){dir.create(treport.dir)}

# Set options
rnb.options(
  analysis.name = "Blueprint Bcell VS Tcell",
  assembly = "hg38",
  import.table.separator = "\t",
  region.aggregation = "sum",
  region.types = c("cpgislands"),
  import.default.data.type = "data.dir",
  import.bed.style = "EPP",
  analyze.sites = FALSE,
  disk.dump.big.matrices = TRUE,
  strand.specific = FALSE,
  filtering.sex.chromosomes.removal = TRUE,
  differential.enrichment.lola = FALSE,
  #differential.enrichment.lola.dbs = "/icbb/projects/share/annotations/lolaDB",
  identifiers.column = "bedFile",
  differential.comparison.columns = "cellTypeShort" #  exclusive cell-types
)

# Differential methylation
rnb.run.differential(rnbset_bcells, breport.dir)
rnb.run.differential(rnbset_tcells, treport.dir)
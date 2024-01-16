#!/usr/bin/env Rscript

#####################################################################
# 03_subset_analysis.R
# created on 2023-11-09 by Irem Gunduz
# Run RnBeads vanilla analysis for Blueprint Methylation Data
# Subset to only include B cells and T cells
#####################################################################

set.seed(42)
muLogR::logger.info("Loading libraries...")
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(logger)
  library(muLogR)
  library(RnBeads)
})

# Directory where your data is located
data.dir <- "/scratch/icbb/mtfr_manuscript/"
bed.dir <- file.path(data.dir, "blueprint_data")
sample.annotation <- file.path(bed.dir, "samples_subset.tsv")
num.cores <- 30

# Directory where the output should be written to
analysis.dir <- "/scratch/icbb/mtfr_manuscript/"
if (!dir.exists(analysis.dir)) dir.create(analysis.dir)
analysis.dir <- file.path(analysis.dir, "BLUEPRINT")
if (!dir.exists(analysis.dir)) dir.create(analysis.dir)

# Multiprocess
parallel.setup(30)

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
  differential.enrichment.lola = TRUE,
  differential.enrichment.lola.dbs = "/icbb/projects/share/annotations/lolaDB",
  identifiers.column = "bedFile",
  differential.comparison.columns = "cellTypeShort" #  exclusive cell-types
)


# Directory where the report files should be written to
report.dir <- file.path(analysis.dir, "reports")
rnb.initialize.reports(report.dir)

# Multiprocess
parallel.setup(num.cores)

# Set up the analysis
rnb.options(
  analysis.name = "Blueprint Healthy Cells",
  assembly = "hg38",
  import.table.separator = "\t",
  region.aggregation = "sum",
  region.types = c("cpgislands", "promoters", "genes", "tiling"),
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
} 

## Exploratory analysis
rnb.run.exploratory(rnbset, report.dir)


#Subset to only include B cells
idx <- rnb_set@pheno[rnb_set@pheno$cellTypeGroup != "Bcell", ]$bedFile
idxt <- rnb_set@pheno[rnb_set@pheno$cellTypeGroup == "Bcell", ]$bedFile
rnb_set_bcells <- remove.samples(rnb_set,idx)
rnb_set_tcells <- remove.samples(rnb_set,idxt)

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
  differential.enrichment.lola = TRUE,
  differential.enrichment.lola.dbs = "/icbb/projects/share/annotations/lolaDB",
  identifiers.column = "bedFile",
  differential.comparison.columns = "cellTypeShort" #  exclusive cell-types
)

# Differential methylation
rnb.run.differential(rnb_set_bcells, breport.dir)
rnb.run.differential(rnb_set_tcells, treport.dir)
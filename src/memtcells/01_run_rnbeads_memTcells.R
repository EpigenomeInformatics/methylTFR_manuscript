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
data.dir <- "/icbb/projects/share/datasets/"
bed.dir <- file.path(data.dir, "memoryTcells")
sample.annotation <- file.path(bed.dir, "samples.tsv")
num.cores <- 30

# Directory where the output should be written to
analysis.dir <- "/icbb/projects/igunduz/methylTFR_manuscript/results"
if (!dir.exists(analysis.dir)) dir.create(analysis.dir)
analysis.dir <- file.path(analysis.dir, "memoryTcells_060525")
if (!dir.exists(analysis.dir)) dir.create(analysis.dir)

tiling1kb <- muRtools::getTilingRegions("hg38", width=1000L, onlyMainChrs=TRUE)%>%
  data.table::as.data.table() %>%
  dplyr::select(seqnames, start, end) %>%
  as.data.frame()
colnames(tiling1kb) <- c("Chromosome", "Start", "End")
rnb.set.annotation(type = "tiling1kb", regions = tiling1kb, assembly = "hg38")

distal <- readRDS("/icbb/projects/share/annotations/methylTFRAnnotationHg38/inst/extdata/distal_regions.RDS")
distal <- as.data.frame(distal) %>%
  dplyr::select(seqnames, start, end)
colnames(distal) <- c("Chromosome", "Start", "End")
rnb.set.annotation(type = "distal", regions = distal, assembly = "hg38")

# Directory where the report files should be written to
report.dir <- file.path(analysis.dir, "reports")
rnb.initialize.reports(report.dir)
rnb.options(
  analysis.name = "CD4 T memory cell subtypes  Methylation Data",
  assembly = "hg38",
  import.table.separator = "\t",
  disk.dump.big.matrices = TRUE,
  strand.specific = FALSE,
  region.types = c("tiling1kb","distal"),
  #c("cpgislands", "promoters", "genes", "tiling",),
  import.default.data.type = "data.dir",
  import.bed.style = "BisSNP",
  filtering.sex.chromosomes.removal = TRUE,
  identifiers.column = "bedFile",
  #differential.comparison.columns = "cellType",
  differential.comparison.columns.all.pairwise = "cellType"#,
  #covariate.adjustment.columns=c("cmp_individual")
)

# Multiprocess
parallel.setup(num.cores)

#data.source <- c(bed.dir, sample.annotation, 1)
#result <- rnb.run.import(data.source = data.source, data.type = "bs.bed.dir", dir.reports = report.dir)
#rnbset <- result$rnb.set

### Quality Control
#rnb.run.qc(rnbset, report.dir)

## Preprocessing
#rnbset <- rnb.run.preprocessing(rnbset, dir.reports = report.dir)$rnb.set

## save the object
#save.rnb.set(rnbset, paste0(report.dir, "/data_import_data/rnb.set_preprocessed"), archive = FALSE)
rnbset <- load.rnb.set(paste0(report.dir, "/data_import_data/rnb.set_preprocessed"))

## Differential methylation
rnbset_noTEMRA <- remove.samples(rnbset, "51_Hf03_BlTR_Ct_WGBS_S_1.MCSv3.20170714.GRCh38.cpg.filtered.CG.bed")
report.dir <- file.path(analysis.dir, "reports_080725")

rnb.run.differential(rnbset_noTEMRA, report.dir)

#####################################################################
# LOLA Analysis
#######################################################################

logger.start("Running LOLA")
lolaDb_path <- "/icbb/projects/share/annotations/lolaDB/hg38/"

logger.info("Loading RnBeads objects")
diffMeth <- load.rnb.diffmeth(paste0(analysis.dir, "/reports_080725/differential_methylation_data/differential_rnbDiffMeth/"))
rnb_set <- load.rnb.set(paste0(analysis.dir, "/reports/data_import_data/rnb.set_preprocessed/"))

# Run LOLA
res <- performLolaEnrichment.diffMeth(rnb_set, diffMeth, lolaDb_path)
logger.info("Saving results")
saveRDS(res, paste0(analysis.dir, "/reports_080725/differential_methylation_data/differential_rnbDiffMeth/lola_results.rds"))
logger.completed()

#######################################################################
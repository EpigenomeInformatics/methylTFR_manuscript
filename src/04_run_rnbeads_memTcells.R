#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(RnBeads))
# Directory where your data is located
data.dir <- "/icbb/projects/skumar/memoryTcells"
bed.dir <- file.path(data.dir, "bed")
sample.annotation <- file.path(bed.dir, "samples.tsv")
# Directory where the output should be written to
analysis.dir <- "/icbb/projects/skumar/memoryTcells/analysis"
# Directory where the report files should be written to
report.dir <- file.path(analysis.dir, "reports")
rnb.initialize.reports(report.dir)
rnb.options(analysis.name = "CD4 T memory cell subtypes  Methylation Data",
            assembly = "hg38",
            import.table.separator = "\t",
            region.types = c("cpgislands", "promoters", "genes", "tiling"),
            import.default.data.type = "data.dir",
            import.bed.style = "BisSNP",
            filtering.sex.chromosomes.removal = FALSE, 
            identifiers.column = "bedFile")

data.source <- c(bed.dir, sample.annotation, 1)
result <- rnb.run.import(data.source=data.source, data.type="bs.bed.dir", dir.reports=report.dir)
rnbset <- result$rnb.set

## Quality Control
rnb.run.qc(rnbset, report.dir)

## Preprocessing
rnbset <- rnb.run.preprocessing(rnbset, dir.reports=report.dir)$rnb.set

rnbsSaveDir <- file.path(report.dir, "rnbs_wgbs_tmem_bt")
save.rnb.set(rnbset, rnbsSaveDir, archive=FALSE)

## Data export
rnb.options(export.to.csv = TRUE)
rnb.run.tnt(rnbset, report.dir)

## Exploratory analysis
rnb.run.exploratory(rnbset, report.dir)

## Differential methylation
rnb.run.differential(rnbset, report.dir)
save.rnb.set(rnbset, rnbsSaveDir, archive=FALSE)

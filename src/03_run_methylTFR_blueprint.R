suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(methylTFR)
  library(methylTFRAnnotationHg38)
})

data.dir <- "/icbb/projects/igunduz/methylTFR_manuscript/data/BLUEPRINT/"
bed_files <- list.files(data.dir,pattern=".bed")

sa <- fread("/icbb/projects/skumar/blueprint/bed/samples.tsv")
head(sa)

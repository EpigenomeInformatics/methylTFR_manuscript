suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(methylTFR)
  library(methylTFRAnnotationHg38)
})

data.dir <- "/icbb/projects/igunduz/methylTFR_manuscript/data/BLUEPRINT/"
bed_files <- list.files(data.dir, pattern = ".bed")


sample_dir <- "/icbb/projects/skumar/memoryTcells/bed"
sample_ann <- "samples.tsv"
deviations <- run_methyltfr(sample_ann,
  sample_dir,
  filetype = "BisSNP",
  threads = 16
)

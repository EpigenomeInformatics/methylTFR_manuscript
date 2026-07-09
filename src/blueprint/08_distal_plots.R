set.seed(42)
suppressPackageStartupMessages({
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ggplot2)
library(matrixStats)
library(SummarizedExperiment)
library(methylTFRAnnotationHg38)
library(GenomicRanges)
})

enhancer <- readRDS("/icbb/projects/igunduz/annotation/methylTFRAnnotationHg38/inst/extdata/distal_regions.RDS")
tf_bindsites <- getTFbindsites(motifSet = "JASPAR2020")
all_tfbs_combined <- unlist(tf_bindsites)

# Annotate all sites vs. Distal sites
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
anno_all <- annotatePeak(all_tfbs_combined, TxDb = txdb, verbose = FALSE)
anno_distal <- annotatePeak(subsetByOverlaps(all_tfbs_combined, enhancer), TxDb = txdb, verbose = FALSE)

# Plot comparison
p <- plotAnnoBar(list(All_TFBS = anno_all, Distal_Only = anno_distal)) +
  theme_classic() 

ggsave("/scratch/icbb/igunduz/methylTFR_manuscript/blueprint/mTFR_devs_121125/distal_vs_all_TFBS_annotation_barplot.png",
       plot = p, width = 6, height = 4)
#!/usr/bin/env Rscript

#####################################################################
# 04.plot_blueprint_pca.R
# created on 21-07-25 by Irem Gunduz
# Plot PCA of  methylTFR results
#####################################################################
set.seed(42)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(methylTFR)
  library(muLogR)
  library(gplots)
  library(factoextra)
  library(ggfortify)
})

plot_dir <- "/scratch/icbb/regina/methylTFR_manuscript/figures/memTcells"
if(!dir.exists(plot_dir)) dir.create(plot_dir)
dev_obj <- readRDS("//icbb/projects/igunduz/methylTFR_manuscript/results/memoryTcells/mtfr_final_221223/jaspar2020_distal_deviations.RDS")
motifset <- "jaspar2020_distal"
deviations <- deviations(deviations)
sannot <- fread("/icbb/projects/skumar/memoryTcells/bed/samples.tsv")

# Get cell types
tdf <- as.data.frame(t(deviations))
pca <- prcomp(tdf, scale. = T)
tdf$cell_type <- sannot$cellType


# Save PCA individual plot
fn_pca_ind <- file.path(plot_dir, paste0("pca_memTcells_", motifset, ".pdf"))
pdf(fn_pca_ind)
fviz_pca_ind(res.pca, repel = TRUE, title = "PCA on all memTcells samples")
dev.off()


#!/usr/bin/env Rscript

#####################################################################
# cluster_analysis.R
# Created by RN, edited by IBG on 27-10-25
# Perform clustering analysis on Blueprint methylTFR and compare with RnBeads
#####################################################################

suppressPackageStartupMessages({
  library(methylTFR)
  library(RnBeads)
  library(mclust)
  library(dplyr)
  library(ComplexHeatmap)
  library(Cairo)
})
set.seed(13)
source("/icbb/projects/nitschre/methylTFR/scripts/other/helpers.R")
fig_dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR_manuscript/figures/blueprint/"

# Get mtfr, distal and 1kbtiling matrix
mtfr <- readRDS("/icbb/projects/nitschre/methylTFR/r_objects/blueprintjaspar2020_distal_deviations.RDS")
mtfr <- deviationZScores(mtfr)

rnbeads_distal <- load.rnb.set("/icbb/projects/nitschre/methylTFR/results/blueprint_distal/reports/data_import_data/rnb.set_preprocessed/")
rnbeads_tiling1kb <- load.rnb.set("/icbb/projects/igunduz/methylTFR_manuscript/results/BLUEPRINT_080725/reports/data_import_data/rnb.set_preprocessed/")

distal <- meth(rnbeads_distal, type = "distal")
tiling <- meth(rnbeads_tiling1kb, type = "tiling1kb")

# K means clustering
mtfr_res <- kmeans(t(mtfr), centers = 15)
distal_res <- kmeans(t(distal), centers = 15)
tiling_res <- kmeans(t(tiling), centers = 15)

## Adjusted Rank Index
# Get real cell types
sannot <- read.csv("/icbb/projects/igunduz/methylTFR_manuscript/results/BLUEPRINT_080725/reports/data_import_data/annotation.csv", stringsAsFactors = FALSE)
sannot$bedFile <- as.character(sannot$bedFile)
sample_names <- colnames(mtfr)

# Match and extract cellTypeGroup
cell_types <- sannot$cellTypeGroup[match(sample_names, sannot$bedFile)]

# Calculate ARI
ari_mtfr <- adjustedRandIndex(mtfr_res$cluster, cell_types)
ari_distal <- adjustedRandIndex(distal_res$cluster, cell_types)
ari_tiling <- adjustedRandIndex(tiling_res$cluster, cell_types)

# construct confusion matrix
results <- list(mtfr = mtfr_res) # distal=distal_res, tiling=tiling_res)
aris <- list(ari_mtfr) # ari_distal, ari_tiling)
for (i in seq_along(results)) {
  df <- data.frame(
    sample = names(results[[i]]$cluster),
    cluster = unname(results[[i]]$cluster),
    cell_type = cell_types
  )
  df <- table(df$cluster, df$cell_type)
  cM <- Jaccard(df)
  cM <- prettyOrderMat(t(cM), clusterCols = TRUE)$mat %>% t()
  whitePurple <- c(
    "#f7fcfd", "#e0ecf4", "#bfd3e6", "#9ebcda", "#8c96c6",
    "#8c6bb1", "#88419d", "#810f7c", "#4d004b"
  )
  ht_opt$simple_anno_size <- unit(0.25, "cm")

  # Plot the heatmap
  hm <- BORHeatmap(
    cM,
    dataColorMidPoint = 0.4,
    labelCols = TRUE, labelRows = TRUE,
    dataColors = whitePurple,
    showColDendrogram = F,
    showRowDendrogram = F,
    row_names_side = "left",
    width = ncol(cM) * unit(0.5, "cm"),
    height = nrow(cM) * unit(0.5, "cm"),
    border_gp = gpar(col = "black"),
    borderColor = "#a9a9a9",
    use_raster = TRUE
  )
  pdf(paste0(fig_dir, names(results)[i], "_jaccard_heatmap.pdf"), width = 6, height = 6)
  draw(hm)
  grid.text(
    paste("adjusted rank index =", round(aris[[i]], 3)),
    x = 0.4, y = 0.98, gp = gpar(fontsize = 12, col = "black", fontface = "bold")
  )
  dev.off()
}

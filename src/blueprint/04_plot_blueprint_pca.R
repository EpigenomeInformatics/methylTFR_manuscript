#!/usr/bin/env Rscript

#####################################################################
# 04.plot_blueprint_pca.R
# created on 21-07-25 by Irem Gunduz
# Plot blueprint methylTFR results
#####################################################################

set.seed(42)
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(methylTFR)
  library(muLogR)
  library(gplots)
  library(muRtools)
  library(ComplexHeatmap)
  library(factoextra)
  library(ggfortify)
  library(RnBeads)
})
source("/icbb/projects/igunduz/methylTFR_manuscript/src/utils.R")

cell_type_colors <- c(
  "B-cells" = "#1f77b4",
  "DC" = "#ff7f0e",
  "Erythrocytes" = "#2ca02c",
  "Granulocytes" = "#d62728",
  "Megakaryocytes" = "#9467bd",
  "Mf" = "#8c564b",
  "Monocytes" = "#e377c2",
  "NK" = "#7f7f7f",
  "Osteoclast" = "#1b9e77",
  "Other" = "#bcbd22",
  "Plasma" = "#17becf",
  "Progenitors" = "#2ca4a2",
  "T-cells" = "#ff7f0e",
  "Thymocyte" = "#6a5acd"
)

# Remap to cleaner group names
group_remap <- c(
  "Bcell" = "B-cells",
  "DC" = "DC",
  "eryt" = "Erythrocytes",
  "gran" = "Granulocytes",
  "megK" = "Megakaryocytes",
  "Mf" = "Mf",
  "mono" = "Monocytes",
  "NK" = "NK",
  "osteoclast" = "Osteoclast",
  "other" = "Other",
  "plasma" = "Plasma",
  "progenitor" = "Progenitors",
  "Tcell" = "T-cells",
  "thymocyte" = "Thymocyte"
)

# Paths
plot_dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR_manuscript/figures/blueprint/"
if (!dir.exists(plot_dir)) dir.create(plot_dir)
deviations <- readRDS("/scratch/icbb/igunduz/methylTFR_manuscript/blueprint/mTFR_devs_071125/jaspar2020_distal_deviations.RDS")
deviations <- deviations(deviations)
sannot <- read.csv("/scratch/icbb/igunduz/methylTFR_manuscript/blueprint/RnBeads_291025/reports/data_import_data/annotation.csv", stringsAsFactors = FALSE)

# Ensure consistent character types
sannot$bedFile <- as.character(sannot$bedFile)
sample_names <- colnames(deviations)

# Match and extract cellTypeGroup
cell_types <- sannot$cellTypeGroup[match(sample_names, sannot$bedFile)]

# Apply PCA
tdf <- as.data.frame(t(deviations))
pca_result <- prcomp(tdf, center = FALSE, scale. = FALSE)
tdf$groups <- cell_types

# Apply remap
tdf$groups <- group_remap[tdf$groups]


fig_path <- paste0(plot_dir, "PCA_deviations_blueprint.pdf")
pdf(fig_path, width = 10, height = 10)
# Plot PCA
autoplot(pca_result,
  data = tdf,
  colour = "groups",
  main = "PCA",
  size = 5
) +
  theme_classic() +
  scale_color_manual(values = cell_type_colors) +
  theme(legend.position = "bottom")
dev.off()

#####################################################################
# tiling and distal PCAs
#####################################################################

# Load RnBeads objects
rnbeads <- load.rnb.set("/scratch/icbb/igunduz/methylTFR_manuscript/blueprint/RnBeads_291025/reports/data_import_data/rnb.set_preprocessed")

# Extract methylation matrices
distal <- meth(rnbeads, type = "distal")
tiling <- meth(rnbeads, type = "tiling1kb")

# Apply PCA
tdf_tiling <- as.data.frame(t(tiling))
pca_result_tiling <- prcomp(tdf_tiling, center = FALSE, scale. = FALSE)

# Map cell types based on samples
sample_names <- colnames(tiling)
cell_types <- sannot$cellTypeGroup[match(sample_names, sannot$bedFile)]
tdf_tiling$groups <- cell_types
# Apply remap
tdf_tiling$groups <- group_remap[tdf_tiling$groups]

fig_path <- paste0(plot_dir, "PCA_tiling_blueprint.pdf")
pdf(fig_path, width = 10, height = 10)
# Plot PCA for tiling1kb
autoplot(pca_result_tiling,
  data = tdf_tiling,
  colour = "groups",
  main = "PCA Tiling 1kb",
  size = 5
) +
  theme_classic() +
  scale_color_manual(values = cell_type_colors) +
  theme(legend.position = "bottom")
dev.off()

# Apply PCA for distal
tdf_distal <- as.data.frame(t(distal))
pca_result_distal <- prcomp(tdf_distal, center = FALSE, scale. = FALSE)
# Map cell types based on samples
sample_names <- colnames(distal)
cell_types <- sannot$cellTypeGroup[match(sample_names, sannot$bedFile)]
tdf_distal$groups <- cell_types
# Apply remap
tdf_distal$groups <- group_remap[tdf_distal$groups]
fig_path <- paste0(plot_dir, "PCA_distal_blueprint.pdf")
pdf(fig_path, width = 10, height = 10)
# Plot PCA for distal
autoplot(pca_result_distal,
  data = tdf_distal,
  colour = "groups",
  main = "PCA Distal",
  size = 5
) +
  theme_classic() +
  scale_color_manual(values = cell_type_colors) +
  theme(legend.position = "bottom")
dev.off()

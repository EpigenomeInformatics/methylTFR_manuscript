#!/usr/bin/env Rscript

#####################################################################
# 04.plot_blueprint.R
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

# Paths
plot_dir <- "/icbb/projects/nitschre/methylTFR/figures/blueprint/"
if(!dir.exists(plot_dir)) dir.create(plot_dir)
deviations <- readRDS("/icbb/projects/nitschre/methylTFR/r_objects/blueprintjaspar2020_distal_deviations.RDS")
deviations <- deviations(deviations)
#rnbeads_path <- "/icbb/projects/igunduz/methylTFR_manuscript/results/BLUEPRINT_080725/reports/data_import_data/rnb.set_preprocessed/"
sannot <- read.csv("/icbb/projects/igunduz/methylTFR_manuscript/results/BLUEPRINT_080725/reports/data_import_data/annotation.csv", stringsAsFactors = FALSE)

# Get the tiling matrix
#rnbset <- RnBeads::load.rnb.set(rnbeads_path)
#tiling_matrix <- meth(rnbset,type="tiling1kb")

# Ensure consistent character types
sannot$bedFile <- as.character(sannot$bedFile)
sample_names <- colnames(deviations)

# Match and extract cellTypeGroup
cell_types <- sannot$cellTypeGroup[match(sample_names, sannot$bedFile)]

#tdf <- as.data.frame(t(tiling_matrix))
#pca_result <- prcomp(tdf, center = FALSE, scale. = FALSE)
#tdf$groups <- cell_types

# Remove cancer samples
keep_samples <- cell_types != "cancer"
deviations_filtered <- deviations[, keep_samples]
cell_types_filtered <- cell_types[keep_samples]

# Start from the filtered objects (after removing cancer)
valid_samples <- !is.na(cell_types_filtered)
deviation_matrix_final <- deviations_filtered[, valid_samples]
cell_types_final <- cell_types_filtered[valid_samples]

# PCA
cell_type_colors <- c(
  "B-cells" = "#1f77b4",
  "DC" = "#ff7f0e",
  "Erythrocytes" = "#2ca02c",
  "Granulocytes" = "#d62728",
  "Megakaryocytes" = "#9467bd",
  "Mf" = "#8c564b",
  "Monocytes" = "#e377c2",
  "NK" = "#7f7f7f",
  "Other" = "#bcbd22",
  "Plasma" = "#17becf",
  "Progenitors" = "#2ca4a2",
  "T-cells" = "#ff7f0e",
  "Thymocyte" = "#6a5acd"
)

# Apply PCA
tdf <- as.data.frame(t(deviation_matrix_final))
pca_result <- prcomp(tdf, center = FALSE, scale. = FALSE)
tdf$groups <- cell_types_filtered

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
  "other" = "Other",
  "plasma" = "Plasma",
  "progenitor" = "Progenitors",
  "Tcell" = "T-cells",
  "thymocyte" = "Thymocyte"
)

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


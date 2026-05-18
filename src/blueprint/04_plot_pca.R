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

cell_type_colors <- c(
  "B-cells" = "#980043",
  "Plasma" = "#CD2990",
  "DC" = "#EED5B7",
  "Granulocytes" = "#ff7f50",
  "Erythrocytes" = "#67000d",
  "Mf" = "#864a38",
  "Monocytes" = "#CD7054",
  "Megakaryocytes" = "#4a0221",
  "NK" = "#bf812d",
  "Osteoclast" = "#DEB887",
  "T-cells" = "#40E0D0",
  "Thymocyte" = "#74c476",
  "Progenitors" = "#F4C3DA"
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
  "plasma" = "Plasma",
  "progenitor" = "Progenitors",
  "Tcell" = "T-cells",
  "thymocyte" = "Thymocyte"
)

source("/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR_manuscript/src/utils.R")

# Paths
plot_dir <- "/icbb/projects/nitschre/methylTFR/figures/blueprint/PCA/"
if (!dir.exists(plot_dir)) dir.create(plot_dir)

# Load data
deviations <- readRDS("/scratch/icbb/igunduz/methylTFR_manuscript/blueprint/mTFR_devs_121125/JASPAR2020_distal_deviations.RDS")
deviations_noCorrection <- readRDS("/icbb/projects/nitschre/methylTFR/r_objects/JASPAR2020_distal_deviations_bp_noCorrection.RDS")
deviations <- deviations(deviations)
deviations_noCorrection <- deviations(deviations_noCorrection)

sannot <- read.csv("/scratch/icbb/igunduz/methylTFR_manuscript/blueprint/RnBeads_291025/reports/data_import_data/annotation.csv", stringsAsFactors = FALSE)

# Ensure consistent character types
sannot$bedFile <- as.character(sannot$bedFile)
sample_names <- colnames(deviations)

# Match and extract cellTypeGroup
cell_types <- sannot$cellTypeGroup[match(sample_names, sannot$bedFile)]

# Remove 'other' cell types
rm <- which(cell_types == "other")
deviations <- deviations[,-rm]
deviations_noCorrection <- deviations_noCorrection[,-rm]
cell_types <- cell_types[-rm]

#####################################################################
# Plot pie chart for cell types
#####################################################################
# Count samples per cell type
sample_counts <- table(cell_types)

# Change sample names
names(sample_counts) <- group_remap[names(sample_counts)]

# Cell order
order <- c("Megakaryocytes", "Erythrocytes", "Granulocytes", "Monocytes", "Mf", "DC", "Osteoclast", "NK", "T-cells", "Thymocyte", "B-cells", "Plasma","Progenitors")
sample_counts <- sample_counts[order]
cell_type_colors <- cell_type_colors[order]

# Calculate the percentage of each cell type
percentages <- round(sample_counts / sum(sample_counts) * 100, 1)

# Create labels
labels <- paste(names(sample_counts), "\n", percentages, "%", sep = "")

# Create the pie chart
pdf(file.path(plot_dir, "blueprint_pie.pdf"), width = 15, height = 15, onefile = FALSE)
pie(sample_counts, main = "Number of samples per cell type", col = cell_type_colors, labels = labels)
dev.off()

#####################################################################
# mTFR scores PCA
#####################################################################
deviations_list <- list(
  "PCA_deviations_blueprint" = deviations,
  "PCA_deviations_blueprint_noCorrection" = deviations_noCorrection
)

# PCA for both corrected and uncorrected deviations
lapply(names(deviations_list), function(name){
  deviations <- deviations_list[[name]]
  
  # Apply PCA
  tdf <- as.data.frame(t(deviations))
  pca_result <- prcomp(tdf, center = FALSE, scale. = FALSE)
  tdf$groups <- cell_types

  # Apply remap
  tdf$groups <- group_remap[tdf$groups]

  # Plot
  fig_path <- paste0(plot_dir, name, ".pdf")
  pdf(fig_path, width = 10, height = 10)
  p <- autoplot(pca_result,
      data = tdf,
      colour = "groups",
      main = "PCA",
      size = 5
    ) +
      theme_classic() +
      scale_color_manual(values = cell_type_colors) +
      theme(legend.position = "bottom")
  print(p)
  dev.off()
})

#####################################################################
# tiling and distal PCAs
#####################################################################

# Load RnBeads objects
rnbeads <- load.rnb.set("/scratch/icbb/igunduz/methylTFR_manuscript/blueprint/RnBeads_291025/reports/data_import_data/rnb.set_preprocessed")

# Extract methylation matrices
distal <- meth(rnbeads, type = "distal")
tiling <- meth(rnbeads, type = "tiling1kb")

# Apply PCA for tiling
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


#####################################################################
# PCA for T and B subtypes
#####################################################################
# New colors
cell_type_colors <- c(
  "Bcell_pre" = "#201960",
  "Bcell_naive" = "#506CB2",
  "TCD8" = "#FE8320",
  "Bcell_gc" = "#CBDDE9",
  "TCD8_effM" = "#CC3F0A",
  "Bcell_mem" = "#80B1CB",
  "TCD8_cenM" = "#F15F10",
  "Neut_mat" = "#7f7f7f",
  "TCD4_effM" = "#FCA737",
  "TCD8_term" = "#972605",
  "TCD4_cenM" = "#F9CF55",
  "T-TCD4" = "#F6EB8B",
  "Treg" = "#5E1603")

# Subset to T and B cell types
deviations_sub <- deviations[,grep("TC|Treg|Bcell", colnames(deviations))]
sample_names <- colnames(deviations_sub)

# Match and extract cellTypeShort 
cell_types <- sannot$cellTypeShort[match(sample_names, sannot$bedFile)]

# Apply PCA
tdf <- as.data.frame(t(deviations_sub))
pca_result <- prcomp(tdf, center = FALSE, scale. = FALSE)
tdf$groups <- cell_types

fig_path <- paste0(plot_dir, "PCA_deviations_subtypes.pdf")
pdf(fig_path, width = 10, height = 10)
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
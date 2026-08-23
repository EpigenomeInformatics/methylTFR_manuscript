#!/usr/bin/env Rscript

#####################################################################
# 05_variabilty_analysis.R
# created on 07.04.2026 
# Calculating heatmap of most variable TFs in blueprint data
#####################################################################
suppressPackageStartupMessages({
  library(methylTFR)
  library(ComplexHeatmap)
  library(dplyr)
  library(circlize)
  library(ggplot2)
})


set.seed(13)

# Set up
plot_dir <- "/scratch/icbb/regina/methylTFR_manuscript/figures/blueprint/"
results_dir <- "/scratch/icbb/regina/data/blueprint/"

# Scripts
source("/icbb/projects/nitschre/methylTFR/scripts/other/variability_analysis.R")
source("/icbb/projects/nitschre/methylTFR/scripts/other/helpers.R")

# Cell type color
cell_type_colors <- c(
  "Bcell" = "#980043",
  "plasma" = "#CD2990",
  "DC" = "#EED5B7",
  "gran" = "#ff7f50",
  "Mf" = "#864a38",
  "mono" = "#CD7054",
  "Tcell" = "#40E0D0",
  "thymocyte" = "#74c476"
)

# Loading deviations scores and remove problematic sample
dev_obj <- readRDS("/scratch/icbb/igunduz/methylTFR_manuscript/blueprint/mTFR_devs_121125/JASPAR2020_distal_deviations.RDS")
dev_obj <- dev_obj[,!colnames(dev_obj) %in% "Bcell_naive_VB_NBC_NC11_83.bed"]
dev_obj <- dev_obj[,colnames(dev_obj) %in% c("plasma", "gran", "Bcell", "DC","Mf","mono","Tcell","thymocyte")] # Keep only samples that have at least 5 samples

deviations <- deviations(dev_obj)

# Change sample names to cellTypeGroup
annot <- read.csv("/scratch/icbb/regina/methylTFR_manuscript/data/blueprint/RnBeads/reports/data_import_data/annotation.csv", stringsAsFactors = FALSE)
annot <- subset(annot, bedFile %in% colnames(deviations))

# Match order of annotation samples with deviation matrix and remove "other" sample
deviations <- deviations[,match(annot$bedFile, colnames(deviations))]
colnames(deviations) <- annot$cellTypeGroup
deviations <- deviations[,!colnames(deviations) %in% c("other")]

# Compute zscores columnwise
zscores_col <- computeColZScore(deviations)

# Perform variability analysis
if(file.exists(paste0(results_dir, "variability_bp.RDS"))) {  
  var <- readRDS(paste0(results_dir, "variability_bp.RDS"))
} else {
  var <-computeZScoreVariability(zscores_col)
  saveRDS(var, file=paste0(results_dir, "variability_bp.RDS"))
}

# Filter for  pval < .05 and keep top 50 with smallest pvalue
top50 <- var %>%
  filter(p_value_adj < 0.05) %>%
  arrange(desc(variability)) %>%
  head(50)

# Filter for TFs that are in the top 50 variables
deviations_filtered <- deviations[rownames(deviations) %in% c(rownames(top50)),]

# Calculate Row Zscores
zscores <-  methylTFR:::computeRowZScore(deviations_filtered)

# Heatmap Annotation at the top of the plot
ha <- HeatmapAnnotation(
  celltypes=colnames(zscores),
  col = list(celltypes = cell_type_colors)
)

# Color scheme
col <- muRtools::colpal.cont(100, "cptcity.arendal_temperature")

# Column order 
column_order <- c("megK", "eryt", "gran", "mono", "Mf", "DC", "osteoclast", "NK", "Tcell", "thymocyte", "Bcell", "plasma")

# Column split factor
group <- colnames(zscores)
column_split_factor <- factor(group, levels=column_order)

# Clipped matrix
zscores_clipped <- pmax(pmin(zscores,2),-2)

# Plot heatmap
ht <- Heatmap(
  zscores_clipped, column_names_gp = gpar(fontsize = 9),
  top_annotation = ha,
  column_title = "Z-Scores of variable TFs in Blueprint",
  heatmap_legend_param = list(title = "methylTFR Z-scores"),
  show_row_names = TRUE,
  show_column_names = FALSE,
  column_split = column_split_factor, # Split columns by the desired order
  cluster_column_slices = FALSE, # Prevent clustering of the groups themselves)
  col = col)
pdf(paste0(plot_dir, "heatmap_variable_TFs_bp.pdf"))
draw(ht)
dev.off()

# Variability Rank Plot
p <- plotVariability(var)
pdf(paste0(plot_dir, "variability_rank_plot.pdf"))
print(p)
dev.off()


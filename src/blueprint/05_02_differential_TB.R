#!/usr/bin/env Rscript
#####################################################################
# differentials_heatmap.R
# Created by RN on 13-11-2025
# Script to plot heatmap of differential TFs between Tcells and Bcells in blueprint data
#####################################################################

suppressPackageStartupMessages({
  library(methylTFR)
  library(ComplexHeatmap)
  library(dplyr)
  library(circlize)
})

set.seed(12)

# Set up
plot_dir <- "/icbb/projects/nitschre/methylTFR/figures/blueprint/Heatmaps_differentials/"
results_dir <- "/icbb/projects/nitschre/methylTFR/r_objects/"
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# Loading deviations scores and remove problematic sample
deviations_raw <- readRDS("/icbb/projects/nitschre/methylTFR/r_objects/jaspar2020_distal_deviations_TvsB.RDS")
deviations_raw <- deviations_raw[,!colnames(deviations_raw) %in% "Bcell_naive_VB_NBC_NC11_83.bed"]
deviations <- deviations(deviations_raw)

# Aggregating cell names 
colnames(deviations) <- case_when(
  grepl("Bcell", colnames(deviations), ignore.case = TRUE) ~ "Bcell",
  grepl("TC", colnames(deviations), ignore.case = TRUE) ~ "Tcell"
)

# Set groups to compare
groups <- colnames(deviations)

# Perform differential analysis
if(file.exists(paste0(results_dir, "diff_BvsT.RDS"))) {  
  diff <- readRDS(paste0(results_dir, "diff_BvsT.RDS"))
} else {
  diff <- differential_deviation_test(deviations, groups=groups, 
                                     alternative="two.sided", parametric = TRUE)
  saveRDS(diff, file=paste0(results_dir, "diff_BvsT.RDS"))
}


# Filter for  pval < .05 and keep top 50 with smallest pvalue
top50 <- diff %>%
  filter(p_value_adjusted < 0.05) %>%
  arrange(p_value_adjusted) %>%
  head(50)

# Get Z-scores
zscores <- deviationZScores(deviations_raw)

# Filter for TFs that are in the top 50 differentials of BvsTcells
zscores_filtered <- zscores[rownames(zscores) %in% c(rownames(top50)),]

# Heatmap Annotation at the top of the plot
ha <- HeatmapAnnotation(
  celltypes=colnames(deviations),
  col = list(celltypes = c("Bcell" = "#2B4B9B", "Tcell"="#D76016"))
)

# Color scheme
col <- muRtools::colpal.cont(100, "cptcity.arendal_temperature")

# Plot heatmap
ht <- Heatmap(
  zscores_filtered, column_names_gp = gpar(fontsize = 9),
  top_annotation = ha,
  column_title = "Z-Scores of differential TFs in Tcells vs Bcells",
  heatmap_legend_param = list(title = "methylTFR Z-scores"),
  show_row_names = TRUE,
  show_column_names = FALSE,
  col = col)
pdf(paste0(plot_dir, "heatmap_differential_BvsTcells.pdf"))
draw(ht)
dev.off()

## Additional columns
# Mean z-score difference column
mean_bcells <- rowMeans(zscores_filtered[,grep("Bcell", colnames(zscores_filtered))])
mean_tcells <- rowMeans(zscores_filtered[,grep("TCD", colnames(zscores_filtered))])
zscores_diff <- as.matrix(mean_bcells - mean_tcells)

row_order <- ht %>%
  draw() %>%
  row_order()

# Draw column
pdf(paste0(plot_dir, "meanZscoreDiff_BvsTcells.pdf"))
Heatmap(
  zscores_diff, column_names_gp = gpar(fontsize = 9),
  cluster_rows = FALSE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  row_order = row_order,
  heatmap_legend_param = list(title = "Mean z-score difference"),
)
dev.off()

# Column for pvalue adj
pval <- as.matrix(diff[rownames(diff) %in% rownames(zscores_filtered), "p_value_adjusted", drop = FALSE])

# Color scheme
col_fun <- colorRamp2(c(min(pval), max(pval)), c("#fed43d", "#ff5959"))

# Draw column
pdf(paste0(plot_dir, "pvalue_column.pdf"), width = 1.5)
Heatmap(
  pval, column_names_gp = gpar(fontsize = 9),
  cluster_rows = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  row_order = row_order,
  heatmap_legend_param = list(title = "adjusted pvalue"),
  col = col_fun)
dev.off()
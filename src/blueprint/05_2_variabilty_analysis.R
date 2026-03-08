#!/usr/bin/env Rscript
### Calculating heatmap of most variable TFs in blueprint data 
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

# Scripts
source("/icbb/projects/nitschre/methylTFR/scripts/other/variability_analysis.R")
source("/icbb/projects/nitschre/methylTFR/scripts/other/helpers.R")

# Cell type color
cell_type_colors <- c(
  "Bcell" = "#980043",
  "plasma" = "#CD2990",
  "DC" = "#EED5B7",
  "other" = "#8B8682",
  "gran" = "#ff7f50",
  "Mf" = "#864a38",
  "mono" = "#CD7054",
  "Tcell" = "#40E0D0",
  "thymocyte" = "#74c476"
)


# Loading deviations scores and remove problematic sample
deviations_raw <- readRDS("/scratch/icbb/igunduz/methylTFR_manuscript/blueprint/mTFR_devs_121125/JASPAR2020_distal_deviations.RDS")
deviations_raw <- deviations_raw[,!colnames(deviations_raw) %in% "Bcell_naive_VB_NBC_NC11_83.bed"]
deviations <- deviations(deviations_raw)

# Change sample names to cellTypeGroup
annot <- read.csv("/icbb/projects/nitschre/methylTFR/sample_annotation/samples.tsv", sep="\t")
annot <- subset(annot, bedFile %in% colnames(deviations))

deviations <- deviations[,match(annot$bedFile, colnames(deviations))]
colnames(deviations) <- annot$cellTypeGroup

deviations_raw <- deviations_raw[,match(annot$bedFile, colnames(deviations_raw))] 
colnames(deviations_raw) <- annot$cellTypeGroup 

# Compute zscores columnwise
zscores <- computeColZScore(deviations)

# Perform variability analysis
if(file.exists(paste0(results_dir, "variability_bp.RDS"))) {  
  var <- readRDS(paste0(results_dir, "variability_bp.RDS"))
} else {
  var <-computeZScoreVariability(zscores)
  saveRDS(var, file=paste0(results_dir, "variability_bp.RDS"))
}


# Filter for  pval < .05 and keep top 50 with smallest pvalue
top50 <- var %>%
  filter(p_value_adj < 0.05) %>%
  arrange(desc(variability)) %>%
  head(50)

# Keep only samples that have at least 5 samples
deviations_raw_filtered <- deviations_raw[,colnames(deviations_raw) %in% c("plasma", "gran", "Bcell", "DC","Mf","mono","other","Tcell","thymocyte")]

# Get Z-scores
zscores <- deviationZScores(deviations_raw_filtered)

# Filter for TFs that are in the top 50 variables
zscores_filtered <- zscores[rownames(zscores) %in% c(rownames(top50)),]

# Heatmap Annotation at the top of the plot
ha <- HeatmapAnnotation(
  celltypes=colnames(deviations_raw_filtered),
  col = list(celltypes = cell_type_colors)
)

# Color scheme
col <- muRtools::colpal.cont(100, "cptcity.arendal_temperature")

# Column split factor
group <- colnames(deviations_raw_filtered)
column_split_factor <- factor(group, levels=unique(group))

# Clipped matrix
zscores_clipped <- pmax(pmin(zscores_filtered,2),-2)

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

### Calculating differenitals for TCells vs Bcells and Heatmap of significant TFs
set.seed(13)

library(methylTFR)
library(ComplexHeatmap)
library(dplyr)

plot_dir <- "/icbb/projects/nitschre/methylTFR/figures/blueprint"

# Loading deviations scores
deviations_raw <- readRDS("/icbb/projects/nitschre/methylTFR/r_objects/blueprintjaspar2020_distal_deviations_subset.RDS")
deviations <- deviations(deviations_raw)

### Bcell vs Tcells 
colnames(deviations) <- case_when(
  grepl("Bcell", colnames(deviations), ignore.case = TRUE) ~ "Bcell",
  grepl("TC", colnames(deviations), ignore.case = TRUE) ~ "Tcell",
)

# Set groups to compare
groups <- colnames(deviations)

# differntial analysis
diff <- differential_deviation_test(deviations, groups=groups,
                                        alternative="two.sided", parametric = TRUE)
saveRDS(diff, file="//icbb/projects/nitschre/methylTFR/r_objects/diff_blueprintSubset.RDS")
                                       

# Filter for pval < .05
diff_filtered <- diff[diff$p_value_adjusted < 0.05,]

#### Heatmap
# Get Z-scores
zscores <- deviationZScores(deviations_raw)

# Filter for diff TFs
zscores_filtered <- zscores[rownames(zscores) %in% c(rownames(diff_filtered)),]

# Heatmap
path <- file.path(plot_dir, "heatmap_diffmotifs_blueprintSubset.pdf")

ha <- HeatmapAnnotation(
  celltypes=colnames(deviations),
  col = list(celltypes = c("Bcell" = "blue", "Tcell"="red"))
)

pdf(path)
Heatmap(
  zscores_filtered, column_names_gp = gpar(fontsize = 9),
  top_annotation = ha,
  column_title = "Z-Scores of differential TFs in Tcells vs Bcells",
  show_row_names = FALSE,
  show_column_names = FALSE)
dev.off()

## top 50
top50 <- head(diff_filtered[order(diff_filtered$p_value_adjusted), ], 50)
zscores_filtered_50 <- zscores[rownames(zscores) %in% c(rownames(top50)),]

# Heatmap
path <- file.path(plot_dir, "heatmap_50diffmotifs_blueprintSubset.pdf")

ha <- HeatmapAnnotation(
  celltypes=colnames(deviations),
  col = list(celltypes = c("Bcell" = "blue", "Tcell"="red"))
)

pdf(path)
Heatmap(
  zscores_filtered_50, column_names_gp = gpar(fontsize = 9),
  top_annotation = ha,
  column_title = "Z-Scores of differential TFs in Tcells vs Bcells",
  show_row_names = TRUE,
  show_column_names = FALSE)
dev.off()
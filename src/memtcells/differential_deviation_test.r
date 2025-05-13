### Calculating differenitals for TN vs EM and CM + Heatmap of significant TFs
set.seed(13)

library(methylTFR)
library(ComplexHeatmap)

plot_dir <- "/icbb/projects/nitschre/methylTFR/figures"

# Loading deviations scores
deviations_raw <- readRDS("//icbb/projects/nitschre/mTFR_testmemTcells_260624/jaspar2020_distal_deviations.RDS")
deviations <- deviations(deviations_raw)

### EM vs TN
# Subsetting deviation matrix to only naive and EM cells
deviations_em <- deviations[, grepl("TN|EM", colnames(deviations))]

# Set groups to compare
groups <- c("EM", "TN", "EM", "TN")

# differntial analysis
diff_em <- differential_deviation_test(deviations_em, groups=groups,
                                        alternative="two.sided", parametric = TRUE)
saveRDS(diff_em, file="//icbb/projects/nitschre/methylTFR/r_objects/diff_em.RDS")
                                       

# Filter for pval < .05
diff_em_filtered <- diff_em[diff_em$p_value_adjusted < 0.05,]

### CM vs TN
# Subsetting deviations matrix to only naive and CM cells
deviations_cm <- deviations[, grepl("TN|CM", colnames(deviations))]

# groups
groups <- c("CM", "TN", "CM", "TN")

# differential analysis
diff_cm <- differential_deviation_test(deviations_cm, groups=groups,
                                        alternative="two.sided", parametric = TRUE)
saveRDS(diff_cm, file="//icbb/projects/nitschre/methylTFR/r_objects/diff_cm.RDS")


# Filter for pval < .05
diff_cm_filtered <- diff_cm[diff_cm$p_value_adjusted < 0.05,]

#### Heatmap
# Get Z-scores
zscores <- deviationZScores(deviations_raw)

# Filter for diff TFs
# TFs that are in diff_em or in diff_cm
zscores_filtered <- zscores[rownames(zscores) %in% c(rownames(diff_em_filtered), rownames(diff_cm_filtered)),]

# Heatmap
colnames(zscores_filtered) <- c("Hf03_CM", "Hf03_EM", "Hf03_TN", "Hf04_CM", "Hf04_EM",
                                "Hf04_TN", "Hf04_TEMRA")
path <- file.path(plot_dir, "heatmap_diffmotifs.pdf")

ha <- HeatmapAnnotation(
  celltypes=c("CM","EM", "TN", "CM", "EM","TN","TEMRA"),
  col = list(celltypes = c("TN" = "#ff0000", "TEMRA" = "#ff6600", "CM"="#008cff", "EM"="#0037ff"))
)

pdf(path)
Heatmap(
  zscores_filtered, column_names_gp = gpar(fontsize = 9),
  top_annotation = ha,
  column_title = "Z-Scores of differential TFs",
  show_row_names = FALSE,
  column_names_side = "top")
dev.off()

## TFs TOP 50 
top50em <- head(diff_em_filtered[order(diff_em_filtered$p_value_adjusted), ], 50)
top50cm <- head(diff_cm_filtered[order(diff_cm_filtered$p_value_adjusted), ], 50)

zscores_top50 <- zscores[rownames(zscores) %in% c(rownames(top50cm), rownames(top50em)),]

# Heatmap
colnames(zscores_top50) <- c("Hf03_CM", "Hf03_EM", "Hf03_TN", "Hf04_CM", "Hf04_EM",
                                "Hf04_TN", "Hf04_TEMRA")

# Remove TEMRA
zscores_top50 <- zscores_top50[, -which(colnames(zscores_top50) == "Hf04_TEMRA")]
path <- file.path(plot_dir, "heatmap_diffmotifs_top50.pdf")

ha <- HeatmapAnnotation(
  celltypes=c("CM","EM", "TN", "CM", "EM", "TN"),
  col = list(celltypes = c("TN" = "#ff0000", "CM"="#008cff", "EM"="#0037ff"))
)

pdf(path)
Heatmap(
  zscores_top50, column_names_gp = gpar(fontsize = 9),
  top_annotation = ha,
  column_title = "Z-Scores of differential TFs",
  show_row_names = TRUE,
  column_names_side = "top")
dev.off()
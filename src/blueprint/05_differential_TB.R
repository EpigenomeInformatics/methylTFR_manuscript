#!/usr/bin/env Rscript
#####################################################################
# differentials_heatmap.R
# Created by RN on 13-11-2025
# Script to plot heatmap of differential TFs from mTFR
#####################################################################
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(methylTFR)
  library(ComplexHeatmap)
  library(dplyr)
  library(circlize)
  library(factoextra)
  library(ggfortify)
  library(ggplot2)
  library(RColorBrewer)
})

set.seed(12)

plot_dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR_manuscript/figures/blueprint/"
results_dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR_manuscript/tables/"
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

deviations_raw <- readRDS("/scratch/icbb/igunduz/methylTFR_manuscript/blueprint/mTFR_devs_121125/JASPAR2020_distal_deviations.RDS")
deviations_raw <- deviations_raw[, !colnames(deviations_raw) %in% "Bcell_naive_VB_NBC_NC11_83.bed"]
deviations <- deviations(deviations_raw)

sannot <- read.csv("/scratch/icbb/igunduz/methylTFR_manuscript/blueprint/RnBeads_291025/reports/data_import_data/annotation.csv", stringsAsFactors = FALSE)
sannot$bedFile <- as.character(sannot$bedFile)
sample_names <- colnames(deviations)
sample_names_basenames <- basename(sample_names)
sannot_bedfiles_basenames <- basename(sannot$bedFile)
match_indices <- match(sample_names_basenames, sannot_bedfiles_basenames)

cell_types <- sannot$cellTypeGroup[match_indices]
cell_subtypes <- sannot$cellTypeShort[match_indices]

subset_idx <- grepl("Bcell|Tcell", cell_types, ignore.case = TRUE) &
  !is.na(cell_types) &
  !grepl("Thym", cell_subtypes, ignore.case = TRUE)

deviations <- deviations[, subset_idx]
sample_names_filtered <- sample_names[subset_idx]
cell_types_filtered <- cell_types[subset_idx]
cell_subtypes_filtered <- cell_subtypes[subset_idx]

tdf <- as.data.frame(t(deviations))
pca_result <- prcomp(tdf, center = FALSE, scale. = FALSE)

pca_annot <- data.frame(
  SampleID = sample_names_filtered,
  Group = cell_types_filtered,
  Subtype = cell_subtypes_filtered,
  stringsAsFactors = FALSE
)
tdf$Subtype <- pca_annot$Subtype[match(rownames(tdf), pca_annot$SampleID)]

b_subtypes <- sort(unique(
  pca_annot$Subtype[grepl("Bcell", pca_annot$Group, ignore.case = TRUE)]
))
t_subtypes <- sort(unique(
  pca_annot$Subtype[grepl("Tcell", pca_annot$Group, ignore.case = TRUE)]
))

n_b_subtypes <- length(b_subtypes)
n_t_subtypes <- length(t_subtypes)

b_palette_vals <- colorRampPalette(
  brewer.pal(max(3, min(n_b_subtypes, 9)), "YlOrBr")
)(n_b_subtypes)
t_palette_vals <- colorRampPalette(
  brewer.pal(max(3, min(n_t_subtypes, 9)), "Blues")
)(n_t_subtypes)

names(b_palette_vals) <- b_subtypes
names(t_palette_vals) <- t_subtypes
combined_palette <- c(b_palette_vals, t_palette_vals)

fig_path <- paste0(plot_dir, "PCA_deviations_B_and_T_subsets_noThym.pdf")
p <- autoplot(pca_result,
  data = tdf,
  colour = "Subtype",
  main = "PCA - B-Cell and T-Cell Subsets (No Thymocytes)",
  size = 5
) +
  theme_classic() +
  scale_color_manual(values = combined_palette, name = "Cell Subtype") +
  theme(legend.position = "bottom")

pdf(fig_path, width = 10, height = 10)
print(p)
dev.off()

new_colnames <- case_when(
  grepl("Bcell", cell_types_filtered, ignore.case = TRUE) ~ "Bcell",
  grepl("Tcell", cell_types_filtered, ignore.case = TRUE) ~ "Tcell"
)
colnames(deviations) <- new_colnames
groups <- colnames(deviations)

if (file.exists(paste0(results_dir, "diff_BvsT_noThym.csv"))) {
  diff <- read.csv(paste0(results_dir, "diff_BvsT_noThym.csv"), row.names = 1)
} else {
  diff <- differential_deviation_test(deviations, groups = groups, alternative = "two.sided", parametric = TRUE)
  write.csv(diff, paste0(results_dir, "diff_BvsT_noThym.csv"), row.names = TRUE)
}

top50 <- diff %>%
  filter(p_value_adjusted < 0.05) %>%
  arrange(p_value_adjusted) %>%
  head(50)

zscores_filtered <- deviations[rownames(deviations) %in% c(rownames(top50)), ]
zscores_filtered <- methylTFR:::computeRowZScore(zscores_filtered)

ha <- HeatmapAnnotation(
  celltypes = groups,
  col = list(celltypes = c("Bcell" = "#2B4B9B", "Tcell" = "#D76016"))
)

ht <- Heatmap(
  zscores_filtered,
  column_names_gp = gpar(fontsize = 9),
  top_annotation = ha,
  column_title = "Z-Scores of differential TFs in Tcells vs Bcells (No Thymocytes)",
  heatmap_legend_param = list(title = "methylTFR Z-scores"),
  show_row_names = TRUE,
  show_column_names = FALSE
)

pdf(paste0(plot_dir, "heatmap_BvsTcells_noThym.pdf"))
ht_drawn <- draw(ht)
dev.off()

mean_bcells <- rowMeans(zscores_filtered[, groups == "Bcell"])
mean_tcells <- rowMeans(zscores_filtered[, groups == "Tcell"])
zscores_diff <- as.matrix(mean_bcells - mean_tcells)

row_order <- row_order(ht_drawn)

col_fun <- colorRamp2(c(min(zscores_diff), max(zscores_diff)), c("#DCEDC8", "#1A237E"))

pdf(paste0(plot_dir, "meanZscoreDiff_BvsTcells_noThym.pdf"), width = 2)
Heatmap(
  zscores_diff,
  column_names_gp = gpar(fontsize = 9),
  cluster_rows = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  row_order = row_order,
  heatmap_legend_param = list(title = "Mean z-score difference"),
  col = col_fun
)
dev.off()

pval <- as.matrix(diff[rownames(diff) %in% rownames(zscores_filtered), "p_value_adjusted", drop = FALSE])

col_fun <- colorRamp2(c(min(pval), max(pval)), c("#fed43d", "#ff5959"))

pdf(paste0(plot_dir, "pvalue_column_noThym.pdf"), width = 1.5)
Heatmap(
  pval,
  column_names_gp = gpar(fontsize = 9),
  cluster_rows = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  row_order = row_order,
  heatmap_legend_param = list(title = "adjusted pvalue"),
  col = col_fun
)
dev.off()
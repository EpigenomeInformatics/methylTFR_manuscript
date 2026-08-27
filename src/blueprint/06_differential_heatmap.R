#!/usr/bin/env Rscript

#####################################################################
# 06_differential_heatmap.R
# created on 23-08-2026 by Irem B Gunduz
# Differential deviation test across all cell type groups, heatmap of
# the Z-scores of the top motifs (Figure C), and a barplot of the cell
# type composition of the cohort (replaces the pie chart of Figure A)
#####################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(grid)
  library(circlize)
  library(ComplexHeatmap)
  library(logger)
  library(methylTFR)
})
set.seed(42)

#####################################################################
# Settings
#####################################################################

# Motif sets whose deviations are tested
motifSets <- c("jaspar2020", "jaspar2020_distal")

# Differential test settings
alternative <- "two.sided"
parametric <- TRUE
padj.method <- "BH"
padj.cutoff <- 0.05

# Number of motifs shown in the heatmap
top.n <- 50

# The stored Z-scores are row-wise across every sample of the run, which is
# what the figure shows. Set this to TRUE to recompute them within the
# samples actually plotted.
recompute.z <- FALSE

# Directories
analysis.dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/blueprint/"
rnb.tag <- "RnBeads_230826"
dev.tag <- "mTFR_devs_230826"
sannot.file <- file.path(analysis.dir, rnb.tag, "reports", "data_import_data", "annotation.csv")

# Figures live in the repository, next to the tables they belong with.
# Only the footprints stay under analysis.dir, they are too many for git.
github.dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR_manuscript"
plot.dir <- file.path(github.dir, "figures", "blueprint", "diff_heatmaps_230826")
if (!dir.exists(plot.dir)) dir.create(plot.dir, recursive = TRUE)

table.dir <- file.path(github.dir, "tables")
if (!dir.exists(table.dir)) dir.create(table.dir, recursive = TRUE)

# Remap to cleaner group names, same mapping as 03
group_remap <- c(
  "Bcell" = "B-cells",
  "DC" = "DC",
  "eryt" = "Erythrocytes",
  "gran" = "Granulocytes",
  "megK" = "Megakaryocytes",
  "Mf" = "Macrophages",
  "mono" = "Monocytes",
  "NK" = "NK-cells",
  "osteoclast" = "Osteoclast",
  "other" = "Other",
  "plasma" = "Plasma",
  "progenitor" = "Progenitors",
  "Tcell" = "T-cells",
  "thymocyte" = "Thymocyte"
)

cell_type_colors <- c(
  "B-cells" = "#C2377C",
  "DC" = "#8C6D3F",
  "Erythrocytes" = "#7E4B2A",
  "Granulocytes" = "#E8A33D",
  "Macrophages" = "#B5A38A",
  "Megakaryocytes" = "#6A3D9A",
  "Monocytes" = "#C2703D",
  "NK-cells" = "#2CA02C",
  "Osteoclast" = "#8C8C8C",
  "Other" = "#BCBD22",
  "Plasma" = "#7B1E3D",
  "Progenitors" = "#17BECF",
  "T-cells" = "#4FC3D9",
  "Thymocyte" = "#5B9BD5"
)

# Blue to yellow to red, capped at +/- 2 as in the figure
zscore.limit <- 2
heatmap_colors <- colorRamp2(
  c(-zscore.limit, -1, 0, 1, zscore.limit),
  c("#2B4B9B", "#4EC3E0", "#FFF200", "#F5A623", "#B21212")
)

#####################################################################
# Helper functions
#####################################################################

# Row-wise Z-score, used only when recompute.z is TRUE
row_zscore <- function(mat) {
  centred <- mat - rowMeans(mat, na.rm = TRUE)
  sds <- apply(mat, 1, sd, na.rm = TRUE)
  out <- centred / sds
  out[!is.finite(out)] <- 0
  out
}

# Heatmap of the Z-scores, samples in columns, motifs in rows.
# Columns are split by cell type and clustered within each block, so the
# dendrogram describes the samples of one cell type rather than mixing
# them. Rows are clustered across the whole matrix.
plot_zscore_heatmap <- function(mat, cell_types, title, file) {
  if (nrow(mat) < 2) {
    log_warn(basename(file), ": only ", nrow(mat), " motif(s), skipping the heatmap")
    return(invisible(NULL))
  }
  present <- intersect(names(cell_type_colors), unique(cell_types))
  ha <- HeatmapAnnotation(
    `Cell type` = cell_types,
    col = list(`Cell type` = cell_type_colors[present]),
    annotation_name_gp = gpar(fontsize = 8),
    show_legend = TRUE
  )

  # The palette order fixes the block order, so it matches every other
  # Blueprint panel rather than following the clustering
  column_split_factor <- factor(cell_types, levels = present)

  ht <- Heatmap(
    mat,
    name = "methylTFR\nZ-scores",
    col = heatmap_colors,
    top_annotation = ha,
    column_title = title,
    column_title_gp = gpar(fontsize = 10),
    column_split = column_split_factor,
    cluster_column_slices = FALSE,
    column_gap = unit(0.8, "mm"),
    show_row_names = TRUE,
    show_column_names = FALSE,
    row_names_side = "right",
    row_names_gp = gpar(fontsize = 7),
    row_dend_side = "left",
    clustering_distance_rows = "euclidean",
    clustering_distance_columns = "euclidean",
    clustering_method_rows = "complete",
    clustering_method_columns = "complete",
    heatmap_legend_param = list(
      direction = "horizontal",
      at = seq(-zscore.limit, zscore.limit, by = 1),
      legend_width = unit(4, "cm"),
      title_position = "lefttop"
    )
  )
  pdf(file, width = 9, height = 11)
  draw(ht,
    heatmap_legend_side = "bottom",
    annotation_legend_side = "bottom",
    merge_legend = TRUE
  )
  dev.off()
  log_info("Wrote ", file)
}

# Cell type composition as a horizontal barplot. A pie chart of fourteen
# groups cannot be read, so the groups are sorted by size and each bar
# carries its own count and share.
plot_composition_barplot <- function(cell_types, title, file) {
  composition <- as.data.frame(table(cell_types), stringsAsFactors = FALSE)
  colnames(composition) <- c("cell_type", "n")
  composition$percent <- 100 * composition$n / sum(composition$n)
  composition <- composition[order(composition$n), ]
  composition$cell_type <- factor(composition$cell_type, levels = composition$cell_type)
  composition$label <- paste0(composition$n, " (", sprintf("%.1f", composition$percent), "%)")

  p <- ggplot(composition, aes(x = cell_type, y = n, fill = cell_type)) +
    geom_col(width = 0.75) +
    geom_text(aes(label = label), hjust = -0.15, size = 3.2) +
    coord_flip(clip = "off") +
    scale_fill_manual(values = cell_type_colors, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    xlab(NULL) +
    ylab("Number of samples") +
    ggtitle(title) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 11),
      axis.text.y = element_text(size = 9),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    )

  ggsave(file, plot = p, width = 7, height = 5)
  log_info("Wrote ", file)
}

#####################################################################
# Sample annotation
#####################################################################

sannot <- read.csv(sannot.file, stringsAsFactors = FALSE)
sannot$bedFile <- as.character(sannot$bedFile)

#####################################################################
# Differential analysis, heatmap and composition, per motif set
#####################################################################

for (motifSet in motifSets) {
  dev.file <- file.path(analysis.dir, dev.tag, paste0(motifSet, "_deviations.RDS"))
  if (!file.exists(dev.file)) {
    log_warn(motifSet, ": ", dev.file, " not found, skipping this motif set")
    next
  }

  log_info("Loading deviations from ", dev.file)
  dev_obj <- readRDS(dev.file)

  # Healthy samples carrying a cell type we can map. Groups come from the
  # annotation rather than from the sample names.
  matched <- sannot[match(colnames(dev_obj), sannot$bedFile), ]
  keep <- which(
    !is.na(matched$DISEASE) &
      matched$DISEASE == "None" &
      matched$cellTypeGroup %in% names(group_remap)
  )
  if (length(keep) == 0) {
    log_warn(motifSet, ": no samples left after filtering, skipping")
    next
  }

  # Subsetting the object keeps the deviations and the Z-scores aligned
  dev_obj <- dev_obj[, keep]
  cell_types <- unname(group_remap[matched$cellTypeGroup[keep]])
  log_info(motifSet, ": ", length(cell_types), " samples across ",
    length(unique(cell_types)), " cell type groups")
  if (length(unique(cell_types)) < 2) {
    log_warn(motifSet, ": only one cell type group present, skipping")
    next
  }

  # Composition of the cohort that goes into the heatmap
  plot_composition_barplot(
    cell_types,
    paste0("Cell type composition (n = ", length(cell_types), ")"),
    file.path(plot.dir, paste0("celltype_composition_", motifSet, ".pdf"))
  )

  # Differential test. With more than two groups this is an ANOVA across
  # all cell types, so the motifs that come out are the ones that separate
  # the lineages rather than any single pair.
  log_info(motifSet, ": running the differential deviation test")
  diff <- differential_deviation_test(
    deviations = deviations(dev_obj),
    groups = cell_types,
    alternative = alternative,
    parametric = parametric,
    padjMethod = padj.method
  )
  diff <- diff[order(diff$p_value_adjusted, -diff$mean_difference), ]

  saveRDS(diff, file.path(table.dir, paste0("diff_", motifSet, "_allcelltypes.RDS")))
  write.csv(diff,
    file.path(table.dir, paste0("diff_", motifSet, "_allcelltypes.csv")),
    row.names = FALSE
  )

  # Significant motifs, which() drops the motifs whose test returned NA
  sig <- diff[which(diff$p_value_adjusted < padj.cutoff), ]
  log_info(motifSet, ": ", nrow(sig), " motifs with adjusted p < ", padj.cutoff)
  if (nrow(sig) == 0) {
    log_warn(motifSet, ": nothing significant, no heatmap")
    next
  }

  # Z-scores of the same samples, in the same order as cell_types
  zscores <- if (recompute.z) row_zscore(deviations(dev_obj)) else deviationZScores(dev_obj)

  top <- head(sig, top.n)
  mat <- zscores[top$motifs, , drop = FALSE]

  # Clustering cannot handle non finite rows
  finite_rows <- is.finite(rowSums(mat))
  if (any(!finite_rows)) {
    log_warn(motifSet, ": dropping ", sum(!finite_rows), " motifs with non finite Z-scores")
    mat <- mat[finite_rows, , drop = FALSE]
  }

  plot_zscore_heatmap(
    mat, cell_types,
    paste0("Top ", nrow(mat), " differential TFs (", motifSet, ")"),
    file.path(plot.dir, paste0("heatmap_top", top.n, "_diffmotifs_", motifSet, ".pdf"))
  )
}

log_success("Finished the differential analysis, heatmaps and composition barplots")

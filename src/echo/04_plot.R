#!/usr/bin/env Rscript

#####################################################################
# 04_plot.R
# created on 25-08-2026 by Irem B Gunduz
# The ECHO integration panels that 01 to 04 do not already produce:
#   - the paired methylTFR / chromVAR heatmap with the row-wise
#     correlation and the MOFA factor annotation
#   - the correlation boxplot split by methyl-SELEX call
#   - the per TF chromVAR against methylTFR scatters
#
# 01 already covers the factor scatter, the modality contribution, the
# variance explained, the top TF dotplot and the factor strip plot; 03
# covers the footprints. Nothing here retrains MOFA, it reads the saved
# factor table.
#####################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(grid)
  library(circlize)
  library(RColorBrewer)
  library(ComplexHeatmap)
  library(patchwork)
  library(scales)
  library(logger)
})
set.seed(12)

#####################################################################
# Settings
#####################################################################

# Pseudobulk TF activity matrices, TFs in rows and samples in columns.
# Same two objects 01 reads.
chromvar.file <- "/icbb/projects/nitschre/methylTFR/scripts/other/cvar_zscores_adj_psuedobulk.R"
mtfr.file <- "/icbb/projects/nitschre/methylTFR/scripts/other/mtfr_zscores_adj_psuedobulk.R"

# ECHO MOFA weights, long: feature, factor, value, view. The view is
# appended to the feature name.
factor.table <- "/icbb/projects/nitschre/methylTFR/tables/echo_factor_table.csv"

# Top features per factor and view that make up the heatmap rows
tfs.per.factor <- 6
factors.keep <- c("Factor1", "Factor2", "Factor3", "Factor4", "Factor7")

# TFs shown as individual scatters
scatter.tfs <- c("BATF", "FOXL1", "SPIB", "POU2F3")

# Modality contribution and factor R2, written by 01. Panel A is redrawn
# from these rather than recomputed, because the fraction needs the full
# weight matrix and the factor table holds only the top features.
contrib.file <- "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR_manuscript/tables/echo_mofa_modality_contribution.csv"
r2.file <- "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR_manuscript/tables/echo_mofa_factor_celltype_R2.csv"

# View colours, matching 01
view_colors <- c("chromVAR" = "#8FA8D4", "methylTFR" = "#E4675C")

selex.file <- "/icbb/projects/igunduz/exposure_atlas_manuscript/sample_annots/Selex_data.csv"
selex.groups <- c("MethylMinus", "MethylPlus")

# Directories
# Figures live in the repository, next to the tables they belong with
github.dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR_manuscript"
plot.dir <- file.path(github.dir, "figures", "echo", "integration_230826")
if (!dir.exists(plot.dir)) dir.create(plot.dir, recursive = TRUE)

table.dir <- file.path(github.dir, "tables")
if (!dir.exists(table.dir)) dir.create(table.dir, recursive = TRUE)

# Palettes
cell_type_colors <- c(
  "B-cell" = "#AE017E",
  "Monocyte" = "#CC4C02",
  "NK-cell" = "#A65628",
  "Th-Mem" = "#41B6C4",
  "Tc-Mem" = "#4292C6",
  "Tc-Naive" = "#888FB5",
  "Th-Naive" = "#C7E9B4"
)

mtfr_colors <- colorRamp2(
  c(-4, -2, 0, 2, 4),
  c("#2A2F9E", "#7570B3", "#FFFFFF", "#F16C43", "#E03426")
)
# Diverging, teal for a depleted motif and gold for an enriched one, on a
# symmetric -5 to 5 scale. Sequential would have read a strong depletion
# and a value near zero as the same pale colour.
chromvar_limit <- 5
chromvar_colors <- colorRamp2(
  seq(-chromvar_limit, chromvar_limit, length.out = 11),
  rev(brewer.pal(11, "BrBG"))
)
cor_colors <- colorRamp2(
  c(-1, -0.5, 0, 0.5, 1),
  c("#8E0152", "#DE77AE", "#FFFFFF", "#7FBC41", "#276419")
)
selex_colors <- c("MethylMinus" = "#8B1A1A", "MethylPlus" = "#1B6B3A")

#####################################################################
# Load the two modalities
#####################################################################

for (f in c(chromvar.file, mtfr.file, factor.table)) {
  if (!file.exists(f)) stop("Missing input: ", f)
}

# These carry a .R extension but hold serialised matrices, as in 01
chromVar <- as.matrix(readRDS(gzfile(chromvar.file)))
mtfr <- as.matrix(readRDS(gzfile(mtfr.file)))
log_info("chromVAR: ", nrow(chromVar), " x ", ncol(chromVar),
  ", methylTFR: ", nrow(mtfr), " x ", ncol(mtfr))

# Same samples, same order, on both sides
common_samples <- intersect(colnames(mtfr), colnames(chromVar))
if (length(common_samples) == 0) stop("The two matrices share no samples")
if (length(common_samples) < ncol(mtfr)) {
  log_warn("Only ", length(common_samples), " of ", ncol(mtfr), " samples are shared")
}
mtfr <- mtfr[, common_samples, drop = FALSE]
chromVar <- chromVar[, common_samples, drop = FALSE]

# Cell type is the part of the sample name before the first underscore,
# the convention 01 uses
cell_types <- sub("_.*", "", common_samples)
unmapped <- setdiff(unique(cell_types), names(cell_type_colors))
if (length(unmapped) > 0) {
  stop("No colour for cell type(s): ", paste(unmapped, collapse = ", "))
}
log_info(length(common_samples), " samples across ",
  length(unique(cell_types)), " cell types")

#####################################################################
# MOFA features
#####################################################################

ft <- read.csv(factor.table, header = TRUE, stringsAsFactors = FALSE)
required <- c("feature", "factor", "value", "view")
missing_cols <- setdiff(required, colnames(ft))
if (length(missing_cols) > 0) {
  stop("Missing column(s) in the factor table: ", paste(missing_cols, collapse = ", "))
}

# The view is appended to the feature name, stripped with the row's own
# view so motif names carrying punctuation survive
ft$feature <- mapply(
  function(feat, vw) sub(paste0("_", vw, "$"), "", feat),
  ft$feature, ft$view,
  USE.NAMES = FALSE
)

if (!is.null(factors.keep)) ft <- ft[ft$factor %in% factors.keep, ]
ft$value_abs <- abs(ft$value)
ft <- ft[order(ft$factor, ft$view, -ft$value_abs), ]
rank_in_group <- stats::ave(
  seq_len(nrow(ft)), paste(ft$factor, ft$view),
  FUN = seq_along
)
ft_top <- ft[rank_in_group <= tfs.per.factor, ]

# One factor per feature, the one where it carries the largest weight
feature_factor <- ft_top %>%
  group_by(feature) %>%
  slice_max(value_abs, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(feature, factor) %>%
  as.data.frame()

mofa_tfs <- feature_factor$feature
common_tfs <- Reduce(intersect, list(mofa_tfs, rownames(mtfr), rownames(chromVar)))
absent <- setdiff(mofa_tfs, common_tfs)
if (length(absent) > 0) {
  log_warn(length(absent), " MOFA features are not in both matrices: ",
    paste(head(absent, 10), collapse = ", "))
}
if (length(common_tfs) < 2) stop("Fewer than two MOFA features are in both matrices")
log_info(length(common_tfs), " features on the heatmap")

#####################################################################
# Panel A, modality contribution with the cell type R2 over it
#####################################################################

p_modality <- NULL
if (!file.exists(contrib.file) || !file.exists(r2.file)) {
  log_warn("Modality contribution tables not found, run 01 first. Panel A is skipped")
} else {
  agg <- read.csv(contrib.file, stringsAsFactors = FALSE)
  r2_tab <- read.csv(r2.file, stringsAsFactors = FALSE)

  # Factors ordered by the variance cell type explains, so the R2 line
  # falls from left to right
  factor_order <- r2_tab$Factor[order(-r2_tab$R2)]
  factor_order <- factor_order[factor_order %in% unique(agg$factor)]
  if (length(factor_order) == 0) {
    log_warn("The two tables share no factor, Panel A is skipped")
  } else {
    agg <- agg[agg$factor %in% factor_order, ]
    agg$factor <- factor(agg$factor, levels = factor_order)

    unmapped_views <- setdiff(unique(agg$view), names(view_colors))
    if (length(unmapped_views) > 0) {
      stop("No colour for view(s): ", paste(unmapped_views, collapse = ", "))
    }
    agg$view <- factor(agg$view, levels = names(view_colors))

    r2_line <- r2_tab[r2_tab$Factor %in% factor_order, ]
    r2_line$Factor <- factor(r2_line$Factor, levels = factor_order)

    p_modality <- ggplot(agg, aes(x = factor, y = frac)) +
      geom_bar(aes(fill = view), stat = "identity", width = 0.7) +
      geom_text(
        aes(group = view, label = ifelse(frac > 0.03, paste0(round(frac * 100, 1), "%"), "")),
        position = position_stack(vjust = 0.5), size = 3
      ) +
      geom_line(data = r2_line, aes(x = Factor, y = R2, group = 1, colour = "Variance Explained (R2)")) +
      geom_point(data = r2_line, aes(x = Factor, y = R2, colour = "Variance Explained (R2)")) +
      scale_colour_manual(values = c("Variance Explained (R2)" = "purple"), name = NULL) +
      scale_y_continuous(
        labels = percent_format(accuracy = 1),
        sec.axis = sec_axis(~., name = expression(Variance ~ Explained ~ (R^2)))
      ) +
      scale_fill_manual(values = view_colors) +
      labs(x = "Factor", y = "Relative contribution (|weights| fraction)", fill = "View") +
      theme_classic(base_size = 12) +
      theme(
        legend.position = "bottom",
        axis.title.y.right = element_text(colour = "purple")
      )

    file <- file.path(plot.dir, "modality_contribution_fraction_R2.pdf")
    ggsave(file, p_modality, width = 8, height = 4.5)
    log_info("Wrote ", file)
  }
}

#####################################################################
# Paired heatmap, methylTFR left, chromVAR right
#####################################################################

mtfr_heat <- mtfr[common_tfs, , drop = FALSE]
chromvar_heat <- chromVar[common_tfs, , drop = FALSE]

row_cor <- vapply(common_tfs, function(tf) {
  suppressWarnings(cor(
    as.numeric(mtfr_heat[tf, ]), as.numeric(chromvar_heat[tf, ]),
    method = "pearson", use = "complete.obs"
  ))
}, numeric(1))
row_cor[!is.finite(row_cor)] <- 0

split_factor <- factor(cell_types, levels = intersect(
  names(cell_type_colors), unique(cell_types)
))
cell_annotation <- function() {
  HeatmapAnnotation(
    celltypes = cell_types,
    col = list(celltypes = cell_type_colors[levels(split_factor)]),
    annotation_name_gp = gpar(fontsize = 8),
    show_legend = FALSE
  )
}

ht_mtfr <- Heatmap(mtfr_heat,
  name = "mTFR\nZ-scores",
  col = mtfr_colors,
  top_annotation = cell_annotation(),
  column_split = split_factor,
  cluster_columns = TRUE,
  cluster_column_slices = FALSE,
  cluster_rows = TRUE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  row_dend_side = "left",
  column_title_gp = gpar(fontsize = 8)
)

ht_chromvar <- Heatmap(chromvar_heat,
  name = "chromVAR\nZ-scores",
  col = chromvar_colors,
  heatmap_legend_param = list(at = seq(-chromvar_limit, chromvar_limit, by = 2)),
  top_annotation = cell_annotation(),
  column_split = split_factor,
  cluster_columns = TRUE,
  cluster_column_slices = FALSE,
  cluster_rows = FALSE, # row order comes from the methylTFR half
  show_row_names = FALSE,
  show_column_names = FALSE,
  column_title_gp = gpar(fontsize = 8)
)

factor_vec <- feature_factor$factor[match(common_tfs, feature_factor$feature)]
factor_levels <- sort(unique(factor_vec))
factor_colors <- setNames(
  colorRampPalette(brewer.pal(max(3, min(9, length(factor_levels))), "Set2"))(
    length(factor_levels)
  ),
  factor_levels
)

ha_right <- rowAnnotation(
  `Row Corr.` = row_cor,
  Factors = factor_vec,
  col = list(`Row Corr.` = cor_colors, Factors = factor_colors),
  TF = anno_text(common_tfs, gp = gpar(fontsize = 6)),
  annotation_name_gp = gpar(fontsize = 7),
  width = unit(28, "mm")
)

# Captured once as a grob, so the same drawing serves the standalone PDF
# and the assembled figure below
ht_grob <- grid.grabExpr(
  draw(ht_mtfr + ht_chromvar + ha_right,
    merge_legend = TRUE, heatmap_legend_side = "right"
  )
)

file <- file.path(plot.dir, "paired_heatmap_mtfr_chromvar.pdf")
pdf(file, width = 14, height = 10)
grid.draw(ht_grob)
dev.off()
log_info("Wrote ", file)

#####################################################################
# Per TF scatters, chromVAR against methylTFR
#####################################################################

plot_tf_scatter <- function(tf) {
  df <- data.frame(
    chromVar = as.numeric(chromVar[tf, ]),
    mtfr = as.numeric(mtfr[tf, ]),
    celltype = cell_types,
    stringsAsFactors = FALSE
  )
  cor_val <- suppressWarnings(cor(df$chromVar, df$mtfr,
    method = "pearson", use = "complete.obs"
  ))

  ggplot(df, aes(x = chromVar, y = mtfr)) +
    geom_point(aes(color = celltype), size = 1.6, alpha = 0.9) +
    geom_smooth(method = "lm", formula = y ~ x, colour = "red", se = FALSE, linewidth = 0.6) +
    scale_color_manual(values = cell_type_colors, name = "Cell type") +
    annotate("text",
      x = Inf, y = Inf, hjust = 1.05, vjust = 1.4, size = 4,
      label = paste0("Cor. Coef.: ", format(round(cor_val, 2), nsmall = 2))
    ) +
    labs(title = tf, x = "chromVar", y = "methylTFR") +
    theme_classic(base_size = 13) +
    theme(plot.title = element_text(hjust = 0, face = "plain"))
}

available_tfs <- intersect(scatter.tfs, intersect(rownames(mtfr), rownames(chromVar)))
missing_scatter <- setdiff(scatter.tfs, available_tfs)
if (length(missing_scatter) > 0) {
  log_warn("Not in both matrices, no scatter: ", paste(missing_scatter, collapse = ", "))
}
# One page for the whole set. The panels are not written individually as
# well, the combined figure is the one that goes into the manuscript.
if (length(available_tfs) > 0) {
  panels <- lapply(available_tfs, plot_tf_scatter)
  combined <- patchwork::wrap_plots(panels, ncol = 2, guides = "collect")
  file <- file.path(plot.dir, "scatter_chromvar_vs_mtfr_panels.pdf")
  ggsave(file, combined,
    width = 11, height = 4.5 * ceiling(length(available_tfs) / 2)
  )
  log_info("Wrote ", file)
}

#####################################################################
# Correlation by methyl-SELEX call, over every shared TF
#####################################################################

p_selex <- NULL
if (!file.exists(selex.file)) {
  log_warn("SELEX table not found, skipping that panel: ", selex.file)
} else {
  selex <- read.csv(selex.file, skip = 20, header = TRUE, sep = ";")
  call_col <- intersect(
    c("methyl.SELEX.call", "Call", "methylSELEXcall"), colnames(selex)
  )
  name_col <- intersect(c("TF.name", "TF", "TFname"), colnames(selex))

  if (length(call_col) == 0 || length(name_col) == 0) {
    log_warn(
      "No SELEX call or TF name column. Present: ",
      paste(colnames(selex), collapse = ", ")
    )
  } else {
    call_col <- call_col[1]
    name_col <- name_col[1]
    log_info("SELEX columns: name = ", name_col, ", call = ", call_col)

    # Correlation across every TF the two modalities share, not only the
    # MOFA features
    all_tfs <- intersect(rownames(mtfr), rownames(chromVar))
    cor_all <- vapply(all_tfs, function(tf) {
      suppressWarnings(cor(
        as.numeric(mtfr[tf, ]), as.numeric(chromVar[tf, ]),
        method = "pearson", use = "complete.obs"
      ))
    }, numeric(1))

    cor_df <- data.frame(
      TF.name = all_tfs, Correlation = cor_all,
      stringsAsFactors = FALSE
    )
    cor_df <- cor_df[is.finite(cor_df$Correlation), ]
    cor_df$call <- selex[[call_col]][match(cor_df$TF.name, selex[[name_col]])]
    cor_df <- cor_df[cor_df$call %in% selex.groups, ]

    if (nrow(cor_df) == 0) {
      log_warn("No TF carries one of ", paste(selex.groups, collapse = " or "))
    } else {
      cor_df$call <- factor(cor_df$call, levels = selex.groups)
      counts <- as.data.frame(table(cor_df$call))
      colnames(counts) <- c("call", "count")
      log_info("SELEX groups: ",
        paste(counts$call, counts$count, sep = " = ", collapse = ", "))

      p_selex <- ggplot(cor_df, aes(x = call, y = Correlation, fill = call)) +
        stat_boxplot(geom = "errorbar", width = 0.3, linewidth = 0.6) +
        geom_boxplot(colour = "black", outlier.shape = NA, width = 0.6) +
        geom_jitter(aes(color = call), width = 0.2, size = 1.2, alpha = 0.7) +
        scale_fill_manual(values = selex_colors) +
        scale_color_manual(values = selex_colors) +
        scale_y_continuous(limits = c(-1, 1.1), breaks = seq(-1, 1, 0.2)) +
        labs(x = "SELEX Group", y = "Correlation") +
        theme_classic(base_size = 13) +
        theme(legend.position = "none", axis.title = element_text(face = "plain")) +
        geom_text(
          data = counts, inherit.aes = FALSE,
          aes(x = call, y = 1.06, label = count), vjust = 0, size = 4.5
        )

      file <- file.path(plot.dir, "correlation_boxplot_by_selex.pdf")
      ggsave(file, p_selex, width = 5.5, height = 6)
      log_info("Wrote ", file)

      write.csv(cor_df, file.path(table.dir, "echo_correlation_by_selex.csv"),
        row.names = FALSE
      )
    }
  }
}

#####################################################################
# The assembled figure: contribution and SELEX over the paired heatmap
#####################################################################

if (is.null(p_modality) || is.null(p_selex)) {
  log_warn("Panel A or the SELEX panel is missing, the combined figure is skipped")
} else {
  top_row <- p_modality + p_selex + plot_layout(widths = c(2, 1))
  combined <- (top_row / wrap_elements(full = ht_grob)) +
    plot_layout(heights = c(1, 1.7)) +
    plot_annotation(tag_levels = "A")

  file <- file.path(plot.dir, "figure_echo_integration.pdf")
  ggsave(file, combined, width = 16, height = 16, limitsize = FALSE)
  log_info("Wrote ", file)
}

log_success("Finished the ECHO integration panels, figures in ", plot.dir)

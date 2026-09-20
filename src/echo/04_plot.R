#!/usr/bin/env Rscript

#####################################################################
# 04_plot.R
# created on 25-08-2026 by Irem B Gunduz
# The ECHO integration panels that 01 to 04 do not already produce:
#   - figure_echo_integration: the modality contribution, the factor
#     strip plot and the paired methylTFR / chromVAR heatmap
#   - sup_figure_echo_integration: the correlation boxplot by methyl-SELEX
#     call, the sample wise agreement of the two modalities and the
#     correlation distribution by cell type, faceted by latent factor
#
# Nothing here retrains MOFA, it reads the model and the tables 01 wrote.
#####################################################################

suppressPackageStartupMessages({
  library(MOFA2)
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

# Rows the heatmap shows for each factor, split evenly over the two views
tfs.per.factor <- 10
factors.keep <- c("Factor1", "Factor2", "Factor3", "Factor4", "Factor7")
# These two scatters always. The rest are filled with the best marker of a
# cell type, taken from the heatmap rows so every scatter is a row the
# reader can find again in the heatmap.
scatter.tfs <- c("SPIB", "POU2F3", "TBX20", "FOXL1")
scatter.panels <- 6
# The remaining panels are filled from this lineage
scatter.lineage <- "T-cell"
# A motif has to pull its cell type at least this far, in pooled standard
# deviations and in both modalities, before it earns a panel
marker.min.score <- 0.75
# Of those, one panel goes to a motif the two modalities agree on
# positively, so the figure is not read as methylTFR always opposing chromVAR
scatter.include.positive <- TRUE

# The four T cell subsets sit on top of each other, so ellipses are drawn
# per lineage and read as one T cell cluster rather than four outlines
lineage_of <- c(
  "B-cell" = "B-cell",
  "Monocyte" = "Monocyte",
  "NK-cell" = "NK-cell",
  "Tc-Mem" = "T-cell",
  "Th-Mem" = "T-cell",
  "Tc-Naive" = "T-cell",
  "Th-Naive" = "T-cell"
)
ellipse.min.samples <- 4
ellipse.level <- 0.9
# A lineage is circled when its centroid sits at least this many pooled
# standard deviations away from the centroid of the rest
ellipse.min.separation <- 1.0

# Modality contribution and factor R2, written by 01. Panel A is redrawn
# from these rather than recomputed, because the fraction needs the full
# weight matrix and the factor table holds only the top features.
contrib.file <- "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR_manuscript/tables/echo_mofa_modality_contribution.csv"
r2.file <- "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR_manuscript/tables/echo_mofa_factor_celltype_R2.csv"

# Trained model from 01, read for the factor values of the strip plot
model.rds <- "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR_manuscript/figures/echo/mofa_integration/MOFAobject_trained.rds"

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
  c(-2, -1, 0, 1, 2),
  c("#2B4B9B", "#4EC3E0", "#FFF200", "#F5A623", "#B21212")
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

# Illustrator sees every opaque panel, plot and legend background as its own
# white rectangle to select and delete, so none of them are drawn
no_bg <- theme(
  plot.background = element_blank(),
  panel.background = element_blank(),
  legend.background = element_blank(),
  legend.box.background = element_blank(),
  legend.key = element_blank(),
  strip.background = element_blank()
)

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

# Factors ordered by the variance cell type explains, as panel A orders them
heat_factors <- unique(ft$factor)
if (file.exists(r2.file)) {
  ranked <- read.csv(r2.file, stringsAsFactors = FALSE)
  ranked <- ranked$Factor[order(-ranked$R2)]
  heat_factors <- c(intersect(ranked, heat_factors), setdiff(heat_factors, ranked))
}

# A flat row carries no signal and would draw as a blank stripe, so it is
# out of the pool on both sides
constant_rows <- function(m, rows) {
  s <- apply(m[rows, , drop = FALSE], 1, sd, na.rm = TRUE)
  rows[!is.finite(s) | s == 0]
}
pool <- intersect(rownames(mtfr), rownames(chromVar))
flat <- union(constant_rows(mtfr, pool), constant_rows(chromVar, pool))
if (length(flat) > 0) {
  log_warn(length(flat), " flat row(s) dropped: ", paste(head(flat, 10), collapse = ", "))
  pool <- setdiff(pool, flat)
}

# The view is kept. MOFA weights are not on a common scale across views, so
# collapsing them with a max ranks by whichever view carries the larger
# numbers and buries the features the other one calls top.
candidates <- ft %>%
  filter(feature %in% pool) %>%
  group_by(factor, view, feature) %>%
  summarise(value_abs = max(value_abs), .groups = "drop") %>%
  as.data.frame()
absent <- setdiff(unique(ft$feature), pool)
if (length(absent) > 0) {
  log_warn(length(absent), " MOFA features are not in both matrices: ",
    paste(head(absent, 10), collapse = ", "))
}
if (nrow(candidates) == 0) stop("No MOFA feature is in both matrices")

views <- unique(candidates$view)
per_view <- ceiling(tfs.per.factor / length(views))

taken <- character(0)
selected <- vector("list", length(heat_factors))
names(selected) <- heat_factors

# One row per factor per view per round, so each factor is described by the
# top features of both modalities and no factor drains the shared pool
for (i in seq_len(per_view)) {
  for (f in heat_factors) {
    for (v in views) {
      if (NROW(selected[[f]]) >= tfs.per.factor) next
      pick <- candidates %>%
        filter(factor == f, view == v, !feature %in% taken) %>%
        arrange(desc(value_abs)) %>%
        head(1)
      if (nrow(pick) == 0) next
      selected[[f]] <- bind_rows(selected[[f]], pick)
      taken <- c(taken, pick$feature)
    }
  }
}

feature_factor <- bind_rows(selected[heat_factors]) %>%
  arrange(match(factor, heat_factors), view, desc(value_abs)) %>%
  as.data.frame()

log_info(
  "Heatmap TFs per factor: ",
  paste(names(table(feature_factor$factor)), table(feature_factor$factor),
    sep = " = ", collapse = ", "
  )
)
log_info(
  "Heatmap TFs per view: ",
  paste(names(table(feature_factor$view)), table(feature_factor$view),
    sep = " = ", collapse = ", "
  )
)

common_tfs <- feature_factor$feature
if (length(common_tfs) < 2) stop("Fewer than two MOFA features are in both matrices")
log_info(length(common_tfs), " features on the heatmap")

#####################################################################
# methyl-SELEX calls
#####################################################################

tf_root <- function(x) {
  parts <- unlist(strsplit(x, "::", fixed = TRUE))
  toupper(trimws(sub("\\s*\\(var\\.[0-9]+\\)$", "", parts)))
}

selex_calls <- NULL
if (!file.exists(selex.file)) {
  log_warn("SELEX table not found, the SELEX column and panel are skipped: ", selex.file)
} else {
  selex_tab <- read.csv(selex.file, skip = 20, header = TRUE, sep = ";")
  call_col <- intersect(
    c("methyl.SELEX.call", "Call", "methylSELEXcall"), colnames(selex_tab)
  )
  name_col <- intersect(c("TF.name", "TF", "TFname"), colnames(selex_tab))
  if (length(call_col) == 0 || length(name_col) == 0) {
    log_warn(
      "No SELEX call or TF name column. Present: ",
      paste(colnames(selex_tab), collapse = ", ")
    )
  } else {
    log_info("SELEX columns: name = ", name_col[1], ", call = ", call_col[1])
    selex_calls <- setNames(
      as.character(selex_tab[[call_col[1]]]),
      toupper(as.character(selex_tab[[name_col[1]]]))
    )
  }
}

# A TF with no entry was never assayed rather than missing, so it reads
# Inconclusive
selex_call_of <- function(tfs) {
  if (is.null(selex_calls)) {
    return(NULL)
  }
  vapply(tfs, function(t) {
    hit <- tf_root(t)
    hit <- hit[hit %in% names(selex_calls)]
    if (length(hit) == 0) {
      return("Inconclusive")
    }
    call <- selex_calls[[hit[1]]]
    if (is.na(call) || !nzchar(call)) "Inconclusive" else call
  }, character(1), USE.NAMES = FALSE)
}

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
    ggsave(file, p_modality & no_bg, width = 8, height = 4.5, bg = "transparent")
    log_info("Wrote ", file)
  }
}

#####################################################################
# Panel B, factor values per sample
#####################################################################

p_strip <- NULL
if (!file.exists(model.rds)) {
  log_warn("Trained model not found, run 01 first. Panel B is skipped: ", model.rds)
} else {
  factors_long <- get_factors(readRDS(model.rds), as.data.frame = TRUE)
  factors_long$sample <- as.character(factors_long$sample)
  factors_long$celltype <- sub("_.*", "", factors_long$sample)
  factors_long <- factors_long[factors_long$factor %in% heat_factors, ]

  unmapped_strip <- setdiff(unique(factors_long$celltype), names(cell_type_colors))
  if (length(unmapped_strip) > 0) {
    stop("No colour for cell type(s): ", paste(unmapped_strip, collapse = ", "))
  }
  factors_long$factor <- factor(factors_long$factor, levels = heat_factors)
  factors_long$celltype <- factor(factors_long$celltype, levels = names(cell_type_colors))

  p_strip <- ggplot(factors_long, aes(x = factor, y = value, color = celltype)) +
    geom_jitter(width = 0.2, height = 0, size = 2, alpha = 0.8) +
    scale_color_manual(values = cell_type_colors, name = "Cell type") +
    labs(x = "Factor", y = "Factor value") +
    theme_classic(base_size = 13) +
    theme(panel.grid.major.x = element_blank())

  file <- file.path(plot.dir, "factors_stripplot_by_celltype.pdf")
  ggsave(file, p_strip & no_bg, width = 8, height = 5, bg = "transparent")
  log_info("Wrote ", file)
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
row_split_factor <- factor(feature_factor$factor,
  levels = heat_factors[heat_factors %in% feature_factor$factor]
)
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
  row_split = row_split_factor,
  cluster_row_slices = FALSE,
  row_gap = unit(1, "mm"),
  row_title_gp = gpar(fontsize = 8),
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

selex_call <- selex_call_of(common_tfs)
selex_col_map <- NULL
if (!is.null(selex_call)) {
  log_info(
    "SELEX calls on the heatmap: ",
    paste(names(table(selex_call)), table(selex_call), sep = " = ", collapse = ", ")
  )
  known <- c(
    "Inconclusive" = "#D9D9D9",
    "Little effect" = "#BDB76B",
    selex_colors
  )
  extra <- setdiff(unique(selex_call), names(known))
  if (length(extra) > 0) {
    log_warn("Unlisted SELEX call(s), coloured grey: ", paste(extra, collapse = ", "))
    known <- c(known, setNames(rep("#8C8C8C", length(extra)), extra))
  }
  selex_col_map <- known[intersect(names(known), unique(selex_call))]
}

# The factor is already the row split, so it needs no colour column here
annotation_args <- list(
  `Row-wise Pearson Correlation` = row_cor,
  col = list(`Row-wise Pearson Correlation` = cor_colors),
  annotation_name_gp = gpar(fontsize = 7),
  width = unit(30, "mm")
)
if (!is.null(selex_call)) {
  annotation_args[["SELEX Call"]] <- selex_call
  annotation_args$col[["SELEX Call"]] <- selex_col_map
}
# which = "row" is explicit because do.call evaluates anno_text() before
# rowAnnotation can infer it
annotation_args$TF <- anno_text(
  common_tfs,
  gp = gpar(fontsize = 6), which = "row"
)

ha_right <- do.call(rowAnnotation, annotation_args)

# Captured once as a grob, so the same drawing serves the standalone PDF
# and the assembled figure below
ht_grob <- grid.grabExpr(
  draw(ht_mtfr + ht_chromvar + ha_right,
    merge_legend = TRUE, heatmap_legend_side = "right"
  )
)

file <- file.path(plot.dir, "paired_heatmap_mtfr_chromvar.pdf")
pdf(file, width = 14, height = 10, bg = "transparent")
grid.draw(ht_grob)
dev.off()
log_info("Wrote ", file)

#####################################################################
# Per TF scatters, chromVAR against methylTFR
#####################################################################

# Circling every group marks the crowd in the middle as though it were a
# cluster. Only a lineage whose centroid stands away from the rest, on the
# scale of the scatter itself, is worth outlining.
separating_groups <- function(df) {
  sds <- c(sd(df$chromVar, na.rm = TRUE), sd(df$mtfr, na.rm = TRUE))
  sds[!is.finite(sds) | sds == 0] <- 1
  keep <- character(0)
  for (g in unique(df$lineage)) {
    inside <- df$lineage == g
    if (sum(inside) < ellipse.min.samples || all(inside)) next
    delta <- c(
      mean(df$chromVar[inside]) - mean(df$chromVar[!inside]),
      mean(df$mtfr[inside]) - mean(df$mtfr[!inside])
    ) / sds
    if (sqrt(sum(delta^2)) >= ellipse.min.separation) keep <- c(keep, g)
  }
  keep
}

plot_tf_scatter <- function(tf) {
  df <- data.frame(
    chromVar = as.numeric(chromVar[tf, ]),
    mtfr = as.numeric(mtfr[tf, ]),
    celltype = cell_types,
    lineage = unname(lineage_of[cell_types]),
    stringsAsFactors = FALSE
  )
  cor_val <- suppressWarnings(cor(df$chromVar, df$mtfr,
    method = "pearson", use = "complete.obs"
  ))

  groups <- separating_groups(df)
  log_info(tf, ": ellipse on ",
    if (length(groups) == 0) "no lineage" else paste(groups, collapse = ", "))
  ell <- df[df$lineage %in% groups, , drop = FALSE]

  p <- ggplot(df, aes(x = chromVar, y = mtfr)) +
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

  if (nrow(ell) > 0) {
    p <- p + stat_ellipse(
      data = ell, aes(group = lineage),
      type = "norm", level = ellipse.level, colour = "grey35",
      linetype = "dashed", linewidth = 0.4, show.legend = FALSE
    )
  }
  p
}

available_tfs <- intersect(scatter.tfs, pool)
missing_scatter <- setdiff(scatter.tfs, available_tfs)
if (length(missing_scatter) > 0) {
  log_warn("Not in both matrices, no scatter: ", paste(missing_scatter, collapse = ", "))
}

# A TF marks a cell type when both modalities pull that cell type away from
# the rest, so the score is the weaker of the two standardised differences.
# Only heatmap rows are scored, so every scatter is a row of the heatmap.
marker_score <- function(tf, ct) {
  inside <- cell_types == ct
  gap <- function(m) {
    v <- as.numeric(m[tf, ])
    sdv <- sd(v, na.rm = TRUE)
    if (!is.finite(sdv) || sdv == 0) {
      return(0)
    }
    abs(mean(v[inside], na.rm = TRUE) - mean(v[!inside], na.rm = TRUE)) / sdv
  }
  min(gap(mtfr), gap(chromVar))
}

marker_tab <- expand.grid(
  tf = common_tfs, celltype = unique(cell_types),
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
)
marker_tab$score <- mapply(marker_score, marker_tab$tf, marker_tab$celltype)
marker_tab <- marker_tab[order(-marker_tab$score), ]

write.csv(marker_tab, file.path(table.dir, "echo_marker_scores.csv"), row.names = FALSE)

# Each motif is scored against the one cell type it separates best. The
# score gates the extra panels, so only a motif that visibly pulls a cell
# type apart in both modalities is drawn.
tf_best <- marker_tab %>%
  group_by(tf) %>%
  slice_max(score, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  as.data.frame()

tf_cor <- vapply(common_tfs, function(tf) {
  suppressWarnings(cor(
    as.numeric(mtfr[tf, ]), as.numeric(chromVar[tf, ]),
    method = "pearson", use = "complete.obs"
  ))
}, numeric(1))
tf_best$correlation <- unname(tf_cor[tf_best$tf])
tf_best$lineage <- unname(lineage_of[tf_best$celltype])

extra_pool <- tf_best %>%
  filter(
    score >= marker.min.score,
    lineage == scatter.lineage,
    !tf %in% available_tfs
  ) %>%
  arrange(desc(score))

need <- scatter.panels - length(available_tfs)
if (nrow(extra_pool) == 0) {
  log_warn(
    "No ", scatter.lineage, " motif reaches a score of ", marker.min.score,
    ", the scatters stay at ", length(available_tfs)
  )
} else {
  # The best positively correlated one first, when nothing chosen so far is
  positive_first <- character(0)
  if (scatter.include.positive && !any(tf_cor[available_tfs] > 0, na.rm = TRUE)) {
    positive_first <- head(extra_pool$tf[extra_pool$correlation > 0], 1)
    if (length(positive_first) == 0) {
      log_warn("No ", scatter.lineage, " motif clearing the score correlates positively")
    }
  }
  picked <- head(unique(c(positive_first, extra_pool$tf)), need)
  log_info(
    scatter.lineage, " panels: ",
    paste(picked, round(tf_best$score[match(picked, tf_best$tf)], 2),
      sep = " score ", collapse = ", ")
  )
  available_tfs <- c(available_tfs, picked)
}

if (length(available_tfs) < scatter.panels) {
  log_warn("Only ", length(available_tfs), " of ", scatter.panels, " scatter panels")
}

p_scatter <- NULL
if (length(available_tfs) > 0) {
  p_scatter <- wrap_plots(
    lapply(available_tfs, function(tf) {
      plot_tf_scatter(tf)
    }),
    ncol = 2, guides = "collect"
  )
  file <- file.path(plot.dir, "scatter_chromvar_vs_mtfr_panels.pdf")
  ggsave(file, p_scatter & no_bg,
    width = 11, height = 4.5 * ceiling(length(available_tfs) / 2),
    bg = "transparent"
  )
  log_info("Wrote ", file)
}

#####################################################################
# mTFR / chromVAR correlation by cell type, faceted by latent factor
#####################################################################

# Correlation within each cell type, so a cell type needs enough samples
min.corr.samples <- 3

tf_top_factor <- ft %>%
  group_by(feature) %>%
  slice_max(abs(value), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(feature, factor) %>%
  as.data.frame()
factor_of <- setNames(tf_top_factor$factor, tf_top_factor$feature)

ct_levels <- intersect(names(cell_type_colors), unique(cell_types))
factor_cor <- do.call(rbind, lapply(ct_levels, function(ct) {
  cols <- which(cell_types == ct)
  if (length(cols) < min.corr.samples) {
    log_warn(ct, ": ", length(cols), " samples, no per cell type correlation")
    return(NULL)
  }
  vals <- vapply(common_tfs, function(tf) {
    suppressWarnings(cor(
      as.numeric(mtfr[tf, cols]), as.numeric(chromVar[tf, cols]),
      method = "pearson", use = "complete.obs"
    ))
  }, numeric(1))
  data.frame(tf = common_tfs, correlation = unname(vals),
    celltype = ct, stringsAsFactors = FALSE)
}))

factor_cor$factor <- factor_of[factor_cor$tf]
factor_cor <- factor_cor[!is.na(factor_cor$factor) & is.finite(factor_cor$correlation), ]
factor_cor$factor <- droplevels(factor(factor_cor$factor, levels = heat_factors))
factor_cor$celltype <- factor(factor_cor$celltype, levels = ct_levels)

p_factor_cor <- ggplot(factor_cor, aes(x = celltype, y = correlation, fill = celltype)) +
  stat_boxplot(geom = "errorbar", width = 0.3, linewidth = 0.5) +
  geom_boxplot(colour = "black", outlier.shape = NA, width = 0.6) +
  scale_fill_manual(values = cell_type_colors) +
  scale_y_continuous(breaks = seq(-1, 1, 0.5)) +
  coord_cartesian(ylim = c(-1, 1)) +
  facet_wrap(~factor) +
  labs(x = "Cell type", y = "mTFR / chromVAR correlation") +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "none",
    axis.title = element_text(face = "plain"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(file.path(plot.dir, "correlation_boxplot_by_celltype_factor.pdf"),
  p_factor_cor & no_bg, width = 9, height = 7, bg = "transparent")
log_info("Wrote the per cell type correlation boxplots by latent factor")

#####################################################################
# Sample wise agreement of the two modalities, by cell type
#####################################################################

# The heatmap and the SELEX panel both read across samples for one TF. This
# reads the other way, across TFs within one sample, so a cell type whose
# two modalities disagree shows up as a low box rather than being averaged
# away.
sample_cor <- vapply(common_samples, function(sm) {
  suppressWarnings(cor(
    as.numeric(mtfr[pool, sm]), as.numeric(chromVar[pool, sm]),
    method = "pearson", use = "complete.obs"
  ))
}, numeric(1))

sample_cor_df <- data.frame(
  sample = common_samples,
  celltype = factor(cell_types, levels = names(cell_type_colors)),
  correlation = sample_cor,
  stringsAsFactors = FALSE
)
sample_cor_df <- sample_cor_df[is.finite(sample_cor_df$correlation), ]

p_sample_cor <- ggplot(sample_cor_df, aes(x = celltype, y = correlation, fill = celltype)) +
  stat_boxplot(geom = "errorbar", width = 0.3, linewidth = 0.6) +
  geom_boxplot(colour = "black", outlier.shape = NA, width = 0.6) +
  geom_jitter(aes(fill = celltype),
    width = 0.2, size = 1.6, shape = 21, colour = "grey15", stroke = 0.3
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4, colour = "grey40") +
  scale_fill_manual(values = cell_type_colors) +
  labs(
    x = "Cell type",
    y = paste0("Correlation across ", length(pool), " TFs")
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

file <- file.path(plot.dir, "sample_correlation_by_celltype.pdf")
ggsave(file, p_sample_cor & no_bg, width = 6, height = 5, bg = "transparent")
log_info("Wrote ", file)

write.csv(sample_cor_df, file.path(table.dir, "echo_sample_correlation_by_celltype.csv"),
  row.names = FALSE
)

#####################################################################
# Correlation by methyl-SELEX call, over every shared TF
#####################################################################

p_selex <- NULL
if (is.null(selex_calls)) {
  log_warn("No SELEX calls, the correlation boxplot is skipped")
} else {
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
  cor_df$call <- selex_call_of(cor_df$TF.name)
  cor_df <- cor_df[cor_df$call %in% selex.groups, ]

  if (nrow(cor_df) == 0) {
    log_warn("No TF carries one of ", paste(selex.groups, collapse = " or "))
  } else {
    cor_df$call <- factor(cor_df$call, levels = selex.groups)
    counts <- as.data.frame(table(cor_df$call))
    colnames(counts) <- c("call", "count")
    log_info("SELEX groups: ",
      paste(counts$call, counts$count, sep = " = ", collapse = ", "))

    # Two sided Wilcoxon rank sum on the correlations, one test per pair of
    # SELEX groups, Benjamini-Hochberg over those tests
    groups_present <- levels(droplevels(cor_df$call))
    tests <- NULL
    if (length(groups_present) >= 2) {
      tests <- bind_rows(lapply(combn(groups_present, 2, simplify = FALSE), function(g) {
        x <- cor_df$Correlation[cor_df$call == g[1]]
        y <- cor_df$Correlation[cor_df$call == g[2]]
        data.frame(
          group1 = g[1], group2 = g[2], n1 = length(x), n2 = length(y),
          median1 = median(x), median2 = median(y),
          p_value = if (length(x) < 3 || length(y) < 3) {
            NA_real_
          } else {
            suppressWarnings(wilcox.test(x, y)$p.value)
          },
          stringsAsFactors = FALSE
        )
      }))
      tests$p_adjusted <- p.adjust(tests$p_value, method = "BH")

      write.csv(tests, file.path(table.dir, "echo_selex_correlation_test.csv"),
        row.names = FALSE
      )
      log_info(
        "Wilcoxon: ",
        paste(tests$group1, " vs ", tests$group2, " p.adj = ",
          signif(tests$p_adjusted, 3),
          collapse = "; "
        )
      )

      tests$x1 <- match(tests$group1, groups_present)
      tests$x2 <- match(tests$group2, groups_present)
      tests$y <- 0.72 + (seq_len(nrow(tests)) - 1) * 0.16
      tests$label <- ifelse(is.na(tests$p_adjusted), "n.a.",
        paste0("p.adj = ", format.pval(tests$p_adjusted, digits = 2, eps = 1e-16))
      )
    }

    p_selex <- ggplot(cor_df, aes(x = call, y = Correlation, fill = call)) +
      stat_boxplot(geom = "errorbar", width = 0.3, linewidth = 0.6) +
      geom_boxplot(colour = "black", outlier.shape = NA, width = 0.6) +
      geom_jitter(aes(color = call), width = 0.2, size = 1.2, alpha = 0.7) +
      scale_fill_manual(values = selex_colors) +
      scale_color_manual(values = selex_colors) +
      scale_y_continuous(breaks = seq(-1, 1, 0.2)) +
      coord_cartesian(ylim = c(-1, 1.12), clip = "off") +
      labs(x = "SELEX Group", y = "Correlation") +
      theme_classic(base_size = 13) +
      theme(legend.position = "none", axis.title = element_text(face = "plain")) +
      geom_text(
        data = counts, inherit.aes = FALSE,
        aes(x = call, y = 1.06, label = count), vjust = 0, size = 4.5
      )

    if (!is.null(tests)) {
      p_selex <- p_selex +
        geom_segment(
          data = tests, inherit.aes = FALSE,
          aes(x = x1, xend = x2, y = y, yend = y), linewidth = 0.5
        ) +
        geom_segment(
          data = tests, inherit.aes = FALSE,
          aes(x = x1, xend = x1, y = y - 0.03, yend = y), linewidth = 0.5
        ) +
        geom_segment(
          data = tests, inherit.aes = FALSE,
          aes(x = x2, xend = x2, y = y - 0.03, yend = y), linewidth = 0.5
        ) +
        geom_text(
          data = tests, inherit.aes = FALSE,
          aes(x = (x1 + x2) / 2, y = y + 0.02, label = label),
          vjust = 0, size = 3.8, fontface = "bold"
        )
    }

    file <- file.path(plot.dir, "correlation_boxplot_by_selex.pdf")
    ggsave(file, p_selex & no_bg, width = 5.5, height = 6, bg = "transparent")
    log_info("Wrote ", file)

    write.csv(cor_df, file.path(table.dir, "echo_correlation_by_selex.csv"),
      row.names = FALSE
    )
  }
}

#####################################################################
# The main figure: contribution, factor values, paired heatmap
#####################################################################

tag_theme <- theme(plot.tag = element_text(face = "bold", size = 16))

if (is.null(p_modality) || is.null(p_strip)) {
  log_warn("Panel A or B is missing, the main figure is skipped")
} else {
  combined <- wrap_plots(
    wrap_plots(
      p_modality + labs(tag = "A"),
      p_strip + labs(tag = "B"),
      nrow = 1, widths = c(1, 1.15)
    ),
    wrap_elements(full = ht_grob) + labs(tag = "C"),
    ncol = 1, heights = c(1, 2.1)
  ) & tag_theme

  file <- file.path(plot.dir, "figure_echo_integration.pdf")
  ggsave(file, combined & no_bg, width = 18, height = 17, bg = "transparent", limitsize = FALSE)
  log_info("Wrote ", file)
}

#####################################################################
# The supplementary figure: SELEX, sample agreement, per TF scatters
#####################################################################

if (is.null(p_selex) || is.null(p_factor_cor)) {
  log_warn("The SELEX panel or the per factor panel is missing, the supplementary figure is skipped")
} else {
  sup <- wrap_plots(
    wrap_plots(
      p_selex + labs(tag = "A"),
      p_sample_cor + labs(tag = "B"),
      nrow = 1, widths = c(1, 1.2)
    ),
    p_factor_cor + labs(tag = "C"),
    ncol = 1, heights = c(1, 1.6)
  ) & tag_theme

  file <- file.path(plot.dir, "sup_figure_echo_integration.pdf")
  ggsave(file, sup & no_bg,
    width = 12, height = 15,
    bg = "transparent", limitsize = FALSE
  )
  log_info("Wrote ", file)
}

log_success("Finished the ECHO integration panels, figures in ", plot.dir)

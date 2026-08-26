#!/usr/bin/env Rscript

#####################################################################
# 08_bp_mofa.R
# created on 25-10-2025 by Irem B Gunduz and Regina Nitsch
# Updated by IBG on 25-08-2026
# MOFA2 integration of the methylTFR deviations and the matching RNA
# counts of the Blueprint samples, the factor and modality figures,
# and the paired heatmaps with the SELEX annotation
#
# The WGBS to RNA sample map comes from 07 rather than being derived
# here a second time.
#####################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(grid)
  library(circlize)
  library(RColorBrewer)
  library(ComplexHeatmap)
  library(logger)
  library(MOFA2)
  library(reticulate)
  library(methylTFR)
})
set.seed(12)

#####################################################################
# Settings
#####################################################################

use_python(Sys.which("python"), required = TRUE)
motifSet <- "jaspar2020_distal"

# Cell type groups excluded from the figures. Only the "other" catch all
# is dropped, plasma is kept and shown.
drop.cell.types <- "Other"

num.factors <- 14
top.factors <- 7 # factors carried into the strip plot
top.tfs <- 5 # top features per factor and view
drop.factor1 <- TRUE # Factor 1 tracks the number of expressed features

# Display names for the raw cellTypeGroup values of the annotation. The
# order here is the legend order and the heatmap column order.
cell_type_labels <- c(
  "megK" = "Megakaryocytes",
  "eryt" = "Erythrocytes",
  "gran" = "Granulocytes",
  "mono" = "Monocytes",
  "Mf" = "Macrophages",
  "DC" = "Dendritic cells",
  "osteoclast" = "Osteoclast",
  "NK" = "NK",
  "Tcell" = "T-cells",
  "thymocyte" = "Thymocyte",
  "Bcell" = "B-cells",
  "plasma" = "Plasma",
  "progenitor" = "Progenitors",
  "other" = "Other"
)
cell_type_levels <- unname(cell_type_labels)

# Keyed on the display names above. The previous palette had "Progenitors"
# where the annotation says "progenitor", so that group came out grey.
cell_type_colors <- c(
  "Megakaryocytes" = "#4a0221",
  "Erythrocytes" = "#67000d",
  "Granulocytes" = "#ff7f50",
  "Monocytes" = "#CD7054",
  "Macrophages" = "#864a38",
  "Dendritic cells" = "#EED5B7",
  "Osteoclast" = "#DEB887",
  "NK" = "#bf812d",
  "T-cells" = "#40E0D0",
  "Thymocyte" = "#74c476",
  "B-cells" = "#980043",
  "Plasma" = "#CD2990",
  "Progenitors" = "#df65b0",
  "Other" = "#8B8682"
)

# Column order of the heatmaps, in display names
column.order <- setdiff(cell_type_levels, c("Other", "Progenitors"))

# Factor pairs of the scatter panel, one figure each. NULL falls back to
# the two factors with the highest R2.
scatter.pairs <- list(
  c("Factor5", "Factor2"),
  c("Factor3", "Factor8")
)

view_colors <- c("mtfr" = "#ED4B4A", "expr" = "#6BC75A")

# Directories
analysis.dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/blueprint"
dev.tag <- "mTFR_devs_230826"
mofa.dir <- file.path(analysis.dir, "mofa_230826")

dev.file <- file.path(analysis.dir, dev.tag, paste0(motifSet, "_deviations.RDS"))
counts.file <- file.path(mofa.dir, "bp_rawRNAcounts.RDS")
map.file <- file.path(mofa.dir, "bp_rna_wgbs_map.tsv")

# run_mofa writes an HDF5, the trained object is a separate RDS. The old
# script pointed both at the same path, so the HDF5 was overwritten.
model.hdf5 <- file.path(mofa.dir, "mtfr_expr_model_bp.hdf5")
model.rds <- file.path(mofa.dir, "mtfr_expr_model_bp.rds")

plot.dir <- file.path(analysis.dir, "mofa_figures_230826")
if (!dir.exists(plot.dir)) dir.create(plot.dir, recursive = TRUE)

table.dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR_manuscript/tables"
if (!dir.exists(table.dir)) dir.create(table.dir, recursive = TRUE)

selex.file <- "/icbb/projects/igunduz/exposure_atlas_manuscript/sample_annots/Selex_data.csv"

#####################################################################
# Helper functions
#####################################################################

row_zscore <- function(mat) {
  mat <- as.matrix(mat)
  centred <- mat - rowMeans(mat, na.rm = TRUE)
  sds <- apply(mat, 1, sd, na.rm = TRUE)
  out <- centred / sds
  out[!is.finite(out)] <- 0
  out
}

#####################################################################
# Load the two modalities and the sample map
#####################################################################

for (f in c(dev.file, counts.file, map.file)) {
  if (!file.exists(f)) stop("Missing input, run 02 and 07 first: ", f)
}

log_info("Loading deviations from ", dev.file)
mtfr <- readRDS(dev.file)
rna <- readRDS(counts.file)
map <- read.delim(map.file, stringsAsFactors = FALSE)
log_info(nrow(map), " paired samples in the map")

# Keep the pairs present in both modalities, in one shared order
map <- map[map$bedFile %in% colnames(mtfr) & map$rna_id %in% colnames(rna), ]
if (nrow(map) == 0) {
  stop(
    "No sample is present in both modalities. Deviations carry e.g. '",
    colnames(mtfr)[1], "', the counts carry '", colnames(rna)[1], "'."
  )
}
log_info(nrow(map), " samples shared by the deviations and the counts")

mtfr_z <- row_zscore(deviations(mtfr))[, map$bedFile, drop = FALSE]
rna_mat <- rna[, map$rna_id, drop = FALSE]

# One sample name for both views
sample_ids <- map$rna_id
colnames(mtfr_z) <- sample_ids
colnames(rna_mat) <- sample_ids

# The RNA side is restricted to the genes that are also motifs
common_features <- intersect(rownames(rna_mat), rownames(mtfr_z))
log_info(length(common_features), " features shared by the motifs and the genes")
if (length(common_features) < 10) {
  stop("Too few shared features to run MOFA")
}
rna_z <- row_zscore(rna_mat[common_features, , drop = FALSE])

stopifnot(identical(colnames(mtfr_z), colnames(rna_z)))

#####################################################################
# MOFA
#####################################################################

data_list <- list(mtfr = as.matrix(mtfr_z), expr = as.matrix(rna_z))
MOFAobject <- create_mofa(data_list)

# Sample names come from the aligned matrices, not from the unfiltered
# deviations object, which held every Blueprint sample
unmapped <- setdiff(unique(map$cellTypeGroup), names(cell_type_labels))
if (length(unmapped) > 0) {
  stop("No display name for cell type(s): ", paste(unmapped, collapse = ", "))
}
metadata <- data.frame(
  sample = sample_ids,
  celltype = factor(
    unname(cell_type_labels[map$cellTypeGroup]),
    levels = cell_type_levels
  ),
  bedFile = map$bedFile,
  stringsAsFactors = FALSE
)
samples_metadata(MOFAobject) <- metadata

model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- min(num.factors, ncol(mtfr_z) - 1)
log_info("Training with up to ", model_opts$num_factors, " factors")

MOFAobject <- prepare_mofa(MOFAobject, model_options = model_opts)

if (file.exists(model.rds)) {
  log_info("Loading the trained model from ", model.rds)
  MOFAobject.trained <- readRDS(model.rds)
} else {
  MOFAobject.trained <- run_mofa(MOFAobject, model.hdf5, use_basilisk = FALSE)
  saveRDS(MOFAobject.trained, model.rds)
  log_success("Wrote ", model.rds)
}

# Drop the groups that are not shown
metadata <- samples_metadata(MOFAobject.trained)
keep_samples <- metadata$sample[!metadata$celltype %in% drop.cell.types]
log_info("Keeping ", length(keep_samples), " of ", nrow(metadata),
  " samples after dropping ", paste(drop.cell.types, collapse = ", "))
MOFAobject.trained <- subset_samples(MOFAobject.trained, keep_samples)
metadata <- samples_metadata(MOFAobject.trained)

#####################################################################
# Variance in the factors explained by cell type
#####################################################################

# The cell type is joined on the sample, not recycled: get_factors returns
# one row per sample AND factor, so assigning a per-sample vector silently
# mislabelled most rows.
factors_long <- get_factors(MOFAobject.trained, as.data.frame = TRUE)
factors_long$sample <- as.character(factors_long$sample)
factors_long <- left_join(
  factors_long, metadata[, c("sample", "celltype")],
  by = "sample"
)
if (any(is.na(factors_long$celltype))) {
  stop("Some factor rows could not be matched to a cell type")
}

# Fraction of the variance in a factor that cell type accounts for.
# Computed one factor at a time in base R: rename() is masked by
# S4Vectors once MOFA2 is attached, and cur_data_all() is deprecated.
factor_r2 <- function(df) {
  if (length(unique(df$celltype)) < 2) {
    return(NA_real_)
  }
  ss <- summary(aov(value ~ celltype, data = df))[[1]]
  ss["celltype", "Sum Sq"] / (ss["celltype", "Sum Sq"] + ss["Residuals", "Sum Sq"])
}

factor_split <- split(factors_long, as.character(factors_long$factor))

r2_df <- data.frame(
  Factor = names(factor_split),
  R2 = vapply(factor_split, factor_r2, numeric(1)),
  stringsAsFactors = FALSE,
  row.names = NULL
)

if (all(is.na(r2_df$R2))) {
  stop("No factor could be tested against cell type, only one group is present")
}

# Decreasing R2, with any untestable factor last
r2_df <- r2_df[order(-r2_df$R2), ]

if (drop.factor1) {
  r2_df <- r2_df[r2_df$Factor != "Factor1", ]
}
write.csv(r2_df, file.path(table.dir, "mofa_bp_factor_celltype_R2.csv"), row.names = FALSE)

top_factors <- head(r2_df$Factor, top.factors)
log_info("Top factors: ", paste(top_factors, collapse = ", "))

#####################################################################
# Factor scatter
#####################################################################

factors_wide <- factors_long %>%
  filter(!celltype %in% drop.cell.types) %>%
  select(sample, celltype, factor, value) %>%
  pivot_wider(names_from = factor, values_from = value)
factors_wide$celltype <- factor(factors_wide$celltype, levels = cell_type_levels)

available <- setdiff(colnames(factors_wide), c("sample", "celltype"))

# One figure per requested pair. NULL falls back to the two factors that
# explain the most cell type variance.
pairs <- scatter.pairs
if (is.null(pairs)) {
  pairs <- list(head(intersect(top_factors, available), 2))
}

scatter_plot <- function(f_x, f_y) {
  ggplot(factors_wide, aes(x = .data[[f_x]], y = .data[[f_y]], color = celltype)) +
    geom_point(alpha = 0.9, size = 2) +
    stat_ellipse(aes(group = celltype), linetype = 2, alpha = 0.5, linewidth = 0.4) +
    scale_color_manual(values = cell_type_colors, drop = FALSE, name = "Cell type") +
    labs(x = f_x, y = f_y) +
    theme_classic(base_size = 13) +
    theme(legend.position = "right")
}

for (pr in pairs) {
  missing_f <- setdiff(pr, available)
  if (length(missing_f) > 0) {
    log_warn("Skipping the scatter for ", paste(pr, collapse = " vs "),
      ", not among the fitted factors: ", paste(missing_f, collapse = ", "))
    next
  }
  p_scatter <- scatter_plot(pr[1], pr[2])
  file <- file.path(plot.dir, paste0("factor_scatter_", pr[1], "_", pr[2], ".pdf"))
  # The object name was misspelled here, so this figure never rendered
  ggsave(file, plot = p_scatter, width = 7, height = 5)
  log_info("Wrote ", file)
}

#####################################################################
# Modality contribution per factor
#####################################################################

weights_df <- get_weights(MOFAobject.trained, as.data.frame = TRUE) %>%
  mutate(value_abs = abs(value))

agg <- weights_df %>%
  group_by(factor, view) %>%
  summarise(sum_abs = sum(value_abs, na.rm = TRUE), .groups = "drop") %>%
  group_by(factor) %>%
  mutate(frac = sum_abs / sum(sum_abs)) %>%
  ungroup() %>%
  filter(factor %in% r2_df$Factor)

# Every retained factor is shown, ordered by the variance the cell type
# explains, so the R2 line falls from left to right
factor_order <- r2_df$Factor
agg$factor <- factor(agg$factor, levels = factor_order)

r2_line <- r2_df
r2_line$Factor <- factor(r2_line$Factor, levels = factor_order)

p_modality_frac <- ggplot(agg, aes(x = factor, y = frac)) +
  geom_bar(aes(fill = view), stat = "identity", width = 0.7) +
  geom_text(aes(group = view, label = ifelse(frac > 0.03, paste0(round(frac * 100, 1), "%"), "")),
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
  labs(
    x = "Factor", y = "Relative contribution (|weights| fraction)", fill = "View"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.y.right = element_text(colour = "purple")
  )

ggsave(file.path(plot.dir, "mtfr_modality_contribution_fraction_R2.pdf"),
  plot = p_modality_frac, width = 11, height = 4.5
)

#####################################################################
# Factor values per sample
#####################################################################

factors_top <- factors_long %>%
  filter(factor %in% top_factors, !celltype %in% drop.cell.types)
# The level vector referred to a column called Factors, which does not
# exist, so every level came out NA and the panel order was arbitrary
factors_top$factor <- factor(factors_top$factor, levels = top_factors)
factors_top$celltype <- factor(factors_top$celltype, levels = cell_type_levels)

p_strip <- ggplot(factors_top, aes(x = factor, y = value, color = celltype)) +
  geom_jitter(width = 0.2, height = 0, size = 2, alpha = 0.8) +
  scale_color_manual(values = cell_type_colors, drop = FALSE) +
  labs(x = "Factor", y = "Factor value", color = "Cell type") +
  theme_classic(base_size = 13) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

ggsave(file.path(plot.dir, "factors_stripplot_by_celltype.pdf"),
  plot = p_strip, width = 8, height = 5
)

#####################################################################
# Top features per factor and view
#####################################################################

# weights_df, not weights: the latter is a base R function, so this loop
# used to fail with "$ operator is invalid for atomic vectors"
n_factors_tf <- min(5, nrow(r2_df))
topTfs <- lapply(r2_df$Factor[seq_len(n_factors_tf)], function(f) {
  lapply(unique(weights_df$view), function(modality) {
    d <- weights_df[weights_df$view == modality & weights_df$factor == f, ]
    d <- d[order(d$value_abs, decreasing = TRUE), ]
    d <- head(d, top.tfs)[, c("feature", "value", "factor", "view")]
    d$feature <- gsub("_.*$", "", d$feature)
    d
  }) %>% bind_rows()
}) %>% bind_rows()

write.csv(topTfs, file.path(table.dir, "mofa_bp_topTFs.csv"), row.names = FALSE)
log_info(length(unique(topTfs$feature)), " unique features across the top factors")

#####################################################################
# Paired heatmaps of the MOFA features
#####################################################################

mofa_tfs <- unique(topTfs$feature)
common_TFs <- Reduce(intersect, list(mofa_tfs, rownames(rna_z), rownames(mtfr_z)))
if (length(common_TFs) < 2) stop("Fewer than two MOFA features are in both matrices")
log_info(length(common_TFs), " MOFA features present in both modalities")

# Restricted to the samples that are still in the model, and to the
# groups that are plotted
heat_samples <- intersect(metadata$sample, colnames(mtfr_z))
heat_groups <- metadata$celltype[match(heat_samples, metadata$sample)]
keep <- !heat_groups %in% drop.cell.types
heat_samples <- heat_samples[keep]
heat_groups <- heat_groups[keep]

mtfr_heatmap <- mtfr_z[common_TFs, heat_samples, drop = FALSE]
rna_heatmap <- rna_z[common_TFs, heat_samples, drop = FALSE]

present <- intersect(names(cell_type_colors), unique(heat_groups))
ha <- HeatmapAnnotation(
  celltypes = heat_groups,
  col = list(celltypes = cell_type_colors[present])
)
column_split_factor <- factor(heat_groups, levels = intersect(column.order, unique(heat_groups)))

heatmap_mtfr <- Heatmap(
  mtfr_heatmap,
  name = "mTFR Z-score",
  row_names_gp = gpar(fontsize = 4),
  top_annotation = ha,
  column_title = "methylTFR Z-scores of the MOFA TFs",
  show_row_names = TRUE,
  show_column_names = FALSE,
  column_split = column_split_factor,
  cluster_column_slices = FALSE
)

file <- file.path(plot.dir, "heatmap_mofa_mtfr_bp.pdf")
pdf(file, width = 8, height = 10)
ht_mtfr <- draw(heatmap_mtfr)
dev.off()
log_info("Wrote ", file)

# Row order taken from the drawn object, so the expression heatmap and the
# annotation columns line up with it
row_order_indices <- row_order(ht_mtfr)
if (is.list(row_order_indices)) row_order_indices <- unlist(row_order_indices, use.names = FALSE)
row_order_names <- rownames(mtfr_heatmap)[row_order_indices]

col_fun_expr <- colorRamp2(seq(-2, 2, length.out = 11), brewer.pal(11, "PRGn"))

file <- file.path(plot.dir, "heatmap_mofa_expr_bp.pdf")
pdf(file, width = 8, height = 10)
draw(Heatmap(
  rna_heatmap,
  name = "RNA Z-score",
  row_names_gp = gpar(fontsize = 4),
  top_annotation = ha,
  column_title = "RNA Z-scores of the MOFA TFs",
  show_row_names = TRUE,
  show_column_names = FALSE,
  column_split = column_split_factor,
  cluster_column_slices = FALSE,
  cluster_rows = FALSE,
  row_order = row_order_indices,
  col = col_fun_expr
))
dev.off()
log_info("Wrote ", file)

#####################################################################
# Row-wise correlation column
#####################################################################

row_correlation <- vapply(seq_len(nrow(mtfr_heatmap)), function(i) {
  suppressWarnings(cor(
    as.numeric(mtfr_heatmap[i, ]), as.numeric(rna_heatmap[i, ]),
    use = "complete.obs"
  ))
}, numeric(1))
row_correlation[!is.finite(row_correlation)] <- 0
correlation_matrix <- matrix(row_correlation,
  ncol = 1, dimnames = list(rownames(mtfr_heatmap), "Correlation")
)

col_fun_cor <- colorRamp2(
  seq(-1, 1, length.out = 11), rev(brewer.pal(11, "PiYG"))
)

file <- file.path(plot.dir, "correlation_heatmap.pdf")
pdf(file, width = 3, height = 10)
draw(Heatmap(
  correlation_matrix,
  name = "Row Correlation",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 4),
  col = col_fun_cor,
  row_order = row_order_indices
))
dev.off()
log_info("Wrote ", file)

#####################################################################
# Factor annotation column
#####################################################################

# One factor per feature, then aligned to the heatmap rows. The old version
# matched against mtfr_filtered, an object that was never created.
mofa_filt <- topTfs[!duplicated(topTfs$feature), ]
mofa_filt <- mofa_filt[match(rownames(mtfr_heatmap), mofa_filt$feature), ]
mofa_filt$factor[is.na(mofa_filt$factor)] <- "none"

factor_levels <- sort(unique(mofa_filt$factor))
level_col <- setNames(
  colorRampPalette(brewer.pal(max(3, min(9, length(factor_levels))), "Set1"))(length(factor_levels)),
  factor_levels
)

file <- file.path(plot.dir, "factor_column.pdf")
pdf(file, width = 2.2, height = 10)
draw(Heatmap(
  as.matrix(mofa_filt$factor),
  name = "Factors",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  row_order = row_order_indices,
  col = level_col
))
dev.off()
log_info("Wrote ", file)

#####################################################################
# SELEX annotation
#####################################################################

if (!file.exists(selex.file)) {
  log_warn("SELEX table not found, skipping the SELEX panels: ", selex.file)
} else {
  selex <- read.csv(selex.file, skip = 20, header = TRUE, sep = ";")

  # The old script used selex$Call in one place and
  # selex$methyl.SELEX.call in another, so one of the two was always NULL
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

    motifs <- rownames(mtfr_heatmap)
    idx <- match(motifs, selex[[name_col]])
    selex_call <- ifelse(is.na(idx), "Inconclusive", as.character(selex[[call_col]][idx]))
    selex_call[is.na(selex_call)] <- "Inconclusive"

    call_levels <- sort(unique(selex_call))
    call_col_map <- setNames(
      colorRampPalette(c("#CCCCCC", "#BDB76B", "#8B0000", "#008080"))(length(call_levels)),
      call_levels
    )

    file <- file.path(plot.dir, "selex_annotation_heatmap.pdf")
    pdf(file, width = 2.5, height = 10)
    draw(Heatmap(
      as.matrix(selex_call),
      name = "SELEX Call",
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      show_row_names = FALSE,
      show_column_names = FALSE,
      row_order = row_order_indices,
      col = call_col_map,
      width = unit(1, "cm")
    ))
    dev.off()
    log_info("Wrote ", file)

    ###################################################################
    # Correlation by SELEX call, over every shared feature
    ###################################################################

    # Computed on the aligned Z-score matrices, not on the deviations
    # object and the raw counts, which have neither the same rows nor the
    # same class
    all_feats <- intersect(rownames(mtfr_z), rownames(rna_z))
    cor_all <- vapply(all_feats, function(g) {
      suppressWarnings(cor(
        as.numeric(mtfr_z[g, heat_samples]),
        as.numeric(rna_z[g, heat_samples]),
        use = "complete.obs"
      ))
    }, numeric(1))

    cor_df <- data.frame(
      TF.name = all_feats, Correlation = cor_all,
      stringsAsFactors = FALSE
    )
    cor_df <- cor_df[is.finite(cor_df$Correlation), ]
    cor_df$call <- selex[[call_col]][match(cor_df$TF.name, selex[[name_col]])]
    cor_df <- cor_df[cor_df$call %in% c("MethylPlus", "MethylMinus"), ]

    if (nrow(cor_df) == 0) {
      log_warn("No feature carries a MethylPlus or MethylMinus SELEX call")
    } else {
      motif_counts <- as.data.frame(table(cor_df$call))
      colnames(motif_counts) <- c("call", "count")

      p_selex <- ggplot(cor_df, aes(x = call, y = Correlation, fill = call)) +
        stat_boxplot(geom = "errorbar", width = 0.3, linewidth = 0.8) +
        geom_boxplot(colour = "black", outlier.shape = NA, width = 0.6) +
        geom_jitter(aes(color = call), width = 0.2, size = 1.5, alpha = 0.7) +
        scale_fill_manual(values = c("MethylPlus" = "darkgreen", "MethylMinus" = "darkred")) +
        scale_color_manual(values = c("MethylPlus" = "darkgreen", "MethylMinus" = "darkred")) +
        scale_y_continuous(limits = c(-1, 1.05), breaks = seq(-1, 1, 0.2)) +
        labs(x = "SELEX Group", y = "Correlation") +
        theme_classic(base_size = 14) +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.title = element_text(face = "bold"),
          legend.position = "none"
        ) +
        geom_text(
          data = motif_counts, inherit.aes = FALSE,
          aes(x = call, y = 1.02, label = count), vjust = 0, size = 4
        )

      file <- file.path(plot.dir, "correlation_boxplot_by_selex.pdf")
      ggsave(file, plot = p_selex, width = 8, height = 6)
      log_info("Wrote ", file)
      write.csv(cor_df, file.path(table.dir, "mofa_bp_correlation_by_selex.csv"),
        row.names = FALSE
      )
    }
  }
}

#####################################################################
# Tables
#####################################################################

write.csv(agg, file.path(table.dir, "mofa_bp_modality_contribution.csv"), row.names = FALSE)
write.csv(weights_df, file.path(table.dir, "mofa_bp_weights_long.csv"), row.names = FALSE)

log_success("Finished the Blueprint MOFA integration, figures in ", plot.dir)

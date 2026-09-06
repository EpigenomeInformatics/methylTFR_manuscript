#!/usr/bin/env Rscript

#####################################################################
# 08_bp_viper_mofa.R
# created on 27-08-2026 by Irem B Gunduz
# Updated on 06-09-2026 by Irem B Gunduz
# MOFA2 integration of the methylTFR deviations and the VIPER
# transcription factor activities of the Blueprint samples
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
  library(patchwork)
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

drop.cell.types <- c("Other", "Thymocyte")

num.factors <- 14
top.factors <- 7
drop.factor1 <- TRUE

heatmap.padj.cutoff <- 0.05
heatmap.max.motifs <- 60

cell_type_labels <- c(
  "megK" = "Megakaryocytes",
  "eryt" = "Erythrocytes",
  "gran" = "Granulocytes",
  "mono" = "Monocytes",
  "Mf" = "Macrophages",
  "DC" = "Dendritic cells",
  "osteoclast" = "Osteoclast",
  "NK" = "NK-cells",
  "Tcell" = "T-cells",
  "thymocyte" = "Thymocyte",
  "Bcell" = "B-cells",
  "plasma" = "Plasma",
  "progenitor" = "Progenitors",
  "other" = "Other"
)
cell_type_levels <- unname(cell_type_labels)

cell_type_colors <- c(
  "B-cells" = "#C2377C",
  "Dendritic cells" = "#8C6D3F",
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

# Columns of the paired heatmap
heatmap.cell.types <- c(
  "Megakaryocytes", "Erythrocytes", "Granulocytes", "Monocytes",
  "Macrophages", "Dendritic cells", "NK-cells", "T-cells", "B-cells"
)
column.order <- intersect(cell_type_levels, heatmap.cell.types)

scatter.pairs <- NULL
ellipse.min.samples <- 4
ellipse.level <- 0.9

view_colors <- c("mtfr" = "#ED4B4A", "viper" = "#6BC75A")

selex.file <- "/icbb/projects/igunduz/exposure_atlas_manuscript/sample_annots/Selex_data.csv"

# Directories
analysis.dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/blueprint"
dev.tag <- "mTFR_devs_230826"
mofa.dir <- file.path(analysis.dir, "mofa_230826")

dev.file <- file.path(analysis.dir, dev.tag, paste0(motifSet, "_deviations.RDS"))
viper.file <- file.path(mofa.dir, "bp_viper_activity.RDS")
map.file <- file.path(mofa.dir, "bp_rna_wgbs_map.tsv")

github.dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR_manuscript"
table.dir <- file.path(github.dir, "tables")
diff.file <- file.path(table.dir, paste0("diff_", motifSet, "_allcelltypes.RDS"))

model.hdf5 <- file.path(mofa.dir, "mtfr_viper_model_bp.hdf5")
model.rds <- file.path(mofa.dir, "mtfr_viper_model_bp.rds")

plot.dir <- file.path(github.dir, "figures", "blueprint", "mofa_viper_230826")
if (!dir.exists(plot.dir)) dir.create(plot.dir, recursive = TRUE)
if (!dir.exists(table.dir)) dir.create(table.dir, recursive = TRUE)

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

# JASPAR writes heterodimers as FOS::JUNB and variants as JUN(var.2), so a
# motif can map to several TFs and a plain intersect would drop both forms.
motif_to_tfs <- function(motif) {
  parts <- unlist(strsplit(motif, "::", fixed = TRUE))
  parts <- sub("\\s*\\(var\\.[0-9]+\\)$", "", parts)
  toupper(trimws(parts))
}

# The VIPER row a motif is shown against: the first component TF with an
# activity, recorded so a dimer row can be traced back to its subunit.
match_motifs_to_viper <- function(motifs, viper_rows) {
  hits <- lapply(motifs, function(m) {
    tfs <- motif_to_tfs(m)
    found <- tfs[tfs %in% toupper(viper_rows)]
    if (length(found) == 0) {
      return(NULL)
    }
    viper_rows[match(found[1], toupper(viper_rows))]
  })
  names(hits) <- motifs
  keep <- !vapply(hits, is.null, logical(1))
  data.frame(
    motif = motifs[keep],
    viper_tf = unlist(hits[keep], use.names = FALSE),
    stringsAsFactors = FALSE
  )
}

#####################################################################
# Load the two modalities and the sample map
#####################################################################

for (f in c(dev.file, viper.file, map.file)) {
  if (!file.exists(f)) {
    stop("Missing input: ", f, ". Run 02, 07 and 07b first")
  }
}

dev_obj <- readRDS(dev.file)
mtfr_z <- deviationZScores(dev_obj)
viper_act <- readRDS(viper.file)
map <- read.delim(map.file, stringsAsFactors = FALSE)

log_info(
  "methylTFR: ", nrow(mtfr_z), " motifs x ", ncol(mtfr_z),
  " samples, VIPER: ", nrow(viper_act), " TFs x ", ncol(viper_act), " samples"
)

map <- map[map$bedFile %in% colnames(mtfr_z) & map$rna_id %in% colnames(viper_act), ]
if (nrow(map) == 0) {
  stop(
    "The sample map matched nothing. The deviations are keyed on '",
    colnames(mtfr_z)[1], "', the activities on '", colnames(viper_act)[1], "'."
  )
}
log_info(nrow(map), " samples shared by the deviations and the activities")

mtfr_mat <- mtfr_z[, map$bedFile, drop = FALSE]
viper_mat <- viper_act[, map$rna_id, drop = FALSE]

sample_ids <- map$rna_id
colnames(mtfr_mat) <- sample_ids
colnames(viper_mat) <- sample_ids

# VIPER scores are centred over the samples of their own run, so both views
# are Z-scored here to put them on a comparable scale.
mtfr_view <- row_zscore(mtfr_mat)
viper_view <- row_zscore(viper_mat)
stopifnot(identical(colnames(mtfr_view), colnames(viper_view)))

metadata <- data.frame(
  sample = sample_ids,
  celltype = unname(cell_type_labels[map$cellTypeGroup]),
  stringsAsFactors = FALSE
)
if (any(is.na(metadata$celltype))) {
  stop(
    "No display name for cellTypeGroup value(s): ",
    paste(unique(map$cellTypeGroup[is.na(metadata$celltype)]), collapse = ", ")
  )
}

#####################################################################
# Train
#####################################################################

data_list <- list(mtfr = as.matrix(mtfr_view), viper = as.matrix(viper_view))

if (file.exists(model.rds)) {
  log_info("Loading the trained model from ", model.rds)
  MOFAobject.trained <- readRDS(model.rds)
} else {
  MOFAobject <- create_mofa(data_list)
  samples_metadata(MOFAobject) <- metadata

  model_opts <- get_default_model_options(MOFAobject)
  model_opts$num_factors <- num.factors

  MOFAobject <- prepare_mofa(MOFAobject, model_options = model_opts)
  MOFAobject.trained <- run_mofa(MOFAobject, model.hdf5, use_basilisk = FALSE)
  saveRDS(MOFAobject.trained, model.rds)
  log_success("Wrote ", model.rds)
}

keep_samples <- metadata$sample[!metadata$celltype %in% drop.cell.types]
keep_samples <- intersect(keep_samples, samples_names(MOFAobject.trained)[[1]])
log_info("Keeping ", length(keep_samples), " of ", nrow(metadata),
  " samples after dropping ", paste(drop.cell.types, collapse = ", "))
MOFAobject.trained <- subset_samples(MOFAobject.trained, keep_samples)
metadata <- metadata[metadata$sample %in% keep_samples, ]

#####################################################################
# Variance in each factor that cell type accounts for
#####################################################################

factors_long <- get_factors(MOFAobject.trained, as.data.frame = TRUE)
factors_long$sample <- as.character(factors_long$sample)
factors_long <- left_join(factors_long, metadata, by = "sample")
if (any(is.na(factors_long$celltype))) {
  stop("Some factor rows could not be matched to a cell type")
}

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
  stringsAsFactors = FALSE, row.names = NULL
)
if (all(is.na(r2_df$R2))) stop("No factor could be tested against cell type")
r2_df <- r2_df[order(-r2_df$R2), ]
if (drop.factor1) r2_df <- r2_df[r2_df$Factor != "Factor1", ]

write.csv(r2_df, file.path(table.dir, "mofa_bp_viper_factor_celltype_R2.csv"),
  row.names = FALSE
)
top_factors <- head(r2_df$Factor, top.factors)
log_info("Top factors: ", paste(top_factors, collapse = ", "))

#####################################################################
# Weights, and the factor each motif loads on
#####################################################################

weights_df <- get_weights(MOFAobject.trained, as.data.frame = TRUE) %>%
  mutate(
    feature = as.character(feature),
    factor = as.character(factor),
    view = as.character(view),
    value_abs = abs(value)
  )

mtfr_top_factor <- weights_df %>%
  filter(view == "mtfr") %>%
  group_by(feature) %>%
  slice_max(value_abs, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(feature, factor)

#####################################################################
# A) Modality contribution per factor, with the R2 over it
#####################################################################

agg <- weights_df %>%
  group_by(factor, view) %>%
  summarise(sum_abs = sum(value_abs, na.rm = TRUE), .groups = "drop") %>%
  group_by(factor) %>%
  mutate(frac = sum_abs / sum(sum_abs)) %>%
  ungroup() %>%
  filter(factor %in% r2_df$Factor)

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
  labs(x = "Factor", y = "Relative contribution (|weights| fraction)", fill = "View") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.y.right = element_text(colour = "purple")
  )

ggsave(file.path(plot.dir, "viper_modality_contribution_fraction_R2.pdf"),
  plot = p_modality_frac, width = 11, height = 4.5
)

#####################################################################
# B) Factor values per sample
#####################################################################

factors_top <- factors_long %>%
  filter(factor %in% top_factors, !celltype %in% drop.cell.types)
factors_top$factor <- factor(factors_top$factor, levels = top_factors)
factors_top$celltype <- factor(factors_top$celltype, levels = cell_type_levels)

p_strip <- ggplot(factors_top, aes(x = factor, y = value, color = celltype)) +
  geom_jitter(width = 0.2, height = 0, size = 2, alpha = 0.8) +
  scale_color_manual(values = cell_type_colors, name = "Cell type", drop = TRUE) +
  labs(x = "Factor", y = "Factor value") +
  theme_classic(base_size = 13) +
  theme(panel.grid.major.x = element_blank())

ggsave(file.path(plot.dir, "viper_factors_stripplot_by_celltype.pdf"),
  plot = p_strip, width = 9, height = 5
)

#####################################################################
# C) Factor scatters
#####################################################################

factors_wide <- factors_long %>%
  select(sample, factor, value, celltype) %>%
  pivot_wider(names_from = factor, values_from = value) %>%
  mutate(celltype = factor(celltype, levels = cell_type_levels)) %>%
  as.data.frame()

available <- setdiff(colnames(factors_wide), c("sample", "celltype"))
pairs <- scatter.pairs
if (is.null(pairs)) {
  pairs <- list(top_factors[1:2], top_factors[3:4])
  pairs <- Filter(function(p) all(!is.na(p)), pairs)
}

scatter_plot <- function(fx, fy) {
  counts <- table(droplevels(factors_wide$celltype))
  ell <- factors_wide[factors_wide$celltype %in%
    names(counts)[counts >= ellipse.min.samples], , drop = FALSE]

  ggplot(factors_wide, aes(x = .data[[fx]], y = .data[[fy]], colour = celltype)) +
    geom_point(alpha = 0.9, size = 2) +
    stat_ellipse(
      data = ell, aes(group = celltype),
      type = "norm", level = ellipse.level,
      linetype = "dashed", linewidth = 0.4, show.legend = FALSE
    ) +
    scale_colour_manual(values = cell_type_colors, name = "Cell type", drop = TRUE) +
    labs(x = fx, y = fy) +
    theme_classic(base_size = 13)
}

scatter_panels <- list()
for (pr in pairs) {
  if (length(setdiff(pr, available)) > 0) {
    log_warn("Skipping the scatter for ", paste(pr, collapse = " vs "),
      ", not among the fitted factors")
    next
  }
  p_scatter <- scatter_plot(pr[1], pr[2])
  scatter_panels[[paste(pr, collapse = "_")]] <- p_scatter
  file <- file.path(plot.dir, paste0("factor_scatter_", pr[1], "_", pr[2], ".pdf"))
  ggsave(file, p_scatter, width = 7, height = 5)
  log_info("Wrote ", file)
}

#####################################################################
# Rows of the paired heatmap
#####################################################################

if (!file.exists(diff.file)) {
  stop("Differential results not found, run 06 first: ", diff.file)
}
diff <- readRDS(diff.file)
motif_col <- intersect(c("motifs", "motif", "feature"), colnames(diff))
if (length(motif_col) == 0) {
  stop(
    "No motif column in ", basename(diff.file),
    ". Present: ", paste(colnames(diff), collapse = ", ")
  )
}
diff$motif <- as.character(diff[[motif_col[1]]])

sig <- diff[which(diff$p_value_adjusted < heatmap.padj.cutoff), ]
log_info(nrow(sig), " motifs differential at adjusted p < ", heatmap.padj.cutoff)
if (nrow(sig) == 0) stop("No differential motif at the chosen cutoff")

sig <- sig[sig$motif %in% rownames(mtfr_view), ]
pairs_tab <- match_motifs_to_viper(sig$motif, rownames(viper_view))
log_info(
  nrow(pairs_tab), " of ", nrow(sig),
  " differential motifs have a VIPER activity"
)

# Rows are restricted to the factors that carry the cell type variance
pairs_tab$factor <- mtfr_top_factor$factor[match(pairs_tab$motif, mtfr_top_factor$feature)]
pairs_tab <- pairs_tab[which(pairs_tab$factor %in% top_factors), ]
log_info(nrow(pairs_tab), " of them load on ", paste(top_factors, collapse = ", "))
if (nrow(pairs_tab) < 2) {
  stop("Fewer than two differential motifs load on the top factors")
}

pairs_tab <- pairs_tab[order(match(pairs_tab$motif, sig$motif)), ]
if (nrow(pairs_tab) > heatmap.max.motifs) {
  log_info("Showing the ", heatmap.max.motifs, " most significant of them")
  pairs_tab <- head(pairs_tab, heatmap.max.motifs)
}

write.csv(pairs_tab, file.path(table.dir, "bp_motif_viper_pairs.csv"), row.names = FALSE)

n_dimer <- sum(grepl("::", pairs_tab$motif, fixed = TRUE))
if (n_dimer > 0) {
  log_info(n_dimer, " heterodimer motif(s) matched on one subunit, see bp_motif_viper_pairs.csv")
}

#####################################################################
# D) Paired heatmap, methylTFR left, VIPER right
#####################################################################

heat_samples <- metadata$sample[metadata$celltype %in% heatmap.cell.types]
heat_samples <- intersect(heat_samples, colnames(mtfr_view))
heat_groups <- metadata$celltype[match(heat_samples, metadata$sample)]
log_info(
  length(heat_samples), " samples in the heatmap across ",
  length(unique(heat_groups)), " cell types"
)

mtfr_heatmap <- mtfr_view[pairs_tab$motif, heat_samples, drop = FALSE]
viper_heatmap <- viper_view[pairs_tab$viper_tf, heat_samples, drop = FALSE]
rownames(viper_heatmap) <- pairs_tab$motif

present <- intersect(names(cell_type_colors), unique(heat_groups))
column_split_factor <- factor(heat_groups, levels = intersect(column.order, unique(heat_groups)))

cell_annotation <- function() {
  HeatmapAnnotation(
    celltypes = heat_groups,
    col = list(celltypes = cell_type_colors[present]),
    annotation_name_gp = gpar(fontsize = 8),
    show_legend = FALSE
  )
}

col_fun_mtfr <- colorRamp2(
  c(-2, -1, 0, 1, 2),
  c("#2B4B9B", "#4EC3E0", "#FFF200", "#F5A623", "#B21212")
)
col_fun_viper <- colorRamp2(seq(-2, 2, length.out = 11), brewer.pal(11, "PRGn"))
col_fun_cor <- colorRamp2(seq(-1, 1, length.out = 11), rev(brewer.pal(11, "PiYG")))

ht_mtfr <- Heatmap(
  mtfr_heatmap,
  name = "mTFR\nZ-score",
  col = col_fun_mtfr,
  top_annotation = cell_annotation(),
  column_split = column_split_factor,
  cluster_column_slices = FALSE,
  column_gap = unit(0.8, "mm"),
  show_row_names = FALSE,
  show_column_names = FALSE,
  row_dend_side = "left",
  column_title = "methylTFR",
  column_title_gp = gpar(fontsize = 8)
)

ht_viper <- Heatmap(
  viper_heatmap,
  name = "VIPER\nZ-score",
  col = col_fun_viper,
  top_annotation = cell_annotation(),
  column_split = column_split_factor,
  cluster_column_slices = FALSE,
  column_gap = unit(0.8, "mm"),
  cluster_rows = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  column_title = "VIPER",
  column_title_gp = gpar(fontsize = 8)
)

row_correlation <- vapply(seq_len(nrow(mtfr_heatmap)), function(i) {
  suppressWarnings(cor(
    as.numeric(mtfr_heatmap[i, ]), as.numeric(viper_heatmap[i, ]),
    use = "complete.obs"
  ))
}, numeric(1))
row_correlation[!is.finite(row_correlation)] <- 0

# A motif with no methyl-SELEX entry was never assayed, so it reads Inconclusive
selex_call <- NULL
if (!file.exists(selex.file)) {
  log_warn("SELEX table not found, the SELEX column is skipped: ", selex.file)
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
    idx <- match(rownames(mtfr_heatmap), selex_tab[[name_col[1]]])
    selex_call <- ifelse(is.na(idx), "Inconclusive",
      as.character(selex_tab[[call_col[1]]][idx])
    )
    selex_call[is.na(selex_call)] <- "Inconclusive"
    log_info(
      "SELEX calls on the heatmap: ",
      paste(names(table(selex_call)), table(selex_call), sep = " = ", collapse = ", ")
    )
  }
}

selex_col_map <- NULL
if (!is.null(selex_call)) {
  known <- c(
    "Inconclusive" = "#D9D9D9",
    "Little effect" = "#BDB76B",
    "MethylMinus" = "#8B0000",
    "MethylPlus" = "#008080"
  )
  extra <- setdiff(unique(selex_call), names(known))
  if (length(extra) > 0) {
    log_warn("Unlisted SELEX call(s), coloured grey: ", paste(extra, collapse = ", "))
    known <- c(known, setNames(rep("#8C8C8C", length(extra)), extra))
  }
  selex_col_map <- known[intersect(names(known), unique(selex_call))]
}

motif_factor <- pairs_tab$factor[match(rownames(mtfr_heatmap), pairs_tab$motif)]
factor_levels_ann <- top_factors[top_factors %in% motif_factor]
factor_col_map <- setNames(
  colorRampPalette(brewer.pal(max(3, min(9, length(factor_levels_ann))), "Set1"))(
    length(factor_levels_ann)
  ),
  factor_levels_ann
)

annotation_args <- list(
  Factors = factor(motif_factor, levels = factor_levels_ann),
  `Row-wise Pearson Correlation` = row_correlation,
  col = list(
    Factors = factor_col_map,
    `Row-wise Pearson Correlation` = col_fun_cor
  ),
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
  rownames(mtfr_heatmap),
  gp = gpar(fontsize = 6), which = "row"
)

ha_right <- do.call(rowAnnotation, annotation_args)

panel_d <- grid.grabExpr(
  draw(ht_mtfr + ht_viper + ha_right,
    merge_legend = TRUE, heatmap_legend_side = "right"
  )
)

file <- file.path(plot.dir, "heatmap_mtfr_viper_differential_bp.pdf")
pdf(file, width = 15, height = 11)
grid.draw(panel_d)
dev.off()
log_info("Wrote ", file)

write.csv(
  data.frame(
    motif = rownames(mtfr_heatmap),
    viper_tf = pairs_tab$viper_tf,
    factor = motif_factor,
    correlation = row_correlation,
    stringsAsFactors = FALSE
  ),
  file.path(table.dir, "bp_mtfr_viper_row_correlation.csv"),
  row.names = FALSE
)

#####################################################################
# The assembled figure
#####################################################################

tag_theme <- theme(plot.tag = element_text(face = "bold", size = 16))

panel_a <- p_modality_frac + labs(tag = "A")
panel_b <- p_strip + labs(tag = "B")

if (length(scatter_panels) == 0) {
  stop("No factor scatter could be drawn, the assembled figure needs panel C")
}
scatter_panels[[1]] <- scatter_panels[[1]] + labs(tag = "C")
panel_c <- wrap_plots(
  lapply(scatter_panels, function(p) p + theme(legend.position = "none")),
  nrow = 1
)

panel_d_wrapped <- wrap_elements(full = panel_d) + labs(tag = "D")

combined <- wrap_plots(
  panel_a, panel_b, panel_c, panel_d_wrapped,
  ncol = 1, heights = c(1, 1.1, 1.2, 2.8)
) & tag_theme

file <- file.path(plot.dir, "figure_bp_viper_mofa_integration.pdf")
ggsave(file, combined, width = 16, height = 24, limitsize = FALSE)
log_info("Wrote ", file)

log_success("Finished the VIPER integration, figures in ", plot.dir)
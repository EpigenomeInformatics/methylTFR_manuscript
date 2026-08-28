#!/usr/bin/env Rscript

#####################################################################
# 09_bp_viper_mofa.R
# created on 27-08-2026 by Irem B Gunduz
# MOFA2 integration of the methylTFR deviations and the VIPER
# transcription factor activities of the Blueprint samples, and the
# paired methylTFR / VIPER heatmap of the differential motifs
#
# The activity matrix comes from 07b, the sample map from 07 and the
# differential motifs from 06. Nothing here recomputes them.
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

drop.cell.types <- c("Other", "Thymocyte")

num.factors <- 14
top.factors <- 7
top.tfs <- 5
drop.factor1 <- TRUE

# Heatmap rows: motifs called differential by 06 that also have a VIPER
# activity. Set to FALSE to fall back to the top MOFA features, which is
# what 08 shows.
heatmap.differential.only <- TRUE
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
  "NK" = "NK",
  "Tcell" = "T-cells",
  "thymocyte" = "Thymocyte",
  "Bcell" = "B-cells",
  "plasma" = "Plasma",
  "progenitor" = "Progenitors",
  "other" = "Other"
)
cell_type_levels <- unname(cell_type_labels)

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

column.order <- setdiff(cell_type_levels, c("Other", "Progenitors"))

scatter.pairs <- NULL

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

# A motif name is not a gene symbol. JASPAR writes heterodimers as
# FOS::JUNB and variants as JUN(var.2), so a motif can correspond to one
# or several TFs, and a plain intersect against the VIPER rows silently
# drops every dimer and every variant.
motif_to_tfs <- function(motif) {
  parts <- unlist(strsplit(motif, "::", fixed = TRUE))
  parts <- sub("\\s*\\(var\\.[0-9]+\\)$", "", parts)
  toupper(trimws(parts))
}

# The VIPER row a motif is shown against: the first of its component TFs
# that has an activity. Which one was used is recorded, so a dimer row can
# be traced back to the subunit it was matched on.
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

# The map carries the WGBS identifier and the RNA identifier of one pair
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

# One identifier on both sides, so MOFA sees a single sample set
sample_ids <- map$rna_id
colnames(mtfr_mat) <- sample_ids
colnames(viper_mat) <- sample_ids

# VIPER already returns a normalised enrichment score, but it is centred
# per TF over the samples that were in the run, not over these. Both views
# are Z-scored here so their scales are comparable.
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
# Factor scatter
#####################################################################

factors_wide <- factors_long %>%
  select(sample, factor, value, celltype) %>%
  pivot_wider(names_from = factor, values_from = value)

available <- setdiff(colnames(factors_wide), c("sample", "celltype"))
pairs <- scatter.pairs
if (is.null(pairs)) {
  pairs <- list(head(r2_df$Factor, 2))
}

scatter_plot <- function(fx, fy) {
  ggplot(factors_wide, aes(x = .data[[fx]], y = .data[[fy]], color = celltype)) +
    geom_point(alpha = 0.9, size = 2) +
    scale_color_manual(values = cell_type_colors, name = "Cell type") +
    labs(x = fx, y = fy) +
    theme_classic(base_size = 13) +
    theme(legend.position = "right")
}

scatter_panels <- list()
for (pr in pairs) {
  missing_f <- setdiff(pr, available)
  if (length(missing_f) > 0) {
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
# Modality contribution per factor, with the R2 over it
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
# Factor values per sample
#####################################################################

factors_top <- factors_long %>%
  filter(factor %in% top_factors, !celltype %in% drop.cell.types)
factors_top$factor <- factor(factors_top$factor, levels = top_factors)
factors_top$celltype <- factor(factors_top$celltype, levels = cell_type_levels)

p_strip <- ggplot(factors_top, aes(x = factor, y = value, color = celltype)) +
  geom_jitter(width = 0.2, height = 0, size = 2, alpha = 0.8) +
  scale_color_manual(values = cell_type_colors, name = "Cell type") +
  labs(x = "Factor", y = "Factor value") +
  theme_classic(base_size = 13) +
  theme(panel.grid.major.x = element_blank())

ggsave(file.path(plot.dir, "viper_factors_stripplot_by_celltype.pdf"),
  plot = p_strip, width = 9, height = 5
)

#####################################################################
# Rows of the paired heatmap
#####################################################################

viper_rows <- rownames(viper_view)

if (heatmap.differential.only) {
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

  # Differential, present in the deviations, and matched to a VIPER TF
  sig <- sig[sig$motif %in% rownames(mtfr_view), ]
  pairs_tab <- match_motifs_to_viper(sig$motif, viper_rows)
  log_info(
    nrow(pairs_tab), " of ", nrow(sig),
    " differential motifs have a VIPER activity"
  )
  if (nrow(pairs_tab) < 2) {
    stop("Fewer than two differential motifs overlap the VIPER TFs")
  }

  # Most significant first, then capped so the row labels stay readable
  pairs_tab <- pairs_tab[order(match(pairs_tab$motif, sig$motif)), ]
  if (nrow(pairs_tab) > heatmap.max.motifs) {
    log_info("Showing the ", heatmap.max.motifs, " most significant of them")
    pairs_tab <- head(pairs_tab, heatmap.max.motifs)
  }
} else {
  topTfs <- weights_df %>%
    filter(factor %in% top_factors) %>%
    group_by(factor, view) %>%
    slice_max(value_abs, n = top.tfs, with_ties = FALSE) %>%
    ungroup()
  pairs_tab <- match_motifs_to_viper(
    intersect(unique(topTfs$feature), rownames(mtfr_view)), viper_rows
  )
  if (nrow(pairs_tab) < 2) stop("Fewer than two MOFA features overlap the VIPER TFs")
}

write.csv(pairs_tab, file.path(table.dir, "bp_motif_viper_pairs.csv"), row.names = FALSE)

n_dimer <- sum(grepl("::", pairs_tab$motif, fixed = TRUE))
if (n_dimer > 0) {
  log_info(n_dimer, " heterodimer motif(s) matched on one subunit, see bp_motif_viper_pairs.csv")
}

#####################################################################
# Paired heatmap, methylTFR left, VIPER right
#####################################################################

heat_samples <- intersect(metadata$sample, colnames(mtfr_view))
heat_groups <- metadata$celltype[match(heat_samples, metadata$sample)]
keep <- !heat_groups %in% drop.cell.types
heat_samples <- heat_samples[keep]
heat_groups <- heat_groups[keep]

mtfr_heatmap <- mtfr_view[pairs_tab$motif, heat_samples, drop = FALSE]
viper_heatmap <- viper_view[pairs_tab$viper_tf, heat_samples, drop = FALSE]

# Rows are labelled by the motif on both halves, so the two line up by eye
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
  cluster_rows = FALSE, # the row order comes from the methylTFR half
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

# The methyl-SELEX call of each motif, as in 08. A motif with no entry was
# never assayed rather than missing, so it is labelled Inconclusive.
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

# The factor each motif carries its largest weight on.
# get_weights returns feature and factor as factors, and indexing a named
# vector with a factor uses its integer codes rather than its labels, so
# both are coerced to character before anything is looked up by name.
motif_factor <- rep("none", nrow(mtfr_heatmap))
names(motif_factor) <- rownames(mtfr_heatmap)

weight_top <- weights_df %>%
  mutate(
    feature = as.character(feature),
    factor = as.character(factor),
    view = as.character(view)
  ) %>%
  filter(view == "mtfr", feature %in% rownames(mtfr_heatmap)) %>%
  group_by(feature) %>%
  slice_max(value_abs, n = 1, with_ties = FALSE) %>%
  ungroup()

if (nrow(weight_top) == 0) {
  log_warn(
    "No mtfr weight matched a heatmap row. MOFA features look like '",
    as.character(weights_df$feature[1]), "', the heatmap rows like '",
    rownames(mtfr_heatmap)[1], "'. The Factors column will read none."
  )
} else {
  motif_factor[weight_top$feature] <- weight_top$factor
}

stopifnot(length(motif_factor) == nrow(mtfr_heatmap))

factor_levels_ann <- sort(unique(motif_factor))
factor_col_map <- setNames(
  colorRampPalette(brewer.pal(max(3, min(9, length(factor_levels_ann))), "Set1"))(
    length(factor_levels_ann)
  ),
  factor_levels_ann
)

annotation_args <- list(
  Factors = unname(motif_factor),
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
# which = "row" is passed explicitly: rowAnnotation normally infers it from
# its own call, but do.call evaluates anno_text() before that happens
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
    correlation = row_correlation,
    stringsAsFactors = FALSE
  ),
  file.path(table.dir, "bp_mtfr_viper_row_correlation.csv"),
  row.names = FALSE
)

#####################################################################
# The assembled figure
#####################################################################

# A the modality contribution, B the factor values, C the factor scatters,
# D the paired heatmap, which travels as the grob captured above
if (!requireNamespace("patchwork", quietly = TRUE)) {
  log_warn("patchwork is not installed, the assembled figure is skipped")
} else {
  panel_c <- if (length(scatter_panels) == 0) {
    NULL
  } else if (length(scatter_panels) == 1) {
    scatter_panels[[1]]
  } else {
    patchwork::wrap_plots(scatter_panels, nrow = 1, guides = "collect") &
      ggplot2::theme(legend.position = "none")
  }

  mid <- if (is.null(panel_c)) {
    p_strip
  } else {
    patchwork::wrap_plots(p_strip, panel_c, ncol = 2, widths = c(1, 1))
  }

  combined <- patchwork::wrap_plots(
    p_modality_frac, mid, patchwork::wrap_elements(full = panel_d),
    ncol = 1, heights = c(1, 1.1, 2.4)
  ) + patchwork::plot_annotation(tag_levels = "A")

  file <- file.path(plot.dir, "figure_bp_viper_mofa_integration.pdf")
  ggsave(file, combined, width = 16, height = 22, limitsize = FALSE)
  log_info("Wrote ", file)
}

log_success("Finished the VIPER integration, figures in ", plot.dir)

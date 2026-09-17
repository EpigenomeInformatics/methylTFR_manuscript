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
top.factors <- 5
drop.factor1 <- TRUE

heatmap.padj.cutoff <- 0.05
heatmap.tfs.per.factor <- 10

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
selex.groups <- c("MethylMinus", "MethylPlus")
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

# Directories
analysis.dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/blueprint"
dev.tag <- "mTFR_devs_230826"
mofa.dir <- file.path(analysis.dir, "mofa_230826")

dev.file <- file.path(analysis.dir, dev.tag, paste0(motifSet, "_deviations.RDS"))
viper.file <- file.path(mofa.dir, "bp_viper_activity.RDS")
map.file <- file.path(mofa.dir, "bp_rna_wgbs_map.tsv")

github.dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR_manuscript"
table.dir <- file.path(github.dir, "tables")
mixed.dir <- file.path(table.dir, "mixed")
diff.file <- file.path(table.dir, paste0("diff_", motifSet, "_allcelltypes.RDS"))

model.hdf5 <- file.path(mofa.dir, "mtfr_viper_model_bp.hdf5")
model.rds <- file.path(mofa.dir, "mtfr_viper_model_bp.rds")

plot.dir <- file.path(github.dir, "figures", "blueprint", "mofa_viper_230826")
if (!dir.exists(plot.dir)) dir.create(plot.dir, recursive = TRUE)
if (!dir.exists(mixed.dir)) dir.create(mixed.dir, recursive = TRUE)

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

write.csv(r2_df, file.path(mixed.dir, "mofa_bp_viper_factor_celltype_R2.csv"),
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
  select(feature, factor, value_abs)

#####################################################################
# A) Modality contribution per factor, with the R2 over it
#####################################################################

agg <- weights_df %>%
  group_by(factor, view) %>%
  summarise(sum_abs = sum(value_abs, na.rm = TRUE), .groups = "drop") %>%
  group_by(factor) %>%
  mutate(frac = sum_abs / sum(sum_abs)) %>%
  ungroup() %>%
  filter(factor %in% top_factors)

agg$factor <- factor(agg$factor, levels = top_factors)
r2_line <- r2_df[r2_df$Factor %in% top_factors, ]
r2_line$Factor <- factor(r2_line$Factor, levels = top_factors)

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
  plot = p_modality_frac & no_bg, width = 11, height = 4.5, bg = "transparent"
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
  plot = p_strip & no_bg, width = 9, height = 5, bg = "transparent"
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
  ggsave(file, p_scatter & no_bg, width = 7, height = 5, bg = "transparent")
  log_info("Wrote ", file)
}

#####################################################################
# methyl-SELEX calls
#####################################################################

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
    selex_calls <- setNames(
      as.character(selex_tab[[call_col[1]]]),
      toupper(as.character(selex_tab[[name_col[1]]]))
    )
  }
}

# A motif with no assayed subunit was never tested rather than missing,
# so it reads Inconclusive
selex_call_of <- function(motifs) {
  if (is.null(selex_calls)) {
    return(NULL)
  }
  vapply(motifs, function(m) {
    hit <- motif_to_tfs(m)
    hit <- hit[hit %in% names(selex_calls)]
    if (length(hit) == 0) {
      return("Inconclusive")
    }
    call <- selex_calls[[hit[1]]]
    if (is.na(call) || !nzchar(call)) "Inconclusive" else call
  }, character(1), USE.NAMES = FALSE)
}

#####################################################################
# Every motif that pairs with a VIPER activity
#####################################################################

viper_pairs <- match_motifs_to_viper(rownames(mtfr_view), rownames(viper_view))
log_info(
  nrow(viper_pairs), " motifs pair with a VIPER activity, ",
  length(unique(viper_pairs$viper_tf)), " distinct TFs"
)
if (nrow(viper_pairs) < 2) stop("Fewer than two motifs pair with a VIPER activity")

sig_motifs <- character(0)
if (!file.exists(diff.file)) {
  log_warn("Differential results not found, the differential flag stays empty: ", diff.file)
} else {
  diff <- readRDS(diff.file)
  motif_col <- intersect(c("motifs", "motif", "feature"), colnames(diff))
  if (length(motif_col) == 0) {
    log_warn(
      "No motif column in ", basename(diff.file),
      ". Present: ", paste(colnames(diff), collapse = ", ")
    )
  } else {
    diff$motif <- as.character(diff[[motif_col[1]]])
    sig_motifs <- diff$motif[which(diff$p_value_adjusted < heatmap.padj.cutoff)]
    log_info(length(sig_motifs), " motifs differential at adjusted p < ", heatmap.padj.cutoff)
  }
}

#####################################################################
# Rows of the paired heatmap
#####################################################################

candidates <- weights_df %>%
  filter(view == "mtfr", factor %in% top_factors, feature %in% viper_pairs$motif) %>%
  left_join(viper_pairs, by = c("feature" = "motif"))

rank_on <- function(f, n, skip_motifs, skip_tfs) {
  candidates %>%
    filter(factor == f, !feature %in% skip_motifs, !viper_tf %in% skip_tfs) %>%
    arrange(desc(value_abs)) %>%
    head(n)
}

taken_motifs <- character(0)
taken_tfs <- character(0)
selected <- vector("list", length(top_factors))
names(selected) <- top_factors

# One row per factor per round, so a scarce pool of VIPER TFs is shared out
# evenly instead of being used up by the factors that come first
for (i in seq_len(heatmap.tfs.per.factor)) {
  for (f in top_factors) {
    pick <- rank_on(f, 1, taken_motifs, taken_tfs)
    if (nrow(pick) == 0) next
    selected[[f]] <- bind_rows(selected[[f]], pick)
    taken_motifs <- c(taken_motifs, pick$feature)
    taken_tfs <- c(taken_tfs, pick$viper_tf)
  }
}

# A factor still short of its quota reuses TFs another factor already carries
for (f in top_factors) {
  short <- heatmap.tfs.per.factor - NROW(selected[[f]])
  if (short <= 0) next
  extra <- rank_on(f, short, taken_motifs, character(0))
  log_info(f, " reused ", nrow(extra), " TF(s) already shown on another factor")
  selected[[f]] <- bind_rows(selected[[f]], extra)
  taken_motifs <- c(taken_motifs, extra$feature)
  if (NROW(selected[[f]]) < heatmap.tfs.per.factor) {
    log_warn(f, " filled only ", NROW(selected[[f]]), " of ", heatmap.tfs.per.factor, " rows")
  }
}

selected <- bind_rows(selected[top_factors]) %>%
  arrange(match(factor, top_factors), desc(value_abs))

log_info(
  "Heatmap motifs per factor: ",
  paste(names(table(selected$factor)), table(selected$factor),
    sep = " = ", collapse = ", "
  )
)
if (nrow(selected) < 2) stop("Fewer than two motifs could be selected for the heatmap")

pairs_tab <- data.frame(
  motif = selected$feature,
  viper_tf = selected$viper_tf,
  factor = selected$factor,
  stringsAsFactors = FALSE
)

n_dimer <- sum(grepl("::", pairs_tab$motif, fixed = TRUE))
if (n_dimer > 0) {
  log_info(n_dimer, " heterodimer motif(s) matched on one subunit")
}

#####################################################################
# Samples and splits shared by the paired panels
#####################################################################

heat_samples <- metadata$sample[metadata$celltype %in% heatmap.cell.types]
heat_samples <- intersect(heat_samples, colnames(mtfr_view))
heat_groups <- metadata$celltype[match(heat_samples, metadata$sample)]
log_info(
  length(heat_samples), " samples in the paired panels across ",
  length(unique(heat_groups)), " cell types"
)

present <- intersect(names(cell_type_colors), unique(heat_groups))
column_split_factor <- factor(heat_groups, levels = intersect(column.order, unique(heat_groups)))
row_split_factor <- factor(pairs_tab$factor, levels = top_factors[top_factors %in% pairs_tab$factor])

#####################################################################
# The integration table: every motif / VIPER pair, over every factor
#####################################################################

cor_tab <- viper_pairs
cor_tab$top_factor <- mtfr_top_factor$factor[match(cor_tab$motif, mtfr_top_factor$feature)]
cor_tab$weight_abs <- mtfr_top_factor$value_abs[match(cor_tab$motif, mtfr_top_factor$feature)]
cor_tab$correlation <- vapply(seq_len(nrow(cor_tab)), function(i) {
  suppressWarnings(cor(
    as.numeric(mtfr_view[cor_tab$motif[i], heat_samples]),
    as.numeric(viper_view[cor_tab$viper_tf[i], heat_samples]),
    method = "pearson", use = "complete.obs"
  ))
}, numeric(1))
cor_tab <- cor_tab[is.finite(cor_tab$correlation), ]

cor_tab$selex_call <- if (is.null(selex_calls)) {
  NA_character_
} else {
  selex_call_of(cor_tab$motif)
}
cor_tab$differential <- cor_tab$motif %in% sig_motifs
cor_tab$in_heatmap <- cor_tab$motif %in% pairs_tab$motif
cor_tab <- cor_tab[order(-cor_tab$weight_abs), ]

write.csv(cor_tab, file.path(mixed.dir, "bp_mtfr_viper_integration.csv"),
  row.names = FALSE
)
log_info(
  nrow(cor_tab), " motif / VIPER pairs in the integration table, ",
  sum(cor_tab$in_heatmap), " of them in the heatmap"
)
if (!is.null(selex_calls)) {
  log_info(
    "SELEX calls over all of them: ",
    paste(names(table(cor_tab$selex_call)), table(cor_tab$selex_call),
      sep = " = ", collapse = ", "
    )
  )
}

#####################################################################
# D) mTFR / VIPER correlation by methyl-SELEX call
#####################################################################

p_selex <- NULL
cor_sel <- cor_tab[which(cor_tab$selex_call %in% selex.groups), ]
if (nrow(cor_sel) == 0) {
  log_warn("No motif carries one of ", paste(selex.groups, collapse = " or "))
} else {
  cor_sel$selex_call <- factor(cor_sel$selex_call, levels = selex.groups)
  selex_counts <- as.data.frame(table(cor_sel$selex_call))
  colnames(selex_counts) <- c("selex_call", "count")
  log_info("SELEX groups in the boxplot: ",
    paste(selex_counts$selex_call, selex_counts$count, sep = " = ", collapse = ", "))

  # Two sided Wilcoxon rank sum on the correlations, one test per pair of
  # SELEX groups, Benjamini-Hochberg over those tests
  groups_present <- levels(droplevels(cor_sel$selex_call))
  tests <- NULL
  if (length(groups_present) >= 2) {
    tests <- bind_rows(lapply(combn(groups_present, 2, simplify = FALSE), function(g) {
      x <- cor_sel$correlation[cor_sel$selex_call == g[1]]
      y <- cor_sel$correlation[cor_sel$selex_call == g[2]]
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

    write.csv(tests, file.path(mixed.dir, "bp_selex_correlation_test.csv"),
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

  p_selex <- ggplot(cor_sel, aes(x = selex_call, y = correlation, fill = selex_call)) +
    stat_boxplot(geom = "errorbar", width = 0.3, linewidth = 0.6) +
    geom_boxplot(colour = "black", outlier.shape = NA, width = 0.6) +
    geom_jitter(aes(color = selex_call), width = 0.2, size = 1.2, alpha = 0.7) +
    scale_fill_manual(values = selex_colors) +
    scale_color_manual(values = selex_colors) +
    scale_y_continuous(breaks = seq(-1, 1, 0.2)) +
    coord_cartesian(ylim = c(-1, 1.12), clip = "off") +
    geom_text(
      data = selex_counts, inherit.aes = FALSE,
      aes(x = selex_call, y = 1.06, label = count), vjust = 0, size = 4.5
    ) +
    labs(x = "SELEX Group", y = "Correlation") +
    theme_classic(base_size = 13) +
    theme(legend.position = "none", axis.title = element_text(face = "plain"))

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

  file <- file.path(plot.dir, "correlation_boxplot_by_selex_bp.pdf")
  ggsave(file, p_selex & no_bg, width = 5.5, height = 6, bg = "transparent")
  log_info("Wrote ", file)
}

#####################################################################
# E) Paired heatmap, methylTFR left, VIPER right
#####################################################################

mtfr_heatmap <- mtfr_view[pairs_tab$motif, heat_samples, drop = FALSE]
viper_heatmap <- viper_view[pairs_tab$viper_tf, heat_samples, drop = FALSE]
rownames(viper_heatmap) <- pairs_tab$motif

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
  row_split = row_split_factor,
  cluster_row_slices = FALSE,
  row_gap = unit(1, "mm"),
  row_title_gp = gpar(fontsize = 8),
  show_row_names = FALSE,
  show_column_names = FALSE,
  row_dend_side = "left",
  column_title_rot = 45,
  column_title_gp = gpar(fontsize = 7)
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
  column_title_rot = 45,
  column_title_gp = gpar(fontsize = 7)
)

row_correlation <- vapply(seq_len(nrow(mtfr_heatmap)), function(i) {
  suppressWarnings(cor(
    as.numeric(mtfr_heatmap[i, ]), as.numeric(viper_heatmap[i, ]),
    method = "pearson", use = "complete.obs"
  ))
}, numeric(1))
row_correlation[!is.finite(row_correlation)] <- 0

selex_call <- selex_call_of(rownames(mtfr_heatmap))
if (!is.null(selex_call)) {
  log_info(
    "SELEX calls on the heatmap: ",
    paste(names(table(selex_call)), table(selex_call), sep = " = ", collapse = ", ")
  )
}

selex_col_map <- NULL
if (!is.null(selex_call)) {
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

panel_heat <- grid.grabExpr(
  draw(ht_mtfr + ht_viper + ha_right,
    merge_legend = TRUE, heatmap_legend_side = "right"
  )
)

file <- file.path(plot.dir, "heatmap_mtfr_viper_bp.pdf")
pdf(file, width = 15, height = 11, bg = "transparent")
grid.draw(panel_heat)
dev.off()
log_info("Wrote ", file)

#####################################################################
# The assembled figure
#####################################################################

row_ab <- wrap_plots(
  p_modality_frac + labs(tag = "A"),
  p_strip + labs(tag = "B"),
  nrow = 1, widths = c(1, 1.15)
)

if (length(scatter_panels) == 0) {
  stop("No factor scatter could be drawn, the assembled figure needs a scatter panel")
}

scatter_panels[[1]] <- scatter_panels[[1]] + labs(tag = "C")
row_cd <- lapply(scatter_panels, function(p) p + theme(legend.position = "none"))
cd_widths <- rep(1, length(row_cd))
if (!is.null(p_selex)) {
  row_cd <- c(row_cd, list(p_selex + labs(tag = "D")))
  cd_widths <- c(cd_widths, 1)
}
row_cd <- wrap_plots(row_cd, nrow = 1, widths = cd_widths)

row_e <- wrap_elements(full = panel_heat) +
  labs(tag = if (is.null(p_selex)) "D" else "E")

combined <- wrap_plots(row_ab, row_cd, row_e, ncol = 1, heights = c(1, 1.1, 2.6)) &
  theme(plot.tag = element_text(face = "bold", size = 16))

file <- file.path(plot.dir, "figure_bp_viper_mofa_integration.pdf")
ggsave(file, combined & no_bg, width = 18, height = 20, bg = "transparent", limitsize = FALSE)
log_info("Wrote ", file)

log_success("Finished the VIPER integration, figures in ", plot.dir)

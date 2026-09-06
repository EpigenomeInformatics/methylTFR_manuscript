#!/usr/bin/env Rscript

#####################################################################
# 10_bp_viper_supplementary.R
# created on 06-09-2026 by Irem B Gunduz
# Supplementary to the VIPER integration figure
#   A  Association with the expectation, before and after the correction
#   B  Deviation scores over deciles of the expectation, before and after,
#      split by whether the factor is expressed at all
#   C  Cell type variance carried by TF mRNA and by VIPER activity
#   D  Agreement of each of the two with the methylTFR deviations
#   E  Histone ChIP agreement with both modalities, split by whether the
#      two modalities run opposite for that motif
#   F  The sign methyl-SELEX predicts
#   G  Three motifs from different families, sample by sample
#
# The expected deviation scores themselves live in the Blueprint
# supplementary, they are a property of the algorithm rather than of the
# integration. Panels E to G need 11_bp_chip_validation.R to have run and
# are dropped with a warning otherwise
#####################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(logger)
  library(methylTFR)
  library(methylTFRAnnotationHg38)
  library(SummarizedExperiment)
})
set.seed(42)

#####################################################################
# Settings
#####################################################################

# Distal regions avoid CpG islands, the genome wide set keeps them, so the
# second set carries a much wider spread of sequence composition and is the
# harder test of the correction
motif.sets <- c(
  "jaspar2020_distal" = "JASPAR2020 distal",
  "jaspar2020" = "JASPAR2020 all regions"
)
motifSet.main <- "jaspar2020_distal"

# Quantiles of TF expression that define the two ends of panel B
expr.low <- 0.25
expr.high <- 0.75

# Panel B summarises the motifs in quantile bins of the expectation rather
# than drawing all of them, so the trend is readable instead of a point cloud
n.bins <- 10

# Set to a character vector of motif names to pin the examples. Left automatic,
# no two examples may share a TF family, otherwise the whole panel is AP-1
example.motifs <- NULL
example.n <- 3L
scatter.label.n <- 6L

# The per mark panel is a sanity check rather than a result, the differences
# are small against their spread. The numbers still reach the log and the table
plot.mark.panel <- FALSE

drop.cell.types <- c("other", "thymocyte")
min.samples.per.celltype <- 3

fig.width <- 13
fig.height <- 17
base.size <- 10

state_colors <- c("Uncorrected" = "#B5A38A", "Corrected" = "#C2377C")
state_lines <- c("Uncorrected" = "dotted", "Corrected" = "solid")
group_colors <- c(
  "All motifs" = "#6A3D9A",
  "TF not expressed" = "#C2703D",
  "TF expressed" = "#2CA02C"
)
source_colors <- c("TF mRNA" = "#B5A38A", "VIPER activity" = "#6BC75A")
set_colors <- c("Integration motifs" = "#C2377C", "Background motifs" = "#B5A38A")
mark_colors <- c(
  "H3K27ac" = "#2B4B9B", "H3K4me1" = "#4EC3E0", "H3K4me3" = "#17BECF",
  "H3K36me3" = "#8C6D3F", "H3K27me3" = "#C2703D", "H3K9me3" = "#7B1E3D"
)
# H3K36me3 sits over transcribed gene bodies, which are methylated, so it is
# kept apart from the two enhancer and promoter marks rather than called active
mark_kind <- c(
  "H3K27ac" = "Enhancer and promoter", "H3K4me1" = "Enhancer and promoter",
  "H3K4me3" = "Enhancer and promoter", "H3K36me3" = "Gene body",
  "H3K27me3" = "Repressive", "H3K9me3" = "Repressive"
)
kind.levels <- c("Enhancer and promoter", "Gene body", "Repressive")

cell_type_labels <- c(
  "megK" = "Megakaryocytes", "eryt" = "Erythrocytes", "gran" = "Granulocytes",
  "mono" = "Monocytes", "Mf" = "Macrophages", "DC" = "Dendritic cells",
  "osteoclast" = "Osteoclast", "NK" = "NK-cells", "Tcell" = "T-cells",
  "thymocyte" = "Thymocyte", "Bcell" = "B-cells", "plasma" = "Plasma",
  "progenitor" = "Progenitors", "other" = "Other"
)
cell_type_levels <- unname(cell_type_labels)
cell_type_colors <- c(
  "B-cells" = "#C2377C", "Dendritic cells" = "#8C6D3F", "Erythrocytes" = "#7E4B2A",
  "Granulocytes" = "#E8A33D", "Macrophages" = "#B5A38A", "Megakaryocytes" = "#6A3D9A",
  "Monocytes" = "#C2703D", "NK-cells" = "#2CA02C", "Osteoclast" = "#8C8C8C",
  "Other" = "#BCBD22", "Plasma" = "#7B1E3D", "Progenitors" = "#17BECF",
  "T-cells" = "#4FC3D9", "Thymocyte" = "#5B9BD5"
)

analysis.dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/blueprint/"
dev.tag <- "mTFR_devs_230826"
mofa.dir <- file.path(analysis.dir, "mofa_230826")

counts.file <- file.path(mofa.dir, "bp_rawRNAcounts.RDS")
viper.file <- file.path(mofa.dir, "bp_viper_activity.RDS")
map.file <- file.path(mofa.dir, "bp_rna_wgbs_map.tsv")
chip.file <- file.path(mofa.dir, "bp_chip_signal.RDS")

github.dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR_manuscript/"
fig.dir <- file.path(github.dir, "figures", "blueprint")
mixed.dir <- file.path(github.dir, "tables", "mixed")
if (!dir.exists(fig.dir)) dir.create(fig.dir, recursive = TRUE)
if (!dir.exists(mixed.dir)) dir.create(mixed.dir, recursive = TRUE)

no_bg <- theme(
  plot.background = element_blank(),
  panel.background = element_blank(),
  legend.background = element_blank(),
  legend.box.background = element_blank(),
  legend.key = element_blank(),
  strip.background = element_blank()
)
tag_theme <- theme(plot.tag = element_text(face = "bold", size = 14))
base_theme <- theme_classic(base_size = base.size) +
  theme(strip.text = element_text(size = base.size - 1))

#####################################################################
# Helpers
#####################################################################

motif_to_tfs <- function(motif) {
  parts <- unlist(strsplit(motif, "::", fixed = TRUE))
  parts <- sub("\\s*\\(var\\.[0-9]+\\)$", "", parts)
  toupper(trimws(parts))
}

# methylTFR stores the observed deviation score minus the expected one, so
# the uncorrected score is recovered by adding the expectation back
load_scores <- function(set) {
  f <- file.path(analysis.dir, dev.tag, paste0(set, "_deviations.RDS"))
  if (!file.exists(f)) {
    log_warn("No deviations for ", set, ", the motif set is skipped: ", f)
    return(NULL)
  }
  obj <- readRDS(f)
  corrected <- rowMeans(deviations(obj), na.rm = TRUE)
  expected <- rowMeans(assay(obj, "expected"), na.rm = TRUE)
  data.table(
    set = unname(motif.sets[set]),
    motif = names(corrected),
    expected = unname(expected),
    corrected = unname(corrected),
    uncorrected = unname(corrected + expected)
  )
}

celltype_r2 <- function(values, groups) {
  ok <- is.finite(values)
  if (sum(ok) < 4 || length(unique(groups[ok])) < 2 || sd(values[ok]) == 0) {
    return(NA_real_)
  }
  ss <- summary(aov(v ~ g, data = data.frame(v = values[ok], g = groups[ok])))[[1]]
  ss["g", "Sum Sq"] / (ss["g", "Sum Sq"] + ss["Residuals", "Sum Sq"])
}

paired_p <- function(a, b) {
  ok <- is.finite(a) & is.finite(b)
  if (sum(ok) < 3) {
    return(NA_real_)
  }
  suppressWarnings(wilcox.test(a[ok], b[ok], paired = TRUE)$p.value)
}

# A single bracket over the two boxes, drawn from the data range so it clears
# the whiskers whatever the panel holds
add_bracket <- function(plot, values, label) {
  rng <- diff(range(values, na.rm = TRUE))
  y <- max(values, na.rm = TRUE) + rng * 0.08
  plot +
    annotate("segment", x = 1, xend = 2, y = y, yend = y, linewidth = 0.4) +
    annotate("segment", x = 1, xend = 1, y = y - rng * 0.02, yend = y, linewidth = 0.4) +
    annotate("segment", x = 2, xend = 2, y = y - rng * 0.02, yend = y, linewidth = 0.4) +
    annotate("text",
      x = 1.5, y = y + rng * 0.02, label = label,
      vjust = 0, size = 3, fontface = "bold"
    ) +
    coord_cartesian(ylim = c(min(values, na.rm = TRUE), y + rng * 0.14), clip = "off")
}

#####################################################################
# TF expression, shared by both halves of the figure
#####################################################################

if (!file.exists(counts.file)) stop("No RNA counts at ", counts.file, ". Run 07 first")
rna_cpm <- as.matrix(readRDS(counts.file))
rna_cpm <- rna_cpm[rowSums(rna_cpm, na.rm = TRUE) > 0, , drop = FALSE]
rna_cpm <- log2(t(t(rna_cpm) / colSums(rna_cpm, na.rm = TRUE)) * 1e6 + 1)
gene_mean <- rowMeans(rna_cpm, na.rm = TRUE)

#####################################################################
# The three scores per motif, in every available motif set
#####################################################################

scores <- rbindlist(lapply(names(motif.sets), load_scores))
if (nrow(scores) == 0) stop("No deviations were found for any motif set")
scores[, set := droplevels(factor(set, levels = unname(motif.sets)))]

scores[, tf_expr := vapply(motif, function(mo) {
  hit <- motif_to_tfs(mo)
  hit <- hit[hit %in% names(gene_mean)]
  if (length(hit) == 0) NA_real_ else max(gene_mean[hit])
}, numeric(1), USE.NAMES = FALSE)]

# The cut points are taken once, over the main motif set, so the two sets are
# split on the same expression level rather than on their own quantiles
cuts <- quantile(scores[set == unname(motif.sets[motifSet.main])]$tf_expr,
  c(expr.low, expr.high),
  na.rm = TRUE
)
scores[, expr_group := fcase(
  is.na(tf_expr), NA_character_,
  tf_expr <= cuts[1], "TF not expressed",
  tf_expr >= cuts[2], "TF expressed",
  default = NA_character_
)]

log_info(
  "Motifs per set: ",
  paste(levels(scores$set), table(scores$set), sep = " = ", collapse = ", ")
)
write.csv(scores, file.path(mixed.dir, "bp_bias_correction_scores.csv"), row.names = FALSE)

group.levels <- names(group_colors)
grouped <- rbindlist(list(
  copy(scores)[, group := "All motifs"],
  copy(scores[expr_group == "TF not expressed"])[, group := "TF not expressed"],
  copy(scores[expr_group == "TF expressed"])[, group := "TF expressed"]
))
grouped[, group := factor(group, levels = group.levels)]

summary_tab <- grouped[, .(
  n = .N,
  slope_uncorrected = unname(coef(lm(uncorrected ~ expected, .SD))[2]),
  slope_corrected = unname(coef(lm(corrected ~ expected, .SD))[2]),
  r_uncorrected = cor(uncorrected, expected, use = "complete.obs"),
  r_corrected = cor(corrected, expected, use = "complete.obs")
), by = .(set, group)]
summary_tab[, removed := 1 - abs(r_corrected) / abs(r_uncorrected)]
log_info("Correction summary")
print(summary_tab)
write.csv(summary_tab, file.path(mixed.dir, "bp_bias_correction_summary.csv"), row.names = FALSE)

#####################################################################
# A) Association with the expectation, before and after the correction
#####################################################################

dumbbell <- summary_tab[, .(set, group,
  Uncorrected = abs(r_uncorrected), Corrected = abs(r_corrected)
)]
dumbbell[, row := factor(paste(set, group, sep = "\n"),
  levels = rev(unique(paste(set, group, sep = "\n")))
)]
dumbbell_long <- melt(dumbbell,
  id.vars = "row", measure.vars = c("Uncorrected", "Corrected"),
  variable.name = "state", value.name = "abs_r"
)
dumbbell_long[, state := factor(state, levels = names(state_colors))]

p_a <- ggplot(dumbbell, aes(y = row)) +
  geom_segment(aes(x = Uncorrected, xend = Corrected, yend = row),
    colour = "grey60", linewidth = 0.5,
    arrow = arrow(length = unit(1.8, "mm"), type = "closed")
  ) +
  geom_point(data = dumbbell_long, aes(x = abs_r, colour = state), size = 2.4) +
  geom_text(aes(x = pmax(Uncorrected, Corrected), label = sprintf("%.2f to %.2f", Uncorrected, Corrected)),
    hjust = -0.15, size = 2.6, colour = "grey30"
  ) +
  scale_colour_manual(values = state_colors, name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.35))) +
  labs(
    x = "|Pearson r| between the deviation score and the expected one\n(lower is less composition dependence)",
    y = NULL
  ) +
  base_theme +
  theme(legend.position = "bottom", axis.text.y = element_text(size = base.size - 2))

#####################################################################
# B) Deviation scores over deciles of the expectation
#####################################################################

long_b <- melt(grouped,
  id.vars = c("set", "group", "motif", "expected"),
  measure.vars = c("uncorrected", "corrected"),
  variable.name = "state", value.name = "score"
)
long_b[, state := factor(fifelse(state == "uncorrected", "Uncorrected", "Corrected"),
  levels = names(state_colors)
)]
# The two states differ by roughly one unit, so each is centred on its own
# mean and only the tilt is compared
long_b[, score_c := score - mean(score, na.rm = TRUE), by = .(set, group, state)]

# Bin edges come from every motif in the set, so the three groups fall into the
# same bins and their lines can be read against one another
long_b[, bin := cut(expected,
  breaks = unique(quantile(expected, probs = seq(0, 1, length.out = n.bins + 1), na.rm = TRUE)),
  include.lowest = TRUE, labels = FALSE
), by = set]

binned <- long_b[, .(
  x = mean(expected, na.rm = TRUE),
  y = mean(score_c, na.rm = TRUE),
  se = sd(score_c, na.rm = TRUE) / sqrt(.N)
), by = .(set, group, state, bin)]

p_b <- ggplot(binned, aes(x = x, y = y, colour = group, linetype = state)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey60") +
  geom_errorbar(aes(ymin = y - se, ymax = y + se),
    width = 0, linewidth = 0.4, linetype = "solid", show.legend = FALSE
  ) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.1, show.legend = FALSE) +
  facet_wrap(~set, scales = "free", nrow = 1) +
  scale_colour_manual(values = group_colors, name = NULL) +
  scale_linetype_manual(values = state_lines, name = NULL) +
  labs(
    x = "Expected deviation score, decile bin mean",
    y = "Deviation score, bin mean +/- s.e.\n(centred on its own mean, flat = no composition dependence)"
  ) +
  base_theme +
  theme(legend.position = "bottom", legend.box = "horizontal")

#####################################################################
# TF mRNA against VIPER activity
#####################################################################

for (f in c(viper.file, map.file)) {
  if (!file.exists(f)) stop("Missing input: ", f, ". Run 07 and 07b first")
}

dev.main <- file.path(analysis.dir, dev.tag, paste0(motifSet.main, "_deviations.RDS"))
mtfr_z <- deviationZScores(readRDS(dev.main))
viper_act <- readRDS(viper.file)
map <- read.delim(map.file, stringsAsFactors = FALSE)

map <- map[
  map$bedFile %in% colnames(mtfr_z) &
    map$rna_id %in% colnames(viper_act) &
    map$rna_id %in% colnames(rna_cpm) &
    !map$cellTypeGroup %in% drop.cell.types,
]
keep_types <- names(which(table(map$cellTypeGroup) >= min.samples.per.celltype))
map <- map[map$cellTypeGroup %in% keep_types, ]
if (nrow(map) < 6) stop("Too few paired samples to compare mRNA with VIPER")
log_info(
  nrow(map), " paired samples across ", length(keep_types), " cell types: ",
  paste(keep_types, collapse = ", ")
)

viper_mat <- viper_act[, map$rna_id, drop = FALSE]
expr_mat <- rna_cpm[, map$rna_id, drop = FALSE]
mtfr_mat <- mtfr_z[, map$bedFile, drop = FALSE]
celltype <- factor(unname(cell_type_labels[map$cellTypeGroup]),
  levels = intersect(cell_type_levels, unname(cell_type_labels[map$cellTypeGroup]))
)

viper_row <- setNames(rownames(viper_mat), toupper(rownames(viper_mat)))
expr_row <- setNames(rownames(expr_mat), toupper(rownames(expr_mat)))

pairs_tab <- data.table(motif = rownames(mtfr_mat))
pairs_tab[, tf := vapply(motif, function(m) {
  hit <- motif_to_tfs(m)
  hit <- hit[hit %in% names(viper_row) & hit %in% names(expr_row)]
  if (length(hit) == 0) NA_character_ else hit[1]
}, character(1), USE.NAMES = FALSE)]
pairs_tab <- pairs_tab[!is.na(tf)]
if (nrow(pairs_tab) < 10) stop("Fewer than ten motifs carry both an mRNA and a VIPER row")
log_info(
  nrow(pairs_tab), " motifs carry both an mRNA and a VIPER measurement, ",
  uniqueN(pairs_tab$tf), " distinct TFs"
)

#####################################################################
# C) Cell type variance carried by each measure
#####################################################################

tf_tab <- data.table(tf = unique(pairs_tab$tf))
tf_tab[, `TF mRNA` := vapply(tf, function(t) {
  celltype_r2(as.numeric(expr_mat[expr_row[[t]], ]), celltype)
}, numeric(1), USE.NAMES = FALSE)]
tf_tab[, `VIPER activity` := vapply(tf, function(t) {
  celltype_r2(as.numeric(viper_mat[viper_row[[t]], ]), celltype)
}, numeric(1), USE.NAMES = FALSE)]
tf_tab <- tf_tab[is.finite(`TF mRNA`) & is.finite(`VIPER activity`)]

#####################################################################
# D) Agreement with the methylTFR deviations
#####################################################################

cor_with <- function(motif, tf, mat, rows) {
  suppressWarnings(cor(
    as.numeric(mtfr_mat[motif, ]), as.numeric(mat[rows[[tf]], ]),
    method = "pearson", use = "complete.obs"
  ))
}

pairs_tab[, `TF mRNA` := vapply(seq_len(.N), function(i) {
  cor_with(motif[i], tf[i], expr_mat, expr_row)
}, numeric(1))]
pairs_tab[, `VIPER activity` := vapply(seq_len(.N), function(i) {
  cor_with(motif[i], tf[i], viper_mat, viper_row)
}, numeric(1))]
pairs_tab <- pairs_tab[is.finite(`TF mRNA`) & is.finite(`VIPER activity`)]

# Two paired tests over the same comparison, so they are adjusted together
p_raw <- c(
  r2 = paired_p(tf_tab$`VIPER activity`, tf_tab$`TF mRNA`),
  cor = paired_p(abs(pairs_tab$`VIPER activity`), abs(pairs_tab$`TF mRNA`))
)
p_adj <- p.adjust(p_raw, method = "BH")
p_lab <- ifelse(is.na(p_adj), "p.adj n.a.",
  paste0("p.adj = ", format.pval(p_adj, digits = 2, eps = 1e-16))
)
log_info(
  "VIPER against mRNA, paired Wilcoxon: R2 p.adj = ", signif(p_adj["r2"], 3),
  ", correlation p.adj = ", signif(p_adj["cor"], 3)
)

r2_long <- melt(tf_tab, id.vars = "tf", variable.name = "source", value.name = "r2")
r2_long[, source := factor(source, levels = names(source_colors))]

p_c <- ggplot(r2_long, aes(x = source, y = r2, fill = source)) +
  stat_boxplot(geom = "errorbar", width = 0.3, linewidth = 0.5) +
  geom_boxplot(colour = "black", outlier.shape = NA, width = 0.55) +
  geom_jitter(width = 0.15, size = 0.7, alpha = 0.45, colour = "grey20") +
  scale_fill_manual(values = source_colors) +
  labs(
    x = sprintf("n = %d transcription factors", nrow(tf_tab)),
    y = expression("Cell type variance explained (" * R^2 * ")")
  ) +
  base_theme +
  theme(legend.position = "none")
p_c <- add_bracket(p_c, r2_long$r2, p_lab[["r2"]])

cor_long <- melt(pairs_tab,
  id.vars = c("motif", "tf"), variable.name = "source", value.name = "correlation"
)
cor_long[, source := factor(source, levels = names(source_colors))]
cor_long[, abs_cor := abs(correlation)]

p_d <- ggplot(cor_long, aes(x = source, y = abs_cor, fill = source)) +
  stat_boxplot(geom = "errorbar", width = 0.3, linewidth = 0.5) +
  geom_boxplot(colour = "black", outlier.shape = NA, width = 0.55) +
  geom_jitter(width = 0.15, size = 0.7, alpha = 0.45, colour = "grey20") +
  scale_fill_manual(values = source_colors) +
  labs(
    x = sprintf("n = %d motifs", nrow(pairs_tab)),
    y = "|Pearson r| with the methylTFR\ndeviation Z-scores, over samples"
  ) +
  base_theme +
  theme(legend.position = "none")
p_d <- add_bracket(p_d, cor_long$abs_cor, p_lab[["cor"]])

write.csv(pairs_tab, file.path(mixed.dir, "bp_viper_vs_expression_correlation.csv"),
  row.names = FALSE
)
write.csv(tf_tab, file.path(mixed.dir, "bp_viper_vs_expression_celltype_R2.csv"),
  row.names = FALSE
)

#####################################################################
# E to G) Histone ChIP, only if 11 has been run
#####################################################################

p_e <- p_f <- p_g <- p_h <- NULL

if (!file.exists(chip.file)) {
  log_warn("No ChIP signal at ", chip.file, ", panels E to H are skipped. Run 11 first")
} else {
  chip <- readRDS(chip.file)
  mark.primary <- chip$mark.primary
  chip_cor <- as.data.table(chip$correlations)
  chip_cor[, set := factor(set, levels = names(set_colors))]
  primary <- chip_cor[mark == mark.primary & is.finite(r_mtfr) & is.finite(r_viper)]
  if (nrow(primary) < 5) stop("Too few motifs carry a ChIP correlation")

  #############################################################
  # E) The geometry the three modalities should form
  #
  # Whether the integration picked a motif turns out to be the wrong split:
  # with four cell types anything that varies by lineage correlates with
  # anything else that does, and integration and background motifs land in
  # the same places. What does separate them is the motif's own relationship
  # between the two modalities. Where methylTFR and VIPER run opposite, the
  # factor loses deviation score as it gains activity, so the active mark
  # should fall on the methylTFR side and rise on the VIPER side. Where they
  # run together, no such expectation holds
  #############################################################

  primary[, coherent := r_mtfr < 0 & r_viper > 0]
  primary[, stratum := fifelse(
    r_mtfr_viper < 0,
    "methylTFR opposite VIPER", "methylTFR with VIPER"
  )]
  primary[, stratum := factor(stratum,
    levels = c("methylTFR opposite VIPER", "methylTFR with VIPER")
  )]

  strat <- table(primary$stratum, factor(primary$coherent, levels = c(FALSE, TRUE)))
  strat_p <- if (all(dim(strat) == c(2, 2)) && min(rowSums(strat)) > 0) {
    suppressWarnings(fisher.test(strat)$p.value)
  } else {
    NA_real_
  }
  ct <- suppressWarnings(cor.test(primary$r_mtfr, primary$r_viper))

  log_info(
    "Coherent quadrant by stratum: ",
    paste(rownames(strat), strat[, "TRUE"], "of", rowSums(strat), collapse = "; "),
    ", Fisher p = ", signif(strat_p, 3)
  )
  log_info(
    "Coupling of the two ChIP correlations: r = ", round(unname(ct$estimate), 3),
    ", p = ", signif(ct$p.value, 3)
  )
  log_info(
    "Integration motifs in the quadrant ",
    nrow(primary[coherent == TRUE & set == "Integration motifs"]), " of ",
    nrow(primary[set == "Integration motifs"]), ", background ",
    nrow(primary[coherent == TRUE & set == "Background motifs"]), " of ",
    nrow(primary[set == "Background motifs"]),
    ", which is the split that does not separate"
  )

  primary[, combined := abs(r_mtfr) + abs(r_viper)]
  labelled <- primary[, .SD[order(-combined)][seq_len(min(.N, scatter.label.n))], by = stratum]

  strat_lab <- data.table(
    stratum = factor(rownames(strat), levels = levels(primary$stratum)),
    label = sprintf(
      "%d of %d motifs (%.0f%%)\nin the shaded quadrant",
      strat[, "TRUE"], rowSums(strat), 100 * strat[, "TRUE"] / rowSums(strat)
    )
  )

  p_e <- ggplot(primary, aes(x = r_mtfr, y = r_viper)) +
    annotate("rect",
      xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf, fill = "grey92", alpha = 0.6
    ) +
    geom_hline(yintercept = 0, linetype = "dotted", colour = "grey50") +
    geom_vline(xintercept = 0, linetype = "dotted", colour = "grey50") +
    geom_smooth(
      method = "lm", formula = y ~ x, se = FALSE,
      colour = "grey35", linewidth = 0.5
    ) +
    geom_point(aes(colour = set), size = 1.7, alpha = 0.85) +
    geom_text_repel(
      data = labelled, aes(label = motif), size = 2.2,
      max.overlaps = Inf, segment.size = 0.2, show.legend = FALSE
    ) +
    geom_text(
      data = strat_lab, inherit.aes = FALSE,
      aes(x = -Inf, y = -Inf, label = label),
      hjust = -0.06, vjust = -0.4, size = 2.4, colour = "grey30"
    ) +
    facet_wrap(~stratum, nrow = 1) +
    scale_colour_manual(values = set_colors, name = NULL) +
    labs(
      x = sprintf(
        "Pearson r, %s against the methylTFR deviation Z-score\nquadrant enrichment %s, overall r = %.2f (%s)",
        mark.primary,
        if (is.na(strat_p)) "p n.a." else paste0("p = ", format.pval(strat_p, digits = 2, eps = 1e-16)),
        unname(ct$estimate),
        paste0("p = ", format.pval(ct$p.value, digits = 2, eps = 1e-16))
      ),
      y = sprintf("Pearson r, %s against\nthe VIPER activity", mark.primary)
    ) +
    base_theme +
    theme(legend.position = "bottom")

  #############################################################
  # F) The sign methyl-SELEX predicts
  #
  # The strongest test available, because it is specific to methylTFR rather
  # than to any cell type signal: a MethylMinus factor binds unmethylated DNA
  # and should lose deviation score where the active mark rises, a MethylPlus
  # factor should gain it. No background matching can produce that contrast
  #############################################################

  selex_calls <- NULL
  selex.candidates <- c(
    file.path(github.dir, "tables", "Selex_data.csv"),
    "/icbb/projects/igunduz/exposure_atlas_manuscript/sample_annots/Selex_data.csv"
  )
  selex.file <- selex.candidates[file.exists(selex.candidates)][1]
  if (is.na(selex.file)) {
    log_warn("No SELEX table found, panel F is skipped")
  } else {
    selex_tab <- read.csv(selex.file, skip = 20, header = TRUE, sep = ";")
    call_col <- intersect(c("methyl.SELEX.call", "Call", "methylSELEXcall"), colnames(selex_tab))
    name_col <- intersect(c("TF.name", "TF", "TFname"), colnames(selex_tab))
    if (length(call_col) == 0 || length(name_col) == 0) {
      log_warn("No SELEX call or TF name column, panel F is skipped")
    } else {
      selex_calls <- setNames(
        as.character(selex_tab[[call_col[1]]]),
        toupper(as.character(selex_tab[[name_col[1]]]))
      )
    }
  }

  if (!is.null(selex_calls)) {
    primary[, selex := vapply(motif, function(m) {
      hit <- motif_to_tfs(m)
      hit <- hit[hit %in% names(selex_calls)]
      if (length(hit) == 0) {
        return(NA_character_)
      }
      call <- selex_calls[[hit[1]]]
      if (is.na(call) || !nzchar(call)) NA_character_ else call
    }, character(1), USE.NAMES = FALSE)]

    selex_dt <- primary[selex %in% c("MethylMinus", "MethylPlus")]
    selex_dt[, selex := factor(selex, levels = c("MethylMinus", "MethylPlus"))]
    if (nrow(selex_dt) < 6 || uniqueN(selex_dt$selex) < 2) {
      log_warn("Too few motifs carry a SELEX call, panel F is skipped")
    } else {
      selex_p <- suppressWarnings(wilcox.test(
        r_mtfr ~ selex,
        data = selex_dt
      )$p.value)
      medians <- selex_dt[, .(m = median(r_mtfr, na.rm = TRUE), n = .N), by = selex]
      log_info(
        "SELEX sign test: ",
        paste(medians$selex, " median r = ", round(medians$m, 3),
          " (n = ", medians$n, ")",
          collapse = "; "
        ),
        ", Wilcoxon p = ", signif(selex_p, 3)
      )

      p_f <- ggplot(selex_dt, aes(x = selex, y = r_mtfr, fill = selex)) +
        geom_hline(yintercept = 0, linetype = "dotted", colour = "grey50") +
        stat_boxplot(geom = "errorbar", width = 0.3, linewidth = 0.5) +
        geom_boxplot(colour = "black", outlier.shape = NA, width = 0.55) +
        geom_jitter(width = 0.15, size = 0.8, alpha = 0.5, colour = "grey20") +
        scale_fill_manual(values = c("MethylMinus" = "#8B1A1A", "MethylPlus" = "#1B6B3A")) +
        labs(
          x = sprintf(
            "methyl-SELEX call (%s)",
            paste(medians$selex, medians$n, sep = " n = ", collapse = ", ")
          ),
          y = sprintf("Pearson r, %s against the\nmethylTFR deviation Z-score", mark.primary)
        ) +
        base_theme +
        theme(legend.position = "none")
      p_f <- add_bracket(
        p_f, selex_dt$r_mtfr,
        paste0("p = ", format.pval(selex_p, digits = 2, eps = 1e-16))
      )
    }
  }

  #############################################################
  # G) Every mark, integration motifs only
  #############################################################

  by_mark <- chip_cor[set == "Integration motifs" & is.finite(r_mtfr)]
  by_mark[, mark := factor(mark, levels = intersect(names(mark_colors), unique(mark)))]
  by_mark[, kind := factor(unname(mark_kind[as.character(mark)]), levels = kind.levels)]
  mark_summary <- by_mark[, .(
    median = median(r_mtfr, na.rm = TRUE),
    lower = quantile(r_mtfr, 0.25, na.rm = TRUE),
    upper = quantile(r_mtfr, 0.75, na.rm = TRUE),
    n = .N
  ), by = .(mark, kind)]

  # Paired per motif, so the contrast does not ride on which motifs happen to
  # carry which mark
  wide <- dcast(droplevels(by_mark[kind != "Gene body"]), motif ~ kind,
    value.var = "r_mtfr", fun.aggregate = mean, na.rm = TRUE
  )
  wide <- wide[is.finite(`Enhancer and promoter`) & is.finite(Repressive)]
  kind_p <- paired_p(wide$`Enhancer and promoter`, wide$Repressive)
  log_info(
    "Enhancer and promoter marks against repressive, paired over ", nrow(wide),
    " motifs: median difference ",
    round(median(wide$`Enhancer and promoter` - wide$Repressive, na.rm = TRUE), 3),
    ", Wilcoxon p = ", signif(kind_p, 3)
  )

  p_g <- ggplot(mark_summary, aes(x = mark, y = median, colour = mark)) +
    geom_hline(yintercept = 0, linetype = "dotted", colour = "grey60") +
    geom_linerange(aes(ymin = lower, ymax = upper), linewidth = 0.8) +
    geom_point(size = 2.6) +
    facet_grid(~kind, scales = "free_x", space = "free_x") +
    scale_colour_manual(values = mark_colors) +
    labs(
      x = sprintf(
        "n = %d integration motifs\nenhancer against repressive, paired %s",
        uniqueN(by_mark$motif),
        if (is.na(kind_p)) "p n.a." else paste0("p = ", format.pval(kind_p, digits = 2, eps = 1e-16))
      ),
      y = "Pearson r with the methylTFR\ndeviation Z-score, median and IQR"
    ) +
    base_theme +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
  if (!plot.mark.panel) p_g <- NULL

  #############################################################
  # H) The motifs where the three modalities agree most closely
  #############################################################

  # JASPAR is full of AP-1 dimers, and ranking on agreement alone returns
  # nothing else, so a family is allowed one example
  tf_family <- function(tf) {
    fam <- sub("[0-9]+$", "", tf)
    fam <- sub("^FOSL$", "FOS", fam)
    fam <- sub("^JUN[BD]$", "JUN", fam)
    fam
  }

  picked <- if (!is.null(example.motifs)) {
    primary[motif %in% example.motifs]
  } else {
    ranked <- primary[coherent == TRUE][order(-combined)]
    taken <- character(0)
    keep <- integer(0)
    for (i in seq_len(nrow(ranked))) {
      fams <- unique(vapply(motif_to_tfs(ranked$motif[i]), tf_family, character(1),
        USE.NAMES = FALSE
      ))
      if (any(fams %in% taken)) next
      taken <- c(taken, fams)
      keep <- c(keep, i)
      if (length(keep) == example.n) break
    }
    ranked[keep]
  }
  if (nrow(picked) == 0) {
    log_warn("No coherent example motif, panel H is skipped")
  } else {
    log_info("Panel H examples: ", paste(picked$motif, collapse = ", "))
    chip_mat <- chip$signals[[mark.primary]]
    chip_samples <- colnames(chip_mat)
    chip_map <- as.data.table(chip$map)
    chip_rna <- chip_map$rna_id[match(chip_samples, chip_map$bedFile)]
    chip_types <- factor(
      unname(cell_type_labels[chip_map$cellTypeGroup[match(chip_samples, chip_map$bedFile)]]),
      levels = cell_type_levels
    )

    examples <- rbindlist(lapply(seq_len(nrow(picked)), function(i) {
      m <- picked$motif[i]
      tf <- picked$viper_tf[i]
      rbindlist(list(
        data.table(
          motif = m, measure = "methylTFR deviation Z-score",
          r = picked$r_mtfr[i], chip = as.numeric(chip_mat[m, ]),
          value = as.numeric(mtfr_z[m, chip_samples]), celltype = chip_types
        ),
        data.table(
          motif = m, measure = "VIPER activity",
          r = picked$r_viper[i], chip = as.numeric(chip_mat[m, ]),
          value = as.numeric(viper_act[tf, chip_rna]), celltype = chip_types
        )
      ))
    }))
    examples[, motif_lab := factor(motif, levels = picked$motif)]
    examples[, measure_lab := factor(measure,
      levels = c("methylTFR deviation Z-score", "VIPER activity")
    )]
    r_lab <- unique(examples[, .(motif_lab, measure_lab, r)])

    p_h <- ggplot(examples, aes(x = chip, y = value)) +
      geom_smooth(
        method = "lm", formula = y ~ x, se = FALSE,
        colour = "grey40", linewidth = 0.5
      ) +
      geom_point(aes(colour = celltype), size = 1.6, alpha = 0.9) +
      geom_text(
        data = r_lab, inherit.aes = FALSE,
        aes(x = -Inf, y = Inf, label = sprintf("r = %.2f", r)),
        hjust = -0.2, vjust = 1.4, size = 2.6, colour = "grey25"
      ) +
      facet_grid(measure_lab ~ motif_lab, scales = "free_y", switch = "y") +
      scale_colour_manual(values = cell_type_colors, name = NULL, drop = TRUE) +
      labs(
        x = sprintf("%s signal at the motif's distal sites, per sample", mark.primary),
        y = NULL
      ) +
      base_theme +
      theme(legend.position = "bottom", strip.placement = "outside")
  }

  #############################################################
  # Reported in the log, an alternative framing if the per motif
  # correlations stay noisy: agreement taken per sample across the
  # motifs instead of per motif across the samples
  #############################################################

  chip_mat <- chip$signals[[mark.primary]]
  shared <- intersect(rownames(chip_mat), rownames(mtfr_z))
  per_sample <- vapply(colnames(chip_mat), function(s) {
    suppressWarnings(cor(chip_mat[shared, s], mtfr_z[shared, s],
      method = "spearman", use = "complete.obs"
    ))
  }, numeric(1))
  log_info(
    "Per sample agreement across ", length(shared), " motifs: median rho = ",
    round(median(per_sample, na.rm = TRUE), 3), ", ",
    sum(per_sample < 0, na.rm = TRUE), " of ", sum(is.finite(per_sample)),
    " samples negative as predicted"
  )
}

#####################################################################
# The assembled figure
#####################################################################

tag_seq <- 0L
tagged <- function(p) {
  tag_seq <<- tag_seq + 1L
  p + labs(tag = LETTERS[tag_seq])
}

build_row <- function(items) {
  keep <- Filter(function(x) !is.null(x[[1]]), items)
  if (length(keep) == 0) {
    return(NULL)
  }
  wrap_plots(lapply(keep, function(x) tagged(x[[1]])),
    nrow = 1, widths = vapply(keep, function(x) x[[2]], numeric(1))
  )
}

row_specs <- list(
  list(items = list(list(p_a, 1), list(p_b, 1.5)), height = 1),
  list(items = list(list(p_c, 1), list(p_d, 1)), height = 1),
  list(items = list(list(p_e, 2), list(p_f, 1), list(p_g, 1)), height = 1.15),
  list(items = list(list(p_h, 1)), height = 1.3)
)

rows <- list()
heights <- numeric(0)
for (spec in row_specs) {
  r <- build_row(spec$items)
  if (is.null(r)) next
  rows <- c(rows, list(r))
  heights <- c(heights, spec$height)
}

supplementary <- wrap_plots(rows, ncol = 1, heights = heights) & tag_theme

file <- file.path(fig.dir, "blueprint_viper_supplementary.pdf")
ggsave(file, supplementary & no_bg,
  width = fig.width, height = fig.height * sum(heights) / 4.45,
  bg = "transparent", limitsize = FALSE
)
log_success("Wrote ", file)

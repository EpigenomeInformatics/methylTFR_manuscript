#!/usr/bin/env Rscript

#####################################################################
# 10_bp_viper_supplementary.R
# created on 06-09-2026 by Irem B Gunduz
# Supplementary to the VIPER integration figure
#   A  Expected deviation scores across motifs, in two motif sets
#   B  Association with the expectation, before and after the correction
#   C  Deviation scores over deciles of the expectation, before and after,
#      split by whether the factor is expressed at all
#   D  Cell type variance carried by TF mRNA and by VIPER activity
#   E  Agreement of each of the two with the methylTFR deviations
#   F  Three factors where the two disagree, over the cell types
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

# Quantiles of TF expression that define the two ends of panel C
expr.low <- 0.25
expr.high <- 0.75
label.n <- 6

# Panel C summarises the motifs in quantile bins of the expectation rather
# than drawing all of them, so the trend is readable instead of a point cloud
n.bins <- 10

# Panel F, set to a character vector of motif names to pin the examples
example.motifs <- NULL
example.n <- 3

drop.cell.types <- c("other", "thymocyte")
min.samples.per.celltype <- 3

fig.width <- 13
fig.height <- 15
base.size <- 10

state_colors <- c("Uncorrected" = "#B5A38A", "Corrected" = "#C2377C")
state_lines <- c("Uncorrected" = "dotted", "Corrected" = "solid")
group_colors <- c(
  "All motifs" = "#6A3D9A",
  "TF not expressed" = "#C2703D",
  "TF expressed" = "#2CA02C"
)
source_colors <- c("TF mRNA" = "#B5A38A", "VIPER activity" = "#6BC75A")
measure_colors <- c(
  "methylTFR deviation" = "#ED4B4A",
  "VIPER activity" = "#6BC75A",
  "TF mRNA" = "#B5A38A"
)

cell_type_labels <- c(
  "megK" = "Megakaryocytes", "eryt" = "Erythrocytes", "gran" = "Granulocytes",
  "mono" = "Monocytes", "Mf" = "Macrophages", "DC" = "Dendritic cells",
  "osteoclast" = "Osteoclast", "NK" = "NK-cells", "Tcell" = "T-cells",
  "thymocyte" = "Thymocyte", "Bcell" = "B-cells", "plasma" = "Plasma",
  "progenitor" = "Progenitors", "other" = "Other"
)
cell_type_levels <- unname(cell_type_labels)

analysis.dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/blueprint/"
dev.tag <- "mTFR_devs_230826"
mofa.dir <- file.path(analysis.dir, "mofa_230826")

counts.file <- file.path(mofa.dir, "bp_rawRNAcounts.RDS")
viper.file <- file.path(mofa.dir, "bp_viper_activity.RDS")
map.file <- file.path(mofa.dir, "bp_rna_wgbs_map.tsv")

github.dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR_manuscript/"
fig.dir <- file.path(github.dir, "figures", "blueprint")
table.dir <- file.path(github.dir, "tables")
mixed.dir <- file.path(table.dir, "mixed")
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
    annotate("text", x = 1.5, y = y + rng * 0.02, label = label,
      vjust = 0, size = 3, fontface = "bold"
    ) +
    coord_cartesian(
      ylim = c(min(values, na.rm = TRUE), y + rng * 0.14), clip = "off"
    )
}

zscore_vec <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
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
# A) The expected deviation scores
#####################################################################

ranked <- scores[order(set, expected)]
ranked[, rank := seq_len(.N), by = set]
ranked[, set_lab := sprintf("%s (%d motifs, range %.3f)", set, .N, max(expected) - min(expected)),
  by = set
]
extremes <- ranked[, .SD[c(seq_len(label.n), seq.int(.N - label.n + 1, .N))], by = set]

p_a <- ggplot(ranked, aes(x = rank, y = expected)) +
  geom_hline(yintercept = 1, linetype = "dotted", colour = "grey50") +
  geom_point(size = 0.6, colour = "grey30") +
  geom_text_repel(
    data = extremes, aes(label = motif), size = 2.2,
    max.overlaps = Inf, colour = "#7B1E3D", segment.size = 0.2
  ) +
  facet_wrap(~set_lab, scales = "free", nrow = 1) +
  labs(
    x = "Motifs, ranked by their expected deviation score",
    y = "Expected deviation score\n(GC predicted centre / flank ratio, 1 = no drift)"
  ) +
  base_theme

#####################################################################
# B) Association with the expectation, before and after the correction
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

p_b <- ggplot(dumbbell, aes(y = row)) +
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
# C) Deviation scores over deciles of the expectation
#####################################################################

long_c <- melt(grouped,
  id.vars = c("set", "group", "motif", "expected"),
  measure.vars = c("uncorrected", "corrected"),
  variable.name = "state", value.name = "score"
)
long_c[, state := factor(fifelse(state == "uncorrected", "Uncorrected", "Corrected"),
  levels = names(state_colors)
)]
# The two states differ by roughly one unit, so each is centred on its own
# mean and only the tilt is compared
long_c[, score_c := score - mean(score, na.rm = TRUE), by = .(set, group, state)]

# Bin edges come from every motif in the set, so the three groups fall into the
# same bins and their lines can be read against one another
long_c[, bin := cut(expected,
  breaks = unique(quantile(expected, probs = seq(0, 1, length.out = n.bins + 1), na.rm = TRUE)),
  include.lowest = TRUE, labels = FALSE
), by = set]

binned <- long_c[, .(
  x = mean(expected, na.rm = TRUE),
  y = mean(score_c, na.rm = TRUE),
  se = sd(score_c, na.rm = TRUE) / sqrt(.N)
), by = .(set, group, state, bin)]

p_c <- ggplot(binned, aes(x = x, y = y, colour = group, linetype = state)) +
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
# D) Cell type variance carried by each measure
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
# E) Agreement with the methylTFR deviations
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

p_d <- ggplot(r2_long, aes(x = source, y = r2, fill = source)) +
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
p_d <- add_bracket(p_d, r2_long$r2, p_lab[["r2"]])

cor_long <- melt(pairs_tab,
  id.vars = c("motif", "tf"), variable.name = "source", value.name = "correlation"
)
cor_long[, source := factor(source, levels = names(source_colors))]
cor_long[, abs_cor := abs(correlation)]

p_e <- ggplot(cor_long, aes(x = source, y = abs_cor, fill = source)) +
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
p_e <- add_bracket(p_e, cor_long$abs_cor, p_lab[["cor"]])

write.csv(pairs_tab, file.path(mixed.dir, "bp_viper_vs_expression_correlation.csv"),
  row.names = FALSE
)
write.csv(tf_tab, file.path(mixed.dir, "bp_viper_vs_expression_celltype_R2.csv"),
  row.names = FALSE
)

#####################################################################
# F) Examples where the deviations follow activity and not mRNA
#####################################################################

pairs_tab[, gain := abs(`VIPER activity`) - abs(`TF mRNA`)]
if (is.null(example.motifs)) {
  picked <- pairs_tab[`VIPER activity` > 0][order(-gain)]
  picked <- head(picked[!duplicated(tf)], example.n)
} else {
  picked <- pairs_tab[motif %in% example.motifs]
}
if (nrow(picked) == 0) stop("No example motif could be chosen for panel F")
log_info("Panel F examples: ", paste(picked$motif, collapse = ", "))

celltype_mean <- function(values) {
  zscore_vec(tapply(values, celltype, mean, na.rm = TRUE))
}

examples <- rbindlist(lapply(seq_len(nrow(picked)), function(i) {
  m <- picked$motif[i]
  t <- picked$tf[i]
  tracks <- list(
    "methylTFR deviation" = celltype_mean(as.numeric(mtfr_mat[m, ])),
    "VIPER activity" = celltype_mean(as.numeric(viper_mat[viper_row[[t]], ])),
    "TF mRNA" = celltype_mean(as.numeric(expr_mat[expr_row[[t]], ]))
  )
  rbindlist(lapply(names(tracks), function(k) {
    data.table(
      label = sprintf(
        "%s\nVIPER r = %.2f, mRNA r = %.2f", m,
        picked$`VIPER activity`[i], picked$`TF mRNA`[i]
      ),
      celltype = factor(levels(celltype), levels = levels(celltype)),
      measure = k,
      value = as.numeric(tracks[[k]])
    )
  }))
}))
examples[, measure := factor(measure, levels = names(measure_colors))]
examples[, label := factor(label, levels = unique(label))]

p_f <- ggplot(examples, aes(x = celltype, y = value, colour = measure, group = measure)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey60") +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.6) +
  facet_wrap(~label, nrow = 1) +
  scale_colour_manual(values = measure_colors, name = NULL) +
  labs(x = NULL, y = "Cell type mean,\nZ-scored across cell types") +
  base_theme +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

write.csv(examples, file.path(mixed.dir, "bp_viper_vs_expression_examples.csv"),
  row.names = FALSE
)

#####################################################################
# The assembled figure
#####################################################################

supplementary <- wrap_plots(
  wrap_plots(p_a + labs(tag = "A"), p_b + labs(tag = "B"), nrow = 1, widths = c(1.5, 1)),
  p_c + labs(tag = "C"),
  wrap_plots(p_d + labs(tag = "D"), p_e + labs(tag = "E"), nrow = 1),
  p_f + labs(tag = "F"),
  ncol = 1, heights = c(1, 1.1, 1, 1.1)
) & tag_theme

file <- file.path(fig.dir, "blueprint_viper_supplementary.pdf")
ggsave(file, supplementary & no_bg,
  width = fig.width, height = fig.height, bg = "transparent", limitsize = FALSE
)
log_success("Wrote ", file)

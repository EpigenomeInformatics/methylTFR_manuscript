#!/usr/bin/env Rscript

#####################################################################
# 10_bp_viper_supplementary.R  —  Irem B Gunduz, 06-09-2026
# Supplementary to the VIPER integration figure (panels A–I).
# Panels E–H need 13_bp_tf_chip_validation.R to have run, else dropped.
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

# distal set avoids CpG islands, the genome-wide set keeps them (harder test)
motif.sets <- c(
  "jaspar2020_distal" = "JASPAR2020 distal",
  "jaspar2020" = "JASPAR2020 all regions"
)
motifSet.main <- "jaspar2020_distal"

# TF expression quantiles defining the two ends of panel B
expr.low <- 0.25
expr.high <- 0.75

# panel B bins the motifs by expectation so the trend is readable
n.bins <- 10

drop.cell.types <- c("other", "thymocyte")
min.samples.per.celltype <- 3

fig.width <- 13
fig.height <- 22
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
fit_colour <- "#D62728"

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

analysis.dir <- "/icbb_triton/scratch/igunduz/methylTFR_manuscript/blueprint/"
dev.tag <- "mTFR_devs_230826"
mofa.dir <- file.path(analysis.dir, "mofa_230826")

counts.file <- file.path(mofa.dir, "bp_rawRNAcounts.RDS")
viper.file <- file.path(mofa.dir, "bp_viper_activity.RDS")
map.file <- file.path(mofa.dir, "bp_rna_wgbs_map.tsv")

github.dir <- "/icbb_triton/scratch/igunduz/methylTFR_manuscript/github/methylTFR_manuscript/"
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

# stored score is observed minus expected; add expected back for uncorrected
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

# bracket over the two boxes, above the whiskers
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

# cut points taken once on the main set so both sets split on the same level
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
# centre each state on its own mean, compare only the tilt
long_b[, score_c := score - mean(score, na.rm = TRUE), by = .(set, group, state)]

# bin edges from all motifs so the three groups share bins
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

# two paired tests over the same comparison, adjusted together
p_raw <- c(
  r2 = paired_p(tf_tab$`VIPER activity`, tf_tab$`TF mRNA`),
  cor = paired_p(abs(pairs_tab$`VIPER activity`), abs(pairs_tab$`TF mRNA`))
)
p_adj <- p.adjust(p_raw, method = "BH")
p_lab <- ifelse(is.na(p_adj), "p.adj n.a.",
  paste0("p.adj = ", format.pval(p_adj, digits = 2, eps = 1e-16))
)
log_info(
  "VIPER against mRNA, paired Wilcoxon: R2 p.adj = ", signif(unname(p_adj["r2"]), 3),
  ", correlation p.adj = ", signif(unname(p_adj["cor"]), 3)
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
# E) Per motif correlation of H3K27ac with each modality, only if 11 has run
#####################################################################

h27cor.file <- file.path(mixed.dir, "bp_chip_validation_correlations.csv")
p_h27cor <- NULL
if (!file.exists(h27cor.file)) {
  log_warn("No ", h27cor.file, ", panel E is skipped. Run 11 first")
} else {
  h27cor <- fread(h27cor.file)
  h27cor <- h27cor[mark == "H3K27ac"]
  h27cor_long <- rbind(
    data.table(modality = "methylTFR deviation", r = h27cor$r_mtfr),
    data.table(modality = "VIPER activity", r = h27cor$r_viper)
  )
  h27cor_long <- h27cor_long[is.finite(r)]
  h27cor_long[, modality := factor(modality,
    levels = c("methylTFR deviation", "VIPER activity"))]

  h27cor_lab <- h27cor_long[, .(
    med = median(r),
    p = tryCatch(wilcox.test(r)$p.value, error = function(e) NA_real_)
  ), by = modality]
  h27cor_lab[, padj := p.adjust(p, method = "BH")]
  h27cor_lab[, label := sprintf("median %.2f\np.adj = %.3f", med, padj)]
  log_info("H3K27ac correlation medians: ",
    paste(h27cor_lab[, paste0(modality, " ", round(med, 3))], collapse = "; "))

  h27cor.col <- c("methylTFR deviation" = "#ED4B4A", "VIPER activity" = "#6BC75A")
  p_h27cor <- ggplot(h27cor_long, aes(x = modality, y = r, fill = modality)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3, colour = "grey50") +
    stat_boxplot(geom = "errorbar", width = 0.3, linewidth = 0.5) +
    geom_boxplot(colour = "black", outlier.shape = NA, width = 0.55) +
    geom_jitter(width = 0.2, size = 0.5, alpha = 0.5, colour = "grey30") +
    geom_text(data = h27cor_lab, aes(x = modality, y = 1.04, label = label),
      inherit.aes = FALSE, vjust = 0, size = 3) +
    scale_fill_manual(values = h27cor.col, guide = "none") +
    scale_x_discrete(labels = c("methylTFR", "VIPER")) +
    coord_cartesian(ylim = c(-0.75, 1.25), clip = "off") +
    labs(x = NULL,
      y = "Pearson r with the H3K27ac signal, across samples within a motif",
      title = paste0("n = ", nrow(h27cor), " motifs")) +
    base_theme +
    theme(plot.title = element_text(hjust = 0, face = "plain", size = base.size))
}

#####################################################################
# F) H3K27ac at the distal motif sites, monocytes vs T-cells, only if 11 has run
#####################################################################

p_h27 <- NULL
prof.file <- file.path(mixed.dir, "bp_chip_profiles.csv")

if (!file.exists(prof.file)) {
  log_warn("No ", prof.file, ", panel F is skipped. Run 11 first")
} else {
  h27.keep <- c(mono = "Monocytes", Tcell = "T-cells")
  h27.col <- c(Monocytes = "#C2703D", `T-cells` = "#4FC3D9")
  h27.motifs <- c("CEBPB", "SPI1", "JUNB", "RELA")

  prof <- fread(prof.file)
  prof <- prof[cellTypeGroup %chin% names(h27.keep)]
  prof[, celltype := factor(h27.keep[cellTypeGroup], levels = h27.keep)]
  prof[, at := fifelse(group == "Random distal positions",
                       "random distal positions", "motif sites")]
  prof[, motif := fifelse(at == "motif sites", group, NA_character_)]

  # the random curve is the shared baseline, repeated in each facet
  rand <- prof[at == "random distal positions"]
  prof <- rbindlist(lapply(h27.motifs, function(m) {
    rbind(prof[motif == m], copy(rand)[, motif := m])
  }))
  prof[, motif := factor(motif, levels = h27.motifs)]

  prof_s <- prof[, .(
    mean = mean(value, na.rm = TRUE),
    se = sd(value, na.rm = TRUE) / sqrt(sum(is.finite(value)))
  ), by = .(motif, celltype, at, pos)]

  log_info(
    "H3K27ac at the motif centre: ",
    paste(prof_s[abs(pos) <= 250,
                 .(m = round(mean(mean), 3)), by = .(motif, celltype, at)][
                   , paste0(motif, " ", celltype, " ", at, " ", m)],
          collapse = "; ")
  )

  p_h27 <- ggplot(prof_s, aes(x = pos, y = mean, colour = celltype,
                              linetype = at, fill = celltype)) +
    geom_hline(yintercept = 1, linewidth = 0.3, colour = "grey85",
               linetype = "dotted") +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey85",
               linetype = "dotted") +
    geom_ribbon(data = prof_s[at == "motif sites"],
                aes(ymin = mean - se, ymax = mean + se),
                alpha = 0.15, colour = NA) +
    geom_line(linewidth = 0.5) +
    facet_wrap(~motif, nrow = 1) +
    scale_colour_manual(values = h27.col) +
    scale_fill_manual(values = h27.col, guide = "none") +
    scale_linetype_manual(values = c(`motif sites` = "solid",
                                     `random distal positions` = "22")) +
    scale_x_continuous(breaks = c(-2000, 0, 2000),
                       labels = c("-2 kb", "motif", "+2 kb")) +
    labs(
      x = "distance from the motif centre", y = "H3K27ac fold change",
      colour = NULL, linetype = NULL
    ) +
    base_theme +
    theme(legend.position = "bottom",
          legend.box = "vertical",
          legend.spacing.y = unit(0.02, "cm"),
          plot.title = element_text(hjust = 0, face = "plain", size = base.size))
}

#####################################################################
# G to I) TF ChIP occupancy, only if 13 has run
# Bound sites are hypomethylated vs matched controls, lost against another lineage.
#####################################################################

p_cov <- p_split <- p_spec <- NULL

chip.dir <- file.path(github.dir, "tables", "chip_validation")
chip.A <- file.path(chip.dir, "A_site_level_split.csv")

if (!file.exists(chip.A)) {
  log_warn("No TF ChIP table at ", chip.A, ", panels G to I are skipped. Run 13 first")
} else {
  col.matched <- "#C2377C"
  col.control <- "#2B4B9B"
  ink.muted <- "grey40"
  ink.mark <- "grey25"

  chipA <- fread(chip.A)

  ct.label <- c(
    GM12878 = "GM12878\n(B cell line)", monocyte = "Monocytes",
    Tcell = "T-cells", Bcell = "B-cells"
  )
  chipA[, ct := factor(ct.label[methylome], levels = ct.label)]
  chipA[, primary := match != "proxy"]
  chipA <- chipA[order(delta)]
  # a motif can appear in more than one methylome, so the row key carries both
  chipA[, row_key := factor(paste(motif, methylome, sep = "@"),
                            levels = paste(motif, methylome, sep = "@"))]
  strip_key <- function(x) sub("@.*$", "", x)

  log_info(
    "TF ChIP split: median ", round(median(chipA$delta, na.rm = TRUE), 3),
    " against its own methylome, ",
    round(median(chipA$delta_control, na.rm = TRUE), 3),
    " against another lineage; ",
    sum(chipA$specificity > 0, na.rm = TRUE), " of ",
    sum(is.finite(chipA$specificity)), " stronger where the factor is bound"
  )

  # G) coverage at motif level: a dot per motif x cell type where ChIP data exists
  ct.order <- c("Monocytes", "T-cells", "B-cells", "GM12878\n(cell line)")
  cov.ct <- c(monocyte = "Monocytes", Tcell = "T-cells", Bcell = "B-cells",
              GM12878 = "GM12878\n(cell line)")

  coverage <- chipA[, .(sites = max(n_supported, na.rm = TRUE)),
                    by = .(motif, factor, methylome)]
  coverage[, cell := factor(cov.ct[methylome], levels = ct.order)]

  # motifs grouped by their TF so variants sit together and line up with H
  m.order <- coverage[order(factor, motif), unique(motif)]
  grid <- CJ(motif = m.order, cell = factor(ct.order, levels = ct.order))
  coverage <- merge(grid, coverage, by = c("motif", "cell"), all.x = TRUE)
  coverage[, motif := factor(motif, levels = rev(m.order))]
  coverage[, tested := !is.na(sites)]

  only.line <- setdiff(
    coverage[cell == "GM12878\n(cell line)" & tested, as.character(motif)],
    coverage[cell != "GM12878\n(cell line)" & tested, as.character(motif)]
  )
  log_info(
    "TF ChIP coverage: ", uniqueN(coverage$motif), " motifs across four cell types, ",
    length(only.line), " tested only in GM12878"
  )

  p_cov <- ggplot(coverage, aes(x = cell, y = motif)) +
    geom_point(data = coverage[tested == FALSE], shape = 4, size = 1,
               colour = "grey82", stroke = 0.4) +
    geom_point(data = coverage[tested == TRUE], size = 2.2,
               colour = ink.mark, alpha = 0.85) +
    scale_x_discrete(position = "top") +
    labs(x = NULL, y = "ChIP motif") +
    base_theme +
    theme(
      panel.grid.major.y = element_line(linewidth = 0.2, colour = "grey94"),
      axis.text.y = element_text(size = base.size - 3),
      axis.text.x = element_text(size = base.size - 2, angle = 45, hjust = 0),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "bottom"
    )

  # H) the two deltas against each other, diagonal = no cell type specificity
  lim <- range(c(chipA$delta, chipA$delta_control), na.rm = TRUE)
  spec_lab <- chipA[primary == TRUE | specificity <= 0]

  p_spec <- ggplot(chipA, aes(x = delta_control, y = delta)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                linewidth = 0.3, colour = "grey60") +
    geom_point(aes(shape = primary), size = 2, colour = ink.mark) +
    geom_text_repel(data = spec_lab, aes(label = motif), size = 2.4,
                    colour = ink.muted, min.segment.length = 0,
                    segment.size = 0.2, box.padding = 0.3, max.overlaps = Inf) +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
                       labels = c(`TRUE` = "primary cells", `FALSE` = "cell line")) +
    coord_equal(xlim = lim, ylim = lim) +
    labs(
      x = "scored against another lineage",
      y = "scored against its own methylome", shape = NULL
    ) +
    base_theme +
    theme(legend.position = "bottom",
          plot.title = element_text(hjust = 0, face = "plain", size = base.size))

  # I) paired shift summarised, cell line kept apart from the primary cells
  chipA[, cells := factor(
    ifelse(primary, "Primary cells", "GM12878 (cell line)"),
    levels = c("GM12878 (cell line)", "Primary cells")
  )]
  chipA[, cells_n := paste0(cells, "\nn = ", .N), by = cells]
  chipA[, cells_n := factor(cells_n, levels = unique(cells_n[order(cells)]))]

  split_long <- melt(
    chipA[, .(row_key, cells_n, primary,
              `own methylome` = delta, `other lineage` = delta_control)],
    id.vars = c("row_key", "cells_n", "primary"),
    variable.name = "scored_against", value.name = "difference"
  )
  split_med <- split_long[, .(m = median(difference, na.rm = TRUE)),
                          by = .(cells_n, scored_against)]

  # boxes rather than paired lines: at 23 comparisons the lines cross into spaghetti
  p_split <- ggplot(split_long, aes(x = scored_against, y = difference)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey75") +
    geom_boxplot(aes(colour = scored_against), fill = NA, width = 0.5,
                 linewidth = 0.4, outlier.shape = NA) +
    geom_point(aes(colour = scored_against), size = 1.5, alpha = 0.7,
               position = position_jitter(width = 0.12, height = 0, seed = 42)) +
    geom_text(
      data = split_med, aes(y = m, label = sprintf("%.2f", m)),
      nudge_x = 0.38, hjust = 0, size = 2.6, colour = ink.muted
    ) +
    scale_x_discrete(expand = expansion(add = c(0.45, 0.8))) +
    facet_wrap(~cells_n, nrow = 1) +
    scale_colour_manual(values = c("own methylome" = col.matched,
                                   "other lineage" = col.control),
                        guide = "none") +
    labs(
      x = NULL, y = "ChIP supported minus\nCpG/GC matched unsupported sites"
    ) +
    base_theme +
    theme(plot.title = element_text(hjust = 0, face = "plain", size = base.size))

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
  list(items = list(list(p_c, 1), list(p_d, 1), list(p_h27cor, 1.5)), height = 1.1),
  list(items = list(list(p_h27, 1)), height = 1.3),
  list(items = list(list(p_cov, 1.2), list(p_spec, 1.2), list(p_split, 1.5)), height = 1.4)
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
  width = fig.width, height = fig.height,
  bg = "transparent", limitsize = FALSE
)
log_success("Wrote ", file)
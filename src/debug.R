#!/usr/bin/env Rscript

#####################################################################
# debug_estimators.R
# created on 07-09-2026 by Irem B Gunduz
# Does subtracting the expected methylation first give a better score?
#
# Three ways to turn one footprint into one number:
#   current   obs_centre/obs_flank - exp_centre/exp_flank   what methylTFR does
#   residual difference   d_centre - d_flank     with d = observed - expected
#   residual ratio        d_centre / d_flank
#
#   A  the recomputed current score against the package's own, which says
#      whether the centre and flank windows below are the right ones
#   B  current against the residual difference
#   C  current against the residual ratio
#   D  bias: how strongly each score still follows the composition
#      prediction across motifs. Lower looks better
#   E  reliability: each score computed twice, once from the left half of
#      the footprint and once from the right, then correlated. A score that
#      cannot reproduce itself cannot be trusted
#   F  bias corrected for reliability. A noisy score correlates weakly with
#      everything, composition included, so a low bias in D only counts if
#      E is high. This is the panel that decides it
#####################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(logger)
  library(GenomicRanges)
  library(methylTFR)
  library(methylTFRAnnotationHg38)
  library(SummarizedExperiment)
})
set.seed(42)

#####################################################################
# Settings
#####################################################################

motifSet <- "jaspar2020_distal"
tfSet <- "jaspar2020"

# Motifs are taken evenly along the expected score, so the composition axis
# panel D needs is covered rather than sampled from one end of it
n.motifs <- 30L
motif.groups <- c("Bcell_mem", "Tcell", "Granulocytes")

# The windows that define the footprint. methylTFR does not expose the ones
# it uses internally, so panel A checks these against its own scores before
# anything else is read into the comparison
centre.window <- 25L
flank.inner <- 150L

fig.width <- 12
fig.height <- 14
base.size <- 10

est_levels <- c("Current", "Residual difference", "Residual ratio")
est_colors <- c(
  "Current" = "#C2377C",
  "Residual difference" = "#2B4B9B",
  "Residual ratio" = "#B5A38A"
)
est_columns <- c(
  "Current" = "current",
  "Residual difference" = "residual_diff",
  "Residual ratio" = "residual_ratio"
)

analysis.dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/blueprint/"
rnb.tag <- "RnBeads_230826"
dev.tag <- "mTFR_devs_230826"
dev.file <- file.path(analysis.dir, dev.tag, paste0(motifSet, "_deviations.RDS"))
sannot.file <- file.path(analysis.dir, rnb.tag, "reports", "data_import_data", "annotation.csv")
cache.file5 <- file.path(
  analysis.dir, "debug", paste0("msites_merged_celltype5G_", rnb.tag, ".rds")
)
distal.file <- "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFRAnnotationHg38_old/inst/extdata/distal_regions.RDS"

github.dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR_manuscript/"
fig.dir <- file.path(github.dir, "figures", "blueprint")
table.dir <- file.path(github.dir, "tables")
if (!dir.exists(fig.dir)) dir.create(fig.dir, recursive = TRUE)

no_bg <- theme(
  plot.background = element_blank(),
  panel.background = element_blank(),
  legend.background = element_blank(),
  legend.box.background = element_blank(),
  legend.key = element_blank(),
  strip.background = element_blank()
)
tag_theme <- theme(plot.tag = element_text(face = "bold", size = 12))
base_theme <- theme_classic(base_size = base.size) +
  theme(plot.title = element_text(hjust = 0, face = "plain", size = base.size))

#####################################################################
# Inputs
#####################################################################

for (f in c(dev.file, cache.file5)) {
  if (!file.exists(f)) {
    stop("Missing input: ", f, ". Run 02 and then 09 so the cache exists")
  }
}

dev_obj <- readRDS(dev.file)
sannot <- read.csv(sannot.file, stringsAsFactors = FALSE)
sannot$bedFile <- as.character(sannot$bedFile)

healthy <- sannot$bedFile[which(!is.na(sannot$DISEASE) & sannot$DISEASE == "None")]
samples <- intersect(healthy, colnames(deviations(dev_obj)))
log_info(length(samples), " healthy samples")

package_dev <- rowMeans(deviations(dev_obj)[, samples, drop = FALSE], na.rm = TRUE)
expected_score <- rowMeans(
  SummarizedExperiment::assay(dev_obj, "expected")[, samples, drop = FALSE],
  na.rm = TRUE
)

tf_bindsites <- getTFbindsites(motifSet = tfSet)
gcfreqs <- getGCfreq(motifSet = motifSet)
gc_dist <- getGenomeGC()
enhancer <- if (motifSet == "jaspar2020_distal") readRDS(distal.file) else NULL

msites <- readRDS(cache.file5)
msites <- as.list(msites)
msites <- msites[intersect(motif.groups, names(msites))]
if (length(msites) == 0) {
  stop("None of ", paste(motif.groups, collapse = ", "), " is in the cache")
}
log_info("Groups: ", paste(names(msites), collapse = ", "))

candidates <- intersect(names(sort(expected_score)), names(tf_bindsites))
picked <- unique(candidates[round(seq(1, length(candidates), length.out = n.motifs))])
log_info(
  length(picked), " motifs spanning expected ",
  round(min(expected_score[picked]), 4), " to ",
  round(max(expected_score[picked]), 4)
)

#####################################################################
# The three scores, from one footprint each
#
# side selects the half of the profile used, which is how the reliability
# in panel E is measured: the same score twice from independent halves
#####################################################################

score_from_wide <- function(wide, side = "all") {
  sel <- switch(side,
    left = wide$x < 0,
    right = wide$x > 0,
    rep(TRUE, nrow(wide))
  )
  w <- wide[sel, , drop = FALSE]
  centre <- abs(w$x) <= centre.window
  flank <- abs(w$x) >= flank.inner
  if (sum(centre) < 3 || sum(flank) < 3) {
    return(NULL)
  }
  obs_c <- mean(w$obs[centre], na.rm = TRUE)
  obs_f <- mean(w$obs[flank], na.rm = TRUE)
  exp_c <- mean(w$exp[centre], na.rm = TRUE)
  exp_f <- mean(w$exp[flank], na.rm = TRUE)
  d <- w$obs - w$exp
  d_c <- mean(d[centre], na.rm = TRUE)
  d_f <- mean(d[flank], na.rm = TRUE)
  data.table(
    current = obs_c / obs_f - exp_c / exp_f,
    residual_diff = d_c - d_f,
    residual_ratio = d_c / d_f,
    flank_meth = obs_f,
    flank_residual = d_f
  )
}

profile_scores <- function(motif, group) {
  df <- tryCatch(
    plotExpectedFootprint(
      motif = motif, tf_bindsites = tf_bindsites, msites = msites[[group]],
      sample_name = group, gc_dist = gc_dist, gcfreqs = gcfreqs,
      enhancer = enhancer, returnPlotData = TRUE
    )$plotDF,
    error = function(e) {
      log_warn(motif, " / ", group, ": ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(df)) {
    return(NULL)
  }
  df <- as.data.table(df)
  wide <- merge(
    df[type == "Observed", .(x, obs = avg_methyl)],
    df[type == "Expected", .(x, exp = avg_methyl)],
    by = "x"
  )
  out <- rbindlist(lapply(c("all", "left", "right"), function(s) {
    r <- score_from_wide(wide, s)
    if (is.null(r)) {
      return(NULL)
    }
    r[, side := s]
    r
  }))
  if (nrow(out) == 0) {
    return(NULL)
  }
  out[, `:=`(motif = motif, group = group)]
  out
}

scores <- rbindlist(lapply(picked, function(m) {
  rbindlist(Filter(Negate(is.null), lapply(names(msites), function(g) {
    profile_scores(m, g)
  })))
}))
if (nrow(scores) == 0) stop("No motif produced a footprint")
log_info(uniqueN(scores$motif), " motifs scored across ", uniqueN(scores$group), " groups")

per_motif <- scores[, lapply(.SD, mean, na.rm = TRUE),
  by = .(motif, side),
  .SDcols = c("current", "residual_diff", "residual_ratio", "flank_meth", "flank_residual")
]
per_motif[, package := package_dev[motif]]
per_motif[, expected := expected_score[motif]]

whole <- per_motif[side == "all"]
write.csv(per_motif, file.path(table.dir, "debug_estimator_scores.csv"), row.names = FALSE)

#####################################################################
# A, are the windows right
#####################################################################

fit_lab <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 4) {
    return("too few points")
  }
  sprintf(
    "Pearson %.2f, Spearman %.2f",
    cor(x[ok], y[ok]), cor(x[ok], y[ok], method = "spearman")
  )
}

p_a <- ggplot(whole, aes(x = package, y = current)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = "grey55") +
  geom_point(size = 1.8, colour = "grey25") +
  labs(
    title = paste0(
      "Recomputed against the package's own score\n",
      fit_lab(whole$package, whole$current)
    ),
    x = "methylTFR deviation score", y = "Recomputed current score"
  ) +
  base_theme

#####################################################################
# B and C, the two alternatives against the current one
#####################################################################

alt_panel <- function(column, label) {
  dt <- whole[is.finite(get(column))]
  ggplot(dt, aes(x = current, y = .data[[column]])) +
    geom_hline(yintercept = 0, linetype = "dotted", colour = "grey60") +
    geom_vline(xintercept = 0, linetype = "dotted", colour = "grey60") +
    geom_point(size = 1.8, colour = est_colors[[label]]) +
    labs(
      title = paste0(label, " against current\n", fit_lab(dt$current, dt[[column]])),
      x = "Current score", y = label
    ) +
    base_theme
}

p_b <- alt_panel("residual_diff", "Residual difference")

ratio_dt <- whole[is.finite(residual_ratio)]
ratio_lim <- quantile(abs(ratio_dt$residual_ratio), 0.9, na.rm = TRUE) * 1.5
outside <- sum(abs(ratio_dt$residual_ratio) > ratio_lim, na.rm = TRUE)
log_info(
  "Residual ratio spans ", signif(min(ratio_dt$residual_ratio), 3), " to ",
  signif(max(ratio_dt$residual_ratio), 3), ", ", outside,
  " motif(s) outside the drawn range"
)

p_c <- alt_panel("residual_ratio", "Residual ratio") +
  coord_cartesian(ylim = c(-ratio_lim, ratio_lim)) +
  labs(caption = sprintf(
    "%d motif(s) outside the drawn range, the flank residual sits near zero",
    outside
  )) +
  theme(plot.caption = element_text(size = base.size - 2, hjust = 0))

#####################################################################
# D to F, bias, reliability, and the one against the other
#####################################################################

left <- per_motif[side == "left"]
right <- per_motif[side == "right"]
setkey(left, motif)
setkey(right, motif)

safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 4) NA_real_ else cor(x[ok], y[ok])
}

summary_tab <- rbindlist(lapply(est_levels, function(lab) {
  col <- est_columns[[lab]]
  rel <- safe_cor(left[[col]][match(right$motif, left$motif)], right[[col]])
  bias <- abs(safe_cor(whole[[col]], whole$expected))
  data.table(
    estimator = factor(lab, levels = est_levels),
    bias = bias,
    reliability = rel,
    # A score can only correlate with anything up to the square root of its
    # own reliability, so dividing by it puts the three on equal footing
    bias_corrected = if (is.na(rel) || rel <= 0) NA_real_ else bias / sqrt(rel)
  )
}))
print(summary_tab)
write.csv(summary_tab, file.path(table.dir, "debug_estimator_summary.csv"), row.names = FALSE)
log_info(
  "bias / reliability / corrected bias: ",
  paste(summary_tab$estimator,
    sprintf(
      "%.2f / %.2f / %s", summary_tab$bias, summary_tab$reliability,
      ifelse(is.na(summary_tab$bias_corrected), "NA",
        sprintf("%.2f", summary_tab$bias_corrected)
      )
    ),
    sep = " = ", collapse = ", "
  )
)

bar_panel <- function(column, title, ylab) {
  dt <- summary_tab[is.finite(get(column))]
  ggplot(dt, aes(x = estimator, y = .data[[column]], fill = estimator)) +
    geom_col(width = 0.6) +
    geom_text(aes(label = sprintf("%.2f", .data[[column]])), vjust = -0.4, size = 3) +
    scale_fill_manual(values = est_colors, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(title = title, x = NULL, y = ylab) +
    base_theme +
    theme(axis.text.x = element_text(angle = 20, hjust = 1))
}

p_d <- bar_panel(
  "bias", "Bias: how much of the composition prediction survives",
  "|r| with the expected deviation score"
)
p_e <- bar_panel(
  "reliability", "Reliability: the left half against the right half",
  "r between the two halves"
)
p_f <- bar_panel(
  "bias_corrected", "Bias once reliability is taken out",
  "|r| divided by the square root of reliability"
)

#####################################################################
# The figure
#####################################################################

fig <- wrap_plots(
  wrap_plots(p_a + labs(tag = "A"), p_b + labs(tag = "B"), nrow = 1),
  wrap_plots(p_c + labs(tag = "C"), p_d + labs(tag = "D"), nrow = 1),
  wrap_plots(p_e + labs(tag = "E"), p_f + labs(tag = "F"), nrow = 1),
  ncol = 1
) & tag_theme

file <- file.path(fig.dir, "debug_estimators.pdf")
ggsave(file, fig & no_bg,
  width = fig.width, height = fig.height, bg = "transparent", limitsize = FALSE
)
log_success("Wrote ", file)
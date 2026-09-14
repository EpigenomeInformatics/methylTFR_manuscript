#!/usr/bin/env Rscript

#####################################################################
# 14_bp_tf_chip_figures.R
# created on 10-09-2026 by Irem B Gunduz
#
# The detail behind panels G and H of the Blueprint supplementary, drawn as a
# supplementary of its own so neither figure has to carry all of it.
#
#   A  every comparison, scored against its own methylome and against a
#      methylome from another lineage, grouped by cell type
#   B  bound fraction against deviation score in GM12878
#   C  the covariate balance the matching achieved
#
# Needs 13_bp_tf_chip_validation.R to have run.
#####################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(logger)
})

base.size <- 9
col.matched <- "#C2377C"
col.control <- "#2B4B9B"
ink.muted <- "grey40"
ink.mark <- "grey25"

github.dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR_manuscript/"
tab.dir <- file.path(github.dir, "tables", "chip_validation")
fig.dir <- file.path(github.dir, "figures", "blueprint")
if (!dir.exists(fig.dir)) dir.create(fig.dir, recursive = TRUE)

base_theme <- theme_classic(base_size = base.size) +
  theme(
    plot.title = element_text(hjust = 0, face = "plain", size = base.size + 1),
    plot.subtitle = element_text(hjust = 0, size = base.size - 1, colour = ink.muted),
    plot.tag = element_text(face = "bold", size = 12),
    axis.line = element_line(linewidth = 0.3, colour = "grey60"),
    axis.ticks = element_line(linewidth = 0.3, colour = "grey60"),
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0, hjust = 0, size = base.size - 1),
    legend.key.size = unit(0.4, "cm")
  )

A <- fread(file.path(tab.dir, "A_site_level_split.csv"))
B <- fread(file.path(tab.dir, "B_motif_level.csv"))
Bcor <- fread(file.path(tab.dir, "B_motif_level_correlation.csv"))

ct.label <- c(GM12878 = "GM12878\n(B cell line)", monocyte = "Monocytes",
              Tcell = "T-cells", Bcell = "B-cells")
A[, ct := factor(ct.label[methylome], levels = ct.label)]
A[, primary := match != "proxy"]
A <- A[order(delta)]
# a motif can appear in more than one methylome, so the row key carries both
# and only the motif is printed
A[, row_key := factor(paste(motif, methylome, sep = "@"),
                      levels = paste(motif, methylome, sep = "@"))]
strip_key <- function(x) sub("@.*$", "", x)

#####################################################################
# A, every comparison
#####################################################################

Al <- melt(
  A[, .(row_key, ct, primary,
        `own methylome` = delta, `other lineage` = delta_control)],
  id.vars = c("row_key", "ct", "primary"),
  variable.name = "scored_against", value.name = "difference"
)

pA <- ggplot(Al, aes(x = difference, y = row_key)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey75") +
  geom_line(aes(group = row_key), colour = "grey80", linewidth = 0.4) +
  geom_point(aes(colour = scored_against, shape = primary), size = 1.9) +
  facet_grid(ct ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_y_discrete(labels = strip_key) +
  scale_colour_manual(values = c("own methylome" = col.matched,
                                 "other lineage" = col.control)) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
                     labels = c(`TRUE` = "primary cells", `FALSE` = "cell line")) +
  labs(
    x = "ChIP supported minus CpG/GC matched unsupported sites",
    y = NULL, colour = "scored against", shape = NULL,
    title = "Every comparison, one row per motif and methylome",
    subtitle = "the gap between the two points is what a shared open chromatin explanation cannot produce"
  ) +
  base_theme + theme(legend.position = "bottom")

#####################################################################
# B, motif level in the only cell type with enough factors
#####################################################################

Bg <- B[methylome == "GM12878"]
cg <- Bcor[methylome == "GM12878"]
lab.B <- unique(rbindlist(list(
  head(Bg[order(deviation)], 4), head(Bg[order(-deviation)], 3),
  head(Bg[order(-frac_supported)], 3)
)))

pB <- ggplot(Bg, aes(x = frac_supported, y = deviation)) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
              colour = col.control, fill = col.control,
              alpha = 0.10, linewidth = 0.5) +
  geom_point(size = 1.8, colour = ink.mark) +
  geom_text_repel(data = lab.B, aes(label = motif), size = 2.3,
                  colour = ink.muted, min.segment.length = 0,
                  segment.size = 0.2, box.padding = 0.35, max.overlaps = Inf) +
  labs(
    x = "fraction of distal motif sites with ChIP support",
    y = "methylTFR deviation score",
    title = "Motif level, GM12878",
    subtitle = sprintf(paste("Spearman rho = %.2f, p = %.3f, n = %d motifs.",
                             "Most are redundant AP-1 dimers scored against the",
                             "same few peak sets,\nso n overstates the evidence.",
                             "Extremes labelled."),
                       cg$rho, cg$p, cg$n_factors)
  ) +
  base_theme

#####################################################################
# C, did the matching work
#####################################################################

bal <- rbindlist(list(
  A[, .(motif, ct, covariate = "CpGs per 500 bp",
        supported = cpg_sup, matched = cpg_mat)],
  A[, .(motif, ct, covariate = "GC fraction",
        supported = gc_sup, matched = gc_mat)]
))

pC <- ggplot(bal, aes(x = matched, y = supported)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              linewidth = 0.3, colour = "grey60") +
  geom_point(size = 1.6, colour = ink.mark, alpha = 0.8) +
  facet_wrap(~covariate, scales = "free") +
  labs(
    x = "matched unsupported sites", y = "ChIP supported sites",
    title = "Covariate balance after matching",
    subtitle = "points on the diagonal mean the two groups are compositional twins"
  ) +
  base_theme

#####################################################################

p <- pA + (pB / pC) +
  plot_layout(widths = c(1.15, 1)) +
  plot_annotation(
    tag_levels = "A",
    caption = paste(
      "Peaks: ReMap2022, UniBind and ENCODE, hg38, consensus of >= 2 datasets where available, trimmed to +/- 250 bp.",
      "Sites: methylTFRAnnotationHg38 jaspar2020, distal only; supported = motif centre inside a peak.",
      "Background: unsupported sites matched 1:1 on CpG count (exact) and GC (+/- 0.02) in +/- 250 bp, without replacement.",
      "Methylomes: ENCODE WGBS, GRCh38, coverage >= 5. GM12878 shares a donor with its ChIP; the three primary methylomes are one donor (ENCDO661BYS).",
      "GM12878 supplies breadth across factors, the primary cells supply relevance. K562 was excluded as a leukaemia line with a distorted methylome.",
      "Expected profiles are rebuilt per site subset, since the stored per-motif table cannot distinguish subsets of one motif.",
      sep = "\n"),
    theme = theme(plot.caption = element_text(size = 6, hjust = 0, colour = ink.muted))
  )

out <- file.path(fig.dir, "tf_chip_validation_supplementary.pdf")
ggsave(out, p, width = 11.5, height = 8.5, useDingbats = FALSE)
log_success("Wrote ", out)

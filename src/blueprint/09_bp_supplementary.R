#!/usr/bin/env Rscript

#####################################################################
# 09_bp_supplementary.R
# created on 06-09-2026 by Irem B Gunduz
# The Blueprint supplementary figure
#   A  PCA of the bias corrected deviations, with the cell type R2
#   B  PCA of the uncorrected deviations, with the cell type R2
#   C  Motif variability across all motifs, top motifs annotated
#   D  Expected deviation scores, with their spread across samples
#   E  Where the variance of the expectation sits, motif against sample
#   F  PCA of the B and T cell subsets
#   G  Observed and expected footprints of two motifs, B and T cells
#   H  B versus T differential, one marker motif per lineage
#   I  PCA of the Altius deviations, with the cell type R2
#   J  Altius AP-1 and CEBP differentials across the cell type groups
#####################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(GenomicRanges)
  library(logger)
  library(RnBeads)
  library(methylTFR)
  library(methylTFRAnnotationHg38)
})
set.seed(42)

#####################################################################
# Settings
#####################################################################

motifSet.distal <- "jaspar2020_distal"
motifSet.genome <- "jaspar2020"
motifSet.altius <- "altius"
tfSet <- "jaspar2020"

drop.cell.types <- "Other"

pca.center <- FALSE
pca.scale <- FALSE
# On the raw uncentred deviations PC1 is the mean methylation profile and
# takes nearly all of the variance, which leaves the B and T subsets sitting
# on top of each other. Panel F standardises each motif across the subset
# first, so the components describe how the samples differ from one another.
subset.zscore <- TRUE
# Panels A and B get the same treatment as each other, so the only thing
# separating them is the correction itself. Standardising each motif also
# removes any part of the expectation that is constant for that motif, so
# what remains between the two panels is the part that varies by sample
pca.zscore <- TRUE
# Each panel names its own motif set and whether the correction is applied,
# so the two can differ. set is "distal" or "genome"
pca.panelA <- list(set = "distal", corrected = FALSE)
pca.panelB <- list(set = "genome", corrected = TRUE)
pca.panelI <- list(set = "altius", corrected = TRUE)
# The share of a component's variance that cell type accounts for, printed
# under each title. The two panels use different motif sets, so the numbers
# describe each panel rather than comparing them. Set FALSE to leave it out
pca.show.r2 <- TRUE

# Canvas. A4 squeezed the three column rows, so the figure is drawn larger
# and scaled down at layout time instead
fig.width <- 16
fig.height <- 24
base.size <- 10
ellipse.level <- 0.9

variability.motifSet <- motifSet.distal
variability.bootstrap <- TRUE
variability.iterations <- 1000L
variability.top <- 8
expected.label.n <- 6

footprint.motifs.EO <- c("FOSL1::JUND", "SPI1")
# Panel H, one motif per lineage
footprint.motif.B <- "BATF"
footprint.motif.T <- "PAX1"
# Panel J, the Altius archetypes drawn across the cell type groups of the
# main figure, with the naive B cells left out
altius.motifs <- c("ap1_1", "ccaat_cebp")
altius.group.order <- c("Monocytes", "Granulocytes", "Bcell_mem", "Tcell")
flank.norm <- 30

distal.file <- "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFRAnnotationHg38_old/inst/extdata/distal_regions.RDS"

# Directories
analysis.dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/blueprint/"
rnb.tag <- "RnBeads_230826"
dev.tag <- "mTFR_devs_230826"

rnb.set.path <- file.path(analysis.dir, rnb.tag, "reports", "data_import_data", "rnb.set_preprocessed")
sannot.file <- file.path(analysis.dir, rnb.tag, "reports", "data_import_data", "annotation.csv")
dev.file.distal <- file.path(analysis.dir, dev.tag, paste0(motifSet.distal, "_deviations.RDS"))
dev.file.genome <- file.path(analysis.dir, dev.tag, paste0(motifSet.genome, "_deviations.RDS"))
dev.file.altius <- file.path(analysis.dir, dev.tag, paste0(motifSet.altius, "_deviations.RDS"))

cache.dir <- file.path(analysis.dir, "debug")
if (!dir.exists(cache.dir)) dir.create(cache.dir, recursive = TRUE)
cache.file5 <- file.path(cache.dir, paste0("msites_merged_celltype5G_", rnb.tag, ".rds"))

github.dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR_manuscript/"
fig.dir <- file.path(github.dir, "figures", "blueprint")
table.dir <- file.path(github.dir, "tables")
if (!dir.exists(fig.dir)) dir.create(fig.dir, recursive = TRUE)
if (!dir.exists(table.dir)) dir.create(table.dir, recursive = TRUE)

# Display names and palette, as in 03 and 06
group_remap <- c(
  "Bcell" = "B-cells",
  "DC" = "Dendritic cells",
  "eryt" = "Erythrocytes",
  "gran" = "Granulocytes",
  "megK" = "Megakaryocytes",
  "Mf" = "Macrophages",
  "mono" = "Monocytes",
  "NK" = "NK-cells",
  "osteoclast" = "Osteoclast",
  "other" = "Other",
  "plasma" = "Plasma",
  "progenitor" = "Progenitors",
  "Tcell" = "T-cells",
  "thymocyte" = "Thymocyte"
)

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
cell_type_levels <- names(cell_type_colors)

# The B and T subsets of panel F, shaded from the lineage colour so the
# panel still reads as B against T at a glance
bcell_subsets <- c("Bcell_pre", "Bcell_naive", "Bcell_gc", "Bcell_mem")
tcell_subsets <- c(
  "TCD4", "TCD4_cenM", "TCD4_effM", "TCD8",
  "TCD8_cenM", "TCD8_effM", "TCD8_term", "Treg"
)
subset_colors <- c(
  setNames(
    colorRampPalette(c("#F2B8D2", "#7B1E3D"))(length(bcell_subsets)),
    bcell_subsets
  ),
  setNames(
    colorRampPalette(c("#BFEAF4", "#1F6F86"))(length(tcell_subsets)),
    tcell_subsets
  )
)

# Footprint groups. Observed keeps the saturated colour, expected a tint
eo_group_order <- c("Bcell_naive", "Bcell_mem", "Tcell")
eo_colors <- c(
  "Expected_Bcell_naive" = "#E5B9C8",
  "Observed_Bcell_naive" = "#C46B8C",
  "Expected_Bcell_mem" = "#C08099",
  "Observed_Bcell_mem" = "#7B1E3D",
  "Expected_Tcell" = "#B3E4EE",
  "Observed_Tcell" = "#4FC3D9"
)
# Panel H draws the same three groups as G, in the observed colours
foot_colors <- c(
  "Bcell_naive" = unname(eo_colors[["Observed_Bcell_naive"]]),
  "Bcell_mem" = unname(eo_colors[["Observed_Bcell_mem"]]),
  "Tcell" = unname(eo_colors[["Observed_Tcell"]])
)

# Panel J, the group colours the main figure uses
altius_colors <- c(
  "Monocytes" = "#C2703D",
  "Granulocytes" = "#E8A33D",
  "Bcell_mem" = "#7B1E3D",
  "Tcell" = "#4FC3D9"
)

# Where the variance of the expectation sits. The first is the part a per
# motif standardisation removes, the other two are what only the correction
# can take out, and the interaction is the larger of them
var_source_colors <- c(
  "Same in every sample" = "#B5A38A",
  "Motif by sample" = "#C2377C",
  "Whole sample shift" = "#E39BBE"
)

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
tag_theme <- theme(plot.tag = element_text(face = "bold", size = 10))

#####################################################################
# Helper functions
#####################################################################

row_zscore <- function(mat) {
  centred <- mat - rowMeans(mat, na.rm = TRUE)
  sds <- apply(mat, 1, sd, na.rm = TRUE)
  out <- centred / sds
  out[!is.finite(out)] <- 0
  out
}

clean_matrix <- function(mat, label) {
  mat <- as.matrix(mat)
  keep <- is.finite(rowSums(mat))
  if (any(!keep)) {
    log_info(label, ": dropping ", sum(!keep), " features with non finite values")
  }
  mat[keep, , drop = FALSE]
}

# Variance of one component that the grouping accounts for. Panels A and B
# differ by little to the eye once the scores are standardised, so the
# separation they achieve is reported as a number rather than left to it
component_r2 <- function(values, groups) {
  if (length(unique(groups)) < 2) {
    return(NA_real_)
  }
  ss <- summary(aov(v ~ g, data = data.frame(v = values, g = groups)))[[1]]
  ss["g", "Sum Sq"] / sum(ss[, "Sum Sq"])
}

pca_panel <- function(mat, groups, palette, title, label, zscore = FALSE,
                      ellipse = NULL, show.r2 = FALSE) {
  mat <- clean_matrix(mat, label)
  if (zscore) mat <- row_zscore(mat)
  log_info(
    label, ": PCA on ", nrow(mat), " features x ", ncol(mat), " samples",
    if (zscore) ", motif wise Z-scores" else ""
  )
  pca <- prcomp(t(mat), center = zscore || pca.center, scale. = pca.scale)
  pct <- 100 * pca$sdev^2 / sum(pca$sdev^2)

  present <- intersect(names(palette), unique(groups))
  df <- data.frame(
    PC1 = pca$x[, 1], PC2 = pca$x[, 2],
    group = factor(groups, levels = present)
  )

  r2 <- c(component_r2(df$PC1, df$group), component_r2(df$PC2, df$group))
  log_info(
    label, ": cell type R2 on PC1 ", round(r2[1], 3), ", PC2 ", round(r2[2], 3)
  )
  if (show.r2) {
    title <- sprintf(
      "%s\ncell type R2: PC1 %.2f, PC2 %.2f", title, r2[1], r2[2]
    )
  }

  p <- ggplot(df, aes(x = PC1, y = PC2, colour = group)) +
    geom_point(size = 1.6, alpha = 0.9) +
    scale_colour_manual(values = palette[present], name = "Cell type") +
    labs(
      title = title,
      x = paste0("PC1 (", format(round(pct[1], 2), nsmall = 2), "%)"),
      y = paste0("PC2 (", format(round(pct[2], 2), nsmall = 2), "%)")
    ) +
    theme_classic(base_size = base.size) +
    theme(
      plot.title = element_text(hjust = 0, face = "plain", size = base.size),
      legend.key.size = unit(3.5, "mm"),
      legend.text = element_text(size = base.size - 2),
      legend.title = element_text(size = base.size - 1)
    )

  if (!is.null(ellipse)) {
    df$lineage <- ellipse
    keep <- names(table(df$lineage))[table(df$lineage) >= 4]
    if (length(keep) > 0) {
      p <- p + stat_ellipse(
        data = df[df$lineage %in% keep, ], aes(group = lineage),
        type = "norm", level = ellipse.level, colour = "grey35",
        linetype = "dashed", linewidth = 0.4, show.legend = FALSE
      )
    }
  }
  p
}

map_cellType5Group <- function(cts) {
  cts <- as.character(cts)
  dplyr::case_when(
    cts %in% c("Bcell_mem", "Bcell_gc") ~ "Bcell_mem",
    cts %in% c("Bcell_naive", "Bcell_pre") ~ "Bcell_naive",
    cts %in% tcell_subsets ~ "Tcell",
    cts == "Mono" ~ "Monocytes",
    cts %in% c(
      "Neut_band", "Neut_mat", "Neut_meta",
      "Neut_myel", "Neut_segm", "Eosi_mat"
    ) ~ "Granulocytes",
    TRUE ~ NA_character_
  )
}

ensure_msites_cols <- function(gr, label) {
  if (!"score" %in% names(mcols(gr))) {
    alt <- intersect(c("beta", "meth", "methylation", "avg_methyl"), names(mcols(gr)))
    if (length(alt) == 0) {
      stop(label, ": no methylation column found, available: ",
        paste(names(mcols(gr)), collapse = ", "))
    }
    mcols(gr)$score <- as.numeric(mcols(gr)[[alt[1]]])
  }
  if (!"coverage" %in% names(mcols(gr))) {
    alt <- intersect(c("covg", "cov", "reads"), names(mcols(gr)))
    mcols(gr)$coverage <- if (length(alt) > 0) {
      as.numeric(mcols(gr)[[alt[1]]])
    } else {
      NA_real_
    }
  }
  gr
}

# Merged methylation per group, cached because it needs the whole RnBeads set
merged_msites <- function(cache.file, group_column, group_fun) {
  if (file.exists(cache.file)) {
    log_info("Reusing merged methylation from ", cache.file)
    msites <- readRDS(cache.file)
  } else {
    log_info("Loading RnBeads object from ", rnb.set.path)
    rnb_set <- RnBeads::load.rnb.set(rnb.set.path)
    grouping <- group_fun(pheno(rnb_set))
    rnb_set@pheno[[group_column]] <- grouping

    drop_idx <- unique(c(
      which(as.character(pheno(rnb_set)$DISEASE) != "None"),
      which(is.na(grouping))
    ))
    if (length(drop_idx) > 0) {
      log_info(length(drop_idx), " samples excluded for ", group_column)
      rnb_set <- RnBeads::remove.samples(rnb_set, drop_idx)
    }
    msites <- rnb.RnBSet.to.GRangesList(mergeSamples(rnb_set, group_column))
    saveRDS(msites, cache.file)
    log_info("Wrote ", cache.file)
  }
  msites <- as.list(msites)
  mapply(ensure_msites_cols, msites, names(msites), SIMPLIFY = FALSE)
}

# Mean deviation of one motif per group, on whichever matrix is passed
motif_group_means <- function(motif, mat, grps, groups) {
  idx <- which(rownames(mat) == motif)
  vapply(groups, function(g) {
    if (length(idx) == 0) {
      return(NA_real_)
    }
    vals <- mat[idx[1], grps == g]
    if (length(vals) == 0 || all(is.na(vals))) NA_real_ else mean(vals, na.rm = TRUE)
  }, numeric(1))
}

fmt <- function(x) ifelse(is.na(x), "NA", format(round(x, 2), nsmall = 2))

# The expected ratio sits within half a percent of 1, so two decimals show
# 1.00 for every group and hide the differences between them
fmt3 <- function(x) ifelse(is.na(x), "NA", format(round(x, 3), nsmall = 3))

footprint_data <- function(motif, msites, bindsites = tf_bindsites,
                          gcfreq = gcfreqs, enh = enhancer) {
  per_group <- lapply(names(msites), function(g) {
    df <- tryCatch(
      plotExpectedFootprint(
        motif = motif,
        tf_bindsites = bindsites,
        msites = msites[[g]],
        sample_name = g,
        gc_dist = gc_dist,
        gcfreqs = gcfreq,
        enhancer = enh,
        returnPlotData = TRUE
      )$plotDF,
      error = function(e) {
        log_warn(motif, " / ", g, ": ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(df)) {
      return(NULL)
    }
    df <- copy(df)
    df[, group := g]
    df
  })
  per_group <- Filter(Negate(is.null), per_group)
  if (length(per_group) == 0) {
    return(NULL)
  }
  rbindlist(per_group)
}

#####################################################################
# Annotation and deviations
#####################################################################

for (f in c(dev.file.distal, dev.file.genome, dev.file.altius)) {
  if (!file.exists(f)) stop("Missing input: ", f, ". Run 02 first")
}

sannot <- read.csv(sannot.file, stringsAsFactors = FALSE)
sannot$bedFile <- as.character(sannot$bedFile)

dev_sets <- list(
  distal = readRDS(dev.file.distal),
  genome = readRDS(dev.file.genome),
  altius = readRDS(dev.file.altius)
)
set_labels <- c(
  distal = "JASPAR2020 distal", genome = "JASPAR2020", altius = "Altius"
)

for (nm in names(dev_sets)) {
  if (!"expected" %in% SummarizedExperiment::assayNames(dev_sets[[nm]])) {
    stop("The ", nm, " deviations carry no expected assay, panel A cannot be drawn")
  }
}
uncorrected_sets <- lapply(dev_sets, function(d) {
  deviations(d) + SummarizedExperiment::assay(d, "expected")
})

dev_distal <- dev_sets$distal
dev_altius <- dev_sets$altius
uncorrected_distal <- uncorrected_sets$distal

# Healthy samples with a mappable cell type
ann_ok <- !is.na(sannot$DISEASE) & sannot$DISEASE == "None" &
  sannot$cellTypeGroup %in% names(group_remap)
mapped <- unname(group_remap[sannot$cellTypeGroup])
ann_ok[which(ann_ok & mapped %in% drop.cell.types)] <- FALSE

healthy <- sannot$bedFile[which(ann_ok)]
samples_sets <- lapply(dev_sets, function(d) {
  intersect(healthy, colnames(deviations(d)))
})
groups_sets <- lapply(samples_sets, function(s) {
  unname(group_remap[sannot$cellTypeGroup[match(s, sannot$bedFile)]])
})
samples_distal <- samples_sets$distal
samples_altius <- samples_sets$altius
log_info(
  "Samples after the annotation filter: ",
  paste(names(samples_sets), lengths(samples_sets), sep = " = ", collapse = ", ")
)

#####################################################################
# A, B and I, PCA of the uncorrected and the corrected deviations
#####################################################################

build_pca <- function(spec, tag) {
  mat <- if (spec$corrected) {
    deviations(dev_sets[[spec$set]])
  } else {
    uncorrected_sets[[spec$set]]
  }
  label <- paste0(
    if (spec$corrected) "Bias corrected deviations (" else "Uncorrected deviations (",
    set_labels[[spec$set]], ")"
  )
  pca_panel(
    mat[, samples_sets[[spec$set]], drop = FALSE], groups_sets[[spec$set]],
    cell_type_colors, label, tag, zscore = pca.zscore, show.r2 = pca.show.r2
  )
}

p_a <- build_pca(pca.panelA, "A")
p_b <- build_pca(pca.panelB, "B")
# I stands on its own row, so it carries its own key
p_i <- build_pca(pca.panelI, "I") +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(nrow = 4, byrow = TRUE))

#####################################################################
# C, motif variability
#####################################################################

var_set <- names(set_labels)[match(variability.motifSet, c(motifSet.distal, motifSet.genome, motifSet.altius))]
var_input <- dev_sets[[var_set]]
var_samples <- samples_sets[[var_set]]

var_res <- computeZScoreVariability(
  deviations(var_input)[, var_samples, drop = FALSE],
  bootstrap = variability.bootstrap,
  niterations = variability.iterations
)
var_res <- var_res[order(-var_res$variability), ]
var_res$rank <- seq_len(nrow(var_res))

write.csv(var_res,
  file.path(table.dir, paste0("variability_", variability.motifSet, "_blueprint.csv")),
  row.names = FALSE
)
log_info(
  nrow(var_res), " motifs ranked, top: ",
  paste(head(var_res$motifs, variability.top), collapse = ", ")
)

p_c <- ggplot(var_res, aes(x = rank, y = variability))
if (all(c("bootstrap_lower_bound", "bootstrap_upper_bound") %in% colnames(var_res))) {
  p_c <- p_c + geom_errorbar(
    aes(ymin = bootstrap_lower_bound, ymax = bootstrap_upper_bound),
    colour = "grey65", linewidth = 0.25, width = 0
  )
}
p_c <- p_c +
  geom_point(size = 0.7, colour = "grey15") +
  # The top motifs sit on top of each other at the left edge of the curve, so
  # the labels are fanned out into the empty right of the panel instead
  geom_text_repel(
    data = head(var_res, variability.top),
    aes(label = motifs), size = 2.8,
    nudge_x = 90, direction = "y", hjust = 0,
    min.segment.length = 0, segment.size = 0.3, segment.colour = "grey55",
    box.padding = 0.3, max.overlaps = Inf
  ) +
  labs(
    title = paste("Motif variability", variability.motifSet),
    x = "Sorted TFs", y = "Variability"
  ) +
  theme_classic(base_size = base.size) +
  theme(plot.title = element_text(hjust = 0, face = "plain", size = base.size))

#####################################################################
# D, the expected deviation scores
#
# What the GC model alone predicts for each motif, before any methylation
# is measured. The correction subtracts this, so the panel shows the size
# and the direction of what is being removed. The vertical range is the
# spread of the same motif's expectation across the samples
#####################################################################

exp_mat <- SummarizedExperiment::assay(dev_distal, "expected")[, samples_distal, drop = FALSE]
exp_v <- rowMeans(exp_mat, na.rm = TRUE)
exp_tab <- data.frame(
  motif = names(exp_v), expected = unname(exp_v), stringsAsFactors = FALSE
)
exp_tab <- exp_tab[order(exp_tab$expected), ]
exp_tab$rank <- seq_len(nrow(exp_tab))
exp_low <- head(exp_tab, expected.label.n)
exp_high <- tail(exp_tab, expected.label.n)
log_info(
  "Expected deviation scores run ", round(min(exp_tab$expected), 4), " to ",
  round(max(exp_tab$expected), 4)
)

p_expected <- ggplot(exp_tab, aes(x = rank, y = expected)) +
  geom_hline(yintercept = 1, linetype = "dotted", colour = "grey50") +
  geom_point(size = 0.7, colour = "grey15") +
  # Both ends of a monotone curve crowd their own corner and leave the middle
  # of the panel empty, so each set of labels is fanned into that space
  geom_text_repel(
    data = exp_low, aes(label = motif), size = 2.6,
    nudge_x = 120, direction = "y", hjust = 0,
    min.segment.length = 0, segment.size = 0.3, segment.colour = "grey55",
    box.padding = 0.3, max.overlaps = Inf
  ) +
  geom_text_repel(
    data = exp_high, aes(label = motif), size = 2.6,
    nudge_x = -120, direction = "y", hjust = 1,
    min.segment.length = 0, segment.size = 0.3, segment.colour = "grey55",
    box.padding = 0.3, max.overlaps = Inf
  ) +
  labs(
    title = paste("Expected deviation scores", motifSet.distal),
    x = "Sorted TFs", y = "Expected deviation score\n(1 = no composition drift)"
  ) +
  theme_classic(base_size = base.size) +
  theme(plot.title = element_text(hjust = 0, face = "plain", size = base.size))

#####################################################################
# E, where the variance of the expectation sits
#
# The expected score stays close to 1 for every motif, which invites the
# question of why it is subtracted at all. It is computed per motif and per
# sample. Standardising the scores removes only the motif part, so the
# sample part is the component the correction is there for
#####################################################################

exp_grand <- mean(exp_mat, na.rm = TRUE)
ss_total <- sum((exp_mat - exp_grand)^2, na.rm = TRUE)
ss_motif <- ncol(exp_mat) * sum((rowMeans(exp_mat, na.rm = TRUE) - exp_grand)^2, na.rm = TRUE)
ss_sample <- nrow(exp_mat) * sum((colMeans(exp_mat, na.rm = TRUE) - exp_grand)^2, na.rm = TRUE)

# Three plain bars. Stacking them was what overprinted the labels, and the
# split is easier to read side by side than as segments of one bar
var_tab <- data.frame(
  source = factor(names(var_source_colors), levels = rev(names(var_source_colors))),
  frac = c(
    ss_motif,
    max(ss_total - ss_motif - ss_sample, 0),
    ss_sample
  ) / ss_total
)

sd_within <- median(apply(exp_mat, 1, sd, na.rm = TRUE), na.rm = TRUE)
sd_dev <- median(
  apply(deviations(dev_distal)[, samples_distal, drop = FALSE], 1, sd, na.rm = TRUE),
  na.rm = TRUE
)
log_info(
  "Expected score variance: ",
  paste(as.character(var_tab$source), paste0(round(100 * var_tab$frac, 1), "%"),
    sep = " = ", collapse = ", "
  ),
  ". Everything but the first survives a per motif standardisation"
)
log_info(
  "Median s.d. across samples within a motif: expected ", signif(sd_within, 3),
  ", corrected deviation ", signif(sd_dev, 3),
  ", ratio ", round(sd_within / sd_dev, 3)
)
write.csv(var_tab, file.path(table.dir, "bp_expected_variance_decomposition.csv"),
  row.names = FALSE
)

p_expvar <- ggplot(var_tab, aes(x = frac, y = source, fill = source)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%", 100 * frac)),
    hjust = -0.15, size = 3
  ) +
  scale_fill_manual(values = var_source_colors, guide = "none") +
  # Neither bar passes 60%, so the scale stops there rather than running to
  # 100 and leaving most of the panel empty
  scale_x_continuous(
    breaks = seq(0, 0.6, 0.2),
    labels = function(x) paste0(round(100 * x), "%"),
    expand = expansion(mult = c(0, 0.2))
  ) +
  labs(title = "Expected score variance", x = NULL, y = NULL) +
  theme_classic(base_size = base.size) +
  theme(
    plot.title = element_text(hjust = 0, face = "plain", size = base.size),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  )

#####################################################################
# F, PCA of the B and T cell subsets
#####################################################################

subset_short <- sannot$cellTypeShort[match(samples_distal, sannot$bedFile)]
keep_bt <- which(subset_short %in% names(subset_colors))
if (length(keep_bt) < 3) {
  stop("Fewer than three B or T cell samples carry a cellTypeShort label")
}
log_info(
  length(keep_bt), " B and T cell samples across ",
  length(unique(subset_short[keep_bt])), " subsets"
)

p_d <- pca_panel(
  deviations(dev_distal)[, samples_distal[keep_bt], drop = FALSE],
  subset_short[keep_bt], subset_colors,
  "B and T cell subsets (JASPAR2020 distal)", "F",
  zscore = subset.zscore,
  ellipse = ifelse(subset_short[keep_bt] %in% bcell_subsets, "B-cells", "T-cells")
) +
  labs(colour = "Subset") +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(nrow = 4, byrow = TRUE))

#####################################################################
# Motif annotation for the footprints
#####################################################################

tf_bindsites <- getTFbindsites(motifSet = tfSet)
gcfreqs <- getGCfreq(motifSet = motifSet.distal)
gc_dist <- getGenomeGC()

# The Altius run is genome wide, so it takes no distal regions
tf_bindsites_altius <- getTFbindsites(motifSet = motifSet.altius)
gcfreqs_altius <- getGCfreq(motifSet = motifSet.altius)
absent_altius <- setdiff(altius.motifs, names(tf_bindsites_altius))
if (length(absent_altius) > 0) {
  stop("Not in ", motifSet.altius, ": ", paste(absent_altius, collapse = ", "))
}

enhancer <- NULL
if (motifSet.distal == "jaspar2020_distal") {
  if (!file.exists(distal.file)) stop("Distal regions file does not exist: ", distal.file)
  enhancer <- readRDS(distal.file)
  log_info("Loaded ", length(enhancer), " distal regions")
}

#####################################################################
# H, the two marker motifs
#####################################################################

diff_motifs <- c(footprint.motif.B, footprint.motif.T)
missing_motifs <- setdiff(diff_motifs, names(tf_bindsites))
if (length(missing_motifs) > 0) {
  stop("Motif has no binding sites: ", paste(missing_motifs, collapse = ", "))
}
log_info("Panel H motifs: ", paste(diff_motifs, collapse = ", "))

#####################################################################
# Merged methylation for the footprints
#####################################################################

msites5 <- merged_msites(cache.file5, "cellType5Group", function(ph) {
  map_cellType5Group(ph$cellTypeShort)
})
msites_eo <- msites5[intersect(eo_group_order, names(msites5))]
if (length(msites_eo) == 0) stop("None of the B or T groups are in the merged methylation")

#####################################################################
# G, observed and expected footprints
#####################################################################

# The uncorrected score is the deviation with the GC expectation added back,
# the corrected one is the deviation itself
short5 <- map_cellType5Group(sannot$cellTypeShort[match(samples_distal, sannot$bedFile)])

eo_panel <- function(motif) {
  df <- footprint_data(motif, msites_eo)
  if (is.null(df)) {
    log_warn(motif, ": no footprint, panel G entry skipped")
    return(NULL)
  }
  df[, key := paste(type, group, sep = "_")]
  present <- intersect(names(eo_colors), unique(df$key))

  groups_here <- names(msites_eo)
  corrected <- motif_group_means(
    motif, deviations(dev_distal)[, samples_distal, drop = FALSE], short5, groups_here
  )
  uncorrected <- motif_group_means(
    motif, uncorrected_distal[, samples_distal, drop = FALSE], short5, groups_here
  )
  expected <- motif_group_means(
    motif,
    SummarizedExperiment::assay(dev_distal, "expected")[, samples_distal, drop = FALSE],
    short5, groups_here
  )

  # The two deviation scores sit on the observed key of each group, so they
  # are read off the colour legend rather than a block in the corner
  key_label <- vapply(present, function(k) {
    g <- sub("^(Observed|Expected)_", "", k)
    if (startsWith(k, "Expected")) {
      return(paste0(g, " exp (", fmt3(expected[[g]]), ")"))
    }
    paste0(g, " obs (", fmt(uncorrected[[g]]), ", dev ", fmt(corrected[[g]]), ")")
  }, character(1))

  df[, key := factor(key_label[match(key, present)], levels = key_label)]
  colours_here <- setNames(unname(eo_colors[present]), key_label)

  ggplot(df, aes(x = x, y = avg_methyl, colour = key)) +
    geom_line(linewidth = 0.4) +
    scale_colour_manual(values = colours_here, name = "Centre / flank methylation ratio") +
    coord_cartesian(xlim = c(-200, 200)) +
    labs(
      title = paste0(motif, " (observed vs expected)"),
      x = "Distance from motif center", y = "Average methylation"
    ) +
    theme_classic(base_size = base.size) +
    theme(
      legend.position = "bottom",
      legend.key.size = unit(3, "mm"),
      legend.text = element_text(size = base.size - 3),
      legend.title = element_text(size = base.size - 3),
      plot.title = element_text(hjust = 0, face = "plain", size = base.size)
    ) +
    guides(colour = guide_legend(ncol = 1))
}

p_e <- Filter(Negate(is.null), lapply(footprint.motifs.EO, eo_panel))

#####################################################################
# H, observed minus expected footprints of the two marker motifs
#####################################################################

diff_panel <- function(motif, msites = msites_eo, colours = foot_colors,
                       dev_mat = deviations(dev_distal)[, samples_distal, drop = FALSE],
                       dev_groups = short5, bindsites = tf_bindsites,
                       gcfreq = gcfreqs, enh = enhancer) {
  df <- footprint_data(motif, msites, bindsites, gcfreq, enh)
  if (is.null(df)) {
    log_warn(motif, ": no footprint, skipped")
    return(NULL)
  }
  df <- df[, .(
    avg_methyl = avg_methyl[type == "Observed"] - avg_methyl[type == "Expected"]
  ), by = .(x, group)]

  # The flanks are the baseline, so the curves start from zero away from the motif
  flank <- max(abs(df$x), na.rm = TRUE)
  df[, avg_methyl := avg_methyl -
    mean(avg_methyl[abs(x) >= flank - flank.norm], na.rm = TRUE), by = group]

  groups_here <- names(msites)
  raw <- motif_group_means(motif, dev_mat, dev_groups, groups_here)

  present <- intersect(names(colours), unique(df$group))
  group_label <- paste0(present, " (dev ", fmt(raw[present]), ")")
  df[, group := factor(group_label[match(group, present)], levels = group_label)]
  colours_here <- setNames(unname(colours[present]), group_label)

  ggplot(df, aes(x = x, y = avg_methyl, colour = group)) +
    geom_line(linewidth = 0.4) +
    scale_colour_manual(values = colours_here, name = "Cell type") +
    coord_cartesian(xlim = c(-200, 200)) +
    labs(
      title = motif,
      x = "Distance from motif center",
      y = "Methylation difference (Observed - Expected)"
    ) +
    theme_classic(base_size = base.size) +
    theme(
      legend.position = "bottom",
      legend.key.size = unit(3, "mm"),
      legend.text = element_text(size = base.size - 2),
      legend.title = element_text(size = base.size - 1),
      plot.title = element_text(hjust = 0, face = "plain", size = base.size)
    )
}

p_f <- Filter(Negate(is.null), lapply(diff_motifs, diff_panel))

#####################################################################
# J, the Altius archetypes across the cell type groups
#####################################################################

msites_altius <- msites5[intersect(altius.group.order, names(msites5))]
if (length(msites_altius) == 0) stop("None of the panel J groups are in the merged methylation")

short5_altius <- map_cellType5Group(
  sannot$cellTypeShort[match(samples_altius, sannot$bedFile)]
)

p_j <- Filter(Negate(is.null), lapply(altius.motifs, function(m) {
  diff_panel(m,
    msites = msites_altius, colours = altius_colors,
    dev_mat = deviations(dev_altius)[, samples_altius, drop = FALSE],
    dev_groups = short5_altius, bindsites = tf_bindsites_altius,
    gcfreq = gcfreqs_altius, enh = NULL
  )
}))

#####################################################################
# The assembled supplementary figure
#####################################################################

if (length(p_e) == 0 || length(p_f) == 0 || length(p_j) == 0) {
  stop("Panel G, H or J produced no footprint, the supplementary figure is not written")
}

p_e[[1]] <- p_e[[1]] + labs(tag = "G")
p_f[[1]] <- p_f[[1]] + labs(tag = "H")
p_j[[1]] <- p_j[[1]] + labs(tag = "J")

# A and B carry the same key, so it is collected once and sits under the two
# of them rather than beside the variability panel
row_ab <- wrap_plots(
  p_a + labs(tag = "A"), p_b + labs(tag = "B"),
  nrow = 1, guides = "collect"
) &
  theme(legend.position = "bottom") &
  guides(colour = guide_legend(nrow = 3, byrow = TRUE))

row1 <- wrap_plots(row_ab, p_c + labs(tag = "C"),
  wrap_plots(p_expected + labs(tag = "D"), p_expvar + labs(tag = "E"),
    ncol = 1, heights = c(2, 1)
  ),
  nrow = 1, widths = c(2, 1, 1.15)
)
row2 <- wrap_plots(c(list(p_d + labs(tag = "F")), p_e), nrow = 1, widths = c(1, 1.2, 1.2))
row3 <- wrap_plots(p_f, nrow = 1)
row4 <- wrap_plots(c(list(p_i + labs(tag = "I")), p_j), nrow = 1, widths = c(1, 1.2, 1.2))

supplementary <- wrap_plots(row1, row2, row3, row4,
  ncol = 1, heights = c(1, 1.5, 1.15, 1.5)
) &
  tag_theme

file <- file.path(fig.dir, "blueprint_supplementary.pdf")
ggsave(file, supplementary & no_bg,
  width = fig.width, height = fig.height, bg = "transparent", limitsize = FALSE
)
log_success("Wrote ", file)
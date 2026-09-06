#!/usr/bin/env Rscript

#####################################################################
# 09_bp_supplementary.R
# created on 06-09-2026 by Irem B Gunduz
# The Blueprint supplementary figure
#   A  PCA of the uncorrected deviations, JASPAR2020
#   B  PCA of the bias corrected deviations, JASPAR2020
#   C  Motif variability across all motifs, top motifs annotated
#   D  PCA of the B and T cell subsets
#   E  Observed and expected footprints of two motifs, B and T cells
#   F  B versus T differential, one marker motif per lineage
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
tfSet <- "jaspar2020"

drop.cell.types <- "Other"

pca.center <- FALSE
pca.scale <- FALSE
# On the raw uncentred deviations PC1 is the mean methylation profile and
# takes nearly all of the variance, which leaves the B and T subsets sitting
# on top of each other. Panel D standardises each motif across the subset
# first, so the components describe how the samples differ from one another.
subset.zscore <- TRUE

# Canvas. A4 squeezed the three column rows, so the figure is drawn larger
# and scaled down at layout time instead
fig.width <- 14
fig.height <- 18
base.size <- 10
ellipse.level <- 0.9

variability.motifSet <- motifSet.distal
variability.bootstrap <- TRUE
variability.iterations <- 1000L
variability.top <- 8

footprint.motifs.EO <- c("FOSL1::JUND", "SPI1")
# Panel F, one motif per lineage
footprint.motif.B <- "BATF"
footprint.motif.T <- "PAX1"
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

# The B and T subsets of panel D, shaded from the lineage colour so the
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
# Panel F draws the same three groups as E, in the observed colours
foot_colors <- c(
  "Bcell_naive" = unname(eo_colors[["Observed_Bcell_naive"]]),
  "Bcell_mem" = unname(eo_colors[["Observed_Bcell_mem"]]),
  "Tcell" = unname(eo_colors[["Observed_Tcell"]])
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

pca_panel <- function(mat, groups, palette, title, label, zscore = FALSE,
                      ellipse = NULL) {
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

footprint_data <- function(motif, msites) {
  per_group <- lapply(names(msites), function(g) {
    df <- tryCatch(
      plotExpectedFootprint(
        motif = motif,
        tf_bindsites = tf_bindsites,
        msites = msites[[g]],
        sample_name = g,
        gc_dist = gc_dist,
        gcfreqs = gcfreqs,
        enhancer = enhancer,
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

for (f in c(dev.file.distal, dev.file.genome)) {
  if (!file.exists(f)) stop("Missing input: ", f, ". Run 02 first")
}

sannot <- read.csv(sannot.file, stringsAsFactors = FALSE)
sannot$bedFile <- as.character(sannot$bedFile)

dev_distal <- readRDS(dev.file.distal)
dev_genome <- readRDS(dev.file.genome)

if (!"expected" %in% SummarizedExperiment::assayNames(dev_genome)) {
  stop("The genome wide deviations carry no expected assay, panel A cannot be drawn")
}
uncorrected_genome <- deviations(dev_genome) +
  SummarizedExperiment::assay(dev_genome, "expected")

# Healthy samples with a mappable cell type
ann_ok <- !is.na(sannot$DISEASE) & sannot$DISEASE == "None" &
  sannot$cellTypeGroup %in% names(group_remap)
mapped <- unname(group_remap[sannot$cellTypeGroup])
ann_ok[which(ann_ok & mapped %in% drop.cell.types)] <- FALSE

healthy <- sannot$bedFile[which(ann_ok)]
samples_genome <- intersect(healthy, colnames(deviations(dev_genome)))
samples_distal <- intersect(healthy, colnames(deviations(dev_distal)))
log_info(
  length(samples_genome), " genome wide and ", length(samples_distal),
  " distal samples after the annotation filter"
)

groups_genome <- unname(group_remap[
  sannot$cellTypeGroup[match(samples_genome, sannot$bedFile)]
])
groups_distal <- unname(group_remap[
  sannot$cellTypeGroup[match(samples_distal, sannot$bedFile)]
])

#####################################################################
# A and B, PCA of the uncorrected and the corrected deviations
#####################################################################

p_a <- pca_panel(
  uncorrected_genome[, samples_genome, drop = FALSE], groups_genome,
  cell_type_colors, "Uncorrected deviations (JASPAR2020)", "A"
)

p_b <- pca_panel(
  deviations(dev_genome)[, samples_genome, drop = FALSE], groups_genome,
  cell_type_colors, "Bias corrected deviations (JASPAR2020)", "B"
)

#####################################################################
# C, motif variability
#####################################################################

var_input <- if (variability.motifSet == motifSet.distal) dev_distal else dev_genome
var_samples <- if (variability.motifSet == motifSet.distal) samples_distal else samples_genome

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
  geom_text_repel(
    data = head(var_res, variability.top),
    aes(label = motifs), size = 2, min.segment.length = 0,
    segment.size = 0.3, box.padding = 0.4, max.overlaps = Inf
  ) +
  labs(
    title = paste("Motif variability", variability.motifSet),
    x = "Sorted TFs", y = "Variability"
  ) +
  theme_classic(base_size = base.size) +
  theme(plot.title = element_text(hjust = 0, face = "plain", size = base.size))

#####################################################################
# D, PCA of the B and T cell subsets
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
  "B and T cell subsets (JASPAR2020 distal)", "D",
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

enhancer <- NULL
if (motifSet.distal == "jaspar2020_distal") {
  if (!file.exists(distal.file)) stop("Distal regions file does not exist: ", distal.file)
  enhancer <- readRDS(distal.file)
  log_info("Loaded ", length(enhancer), " distal regions")
}

#####################################################################
# F, the two marker motifs
#####################################################################

diff_motifs <- c(footprint.motif.B, footprint.motif.T)
missing_motifs <- setdiff(diff_motifs, names(tf_bindsites))
if (length(missing_motifs) > 0) {
  stop("Motif has no binding sites: ", paste(missing_motifs, collapse = ", "))
}
log_info("Panel F motifs: ", paste(diff_motifs, collapse = ", "))

#####################################################################
# Merged methylation for the footprints
#####################################################################

msites5 <- merged_msites(cache.file5, "cellType5Group", function(ph) {
  map_cellType5Group(ph$cellTypeShort)
})
msites_eo <- msites5[intersect(eo_group_order, names(msites5))]
if (length(msites_eo) == 0) stop("None of the B or T groups are in the merged methylation")

#####################################################################
# E, observed and expected footprints
#####################################################################

# The uncorrected score is the deviation with the GC expectation added back,
# the corrected one is the deviation itself
short5 <- map_cellType5Group(sannot$cellTypeShort[match(samples_distal, sannot$bedFile)])
uncorrected_distal <- deviations(dev_distal) +
  SummarizedExperiment::assay(dev_distal, "expected")

eo_panel <- function(motif) {
  df <- footprint_data(motif, msites_eo)
  if (is.null(df)) {
    log_warn(motif, ": no footprint, panel E entry skipped")
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

  # The two deviation scores sit on the observed key of each group, so they
  # are read off the colour legend rather than a block in the corner
  key_label <- vapply(present, function(k) {
    g <- sub("^(Observed|Expected)_", "", k)
    if (startsWith(k, "Expected")) {
      return(paste0(g, " exp"))
    }
    paste0(g, " obs (", fmt(uncorrected[[g]]), " / ", fmt(corrected[[g]]), ")")
  }, character(1))

  df[, key := factor(key_label[match(key, present)], levels = key_label)]
  colours_here <- setNames(unname(eo_colors[present]), key_label)

  ggplot(df, aes(x = x, y = avg_methyl, colour = key)) +
    geom_line(linewidth = 0.4) +
    scale_colour_manual(values = colours_here, name = "obs (uncorrected / corrected)") +
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
# F, observed minus expected footprints of the two marker motifs
#####################################################################

diff_panel <- function(motif) {
  df <- footprint_data(motif, msites_eo)
  if (is.null(df)) {
    log_warn(motif, ": no footprint, panel F entry skipped")
    return(NULL)
  }
  df <- df[, .(
    avg_methyl = avg_methyl[type == "Observed"] - avg_methyl[type == "Expected"]
  ), by = .(x, group)]

  # The flanks are the baseline, so the curves start from zero away from the motif
  flank <- max(abs(df$x), na.rm = TRUE)
  df[, avg_methyl := avg_methyl -
    mean(avg_methyl[abs(x) >= flank - flank.norm], na.rm = TRUE), by = group]

  groups_here <- names(msites_eo)
  raw <- motif_group_means(
    motif, deviations(dev_distal)[, samples_distal, drop = FALSE], short5, groups_here
  )

  present <- intersect(names(foot_colors), unique(df$group))
  group_label <- paste0(present, " (dev ", fmt(raw[present]), ")")
  df[, group := factor(group_label[match(group, present)], levels = group_label)]
  colours_here <- setNames(unname(foot_colors[present]), group_label)

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
# The assembled supplementary figure
#####################################################################

if (length(p_e) == 0 || length(p_f) == 0) {
  stop("Panel E or F produced no footprint, the supplementary figure is not written")
}

p_e[[1]] <- p_e[[1]] + labs(tag = "E")
p_f[[1]] <- p_f[[1]] + labs(tag = "F")

# A and B carry the same key, so it is collected once and sits under the two
# of them rather than beside the variability panel
row_ab <- wrap_plots(
  p_a + labs(tag = "A"), p_b + labs(tag = "B"),
  nrow = 1, guides = "collect"
) &
  theme(legend.position = "bottom") &
  guides(colour = guide_legend(nrow = 3, byrow = TRUE))

row1 <- wrap_plots(row_ab, p_c + labs(tag = "C"), nrow = 1, widths = c(2, 1))
row2 <- wrap_plots(c(list(p_d + labs(tag = "D")), p_e), nrow = 1, widths = c(1, 1.2, 1.2))
row3 <- wrap_plots(p_f, nrow = 1)

supplementary <- wrap_plots(row1, row2, row3, ncol = 1, heights = c(1, 1.5, 1.15)) &
  tag_theme

file <- file.path(fig.dir, "blueprint_supplementary.pdf")
ggsave(file, supplementary & no_bg,
  width = fig.width, height = fig.height, bg = "transparent", limitsize = FALSE
)
log_success("Wrote ", file)

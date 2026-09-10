#!/usr/bin/env Rscript

#####################################################################
# 07_memT_supplementary.R
# created on 06-09-2026 by Irem B Gunduz
# The memory T cell supplementary figure
#   A  LOLA motif enrichment volcanoes, one per contrast
#   B  Observed minus expected footprints of four motifs
#
# No MA plots: they need the RnBDiffMeth object relinked, which takes minutes
# and a great deal of memory for a panel that says little.
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
set.seed(13)

src.dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR_manuscript/src"
source(file.path(src.dir, "utils.R"))
source(file.path(src.dir, "lola_utils.R"))

#####################################################################
# Settings
#####################################################################

motifSet <- "jaspar2020_distal"
tfSet <- "jaspar2020"

comparisons <- list(
  EM = list(test = "TEM", ref = "TN"),
  CM = list(test = "TCM", ref = "TN")
)
drop.cell.types <- "TEMRA"
group_order <- c("TN", "TCM", "TEM")

# Panel C
footprint.motifs <- c("JUN", "FOSL2", "BATF", "SPIB")
flank.norm <- 50

# Panel A. The rank cut is chosen per comparison, as in 06
diffmeth.rank.cut <- 1000
diffmeth.auto.rank.cut <- TRUE

# Panel B, matching 06 so the two figures name the same motifs
top.label.lola <- 10
lola.label.motifs <- c(
  "BATF::JUN", "FOS", "FOS::JUN", "FOS::JUNB", "FOS::JUND",
  "FOSB::JUNB", "FOSL1::JUN", "FOSL1::JUNB", "FOSL1::JUND",
  "FOSL2::JUND", "JUND", "JUN(var.2)",
  "ELK4", "ERF", "ERG", "ETS1", "ETV1", "ETV2", "ETV6",
  "EWSR1-FLI1", "FLI1", "GATA1::TAL1", "LEF1", "RUNX2", "RUNX3",
  "SPI1", "SPIB", "TCF7L1", "TCF7L2", "ZBTB7A"
)
lola.comparison <- c(EM = NA, CM = NA)
lola.region <- c("tiling1kb", "distal", "tiling", "sites")
lola.userSets <- c("rankCut_1000_hyper", "rankCut_1000_hypo")

cell_type_colors <- CELL_TYPE_COLORS

# Canvas
fig.width <- 14
fig.height <- 16
base.size <- 10

distal.file <- "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFRAnnotationHg38_old/inst/extdata/distal_regions.RDS"

# Directories
analysis.dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/memoryTcells"
dev.tag <- "mTFR_devs_230826"
dev.file <- file.path(analysis.dir, dev.tag, paste0(motifSet, "_deviations.RDS"))
rnb.set.path <- file.path(analysis.dir, "reports", "data_import_data", "rnb.set_preprocessed")
diffmeth.dir <- file.path(
  analysis.dir, "reports", "differential_methylation_data",
  "differential_rnbDiffMeth"
)
lola.file <- file.path(diffmeth.dir, "TF_motifs_lola.rds")

cache.dir <- file.path(analysis.dir, "debug")
if (!dir.exists(cache.dir)) dir.create(cache.dir, recursive = TRUE)
cache.file <- file.path(cache.dir, "msites_merged_cellType_230826.rds")

github.dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR_manuscript"
plot.dir <- file.path(github.dir, "figures", "memtcells", "diff_TFs_230826")
if (!dir.exists(plot.dir)) dir.create(plot.dir, recursive = TRUE)

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
tag_theme <- theme(plot.tag = element_text(face = "bold", size = 14))

#####################################################################
# Helper functions
#####################################################################

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

# Deviations split by the cellType column of the deviations object
load_cell_type_deviations <- function(file) {
  if (!file.exists(file)) {
    log_warn("No deviations at ", file, ", the legend will omit the scores")
    return(NULL)
  }
  obj <- readRDS(file)
  cd <- as.data.frame(colData(obj), stringsAsFactors = FALSE)
  if (!"cellType" %in% colnames(cd)) {
    log_warn("The deviations object carries no cellType column")
    return(NULL)
  }
  mat <- deviations(obj)
  split(seq_len(ncol(mat)), as.character(cd$cellType)) |>
    lapply(function(idx) mat[, idx, drop = FALSE])
}

# Exact match first, then a fixed substring, because JUN(var.2) would
# otherwise be read as a regular expression
motif_mean_deviations <- function(motif, dev_list, groups) {
  vapply(groups, function(g) {
    mat <- dev_list[[g]]
    if (is.null(mat)) {
      return(NA_real_)
    }
    idx <- which(rownames(mat) == motif)
    if (length(idx) == 0) {
      log_warn(motif, " is not a row of the deviations, no score for ", g)
      return(NA_real_)
    }
    vals <- mat[idx[1], ]
    if (all(is.na(vals))) NA_real_ else mean(vals, na.rm = TRUE)
  }, numeric(1))
}

make_label <- function(group, dev_val) {
  if (is.na(dev_val)) {
    return(group)
  }
  paste0(group, ": ", format(round(dev_val, 2), nsmall = 2))
}

#####################################################################
# A, LOLA enrichment volcanoes
#####################################################################

p_lola <- list()
if (!file.exists(lola.file)) {
  log_warn("LOLA results not found: ", lola.file)
} else {
  res_lola <- readRDS(lola.file)
  lola_idx <- integer(0)
  for (nm in names(comparisons)) {
    hit <- tryCatch(
      {
        h <- resolve_lola_comparison(res_lola,
          test = comparisons[[nm]]$test, ref = comparisons[[nm]]$ref,
          override = lola.comparison[[nm]]
        )
        h$region <- resolve_lola_region(res_lola, h$index, lola.region)
        h
      },
      error = function(e) {
        log_warn(nm, ": ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(hit)) next
    lola_idx[nm] <- hit$index

    pl <- tryCatch(
      lolaVolcanoPlot(
        lolaRes = res_lola, outputDir = plot.dir,
        comparison = hit$index, region = hit$region,
        label = paste0(nm, "_vs_", comparisons[[nm]]$ref),
        grp1 = hit$grp1, grp2 = hit$grp2,
        n.label = top.label.lola,
        motifs = lola.label.motifs,
        userSets = lola.userSets,
        cell_colors = cell_type_colors
      ),
      error = function(e) {
        log_warn(nm, " volcano: ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(pl)) next
    p_lola[[nm]] <- pl
    log_info("LOLA volcano for ", nm, ", region ", hit$region)
  }
  if (anyDuplicated(lola_idx)) {
    stop(
      "Contrasts resolved to the same LOLA comparison: ",
      paste(names(lola_idx), lola_idx, sep = " -> ", collapse = ", "),
      ". Their panels would be identical."
    )
  }
}

#####################################################################
# B, observed minus expected footprints
#####################################################################

if (!file.exists(cache.file)) {
  log_info("Loading RnBeads object from ", rnb.set.path)
  rnb_set <- RnBeads::load.rnb.set(rnb.set.path)
  msites <- rnb.RnBSet.to.GRangesList(mergeSamples(rnb_set, "cellType"))
  saveRDS(msites, cache.file)
  log_info("Wrote ", cache.file)
} else {
  log_info("Reusing merged methylation from ", cache.file)
  msites <- readRDS(cache.file)
}

msites <- as.list(msites)
msites <- mapply(ensure_msites_cols, msites, names(msites), SIMPLIFY = FALSE)
msites <- msites[intersect(group_order, setdiff(names(msites), drop.cell.types))]
if (length(msites) == 0) stop("None of the subtypes are in the merged methylation")
log_info("Footprint groups: ", paste(names(msites), collapse = ", "))

dev_list <- load_cell_type_deviations(dev.file)

tf_bindsites <- getTFbindsites(motifSet = tfSet)
gcfreqs <- getGCfreq(motifSet = motifSet)
gc_dist <- getGenomeGC()

enhancer <- NULL
if (motifSet == "jaspar2020_distal") {
  if (!file.exists(distal.file)) stop("Distal regions file does not exist: ", distal.file)
  enhancer <- readRDS(distal.file)
}

footprint_panel <- function(motif) {
  if (!motif %in% names(tf_bindsites)) {
    log_warn(motif, ": no binding sites, panel skipped")
    return(NULL)
  }
  groups_here <- names(msites)
  mean_devs <- motif_mean_deviations(motif, dev_list, groups_here)
  labels <- vapply(groups_here, function(g) make_label(g, mean_devs[[g]]), character(1))

  per_group <- lapply(groups_here, function(g) {
    df <- tryCatch(
      plotExpectedFootprint(
        motif = motif, tf_bindsites = tf_bindsites, msites = msites[[g]],
        sample_name = g, gc_dist = gc_dist, gcfreqs = gcfreqs,
        enhancer = enhancer, returnPlotData = TRUE
      )$plotDF,
      error = function(e) {
        log_warn(motif, " / ", g, ": ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(df)) {
      return(NULL)
    }
    diff <- df[, .(
      avg_methyl = avg_methyl[type == "Observed"] - avg_methyl[type == "Expected"]
    ), by = x]

    # The flanks are the baseline, so the curves start from zero away from
    # the motif and the three groups are comparable
    flank <- max(abs(diff$x), na.rm = TRUE)
    diff[, avg_methyl := avg_methyl -
      mean(avg_methyl[abs(x) >= flank - flank.norm], na.rm = TRUE)]
    diff[, group := g]
    diff
  })
  per_group <- Filter(Negate(is.null), per_group)
  if (length(per_group) == 0) {
    log_warn(motif, ": no group produced a footprint")
    return(NULL)
  }
  df <- rbindlist(per_group)

  present <- intersect(group_order, unique(df$group))
  colours_here <- setNames(unname(cell_type_colors[present]), labels[present])
  df[, group := factor(labels[match(group, present)], levels = labels[present])]

  ggplot(df, aes(x = x, y = avg_methyl, colour = group)) +
    geom_line(linewidth = 0.4) +
    scale_colour_manual(values = colours_here, name = NULL) +
    coord_cartesian(xlim = c(-200, 200)) +
    labs(
      title = motif,
      x = "Distance from motif center",
      y = "Methylation difference (Observed - Expected)"
    ) +
    theme_classic(base_size = base.size) +
    theme(
      legend.position = "bottom",
      legend.key.size = unit(3.5, "mm"),
      plot.title = element_text(hjust = 0, face = "plain", size = base.size)
    )
}

p_foot <- Filter(Negate(is.null), lapply(footprint.motifs, footprint_panel))

#####################################################################
# The assembled supplementary figure
#####################################################################

if (length(p_lola) == 0 || length(p_foot) == 0) {
  stop("One of the two panels is empty, the supplementary figure is not written")
}

p_lola[[1]] <- p_lola[[1]] + labs(tag = "A")
p_foot[[1]] <- p_foot[[1]] + labs(tag = "B")

supplementary <- wrap_plots(
  wrap_plots(p_lola, ncol = 1),
  wrap_plots(p_foot, nrow = 2),
  ncol = 1, heights = c(1.5, 1.8)
) & tag_theme

file <- file.path(plot.dir, "memT_supplementary.pdf")
ggsave(file, supplementary & no_bg,
  width = fig.width, height = fig.height, bg = "transparent", limitsize = FALSE
)
log_success("Wrote ", file)

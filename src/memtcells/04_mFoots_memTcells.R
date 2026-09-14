#!/usr/bin/env Rscript

#####################################################################
# 04_mFoots_memTcells.R
# created on 08-11-2025 by Irem B Gunduz
# Updated by IBG on 23-08-2026
# Methylation footprints of selected TFs for the CD4 memory T cell
# subtypes. Metric: Observed - Expected, normalised by the flanks
#####################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(GenomicRanges)
  library(logger)
  library(RnBeads)
  library(methylTFR)
  library(methylTFRAnnotationHg38)
  library(SummarizedExperiment)
})
set.seed(42)

src.dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR_manuscript/src"
source(file.path(src.dir, "utils.R"))

#####################################################################
# Settings
#####################################################################

motifSet <- "jaspar2020_distal"
tfSet <- "jaspar2020" # binding sites are shared with the genome wide set

# Distance from the outer edge used to normalise the flanking baseline
flank.norm <- 50

# TEMRA is unreplicated and excluded everywhere
drop.cell.types <- "TEMRA"

# Legend order of the subtypes
group_order <- c("TN", "TCM", "TEM")

# Green scheme from utils.R, shared with every figure that splits by
# subtype. Both spellings are listed because the cellType column of the
# sample annotation and the cleaned RNA sample names disagree on the
# T prefix.
base_colors <- c(
  CELL_TYPE_COLORS,
  "CM" = unname(CELL_TYPE_COLORS[["TCM"]]),
  "EM" = unname(CELL_TYPE_COLORS[["TEM"]])
)

# Distal regulatory regions, kept outside the annotation package
distal.file <- "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFRAnnotationHg38_old/inst/extdata/distal_regions.RDS"

# Directories
analysis.dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/memoryTcells"
rnb.set.path <- file.path(analysis.dir, "reports", "data_import_data", "rnb.set_preprocessed")

# Deviations from 02, averaged per subtype and shown in the legend
dev.tag <- "mTFR_devs_230826"
dev.file <- file.path(analysis.dir, dev.tag, paste0(motifSet, "_deviations.RDS"))

# Cache of the merged per cell type methylation, shared with 05
cache.dir <- file.path(analysis.dir, "debug")
if (!dir.exists(cache.dir)) dir.create(cache.dir, recursive = TRUE)
cache.file <- file.path(cache.dir, "msites_merged_cellType_230826.rds")

plot.dir <- file.path(analysis.dir, "mFoot_230826")
if (!dir.exists(plot.dir)) dir.create(plot.dir, recursive = TRUE)

# Transcription factors of interest
tfs <- unique(c(
  "JUN", "JUN(var.2)", "RELB", "FOS", "FOXP2",
  "BATF", "IRF4", "SP1", "SPIB", "FOSL2"
))

#####################################################################
# Helper functions
#####################################################################

# computeFootprint expects score and coverage metadata columns
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

# One deviation matrix per subtype, taken from the cellType column of the
# deviations object. NULL when 02 has not been run, in which case the
# legend falls back to the plain subtype names.
load_cell_type_deviations <- function(file) {
  if (!file.exists(file)) {
    log_warn("No deviations at ", file, ", the legend will omit them")
    return(NULL)
  }
  log_info("Loading deviations from ", file)
  obj <- readRDS(file)

  cd <- as.data.frame(colData(obj), stringsAsFactors = FALSE)
  if (!"cellType" %in% colnames(cd)) {
    log_warn(
      "The deviations object carries no cellType column. Available: ",
      paste(colnames(cd), collapse = ", ")
    )
    return(NULL)
  }

  cell_types <- as.character(cd$cellType)
  mat <- deviations(obj)
  out <- lapply(split(seq_along(cell_types), cell_types), function(idx) {
    mat[, idx, drop = FALSE]
  })
  log_info(
    "Deviations for: ",
    paste(names(out), lengths(lapply(out, colnames)), sep = " n=", collapse = ", ")
  )
  out
}

# Mean deviation of a motif within a subtype. Matching is exact first, then
# as a fixed substring, because motif names such as JUN(var.2) would
# otherwise be read as a regular expression.
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
    if (length(idx) > 1) {
      log_warn(motif, ": ", length(idx), " matching rows in ", g, ", using the first")
      idx <- idx[1]
    }
    vals <- mat[idx, ]
    if (all(is.na(vals))) NA_real_ else mean(vals, na.rm = TRUE)
  }, numeric(1))
}

# Legend label, "TN: 1.29" when a score is available
make_label <- function(group, dev_val) {
  if (is.na(dev_val)) {
    return(group)
  }
  paste0(group, ": ", format(round(dev_val, 2), nsmall = 2))
}

#####################################################################
# Load and merge the RnBeads methylation data
#####################################################################

if (!file.exists(cache.file)) {
  log_info("Loading RnBeads object from ", rnb.set.path)
  rnb_set <- RnBeads::load.rnb.set(rnb.set.path)

  # Merging the replicates of each cell type
  rnbset_merged <- mergeSamples(rnb_set, "cellType")
  msites <- rnb.RnBSet.to.GRangesList(rnbset_merged)
  saveRDS(msites, cache.file)
  log_info("Wrote ", cache.file)
} else {
  log_info("Reusing merged methylation from ", cache.file)
  msites <- readRDS(cache.file)
}

msites <- as.list(msites)
msites <- mapply(ensure_msites_cols, msites, names(msites), SIMPLIFY = FALSE)

# Drop the unreplicated subtypes, then keep a fixed legend order
msites <- msites[!names(msites) %in% drop.cell.types]
groups <- c(intersect(group_order, names(msites)), setdiff(names(msites), group_order))
if (length(groups) == 0) {
  stop("No cell types left in the merged methylation data")
}
msites <- msites[groups]

missing_colors <- setdiff(groups, names(base_colors))
if (length(missing_colors) > 0) {
  stop("No colour defined for cell type(s): ", paste(missing_colors, collapse = ", "))
}
log_info("Plotting groups: ", paste(groups, collapse = ", "))

#####################################################################
# Load the motif annotation
#####################################################################

tf_bindsites <- getTFbindsites(motifSet = tfSet)
gcfreqs <- getGCfreq(motifSet = motifSet)
gc_dist <- getGenomeGC()

missing_tfs <- setdiff(tfs, names(tf_bindsites))
if (length(missing_tfs) > 0) {
  log_warn("Motifs not found in the binding sites: ", paste(missing_tfs, collapse = ", "))
}
tf_bindsites <- tf_bindsites[names(tf_bindsites) %in% tfs]

# The distal run restricts both the binding sites and the GC background
enhancer <- NULL
if (motifSet == "jaspar2020_distal") {
  if (!file.exists(distal.file)) {
    stop("Distal regions file does not exist: ", distal.file)
  }
  enhancer <- readRDS(distal.file)
  log_info("Loaded ", length(enhancer), " distal regions from ", distal.file)
}

#####################################################################
# Plot the observed minus expected footprints
#####################################################################

plot_and_save_difference <- function(samples, save_dir, dev_list) {
  for (motif in names(tf_bindsites)) {
    out.file <- file.path(save_dir, paste0("TF_footprint_diff_", motif, ".pdf"))
    if (file.exists(out.file)) next

    log_info("Processing motif ", motif, " for ", length(samples), " subtypes")

    mean_devs <- if (is.null(dev_list)) {
      setNames(rep(NA_real_, length(samples)), names(samples))
    } else {
      motif_mean_deviations(motif, dev_list, names(samples))
    }
    labels <- vapply(
      names(samples), function(g) make_label(g, mean_devs[[g]]), character(1)
    )

    per_group <- lapply(names(samples), function(cell_type) {
      plot_data <- tryCatch(
        plotExpectedFootprint(
          motif = motif,
          tf_bindsites = tf_bindsites,
          msites = samples[[cell_type]],
          sample_name = cell_type,
          gc_dist = gc_dist,
          gcfreqs = gcfreqs,
          enhancer = enhancer,
          returnPlotData = TRUE
        )$plotDF,
        error = function(e) {
          log_warn(motif, " / ", cell_type, ": ", conditionMessage(e))
          NULL
        }
      )
      if (is.null(plot_data)) {
        return(NULL)
      }

      # Observed minus expected methylation
      difference_data <- plot_data[, .(
        avg_methyl = avg_methyl[type == "Observed"] - avg_methyl[type == "Expected"]
      ), by = x]

      # Normalise by subtracting the flanking region baseline, which is on
      # the same scale as the difference itself. Dividing by that baseline
      # would rescale the curve by an arbitrary factor, because the mean of
      # a difference over the flanks sits near zero.
      flank <- max(abs(difference_data$x), na.rm = TRUE)
      idx <- abs(difference_data$x) >= flank - flank.norm
      norm_baseline <- mean(difference_data$avg_methyl[idx], na.rm = TRUE)
      difference_data[, avg_methyl := avg_methyl - norm_baseline]

      difference_data[, type := labels[[cell_type]]]
      difference_data[, original_group := cell_type]
      difference_data
    })
    per_group <- Filter(Negate(is.null), per_group)
    if (length(per_group) == 0) {
      log_warn(motif, ": no subtype produced a footprint, skipping")
      next
    }
    combined_data <- rbindlist(per_group)

    # Keep the legend order and the colours tied to the group, not the label
    present <- unique(combined_data$original_group)
    dynamic_colors <- setNames(base_colors[present], labels[present])
    combined_data[, type := factor(type, levels = names(dynamic_colors))]

    p_combined <- ggplot(combined_data, aes(x = x, y = avg_methyl, color = type)) +
      geom_line() +
      xlab("Distance from motif center") +
      ylab("Methylation difference (Observed - Expected)") +
      theme_classic() +
      ggtitle(paste("TF footprint difference for", motif)) +
      scale_color_manual(values = dynamic_colors, name = NULL) +
      theme(legend.position = "bottom") +
      xlim(-200, 200)

    ggsave(filename = out.file, plot = p_combined, width = 12, height = 8)
    log_info("Wrote ", out.file)
  }
}

dev_list <- load_cell_type_deviations(dev.file)

plot_and_save_difference(msites, plot.dir, dev_list)
log_success("Finished the observed minus expected footprints")

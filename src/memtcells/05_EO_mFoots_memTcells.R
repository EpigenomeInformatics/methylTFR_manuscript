#!/usr/bin/env Rscript

#####################################################################
# 05_EO_mFoots_memTcells.R
# created on 14-05-2025 by Irem B Gunduz
# Updated by IBG on 23-08-2026
# Expected vs observed methylation footprints of selected TFs for the
# CD4 memory T cell subtypes. Metric: raw Expected and Observed curves
#
# Counterpart of 04_mFoots_memTcells.R (Observed / Expected). It shares
# the same merged GRangesList cache, so run 04 first or let this script
# build it.
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
})
set.seed(42)

src.dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR_manuscript/src"
source(file.path(src.dir, "utils.R"))

#####################################################################
# Settings
#####################################################################

motifSet <- "jaspar2020_distal"
tfSet <- "jaspar2020" # binding sites are shared with the genome wide set

# TEMRA is unreplicated and excluded everywhere
drop.cell.types <- "TEMRA"

# Legend order of the subtypes
group_order <- c("TN", "TCM", "TEM")

# Observed uses the subtype colour from utils.R, matching 04. Expected
# uses a lighter tint of the same hue, mixed halfway towards white so the
# pairing survives a change to CELL_TYPE_COLORS.
tint <- function(hex, amount = 0.55) {
  rgb_vals <- grDevices::col2rgb(hex)[, 1]
  mixed <- rgb_vals + (255 - rgb_vals) * amount
  grDevices::rgb(mixed[1], mixed[2], mixed[3], maxColorValue = 255)
}

base_colors <- unlist(lapply(names(CELL_TYPE_COLORS), function(ct) {
  setNames(
    c(tint(CELL_TYPE_COLORS[[ct]]), unname(CELL_TYPE_COLORS[[ct]])),
    paste0(c("Expected_", "Observed_"), ct)
  )
}))

# Distal regulatory regions, kept outside the annotation package
distal.file <- "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFRAnnotationHg38_old/inst/extdata/distal_regions.RDS"

# Directories
analysis.dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/memoryTcells"
rnb.set.path <- file.path(analysis.dir, "reports", "data_import_data", "rnb.set_preprocessed")

# Cache of the merged per cell type methylation, shared with 04
cache.dir <- file.path(analysis.dir, "debug")
if (!dir.exists(cache.dir)) dir.create(cache.dir, recursive = TRUE)
cache.file <- file.path(cache.dir, "msites_merged_cellType_230826.rds")

plot.dir <- file.path(analysis.dir, "EO_mFoot_230826")
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

missing_colors <- setdiff(
  c(paste0("Expected_", groups), paste0("Observed_", groups)),
  names(base_colors)
)
if (length(missing_colors) > 0) {
  stop("No colour defined for: ", paste(missing_colors, collapse = ", "))
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
# Plot the raw expected and observed footprints
#####################################################################

plot_and_save_eo <- function(samples, save_dir) {
  for (motif in names(tf_bindsites)) {
    out.file <- file.path(save_dir, paste0("TF_footprint_EO_", motif, ".pdf"))
    if (file.exists(out.file)) next

    log_info("Processing motif ", motif, " for ", length(samples),
      " subtypes (expected vs observed)")

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

      # type is either "Observed" or "Expected", both curves are kept
      plot_data <- copy(plot_data)
      plot_data[, type := paste(type, cell_type, sep = "_")]
      plot_data
    })
    per_group <- Filter(Negate(is.null), per_group)
    if (length(per_group) == 0) {
      log_warn(motif, ": no subtype produced a footprint, skipping")
      next
    }
    combined_data <- rbindlist(per_group)

    # Keep the legend order defined by base_colors
    present <- intersect(names(base_colors), unique(combined_data$type))
    combined_data[, type := factor(type, levels = present)]

    p_combined <- ggplot(combined_data, aes(x = x, y = avg_methyl, color = type)) +
      geom_line() +
      xlab("Distance from motif center") +
      ylab("Average methylation") +
      theme_classic() +
      ggtitle(paste("Expected vs observed footprints for", motif)) +
      scale_color_manual(values = base_colors[present]) +
      theme(legend.position = "bottom") +
      xlim(-200, 200)

    ggsave(filename = out.file, plot = p_combined, width = 12, height = 8)
    log_info("Wrote ", out.file)
  }
}

plot_and_save_eo(msites, plot.dir)
log_success("Finished the expected vs observed footprints")

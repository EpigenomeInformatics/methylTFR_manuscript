#!/usr/bin/env Rscript

#####################################################################
# 04_mFoots_allcelltypes.R
# created on 07-09-2026 by Irem B Gunduz
# Updated by IBG on 23-08-2026
# Combined methylation footprint plots for all cell type groups:
#   Monocytes / Granulocytes / Bcell_mem / Bcell_naive / Tcell
# Metric: Observed - Expected (subtraction, not ratio)
# Mean deviation scores are added to the legend, per cell type
#
# Group mapping (from table(cts)):
#   Bcell_mem    = Bcell_mem, Bcell_gc
#   Bcell_naive  = Bcell_naive, Bcell_pre
#   Tcell        = TCD4, TCD4_cenM, TCD4_effM, TCD8, TCD8_cenM,
#                  TCD8_effM, TCD8_term, Treg
#   Monocytes    = Mono
#   Granulocytes = Neut_band, Neut_mat, Neut_meta, Neut_myel,
#                  Neut_segm, Eosi_mat
# All other cellTypeShort levels (DC*, Endo*, Eryth, Mac_f0/f1/f2, MK,
# MPP, MSC_BM, NK, Ost, Plas, Thym*) are excluded from this analysis.
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

#####################################################################
# Settings
#####################################################################

# Motif set used for the deviations and the expected footprints
motifSet <- "jaspar2020_distal"
tfSet <- "jaspar2020" # binding sites are shared with the genome wide set

# Distance from the outer edge used to normalise the flanking baseline
flank.norm <- 30

# Directory where the annotation resources are stored locally
annotation.dir <- "/icbb/projects/share/annotations/methylTFRAnnotationHg38/inst/extdata"
if (!dir.exists(annotation.dir)) {
  stop("Annotation directory does not exist: ", annotation.dir)
}
options(methylTFRAnnotationHg38.datadir = annotation.dir)

# Distal regulatory regions, kept outside the annotation package
distal.file <- "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFRAnnotationHg38_old/inst/extdata/distal_regions.RDS"

# Directories
analysis.dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/blueprint/"
rnb.tag <- "RnBeads_230826"
dev.tag <- "mTFR_devs_230826"

rnb.set.path <- file.path(analysis.dir, rnb.tag, "reports", "data_import_data", "rnb.set_preprocessed")
sannot.file <- file.path(analysis.dir, rnb.tag, "reports", "data_import_data", "annotation.csv")
dev.file <- file.path(analysis.dir, dev.tag, paste0(motifSet, "_deviations.RDS"))

# Cache of the merged per group methylation, shared with 05
cache.dir <- file.path(analysis.dir, "debug")
if (!dir.exists(cache.dir)) dir.create(cache.dir, recursive = TRUE)
cache.file <- file.path(cache.dir, paste0("msites_merged_celltype5G_", rnb.tag, ".rds"))

plot.dir <- file.path(analysis.dir, "mFoot_allcelltypes_230826")
if (!dir.exists(plot.dir)) dir.create(plot.dir, recursive = TRUE)

# Transcription factors of interest
tfs <- unique(c(
  "TFAP2B", "PAX5", "PAX9", "PAX1", "NHLH2", "ASCL1", "NHLH1", "BHLHE22", "FERD3L", "PAX6",
  "EMX1", "PAX4", "EN1", "LHX1", "ELF4", "ELF2", "ETV5", "ETV6", "ELF5", "SPIC",
  "SPIB", "SPI1", "EHF", "ELF3", "IKZF1", "TFAP2C", "TFAP2A", "TFAP2E",
  "TGIF1", "CREB3L4", "PBX3", "TAL1::TCF3", "MYOG", "ATOH1", "MYF5", "BHLHA15", "ZBTB18", "EBF3",
  "EBF1", "TFAP4", "NEUROD1", "VSX2", "EGR4", "DPRX", "SOX8", "CUX1", "CUX2", "TEAD3", "CEBPB",
  "FOSL1::JUND"
))

# Legend order of the five groups
group_order <- c("Monocytes", "Granulocytes", "Bcell_mem", "Bcell_naive", "Tcell")

base_colors <- c(
  "Monocytes" = "#C2703D",
  "Granulocytes" = "#E8A33D",
  "Bcell_mem" = "#7B1E3D",
  "Bcell_naive" = "#C46B8C",
  "Tcell" = "#4FC3D9"
)

#####################################################################
# Helper functions
#####################################################################

# Collapse cellTypeShort into the five groups used here
map_cellType5Group <- function(cts) {
  cts <- as.character(cts)
  dplyr::case_when(
    cts %in% c("Bcell_mem", "Bcell_gc") ~ "Bcell_mem",
    cts %in% c("Bcell_naive", "Bcell_pre") ~ "Bcell_naive",
    cts %in% c(
      "TCD4", "TCD4_cenM", "TCD4_effM", "TCD8",
      "TCD8_cenM", "TCD8_effM", "TCD8_term", "Treg"
    ) ~ "Tcell",
    cts == "Mono" ~ "Monocytes",
    cts %in% c(
      "Neut_band", "Neut_mat", "Neut_meta",
      "Neut_myel", "Neut_segm", "Eosi_mat"
    ) ~ "Granulocytes",
    TRUE ~ NA_character_
  )
}

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

# Mean deviation of a motif within each group, NA when the motif is absent
motif_mean_deviations <- function(motif, dev_mat, dev_grps, groups) {
  idx <- which(rownames(dev_mat) == motif)
  if (length(idx) == 0) {
    idx <- grep(paste0("\\b", motif, "\\b"), rownames(dev_mat), ignore.case = TRUE)
    if (length(idx) > 1) {
      log_warn(motif, ": ", length(idx), " matching rows in the deviation matrix, using the first")
      idx <- idx[1]
    }
  }
  vapply(groups, function(g) {
    if (length(idx) == 0) {
      return(NA_real_)
    }
    vals <- dev_mat[idx, dev_grps == g]
    if (length(vals) == 0 || all(is.na(vals))) NA_real_ else mean(vals, na.rm = TRUE)
  }, numeric(1))
}

# Legend label carrying the mean deviation score
make_label <- function(group, dev_val) {
  dev_str <- if (is.na(dev_val)) "N/A" else format(round(dev_val, 2), nsmall = 2)
  paste0(group, " (Mean Dev: ", dev_str, ")")
}

#####################################################################
# Load and process the deviation scores
#####################################################################

log_info("Loading deviations from ", dev.file)
dev_obj <- readRDS(dev.file)
dev_matrix <- deviations(dev_obj)

sannot <- read.csv(sannot.file, stringsAsFactors = FALSE)
sannot$bedFile <- as.character(sannot$bedFile)

# Match the annotation rows to the deviation matrix columns
matched_sannot <- sannot[match(colnames(dev_matrix), sannot$bedFile), ]
matched_sannot$cellType5Group <- map_cellType5Group(matched_sannot$cellTypeShort)

# Keep healthy samples that fall into one of the five groups,
# which() drops the samples that are missing from the annotation
keep_dev <- which(
  !is.na(matched_sannot$DISEASE) &
    matched_sannot$DISEASE == "None" &
    !is.na(matched_sannot$cellType5Group)
)
dev_matrix <- dev_matrix[, keep_dev, drop = FALSE]
dev_groups <- matched_sannot$cellType5Group[keep_dev]
log_info("Kept ", length(keep_dev), " samples for the deviation scores")

#####################################################################
# Load and merge the RnBeads methylation data
#####################################################################

if (!file.exists(cache.file)) {
  log_info("Loading RnBeads object from ", rnb.set.path)
  rnb_set <- RnBeads::load.rnb.set(rnb.set.path)

  cellType5Group <- map_cellType5Group(pheno(rnb_set)$cellTypeShort)
  rnb_set@pheno$cellType5Group <- cellType5Group

  # Drop disease samples and everything outside the five target groups
  drop_idx <- unique(c(
    which(as.character(pheno(rnb_set)$DISEASE) != "None"),
    which(is.na(cellType5Group))
  ))
  if (length(drop_idx) > 0) {
    log_info(length(drop_idx), " samples excluded (disease or outside the 5 target groups)")
    rnb_set <- RnBeads::remove.samples(rnb_set, drop_idx)
  }

  rnbset_merged <- mergeSamples(rnb_set, "cellType5Group")
  msites <- rnb.RnBSet.to.GRangesList(rnbset_merged)
  saveRDS(msites, cache.file)
  log_info("Wrote ", cache.file)
} else {
  log_info("Reusing merged methylation from ", cache.file)
  msites <- readRDS(cache.file)
}

msites <- as.list(msites)
msites <- mapply(ensure_msites_cols, msites, names(msites), SIMPLIFY = FALSE)

# Only keep the five target groups, in a fixed order
groups <- intersect(group_order, names(msites))
if (length(groups) == 0) {
  stop("None of the five target groups are present in the merged methylation data")
}
msites <- msites[groups]
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

# The distal run restricts both the binding sites and the GC background.
# Ensembl Regulatory Build v104, hg38, filtered to distal. Same region set
# that jaspar2020_distal_motif_gcfreq.rds was built against, and the one 01
# hands to RnBeads, so all three views of the data stay comparable.
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

plot_and_save_difference <- function(samples, save_dir, dev_mat, dev_grps) {
  for (motif in names(tf_bindsites)) {
    out.file <- file.path(save_dir, paste0("TF_footprint_diff_", motif, ".pdf"))
    if (file.exists(out.file)) next

    log_info("Processing motif ", motif, " for ", length(samples), " groups")
    mean_devs <- motif_mean_deviations(motif, dev_mat, dev_grps, names(samples))
    labels <- vapply(names(samples), function(g) make_label(g, mean_devs[[g]]), character(1))

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

      # Normalise by subtracting the flanking region baseline
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
      log_warn(motif, ": no group produced a footprint, skipping")
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
      scale_color_manual(values = dynamic_colors) +
      theme(legend.position = "bottom") +
      xlim(-200, 200)

    ggsave(filename = out.file, plot = p_combined, width = 12, height = 8)
    log_info("Wrote ", out.file)
  }
}

plot_and_save_difference(msites, plot.dir, dev_matrix, dev_groups)
log_success("Finished the observed minus expected footprints")

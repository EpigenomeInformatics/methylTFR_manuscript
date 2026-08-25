#!/usr/bin/env Rscript

#####################################################################
# 03_pca_memTcells.R
# created on 14-05-2025 by Irem B Gunduz
# Updated by IBG on 23-08-2026
# PCA of the methylTFR deviations of the CD4 memory T cell subtypes.
# Replaces 03_plot_memTcells.R, whose differential test is done
# properly in 06_differential_TFs.R
#####################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggfortify)
  library(logger)
  library(SummarizedExperiment)
  library(methylTFR)
})
set.seed(42)

#####################################################################
# Settings
#####################################################################

motifSets <- c("jaspar2020", "jaspar2020_distal")

# PCA on the raw deviation matrix, as in the blueprint scripts
pca.center <- FALSE
pca.scale <- FALSE

# TEMRA is unreplicated and excluded everywhere
drop.cell.types <- "TEMRA"

# Both spellings are listed because the cellType column of the sample
# annotation and the cleaned RNA sample names disagree on the T prefix
cell_type_colors <- c(
  "TN" = "#C8E0B4",
  "TCM" = "#4492C6",
  "CM" = "#4492C6",
  "TEM" = "#43B6C4",
  "EM" = "#43B6C4",
  "TEMRA" = "#898FB5"
)

# Directories
analysis.dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/memoryTcells"
dev.tag <- "mTFR_devs_230826"

fig.dir <- file.path(analysis.dir, "figures_230826")
if (!dir.exists(fig.dir)) dir.create(fig.dir, recursive = TRUE)

#####################################################################
# Helper functions
#####################################################################

# Cell type of each column, taken from the annotation carried by the
# object rather than from a grepl on the file names
sample_cell_types <- function(dev_obj) {
  cd <- as.data.frame(colData(dev_obj), stringsAsFactors = FALSE)
  if (!"cellType" %in% colnames(cd)) {
    stop(
      "The deviations object carries no cellType column. Available: ",
      paste(colnames(cd), collapse = ", ")
    )
  }
  as.character(cd$cellType)
}

# Donor identifier, Hf03 or Hf04, parsed from the sample name
sample_donors <- function(ids) {
  donors <- rep(NA_character_, length(ids))
  hit <- grepl("Hf[0-9]+", ids)
  donors[hit] <- sub(".*?(Hf[0-9]+).*", "\\1", ids[hit])
  donors
}

#####################################################################
# PCA per motif set
#####################################################################

for (motifSet in motifSets) {
  dev.file <- file.path(analysis.dir, dev.tag, paste0(motifSet, "_deviations.RDS"))
  if (!file.exists(dev.file)) {
    log_warn(motifSet, ": ", dev.file, " not found, skipping")
    next
  }

  log_info("Loading deviations from ", dev.file)
  dev_obj <- readRDS(dev.file)

  cell_types <- sample_cell_types(dev_obj)
  keep <- which(!is.na(cell_types) & !cell_types %in% drop.cell.types)
  if (length(keep) < 3) {
    log_warn(motifSet, ": only ", length(keep), " samples left, skipping")
    next
  }
  dev_obj <- dev_obj[, keep]
  cell_types <- cell_types[keep]

  missing_colors <- setdiff(unique(cell_types), names(cell_type_colors))
  if (length(missing_colors) > 0) {
    stop("No colour defined for cell type(s): ", paste(missing_colors, collapse = ", "))
  }

  mat <- deviations(dev_obj)
  mat <- mat[is.finite(rowSums(mat)), , drop = FALSE]
  log_info(motifSet, ": ", nrow(mat), " motifs x ", ncol(mat), " samples")

  pca <- prcomp(t(mat), center = pca.center, scale. = pca.scale)

  # Only the grouping columns are passed, ggfortify binds the scores itself
  plot_data <- data.frame(
    cell_type = cell_types,
    donor = sample_donors(colnames(mat)),
    stringsAsFactors = FALSE
  )

  file <- file.path(fig.dir, paste0("PCA_deviations_", motifSet, "_memTcells.pdf"))
  pdf(file, width = 8, height = 8)
  print(
    autoplot(pca,
      data = plot_data,
      colour = "cell_type",
      shape = "donor",
      size = 5,
      main = paste("PCA of methylTFR deviations,", motifSet)
    ) +
      theme_classic() +
      scale_color_manual(values = cell_type_colors) +
      theme(legend.position = "bottom")
  )
  dev.off()
  log_info("Wrote ", file)
}

log_success("Finished the memory T cell PCA")

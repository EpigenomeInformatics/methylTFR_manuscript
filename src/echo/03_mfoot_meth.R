#!/usr/bin/env Rscript

#####################################################################
# created on 14-05-2025 by Irem Gunduz
# Plotting the motif footprints for ECHO data
#####################################################################

suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(ggplot2)
    library(GenomicRanges)
    library(tools)
    library(methylTFR)
    library(RnBeads)
    library(methylTFRAnnotationHg38)
})
set.seed(13)

plot_dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/echo/mFoot_230826/"
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# Pseudobulks from ECHO methylation
# Define the paths to the methylome data for Th and Tc samples
pb.dir <- "/icbb/projects/igunduz/DARPA_analysis/ECHO_Cell_PB_261025/"
if (!dir.exists(pb.dir)) stop("Pseudobulk directory not found: ", pb.dir)
paths <- list.files(pb.dir, full.names = TRUE)
if (length(paths) == 0) stop("No pseudobulk files in ", pb.dir)
message("Reading ", length(paths), " pseudobulks from ", pb.dir)

# Create mSite object
msites <- GRangesList()
for(path in paths){
    msite <- read_methylome(path, "bismarkcytosine")
    name <- tools::file_path_sans_ext(basename(path))
    msites[[name]] <- msite
}

# Prepare motif data. The distal GC frequency table with the genome wide
# binding sites, as everywhere else in the manuscript.
motifSet <- "jaspar2020"
gcfreqs <- getGCfreq("jaspar2020_distal")
tf_bindsites <- getTFbindsites(motifSet)
gc_dist <- getGenomeGC()

# TFs to plot, from ECHO's own MOFA factor table. The original path here
# was topTFs_bp_pseudob.CSV, which is the BLUEPRINT list.
#
# The table is long: one row per feature, factor and view, with the view
# appended to the feature name (FOXF2_chromVar). Copy it into tables/ if
# you want the run to stop depending on another user's directory.
tfs.file <- "/icbb/projects/nitschre/methylTFR/tables/echo_factor_table.csv"

# Views and factors to draw the TFs from, NULL for all of them
tfs.views <- NULL
tfs.factors <- NULL

# Top features per factor and view, ranked by absolute weight
tfs.per.factor <- 6

if (!file.exists(tfs.file)) stop("ECHO factor table not found: ", tfs.file)
tfs_tab <- read.csv(tfs.file, header = TRUE, stringsAsFactors = FALSE)

required <- c("feature", "factor", "value", "view")
missing_cols <- setdiff(required, colnames(tfs_tab))
if (length(missing_cols) > 0) {
  stop(
    "Missing column(s) in ", basename(tfs.file), ": ",
    paste(missing_cols, collapse = ", "), ". Present: ",
    paste(colnames(tfs_tab), collapse = ", ")
  )
}
message(
  nrow(tfs_tab), " rows, views: ", paste(unique(tfs_tab$view), collapse = ", "),
  ", factors: ", length(unique(tfs_tab$factor))
)

if (!is.null(tfs.views)) tfs_tab <- tfs_tab[tfs_tab$view %in% tfs.views, ]
if (!is.null(tfs.factors)) tfs_tab <- tfs_tab[tfs_tab$factor %in% tfs.factors, ]
if (nrow(tfs_tab) == 0) stop("No rows left after filtering by view or factor")

# The view is appended to the feature name. Stripping it with the row's own
# view value keeps motif names containing punctuation intact, which
# gsub("_.*$", "") would not: MZF1(var.2)_chromVar, MAX::MYC_chromVar.
tfs_tab$feature <- mapply(
  function(feat, vw) sub(paste0("_", vw, "$"), "", feat),
  tfs_tab$feature, tfs_tab$view,
  USE.NAMES = FALSE
)

# Top features per factor and view, by absolute weight
tfs_tab$value_abs <- abs(tfs_tab$value)
tfs_tab <- tfs_tab[order(tfs_tab$factor, tfs_tab$view, -tfs_tab$value_abs), ]
rank_in_group <- stats::ave(
  seq_len(nrow(tfs_tab)),
  paste(tfs_tab$factor, tfs_tab$view),
  FUN = seq_along
)
tfs <- unique(tfs_tab$feature[rank_in_group <= tfs.per.factor])
message("Top ", tfs.per.factor, " per factor and view gives ", length(tfs), " unique TFs")

# Indexing a GRangesList by a missing name errors, so the absent ones are
# reported and dropped instead
missing_tfs <- setdiff(tfs, names(tf_bindsites))
if (length(missing_tfs) > 0) {
  message(
    length(missing_tfs), " TFs are not in the binding sites and are skipped: ",
    paste(head(missing_tfs, 10), collapse = ", ")
  )
}
tf_bindsites <- tf_bindsites[names(tf_bindsites) %in% tfs]
if (length(tf_bindsites) == 0) stop("None of the requested TFs have binding sites")
message("Plotting ", length(tf_bindsites), " TFs")

# Ensembl Regulatory Build v104 distal regions. The share path no longer
# holds this file, it lives with the old annotation package.
distal.file <- "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFRAnnotationHg38_old/inst/extdata/distal_regions.RDS"
if (!file.exists(distal.file)) stop("Distal regions not found: ", distal.file)
enhancer <- readRDS(distal.file)

# Function to generate and save a plot for all samples
plot_and_save_difference_covid <- function(samples, save_dir, obs_colors) {

  for (motif in names(tf_bindsites)) {
    if (!file.exists(file.path(save_dir, paste0("TF_footprint_diff_", motif, ".pdf")))) {
      logger.start(paste("Processing motif", motif, "for Blueprint cell samples"))

      # Generate plot data and calculate the difference
      combined_data <- rbindlist(lapply(names(samples), function(cell_type) {
        plot_data <- plotExpectedFootprint(
          motif = motif,
          tf_bindsites = tf_bindsites,
          msites = samples[[cell_type]],
          sample_name = cell_type,
          gc_dist = gc_dist,
          gcfreqs = gcfreqs,
          enhancer = enhancer,
          returnPlotData = TRUE
        )
        # Calculate observed minus expected methylation
        difference_data <- plot_data$plotDF[, .(avg_methyl = avg_methyl[type == "Observed"] - avg_methyl[type == "Expected"]), by = x]
        difference_data[, type := paste("Observed minus Expected", cell_type)]

        # Normalise by subtracting the flanking baseline, which is on the same
        # scale as the difference. Dividing by it would rescale the curve by an
        # arbitrary factor, because the mean of a difference over the flanks
        # sits near zero rather than near one.
        flankNorm <- 50
        flank <- max(abs(difference_data$x), na.rm = TRUE)
        idx <- abs(difference_data$x) >= flank - flankNorm
        norm_baseline <- mean(difference_data$avg_methyl[idx], na.rm = TRUE)
        difference_data[, avg_methyl := avg_methyl - norm_baseline]

        return(difference_data)
      }))
      logger.completed()

      combined_data[, type := factor(type, levels = names(obs_colors))]

      # Plot the combined difference data for all samples
      logger.info("Plotting combined difference data for all samples")
      p_combined <- ggplot(combined_data, aes(x = x, y = avg_methyl, color = type)) +
        geom_line() + # Add lines for each type
        xlab("Distance from motif center") +
        ylab("Methylation difference (Observed - Expected)") +
        theme_classic() +
        ggtitle(paste("TF footprint difference for", motif)) +
        scale_color_manual(values = obs_colors) + # Use the observed color vector
        theme(legend.position = "bottom") +
        xlim(-200, 200) # Set the x-axis limits

      # Save the plot as a PDF in the specified directory
      ggsave(
        filename = file.path(save_dir, paste0("TF_footprint_diff_", motif, ".pdf")),
        plot = p_combined, width = 12, height = 8
      )
    }
  }
}

# Define colors
cell_type_colors <- c(
  "Observed minus Expected B-cell" = "#AE017E",
  "Observed minus Expected Monocyte" = "#CC4C02",
  "Observed minus Expected NK-cell" = "#A65628",
  "Observed minus Expected Th-Mem" = "#41B6C4",
  "Observed minus Expected Tc-Mem" = "#4292C6",
  "Observed minus Expected Tc-Naive" = "#888FB5",
  "Observed minus Expected Th-Naive" = "#C7E9B4"
)

plot_and_save_difference_covid(msites, plot_dir, cell_type_colors)
#####################################################################################
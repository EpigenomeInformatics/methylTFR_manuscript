#!/usr/bin/env Rscript

#####################################################################
# created on 14-05-2025 by Irem Gunduz
# Plotting the motif footprints for ECHO data
#####################################################################

suppressPackageStartupMessages({
    library(methylTFR)
    library(RnBeads)
    library(methylTFRAnnotationHg38)
    library(dplyr)
})
set.seed(13)

plot_dir <- "/scratch/icbb/regina/methylTFR_manuscript/figures/ECHO/mFoot/"

# Pseudobulks from ECHO methylation
# Define the paths to the methylome data for Th and Tc samples
paths <- list.files("/icbb/projects/igunduz/DARPA_analysis/ECHO_Cell_PB_261025/", full.names=TRUE)

# Create mSite object
msites <- GRangesList()
for(path in paths){
    msite <- read_methylome(path, "bismarkcytosine")
    name <- file_path_sans_ext(basename(path))
    msites[[name]] <- msite
}

# Prepare motif data
motifSet <- "JASPAR2020"
gcfreqs <- getGCfreq("JASPAR2020_distal")
tf_bindsites <- getTFbindsites(motifSet)
gc_dist <- getGenomeGC()

# TFs to plot
tfs <- read.csv("/icbb/projects/nitschre/methylTFR/figures/figure4/topTFs_bp_pseudob.CSV", header=TRUE)
tfs <- tfs$feature
tf_bindsites<- tf_bindsites[tfs]

enhancer <- readRDS("/icbb/projects/share/annotations/methylTFRAnnotationHg38/inst/extdata/distal_regions.RDS")

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
          enhancer = distal,
          returnPlotData = TRUE
        )
        # Calculate observed/expected methylation
        difference_data <- plot_data$plotDF[, .(avg_methyl = avg_methyl[type == "Observed"] / avg_methyl[type == "Expected"]), by = x]
        difference_data[, type := paste("Observed divided Expected", cell_type)]

        # Now normalize the ratio by flanking region
        flankNorm= 50
        flank <- max(abs(difference_data$x), na.rm = TRUE)
        idx <- abs(difference_data$x) >= flank - flankNorm
        norm_factor <- mean(difference_data$avg_methyl[idx], na.rm = TRUE)
        difference_data[, avg_methyl := avg_methyl / norm_factor]

        return(difference_data)
      }))
      logger.completed()

      combined_data[, type := factor(type, levels = names(obs_colors))]

      # Plot the combined difference data for all samples
      logger.info("Plotting combined difference data for all samples")
      p_combined <- ggplot(combined_data, aes(x = x, y = avg_methyl, color = type)) +
        geom_line() + # Add lines for each type
        xlab("Distance from motif center") +
        ylab("Methylation difference (Observed / Expected)") +
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
  "Observed divided Expected B-cell" = "#AE017E",
  "Observed divided Expected Monocyte" = "#CC4C02",
  "Observed divided Expected NK-cell" = "#A65628",
  "Observed divided Expected Th-Mem" = "#41B6C4",
  "Observed divided Expected Tc-Mem" = "#4292C6",
  "Observed divided Expected Tc-Naive" = "#888FB5",
  "Observed divided Expected Th-Naive" = "#C7E9B4"
)

plot_and_save_difference_covid(msites, plot_dir, cell_type_colors)
#####################################################################################
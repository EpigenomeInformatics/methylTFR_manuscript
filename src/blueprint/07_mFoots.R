#!/usr/bin/env Rscript

#####################################################################
# 04_plot_mfoot.R
# created on 14-05-2025 by Irem Gunduz
# Plotting the motif footprints for blueprint data
#####################################################################

suppressPackageStartupMessages({
    library(methylTFR)
    library(RnBeads)
    library(methylTFRAnnotationHg38)
    library(dplyr)
})
set.seed(13)

# Set the paths
plot_dir <- "/scratch/icbb/regina/methylTFR_manuscript/figures/blueprint/TF_footprint_diffmotifs/"
if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE)
}
data_dir <- "/scratch/icbb/regina/data/blueprint/"
if (!dir.exists(data_dir)) {
    dir.create(data_dir, recursive = TRUE)
}

if(!file.exists(paste0(data_dir,"bp_msites.RDS"))){
  # Import Rnbeads object
  rnbset <- load.rnb.set("/scratch/icbb/igunduz/methylTFR_manuscript/blueprint/RnBeads_291025/reports/data_import_data/rnb.set_preprocessed")
  # Merging cell type replicates
  rnbset_merged <- mergeSamples(rnbset, "cellTypeGroup")

  # Get the methylation sites as GRangesList
  msites <- rnb.RnBSet.to.GRangesList(rnbset_merged)
  saveRDS(msites, file =paste0(data_dir,"bp_msites.RDS"))
}else{
  msites <- readRDS(paste0(data_dir,"bp_msites.RDS"))
}

# Plot only NK, Tcells, mono and Bcell
msites <- msites[c("Tcell", "NK","Bcell", "mono")]

# Prepare motif data
motifSet <- "JASPAR2020"
gcfreqs <- getGCfreq("JASPAR2020_distal")
tf_bindsites <- getTFbindsites(motifSet)
gc_dist <- getGenomeGC()

# Define the TFs of interest
tfs <-  tfs <- c(
  "CEBPA", "CEBPB", "IRF8", "KLF4", "RUNX1",    # Monocyte Specific
  "EBF1", "PAX5", "TCF3", "POU2F2",             # B Cell Specific
  "TBX21", "GATA3", "RORC", "FOXP3", "RUNX3",   # T Cell Specific
  "LEF1", "TCF7",                               # T Cell Specific (cont.)
  "SPI1", "BCL6"                                # Shared / Lineage
)

tf_bindsites <- tf_bindsites[names(tf_bindsites) %in% tfs]

distal <- readRDS("/icbb/projects/share/annotations/methylTFRAnnotationHg38/inst/extdata/distal_regions.RDS")
distal <- if(motifSet == "JASPAR2020_distal"){distal}else{NULL}

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

# Define the observed colors for CD4T samples
obs_colors <- c(
  "Observed divided Expected Bcell" = "#980043",
  "Observed divided Expected mono" = "#CD7054",
  "Observed divided Expected NK" = "#bf812d",
  "Observed divided Expected Tcell" = "#40E0D0"
)

# Generate and save plots for all samples with observed divided expected
plot_and_save_difference_covid(msites, plot_dir, obs_colors)

#####################################################################
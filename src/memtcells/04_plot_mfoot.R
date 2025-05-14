#!/usr/bin/env Rscript

#####################################################################
# 04_plot_mfoot.R
# created on 14-05-2025 by Irem Gunduz
# Plotting the motif footprints for memory T cells
#####################################################################


suppressPackageStartupMessages({
    library(methylTFR)
    library(RnBeads)
    library(methylTFRAnnotationHg38)
    library(dplyr)
})
set.seed(42)

# Set the paths
plot_dir <- "/icbb/projects/igunduz/irem_github/methylTFR_manuscript/footprints"
if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE)
}
data_dir <- "/icbb/projects/igunduz/irem_github/methylTFR_manuscript/data/"
if (!dir.exists(data_dir)) {
    dir.create(data_dir, recursive = TRUE)
}

if(!file.exists(paste0(data_dir,"cd4t_msites.RDS")){
    # Import Rnbeads object
rnbset <- load.rnb.set("/icbb/projects/nitschre/methylTFR/results/memoryTcells/reports/data_import_data/rnb.set_preprocessed")

# Merging cell type replicates
rnbset_merged <- mergeSamples(rnbset, "cellType")

# Get the methylation sites as GRangesList
msites <- rnb.RnBSet.to.GRangesList(rnbset_merged)
saveRDS(msites, file =paste0(data_dir,"cd4t_msites.RDS"))
})else{
msites <- readRDS(paste0(data_dir,"cd4t_msites.RDS"))
}

# Load annotation files
gcfreqs <- getGCfreq(motifSet = "jaspar2020_distal")
tf_bindsites <- getTFbindsites(motifSet = "jaspar2020")
gc_dist <- getGenomeGC()
enhancer <- readRDS("/icbb/projects/share/annotations/methylTFRAnnotationHg38/inst/extdata/distal_regions.RDS")

samples <- msites
cell_type <- names(msites)[1]

# Function to generate and save a plot for all samples
plot_and_save_difference_covid <- function(samples, save_dir, obs_colors) {
  for (motif in names(tf_bindsites)) {
    if (!file.exists(file.path(save_dir, paste0("TF_footprint_diff_", motif, ".pdf")))) {
      logger.start(paste("Processing motif", motif, "for CD4T samples"))

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
        # Calculate observed/expected methylation
        difference_data <- plot_data$plotDF[, .(avg_methyl = avg_methyl[type == "Observed"] / avg_methyl[type == "Expected"]), by = x]
        difference_data[, type := paste("Observed divided Expected", cell_type)]

        return(difference_data)
      }))
      logger.completed()

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
  "Observed divided Expected TN" = "#C8E0B4",
  "Observed divided Expected TCM" = "#4492C6",
  "Observed divided Expected TEM" = "#43B6C4",
  "Observed divided Expected TEMRA" = "#898FB5"
)


# Generate and save plots for all samples with observed divided expected
plot_and_save_difference_covid(msites, plot_dir, obs_colors)

#####################################################################
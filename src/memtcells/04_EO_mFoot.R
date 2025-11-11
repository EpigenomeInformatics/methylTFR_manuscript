#!/usr/bin/env Rscript

#####################################################################
# 04_plot_mfoot.R
# created on 14-05-2025 by Irem Gunduz
# Plotting the motif footprints for b and t cells
#####################################################################


suppressPackageStartupMessages({
    library(methylTFR)
    library(RnBeads)
    library(methylTFRAnnotationHg38)
    library(dplyr)
})
set.seed(42)

source("/icbb/projects/nitschre/methylTFR/scripts/memoryTcells/plots.R")
source("/icbb/projects/nitschre/methylTFR/R/plot_helpers.R")
source("/icbb/projects/nitschre/methylTFR/R/expected_deviations.R")

# Set the paths
main.dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/memoryTcells/"
sample_dir <- paste0(main.dir, "reports/")
debug <- paste0(main.dir, "debug/")
if(!dir.exists(debug)){dir.create(debug)}
plot_dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR_manuscript/figures/memtcells/"
plot_dir <- paste0(plot_dir, "EO_mFoot/")
if(!dir.exists(plot_dir)){dir.create(plot_dir)}

if(!file.exists(paste0(debug,"memoryTcells_msites.RDS"))){
    # Import Rnbeads object
rnb_set <- RnBeads::load.rnb.set(paste0(sample_dir, "data_import_data/rnb.set_preprocessed"))

# Merging cell type replicates
rnbset_merged <- mergeSamples(rnb_set, "cellType")

# Get the methylation sites as GRangesList
msites <- rnb.RnBSet.to.GRangesList(rnbset_merged)
saveRDS(msites, file =paste0(debug,"memoryTcells_msites.RDS"))
}else{
msites <- readRDS(paste0(debug,"memoryTcells_msites.RDS"))
}

# Load annotation files
motifSet <- "jaspar2020"
gcfreqs <- readRDS("/scratch/icbb/igunduz/methylTFR_manuscript/blueprint/debug/methylTFRAnnotationHg38/inst/extdata/JASPAR2020_distal_motif_gcfreq.rds")
tf_bindsites <- readRDS("/scratch/icbb/igunduz/methylTFR_manuscript/blueprint/debug/methylTFRAnnotationHg38/inst/extdata/JASPAR2020_tf_bindsites.rds")
gc_dist <- readRDS("/scratch/icbb/igunduz/methylTFR_manuscript/blueprint/debug/methylTFRAnnotationHg38/inst/extdata/genomewide_GC_hg38.rds")
enhancer <- readRDS("/icbb/projects/share/annotations/methylTFRAnnotationHg38/inst/extdata/distal_regions.RDS")

# Define the TFs of interest
tfs <- c("JUN","RELB", "FOS","FOXP2", "BATF", "IRF4", "SP1", "FOSL2")
tf_bindsites <- tf_bindsites[tfs]

# Function to generate and save a plot for all samples
plot_and_save_difference<- function(samples, save_dir, obs_colors) {
  for (motif in names(tf_bindsites)) {
    if (!file.exists(file.path(save_dir, paste0("TF_footprint_diff_", motif, ".pdf")))) {
      logger.start(paste("Processing motif", motif, "for memoryTcell samples"))

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
        # Calculate observed methylation
        difference_data <- plot_data$plotDF
        difference_data[, type := paste(type, cell_type, sep= "_")]

        return(difference_data)
      })) 
      logger.completed()

      # Plot the combined difference data for all samples
      logger.info("Plotting observed and expected methylation in some samples")
      p_combined <- ggplot(combined_data, aes(x = x, y = avg_methyl, color = type)) +
        geom_line() + # Add lines for each type
        xlab("Distance from motif center") +
        ylab("Average methylation") +
        theme_classic() +
        ggtitle(paste("Expected vs observed footprints for ", motif)) +
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
  "Expected_TCM" = "#A9CFA4",
  "Observed_TCM" = "#1B7837",

  "Expected_TEM" = "#B0C4DE",
  "Observed_TEM" = "#32527B",

  "Expected_TN" = "#FDD9A0",
  "Observed_TN" = "#E66100",

  "Expected_TEMRA" = "#F4B6C2",
  "Observed_TEMRA" = "#AE017E"
)


# Generate and save plots for all samples with observed divided expected
plot_and_save_difference(msites, plot_dir, obs_colors)

#####################################################################
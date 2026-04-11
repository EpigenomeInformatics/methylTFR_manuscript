#!/usr/bin/env Rscript

#####################################################################
# 05_mFoots.R
# Created by IBG on 08-11-2025
# Create methylation footprint plots for selected TFs for CD4 Tcells
#####################################################################

set.seed(42)
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(methylTFR)
  library(methylTFRAnnotationHg38)
  library(logger)
  library(GenomicRanges)
  library(muLogR)
  library(stringr)
  library(RnBeads)
})

source("/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR_manuscript/src/run_mTFR_RnBeads.R", chdir = TRUE)

# Set the paths
main.dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/memoryTcells/"
sample_dir <- paste0(main.dir, "reports/")
debug <- paste0(main.dir, "debug/")
if(!dir.exists(debug)){dir.create(debug)}
plot_dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR_manuscript/figures/memtcells/"
plot_dir <- paste0(plot_dir, "mFoot/")
if(!dir.exists(plot_dir)){dir.create(plot_dir)}

# Load RnBeads preprocessed data
rnb_set <- RnBeads::load.rnb.set(paste0(sample_dir, "data_import_data/rnb.set_preprocessed"))

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
motifSet <- "JASPAR2020"
tf_bindsites <- getTFbindsites(motifSet)
motifSet <- "JASPAR2020_distal"
gcfreqs <- getGCfreq(motifSet)
gc_dist <- getGenomeGC("hg38")
distal <- readRDS("/icbb/projects/share/annotations/methylTFRAnnotationHg38/inst/extdata/distal_regions.RDS")
distal <- if(motifSet == "JASPAR2020_distal"){distal}else{NULL}


# Define the TFs of interest
tfs <- c("JUN","RELB", "FOS","FOXP2", "BATF", "IRF4", "SP1", "FOSL2", "RUNX3", "ETV1", "ETV2", "SPIB", "CENBP")
tf_bindsites <- tf_bindsites[names(tf_bindsites) %in% tfs]


# Function to generate and save a plot for all samples
plot_and_save_difference <- function(samples, save_dir, obs_colors) {

  for (motif in names(tf_bindsites)) {
    if (!file.exists(file.path(save_dir, paste0("TF_footprint_diff_", motif, ".pdf")))) {
      logger.start(paste("Processing motif", motif, "for memTcell samples"))

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
        flankNorm <- 50
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
  "Observed divided Expected TN" = "#C8E0B4",
  "Observed divided Expected TCM" = "#4492C6",
  "Observed divided Expected TEM" = "#43B6C4",
  "Observed divided Expected TEMRA" = "#898FB5"
)


# Generate and save plots for all samples with observed divided expected
plot_and_save_difference(msites, plot_dir, obs_colors)

#####################################################################
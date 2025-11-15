#!/usr/bin/env Rscript

#####################################################################
# 06_EO_mfoots.R
# Created by IBG on 07-11-2025
# Created for generating expected vs observed methylation footprint plots
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
main.dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/blueprint/"
sample_dir <- paste0(main.dir, "TB_RnBeads_271025/")
debug <- paste0(main.dir, "debug/")
if (!dir.exists(debug)) {
  dir.create(debug)
}
plot_dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR_manuscript/figures/blueprint/"
plot_dir <- paste0(plot_dir, "EO_mFoot/")
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir)
}

if (!file.exists(paste0(debug, "methylation_sites_merged_GRangesList.rds"))) {
  # Load RnBeads preprocessed data
  rnb_set <- RnBeads::load.rnb.set(paste0(sample_dir, "reports/data_import_data/rnb.set_preprocessed"))

  # Merging cell type replicates
  rnbset_merged <- mergeSamples(rnb_set, "cellTypeGroup")

  # Get the methylation sites as GRangesList
  msites <- rnb.RnBSet.to.GRangesList(rnbset_merged)
  saveRDS(msites, paste0(debug, "methylation_sites_merged_GRangesList.rds"))
} else {
  msites <- readRDS(paste0(debug, "methylation_sites_merged_GRangesList.rds"))
}

# Prepare motif data
motifSet <- "JASPAR2020"
tf_bindsites <- getTFbindsites(motifSet)
motifSet <- "JASPAR2020_distal"
gcfreqs <- getGCfreq(motifSet)
gc_dist <- getGenomeGC("hg38")

# Define the TFs of interest
tfs <- c(
  "SPIB", "SPI1", "TFAP2B", "PAX5", "PAX9", "PAX1", "NHLH2", "ASCL1", "NHLH1", "BHLHE22", "FERD3L", "PAX6",
  "EMX1", "PAX4", "EN1", "LHX1", "ELF4", "ELF2", "ETV5", "ETV6", "ELF5", "SPIC",
  "EHF", "ELF3", "IKZF1", "TFAP2C", "TFAP2B", "TFAP2A", "TFAP2E", "TFAP2C",
  "TGIF1", "CREB3L4", "PBX3", "TAL1::TCF3", "MYOG", "ATOH1", "MYF5", "BHLHA15", "ZBTB18", "EBF3",
  "EBF1", "TFAP4", "NEUROD1", "VSX2", "EGR4", "DPRX", "SOX8", "CUX1", "CUX2", "TEAD3"
)
tf_bindsites <- tf_bindsites[names(tf_bindsites) %in% tfs]
tf_bindsites <- tf_bindsites[sort(names(tf_bindsites))]
distal <- readRDS("/icbb/projects/share/annotations/methylTFRAnnotationHg38/inst/extdata/distal_regions.RDS")
distal <- if (motifSet == "JASPAR2020_distal") {
  distal
} else {
  NULL
}

# Function to generate and save a plot for all samples
plot_and_save_difference <- function(samples, save_dir, obs_colors) {
  for (motif in names(tf_bindsites)) {
    if (!file.exists(file.path(save_dir, paste0("TF_footprint_diff_", motif, "2.pdf")))) {
      logger.start(paste("Processing motif", motif, "for B and Tcell samples"))

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
        # Calculate observed methylation
        difference_data <- plot_data$plotDF
        difference_data[, type := paste(type, cell_type, sep = "_")]

        return(difference_data)
      }))
      # fwrite(combined_data, "combined_data2.csv")
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


# Define the observed colors
obs_colors <- c(
  "Expected_Bcell" = "#4492C6",
  "Observed_Bcell" = "#2C4C9B",
  "Expected_Tcell" = "#898FB5",
  "Observed_Tcell" = "#D66117"
)

# Generate and save plots for all samples with observed divided expected
plot_and_save_difference(msites, plot_dir, obs_colors)

#####################################################################

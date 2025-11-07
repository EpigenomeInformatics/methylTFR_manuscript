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
if(!dir.exists(debug)){dir.create(debug)}
plot_dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR_manuscript/figures/blueprint/"
plot_dir <- paste0(plot_dir, "EO_mFoot/")
if(!dir.exists(plot_dir)){dir.create(plot_dir)}

# Load RnBeads preprocessed data
rnb_set <- RnBeads::load.rnb.set(paste0(sample_dir, "reports/data_import_data/rnb.set_preprocessed"))

if(!file.exists(paste0(debug, "methylation_sites_merged_GRangesList.rds"))){
# Merging cell type replicates
rnbset_merged <- mergeSamples(rnb_set, "cellTypeGroup")

# Get the methylation sites as GRangesList
msites <- rnb.RnBSet.to.GRangesList(rnbset_merged)
saveRDS(msites, paste0(debug, "methylation_sites_merged_GRangesList.rds"))

}else{
    msites <- readRDS(paste0(debug, "methylation_sites_merged_GRangesList.rds"))
}

# Prepare motif data
motifSet <- "jaspar2020"
gcfreqs <- readRDS("/scratch/icbb/igunduz/methylTFR_manuscript/blueprint/debug/methylTFRAnnotationHg38/inst/extdata/JASPAR2020_distal_motif_gcfreq.rds")
tf_bindsites <- readRDS("/scratch/icbb/igunduz/methylTFR_manuscript/blueprint/debug/methylTFRAnnotationHg38/inst/extdata/JASPAR2020_tf_bindsites.rds")
gc_dist <- readRDS("/scratch/icbb/igunduz/methylTFR_manuscript/blueprint/debug/methylTFRAnnotationHg38/inst/extdata/genomewide_GC_hg38.rds")

# Define the TFs of interest
tfs <- c("FOXP2", "FOS", "BATF", "IRF4", "GCM1", "FOSL2")
tf_bindsites <- tf_bindsites[tfs]

if(!file.exists("/icbb/projects/share/annotations/methylTFRAnnotationHg38/inst/extdata/distal_regions.RDS")){
    distal <- fread("/icbb/projects/share/annotations/lolaDB/hg38/EnsemblRegBuildBP/regions/regionSet_1.bed", header = FALSE) 
    distal$V6 <- str_replace(distal$V6, ".", "*")
    distal <- GRanges(seqnames = distal$V1,
                  ranges = IRanges(start = distal$V2, end = distal$V3), 
                  strand = distal$V6)
    saveRDS(distal, "/icbb/projects/share/annotations/methylTFRAnnotationHg38/inst/extdata/distal_regions.RDS")
}else{
    distal <- readRDS("/icbb/projects/share/annotations/methylTFRAnnotationHg38/inst/extdata/distal_regions.RDS")
}

# Function to generate and save a plot for all samples
plot_and_save_difference <- function(samples, save_dir, obs_colors) {
  for (motif in names(tf_bindsites)) {
    if (!file.exists(file.path(save_dir, paste0("TF_footprint_diff_", motif, ".pdf")))) {
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
        difference_data[, type := paste(type, cell_type, sep= "_")]

        return(difference_data)
      })) 
      #fwrite(combined_data, "combined_data2.csv")
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

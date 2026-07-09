#!/usr/bin/env Rscript

#####################################################################
# 06_EO_mfoots.R
# Created by IBG on 07-11-2025
# Updated on 07-09-2026: switched from 2-group (Bcell/Tcell) to
#   4-group merging (Bcell_mem / Bcell_naive / TCD4 / TCD8)
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
plot_dir <- paste0(plot_dir, "EO_mFoot090626/")
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir)
}

if (!file.exists(paste0(debug, "methylation_sites_merged_GRangesList_celltype4G.rds"))) {
  # Load RnBeads preprocessed data
  rnb_set <- RnBeads::load.rnb.set(paste0(sample_dir, "reports/data_import_data/rnb.set_preprocessed"))

  # Create a new column in the phenotype data for 4-group cell types.
  # Treg is canonically a CD4+ subset, so it is folded into TCD4.
  cts <- as.character(rnb_set@pheno$cellTypeShort)
  cellType4Group <- dplyr::case_when(
    cts %in% c("Bcell_mem", "Bcell_gc")                       ~ "Bcell_mem",
    cts %in% c("Bcell_naive", "Bcell_pre")                    ~ "Bcell_naive",
    cts %in% c("TCD4", "TCD4_cenM", "TCD4_effM", "Treg")      ~ "TCD4",
    cts %in% c("TCD8", "TCD8_cenM", "TCD8_effM", "TCD8_term") ~ "TCD8",
    TRUE                                                      ~ NA_character_
  )
  stopifnot(!anyNA(cellType4Group)) # catch any unmapped cellTypeShort level
  rnb_set@pheno$cellType4Group <- cellType4Group

  # Merging cell type replicates into the 4 groups
  rnbset_merged <- mergeSamples(rnb_set, "cellType4Group")

  # Get the methylation sites as GRangesList
  msites <- rnb.RnBSet.to.GRangesList(rnbset_merged)
  saveRDS(msites, paste0(debug, "methylation_sites_merged_GRangesList_celltype4G.rds"))
} else {
  msites <- readRDS(paste0(debug, "methylation_sites_merged_GRangesList_celltype4G.rds"))
}

# Prepare motif data
motifSet <- "JASPAR2020"
tf_bindsites <- readRDS("/scratch/icbb/igunduz/methylTFR_manuscript/blueprint/debug/methylTFRAnnotationHg38/inst/extdata/JASPAR2020_tf_bindsites.rds")
motifSet <- "JASPAR2020_distal"
gcfreqs <- readRDS("/scratch/icbb/igunduz/methylTFR_manuscript/blueprint/debug/methylTFRAnnotationHg38/inst/extdata/JASPAR2020_distal_gc_freqs.rds")
gc_dist <- readRDS("/scratch/icbb/igunduz/methylTFR_manuscript/blueprint/debug/methylTFRAnnotationHg38/inst/extdata/genomewide_GC_hg38.rds")


# Define the TFs of interest
tfs <- c(
  "TEAD3","SPI1","SPIB", "TFAP2B", "PAX5", "PAX9", "PAX1", "NHLH2", "ASCL1", "NHLH1", "BHLHE22", "FERD3L", "PAX6",
  "EMX1", "PAX4", "EN1", "LHX1", "ELF4", "ELF2", "ETV5", "ETV6", "ELF5", "SPIC",
  "EHF", "ELF3", "IKZF1", "TFAP2C", "TFAP2B", "TFAP2A", "TFAP2E", "TFAP2C",
  "TGIF1", "CREB3L4", "PBX3", "TAL1::TCF3", "MYOG", "ATOH1", "MYF5", "BHLHA15", "ZBTB18", "EBF3",
  "EBF1", "TFAP4", "NEUROD1", "VSX2", "EGR4", "DPRX", "SOX8", "CUX1", "CUX2"
)[1:2]
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
      logger.start(paste("Processing motif", motif, "for Bcell_mem, Bcell_naive, TCD4 and TCD8 samples"))

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
# Bcell_mem/Bcell_naive stay in a shared blue family (naive lighter, mem
# darker); TCD4/TCD8 get their own distinct families (orange vs green)
# so the two T-cell subsets are readable apart from each other and from
# the B-cell groups.
obs_colors <- c(
  "Expected_Bcell_mem"    = "#4492C6", # medium blue
  "Observed_Bcell_mem"    = "#2C4C9B", # dark blue
  "Expected_Bcell_naive"  = "#7FB8E0", # light blue (close to mem, lighter)
  "Observed_Bcell_naive"  = "#1B3570", # navy blue (close to mem, deeper)
  "Expected_TCD4"         = "#F2B880", # light orange
  "Observed_TCD4"         = "#D66117", # dark orange
  "Expected_TCD8"         = "#A8D5A2", # light green
  "Observed_TCD8"         = "#2E7D32"  # dark green
)

# Generate and save plots for all samples with observed divided expected
plot_and_save_difference(msites, plot_dir, obs_colors)

#####################################################################
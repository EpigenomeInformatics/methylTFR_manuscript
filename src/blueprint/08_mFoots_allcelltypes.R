#!/usr/bin/env Rscript

#####################################################################
# 08_mFoots_allcelltypes.R
# Created by IBG on 07-09-2026
# Combined methylation footprint plots for all cell type groups:
#   Monocytes / Granulocytes / Bcell_mem / Bcell_naive / Tcell
# Metric: Observed - Expected (subtraction, not ratio)
#
# Group mapping (from table(cts)):
#   Bcell_mem    = Bcell_mem, Bcell_gc
#   Bcell_naive  = Bcell_naive, Bcell_pre
#   Tcell        = TCD4, TCD4_cenM, TCD4_effM, TCD8, TCD8_cenM,
#                  TCD8_effM, TCD8_term, Treg
#   Monocytes    = Mono
#   Granulocytes = Neut_band, Neut_mat, Neut_meta, Neut_myel,
#                  Neut_segm, Eosi_mat
# All other cellTypeShort levels (DC*, Endo*, Eryth, Mac_f0/f1/f2, MK,
# MPP, MSC_BM, NK, Ost, Plas, Thym*) are excluded from this analysis --
# macrophages are a distinct differentiated lineage from monocytes, and
# thymocytes are immature precursors rather than mature peripheral T
# cells, so both are dropped rather than folded into Monocytes/Tcell.
#
# TODO before running:
#   1. Set `sample_dir` below to the correct RnBeads report directory
#      for this dataset.
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
# TODO: replace with the sample_dir for the dataset that has Monocytes/Granulocytes
sample_dir <- paste0(main.dir, "TB_RnBeads_271025/")
debug <- paste0(main.dir, "debug/")
if (!dir.exists(debug)) {
  dir.create(debug)
}
plot_dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR_manuscript/figures/blueprint/"
plot_dir <- paste0(plot_dir, "mFoot_allcelltypes_090626/")
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir)
}

if (!file.exists(paste0(debug, "methylation_sites_merged_GRangesList_celltype5G.rds"))) {
  # Load RnBeads preprocessed data
  rnb_set <- RnBeads::load.rnb.set(paste0(sample_dir, "reports/data_import_data/rnb.set_preprocessed"))

  # Create a new column in the phenotype data for the 5-group cell types.
  # Bcell_gc (germinal center, memory-lineage) is folded into Bcell_mem;
  # Bcell_pre (developmentally upstream of naive) is folded into
  # Bcell_naive; all T-cell subtypes collapse into a single Tcell group.
  # Everything not covered by these 5 groups (see mapping note above)
  # is left as NA and dropped before merging.
  cts <- as.character(rnb_set@pheno$cellTypeShort)
  cellType5Group <- dplyr::case_when(
    cts %in% c("Bcell_mem", "Bcell_gc")                       ~ "Bcell_mem",
    cts %in% c("Bcell_naive", "Bcell_pre")                    ~ "Bcell_naive",
    cts %in% c("TCD4", "TCD4_cenM", "TCD4_effM", "TCD8",
               "TCD8_cenM", "TCD8_effM", "TCD8_term", "Treg") ~ "Tcell",
    cts == "Mono"                                             ~ "Monocytes",
    cts %in% c("Neut_band", "Neut_mat", "Neut_meta",
               "Neut_myel", "Neut_segm", "Eosi_mat")          ~ "Granulocytes",
    TRUE                                                      ~ NA_character_
  )
  rnb_set@pheno$cellType5Group <- cellType5Group

  # Drop samples that fall outside the 5 target groups (e.g. DC, Endo,
  # Eryth, Mac_f0/f1/f2, MK, MPP, MSC_BM, NK, Ost, Plas, Thym*)
  drop_idx <- which(is.na(cellType5Group))
  if (length(drop_idx) > 0) {
    logger.info(paste0(
      length(drop_idx), " samples excluded (cell types outside the 5 target groups): ",
      paste(sort(unique(cts[drop_idx])), collapse = ", ")
    ))
    rnb_set <- RnBeads::remove.samples(rnb_set, drop_idx)
  }

  # Merging cell type replicates into the 5 groups
  rnbset_merged <- mergeSamples(rnb_set, "cellType5Group")

  # Get the methylation sites as GRangesList
  msites <- rnb.RnBSet.to.GRangesList(rnbset_merged)
  saveRDS(msites, paste0(debug, "methylation_sites_merged_GRangesList_celltype5G.rds"))
} else {
  msites <- readRDS(paste0(debug, "methylation_sites_merged_GRangesList_celltype5G.rds"))
}

# Prepare motif data
motifSet <- "JASPAR2020"
tf_bindsites <- readRDS("/scratch/icbb/igunduz/methylTFR_manuscript/blueprint/debug/methylTFRAnnotationHg38/inst/extdata/JASPAR2020_tf_bindsites.rds")
motifSet <- "JASPAR2020_distal"
gcfreqs <- readRDS("/scratch/icbb/igunduz/methylTFR_manuscript/blueprint/debug/methylTFRAnnotationHg38/inst/extdata/JASPAR2020_distal_gc_freqs.rds")
gc_dist <- readRDS("/scratch/icbb/igunduz/methylTFR_manuscript/blueprint/debug/methylTFRAnnotationHg38/inst/extdata/genomewide_GC_hg38.rds")


# Define the TFs of interest
tfs <- c(
  "TFAP2B", "PAX5", "PAX9", "PAX1", "NHLH2", "ASCL1", "NHLH1", "BHLHE22", "FERD3L", "PAX6",
  "EMX1", "PAX4", "EN1", "LHX1", "ELF4", "ELF2", "ETV5", "ETV6", "ELF5", "SPIC",
  "SPIB", "SPI1", "EHF", "ELF3", "IKZF1", "TFAP2C", "TFAP2B", "TFAP2A", "TFAP2E", "TFAP2C",
  "TGIF1", "CREB3L4", "PBX3", "TAL1::TCF3", "MYOG", "ATOH1", "MYF5", "BHLHA15", "ZBTB18", "EBF3",
  "EBF1", "TFAP4", "NEUROD1", "VSX2", "EGR4", "DPRX", "SOX8", "CUX1", "CUX2", "TEAD3", "CEBPB"
)
tf_bindsites <- tf_bindsites[names(tf_bindsites) %in% tfs]

distal <- readRDS("/icbb/projects/share/annotations/methylTFRAnnotationHg38/inst/extdata/distal_regions.RDS")
distal <- if (motifSet == "JASPAR2020_distal") {
  distal
} else {
  NULL
}

# Function to generate and save a plot for all samples
plot_and_save_difference <- function(samples, save_dir, obs_colors) {
  for (motif in names(tf_bindsites)) {
    if (!file.exists(file.path(save_dir, paste0("TF_footprint_diff_", motif, ".pdf")))) {
      logger.start(paste("Processing motif", motif, "for Monocytes, Granulocytes, Bcell_mem, Bcell_naive and Tcell samples"))

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
        # Calculate observed minus expected methylation
        difference_data <- plot_data$plotDF[, .(avg_methyl = avg_methyl[type == "Observed"] - avg_methyl[type == "Expected"]), by = x]
        difference_data[, type := cell_type]

        # Normalize by subtracting the flanking-region baseline, so the
        # flanks sit at 0 (analogous to the ratio version centering on 1)
        flankNorm <- 30
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
        geom_line() + # Add lines for each group
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

# Define the colors for all 5 groups.
# Monocytes/Granulocytes/Tcell colors follow the reference figure;
# Bcell_mem keeps the reference figure's B-cell maroon, and Bcell_naive
# is a separate, lighter rose shade in the same family so the two
# B-cell states read as related but distinguishable.
obs_colors <- c(
  "Monocytes"    = "#C2703D", # burnt orange
  "Granulocytes" = "#E8A33D", # gold/amber
  "Bcell_mem"    = "#7B1E3D", # dark maroon (reference B-cell color)
  "Bcell_naive"  = "#C46B8C", # lighter rose, same family as Bcell_mem
  "Tcell"        = "#4FC3D9"  # cyan/turquoise
)

# Generate and save plots for all samples with observed minus expected
plot_and_save_difference(msites, plot_dir, obs_colors)

#####################################################################
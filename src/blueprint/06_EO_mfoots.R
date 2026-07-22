#!/usr/bin/env Rscript

#####################################################################
# 06_EO_mfoots.R
# Created by IBG on 07-11-2025
# Updated on 07-09-2026: switched from 2-group (Bcell/Tcell) to
#   4-group merging (Bcell_mem / Bcell_naive / TCD4 / TCD8)
# Created for generating expected vs observed methylation footprint plots
# Added mean deviation scores to legends (separated by cell type)
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
plot_dir <- paste0(plot_dir, "EO_mFoot_210726/")
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir)
}

# -------------------------------------------------------------------
# Load and Process Deviation Scores
# -------------------------------------------------------------------
logger.info("Loading and filtering deviation scores...")
dev_obj <- readRDS("/scratch/icbb/igunduz/methylTFR_manuscript/blueprint/mTFR_devs_121125/JASPAR2020_distal_deviations.RDS")
dev_matrix <- deviations(dev_obj)

# Load annotation (using the sample_dir annotation to ensure phenotypes match)
sannot <- read.csv(paste0(sample_dir, "reports/data_import_data/annotation.csv"), stringsAsFactors = FALSE)
sannot$bedFile <- as.character(sannot$bedFile)
sample_names <- colnames(dev_matrix)

# Match sannot rows to deviation matrix columns
matched_sannot <- sannot[match(sample_names, sannot$bedFile), ]

# Remove disease samples from deviations
disease_groups <- c("Multiple Myeloma", "Acute Lymphocytic Leukemia")
keep_dev <- !matched_sannot$cellTypeGroup %in% disease_groups

dev_matrix <- dev_matrix[, keep_dev, drop = FALSE]
matched_sannot <- matched_sannot[keep_dev, ]

# Map the kept deviation samples to the 4-group classification
matched_sannot$cellType4Group <- dplyr::case_when(
  matched_sannot$cellTypeShort %in% c("Bcell_mem", "Bcell_gc")                       ~ "Bcell_mem",
  matched_sannot$cellTypeShort %in% c("Bcell_naive", "Bcell_pre")                    ~ "Bcell_naive",
  matched_sannot$cellTypeShort %in% c("TCD4", "TCD4_cenM", "TCD4_effM", "Treg")      ~ "TCD4",
  matched_sannot$cellTypeShort %in% c("TCD8", "TCD8_cenM", "TCD8_effM", "TCD8_term") ~ "TCD8",
  TRUE                                                                               ~ NA_character_
)
dev_groups <- matched_sannot$cellType4Group

# -------------------------------------------------------------------
# Load and Process RnBeads Methylation Data
# -------------------------------------------------------------------
if (!file.exists(paste0(debug, "methylation_sites_merged_GRangesList_celltype_220726.rds"))) {
  # Load RnBeads preprocessed data
  rnb_set <- RnBeads::load.rnb.set(paste0(sample_dir, "reports/data_import_data/rnb.set_preprocessed"))

  # Drop disease samples BEFORE grouping/merging
  disease_idx <- which(as.character(rnb_set@pheno$cellTypeGroup) %in% disease_groups)
  if (length(disease_idx) > 0) {
    rnb_set <- remove.samples(rnb_set, disease_idx)
  }

  # Create a new column in the phenotype data for 4-group cell types.
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
  saveRDS(msites, paste0(debug, "methylation_sites_merged_GRangesList_celltype_220726.rds"))
} else {
  msites <- readRDS(paste0(debug, "methylation_sites_merged_GRangesList_celltype_220726.rds"))
}

# -------------------------------------------------------------------
# Prepare Motifs and TFs
# -------------------------------------------------------------------
motifSet <- "JASPAR2020"
tf_bindsites <- readRDS("/scratch/icbb/igunduz/methylTFR_manuscript/blueprint/debug/methylTFRAnnotationHg38/inst/extdata/JASPAR2020_tf_bindsites.rds")
motifSet <- "JASPAR2020_distal"
gcfreqs <- readRDS("/scratch/icbb/igunduz/methylTFR_manuscript/blueprint/debug/methylTFRAnnotationHg38/inst/extdata/JASPAR2020_DISTAL_motif_gcfreq.rds")
gc_dist <- readRDS("/scratch/icbb/igunduz/methylTFR_manuscript/blueprint/debug/methylTFRAnnotationHg38/inst/extdata/genomewide_GC_hg38.rds")

# Define the TFs of interest
tfs <- c(
  "TEAD3","SPI1","SPIB", "TFAP2B", "PAX5", "PAX9", "PAX1", "NHLH2", "ASCL1", "NHLH1", "BHLHE22", "FERD3L", "PAX6",
  "EMX1", "PAX4", "EN1", "LHX1", "ELF4", "ELF2", "ETV5", "ETV6", "ELF5", "SPIC",
  "EHF", "ELF3", "IKZF1", "TFAP2C", "TFAP2B", "TFAP2A", "TFAP2E", "TFAP2C",
  "TGIF1", "CREB3L4", "PBX3", "TAL1::TCF3", "MYOG", "ATOH1", "MYF5", "BHLHA15", "ZBTB18", "EBF3",
  "EBF1", "TFAP4", "NEUROD1", "VSX2", "EGR4", "DPRX", "SOX8", "CUX1", "CUX2"
)
tf_bindsites <- tf_bindsites[names(tf_bindsites) %in% tfs]
tf_bindsites <- tf_bindsites[sort(names(tf_bindsites))]

distal <- readRDS("/icbb/projects/share/annotations/methylTFRAnnotationHg38/inst/extdata/distal_regions.RDS")
distal <- if (motifSet == "JASPAR2020_distal") {
  distal
} else {
  NULL
}

# -------------------------------------------------------------------
# Plotting Function
# -------------------------------------------------------------------
base_colors <- c(
  "Expected_Bcell_mem"    = "#4492C6", 
  "Observed_Bcell_mem"    = "#2C4C9B", 
  "Expected_Bcell_naive"  = "#7FB8E0", 
  "Observed_Bcell_naive"  = "#1B3570", 
  "Expected_TCD4"         = "#F2B880", 
  "Observed_TCD4"         = "#D66117", 
  "Expected_TCD8"         = "#A8D5A2", 
  "Observed_TCD8"         = "#2E7D32"  
)

plot_and_save_difference <- function(samples, save_dir, base_colors, dev_mat, dev_grps) {
  for (motif in names(tf_bindsites)) {
    # Fixed the check to standard ".pdf" so it correctly skips existing plots
    if (!file.exists(file.path(save_dir, paste0("TF_footprint_diff_", motif, ".pdf")))) {
      logger.start(paste("Processing motif", motif, "for Bcell_mem, Bcell_naive, TCD4 and TCD8 samples"))

      # Find matching row in deviation matrix
      motif_idx <- which(rownames(dev_mat) == motif)
      if (length(motif_idx) == 0) {
        motif_idx <- grep(paste0("\\b", motif, "\\b"), rownames(dev_mat), ignore.case = TRUE)
      }
      
      # Calculate mean deviation for each cell type
      mean_devs <- list()
      for (g in names(samples)) {
        if (length(motif_idx) > 0) {
          mean_devs[[g]] <- mean(dev_mat[motif_idx, dev_grps == g], na.rm = TRUE)
        } else {
          mean_devs[[g]] <- NA
        }
      }

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
        
        difference_data <- plot_data$plotDF
        
        # Format the legend label to include the mean deviation score
        dev_val <- round(mean_devs[[cell_type]], 2)
        dev_str <- ifelse(is.na(dev_val), "N/A", dev_val)
        
        # original type is "Observed" or "Expected"
        difference_data[, original_group := paste(type, cell_type, sep = "_")]
        difference_data[, type := paste0(type, "_", cell_type, " (Mean Dev: ", dev_str, ")")]

        return(difference_data)
      }))
      logger.completed()

      # Map dynamic labels to base colors
      dynamic_colors <- setNames(base_colors[combined_data$original_group], combined_data$type)
      dynamic_colors <- dynamic_colors[!duplicated(names(dynamic_colors))]
      
      # Ensure factor levels keep consistent legend order
      combined_data[, type := factor(type, levels = names(dynamic_colors))]

      # Plot the combined difference data for all samples
      logger.info("Plotting observed and expected methylation in some samples")
      p_combined <- ggplot(combined_data, aes(x = x, y = avg_methyl, color = type)) +
        geom_line() + 
        xlab("Distance from motif center") +
        ylab("Average methylation") +
        theme_classic() +
        ggtitle(paste("Expected vs observed footprints for ", motif)) +
        scale_color_manual(values = dynamic_colors) + 
        theme(legend.position = "bottom") +
        xlim(-200, 200) 

      # Save the plot as a PDF
      ggsave(
        filename = file.path(save_dir, paste0("TF_footprint_diff_", motif, ".pdf")),
        plot = p_combined, width = 12, height = 8
      )
    }
  }
}

# -------------------------------------------------------------------
# Generate and Save Plots
# -------------------------------------------------------------------
plot_and_save_difference(msites, plot_dir, base_colors, dev_matrix, dev_groups)
#!/usr/bin/env Rscript

#####################################################################
# 07_mFoots.R
# Created by IBG on 07-11-2025
# Updated on 07-09-2026: switched from 2-group (Bcell/Tcell) to
#   3-group merging (Bcell_mem / Bcell_naive / Tcell)
# Create methylation footprint plots for selected TFs
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
sample_dir <- paste0(main.dir, "RnBeads_291025/")
debug <- paste0(main.dir, "debug/")
if (!dir.exists(debug)) {
  dir.create(debug)
}
plot_dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR_manuscript/figures/blueprint/"
plot_dir <- paste0(plot_dir, "mFoot_allcells_220726/")
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir)
}

# -------------------------------------------------------------------
# Load and Process Deviation Scores
# -------------------------------------------------------------------
logger.info("Loading and filtering deviation scores...")
dev_obj <- readRDS("/scratch/icbb/igunduz/methylTFR_manuscript/blueprint/mTFR_devs_121125/JASPAR2020_distal_deviations.RDS")
dev_matrix <- deviations(dev_obj)

sannot <- read.csv(paste0(sample_dir, "reports/data_import_data/annotation.csv"), stringsAsFactors = FALSE)
sannot$bedFile <- as.character(sannot$bedFile)
sample_names <- colnames(dev_matrix)

# Match sannot rows to deviation matrix columns
matched_sannot <- sannot[match(sample_names, sannot$bedFile), ]

# Remove disease samples from deviations (keep only non-diseased "None" samples)
keep_dev <- matched_sannot$DISEASE == "None"

dev_matrix <- dev_matrix[, keep_dev, drop = FALSE]
matched_sannot <- matched_sannot[keep_dev, ]

# Map the kept deviation samples to the 3-group classification
matched_sannot$cellType3Group <- dplyr::case_when(
  matched_sannot$cellTypeShort %in% c("Bcell_mem", "Bcell_gc")    ~ "Bcell_mem",
  matched_sannot$cellTypeShort %in% c("Bcell_naive", "Bcell_pre") ~ "Bcell_naive",
  TRUE                                                            ~ "Tcell"
)
dev_groups <- matched_sannot$cellType3Group

# -------------------------------------------------------------------
# Load and Process RnBeads Methylation Data
# -------------------------------------------------------------------
if (!file.exists(paste0(debug, "methylation_sites_merged_GRangesList_celltypes_220726.rds"))) {
  logger.info("Loading RnBeads preprocessed data...")
  rnb_set <- RnBeads::load.rnb.set(paste0(sample_dir, "reports/data_import_data/rnb.set_preprocessed"))
  
  # Drop disease samples BEFORE grouping/merging (keep only DISEASE == "None")
  disease_idx <- which(as.character(rnb_set@pheno$DISEASE) != "None")
  if (length(disease_idx) > 0) {
    rnb_set <- remove.samples(rnb_set, disease_idx)
  }
  
  cts <- as.character(rnb_set@pheno$cellTypeShort)
  cellType3Group <- dplyr::case_when(
    cts %in% c("Bcell_mem", "Bcell_gc")    ~ "Bcell_mem",
    cts %in% c("Bcell_naive", "Bcell_pre") ~ "Bcell_naive",
    TRUE                                    ~ "Tcell"
  )
  rnb_set@pheno$cellType3Group <- cellType3Group

  logger.info("Merging samples into 3 groups...")
  rnbset_merged <- mergeSamples(rnb_set, "cellType3Group")

  msites <- rnb.RnBSet.to.GRangesList(rnbset_merged)
  saveRDS(msites, paste0(debug, "methylation_sites_merged_GRangesList_celltypes_220726.rds"))
} else {
  logger.info("Loading previously merged methylation GRangesList...")
  msites <- readRDS(paste0(debug, "methylation_sites_merged_GRangesList_celltypes_220726.rds"))
}

# -------------------------------------------------------------------
# Prepare Motifs and TFs
# -------------------------------------------------------------------
motifSet <- "JASPAR2020"
tf_bindsites <- getTFbindsites(motifSet)
motifSet_distal <- "JASPAR2020_distal"
gcfreqs <- getGCfreq(motifSet_distal)
gc_dist <- getGenomeGC("hg38")

tfs <- c(
  "TFAP2B", "PAX5", "PAX9", "PAX1", "NHLH2", "ASCL1", "NHLH1", "BHLHE22", "FERD3L", "PAX6",
  "EMX1", "PAX4", "EN1", "LHX1", "ELF4", "ELF2", "ETV5", "ETV6", "ELF5", "SPIC",
  "SPIB", "SPI1", "EHF", "ELF3", "IKZF1", "TFAP2C", "TFAP2B", "TFAP2A", "TFAP2E", "TFAP2C",
  "TGIF1", "CREB3L4", "PBX3", "TAL1::TCF3", "MYOG", "ATOH1", "MYF5", "BHLHA15", "ZBTB18", "EBF3",
  "EBF1", "TFAP4", "NEUROD1", "VSX2", "EGR4", "DPRX", "SOX8", "CUX1", "CUX2", "TEAD3"
)
tf_bindsites <- tf_bindsites[names(tf_bindsites) %in% tfs]

distal <- readRDS("/icbb/projects/share/annotations/methylTFRAnnotationHg38/inst/extdata/distal_regions.RDS")

# -------------------------------------------------------------------
# Plotting Function
# -------------------------------------------------------------------
# Defines base colors assigned to the groups
base_colors <- c(
  "Bcell_mem"   = "#2C4C9B",
  "Bcell_naive" = "#5B9BD5",
  "Tcell"       = "#D66117"
)

plot_and_save_difference <- function(samples, save_dir, base_colors, dev_mat, dev_grps) {
  for (motif in names(tf_bindsites)) {
    if (!file.exists(file.path(save_dir, paste0("TF_footprint_diff_", motif, ".pdf")))) {
      logger.start(paste("Processing motif", motif, "for Bcell_mem, Bcell_naive and Tcell samples"))

      # Find matching row in deviation matrix (exact match or word-boundary partial match)
      motif_idx <- which(rownames(dev_mat) == motif)
      if (length(motif_idx) == 0) {
        motif_idx <- grep(paste0("\\b", motif, "\\b"), rownames(dev_mat), ignore.case = TRUE)
      }
      
      # Calculate mean deviation for each cell type
      mean_devs <- list()
      for (g in names(samples)) {
        if (length(motif_idx) > 0) {
          # Averages across matched rows (if multiple) and group columns
          mean_devs[[g]] <- mean(dev_mat[motif_idx, dev_grps == g], na.rm = TRUE)
        } else {
          mean_devs[[g]] <- NA
        }
      }

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
        
        difference_data <- plot_data$plotDF[, .(avg_methyl = avg_methyl[type == "Observed"] / avg_methyl[type == "Expected"]), by = x]
        
        # Format the legend label to include the mean deviation score
        dev_val <- round(mean_devs[[cell_type]], 2)
        dev_str <- ifelse(is.na(dev_val), "N/A", dev_val)
        label_str <- paste0("Obs/Exp ", cell_type, " (Mean Dev: ", dev_str, ")")
        
        difference_data[, type := label_str]
        difference_data[, original_group := cell_type]

        flankNorm <- 30
        flank <- max(abs(difference_data$x), na.rm = TRUE)
        idx <- abs(difference_data$x) >= flank - flankNorm
        norm_factor <- mean(difference_data$avg_methyl[idx], na.rm = TRUE)
        difference_data[, avg_methyl := avg_methyl / norm_factor]

        return(difference_data)
      }))
      logger.completed()

      # Map dynamic labels to base colors
      dynamic_colors <- setNames(base_colors[combined_data$original_group], combined_data$type)
      dynamic_colors <- dynamic_colors[!duplicated(names(dynamic_colors))]
      
      # Ensure factor levels keep consistent legend order
      combined_data[, type := factor(type, levels = names(dynamic_colors))]

      logger.info("Plotting combined difference data for all samples")
      p_combined <- ggplot(combined_data, aes(x = x, y = avg_methyl, color = type)) +
        geom_line() +
        xlab("Distance from motif center") +
        ylab("Methylation difference (Observed / Expected)") +
        theme_classic() +
        ggtitle(paste("TF footprint difference for", motif)) +
        scale_color_manual(values = dynamic_colors) +
        theme(legend.position = "bottom") +
        xlim(-200, 200)

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
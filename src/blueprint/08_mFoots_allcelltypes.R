#!/usr/bin/env Rscript

#####################################################################
# 08_mFoots_allcelltypes.R
# Created by IBG on 07-09-2026
# Combined methylation footprint plots for all cell type groups:
#   Monocytes / Granulocytes / Bcell_mem / Bcell_naive / Tcell
# Metric: Observed - Expected (subtraction, not ratio)
# Added mean deviation scores to legends (separated by cell type)
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
# MPP, MSC_BM, NK, Ost, Plas, Thym*) are excluded from this analysis.
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

# Updated output folder date
plot_dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR_manuscript/figures/blueprint/"
plot_dir <- paste0(plot_dir, "mFoot_allcelltypes_220726/")
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

# Map the deviation samples to the 5-group classification
matched_sannot$cellType5Group <- dplyr::case_when(
  matched_sannot$cellTypeShort %in% c("Bcell_mem", "Bcell_gc")                       ~ "Bcell_mem",
  matched_sannot$cellTypeShort %in% c("Bcell_naive", "Bcell_pre")                    ~ "Bcell_naive",
  matched_sannot$cellTypeShort %in% c("TCD4", "TCD4_cenM", "TCD4_effM", "TCD8",
                                      "TCD8_cenM", "TCD8_effM", "TCD8_term", "Treg") ~ "Tcell",
  matched_sannot$cellTypeShort == "Mono"                                             ~ "Monocytes",
  matched_sannot$cellTypeShort %in% c("Neut_band", "Neut_mat", "Neut_meta",
                                      "Neut_myel", "Neut_segm", "Eosi_mat")          ~ "Granulocytes",
  TRUE                                                                               ~ NA_character_
)

# Remove disease samples using the DISEASE column
keep_dev <- (matched_sannot$DISEASE == "None") & (!is.na(matched_sannot$cellType5Group))

dev_matrix <- dev_matrix[, keep_dev, drop = FALSE]
matched_sannot <- matched_sannot[keep_dev, ]
dev_groups <- matched_sannot$cellType5Group

# -------------------------------------------------------------------
# Load and Process RnBeads Methylation Data
# -------------------------------------------------------------------
if (!file.exists(paste0(debug, "methylation_sites_merged_GRangesList_celltype5G_220726.rds"))) {
  logger.info("Loading RnBeads preprocessed data...")
  rnb_set <- RnBeads::load.rnb.set(paste0(sample_dir, "reports/data_import_data/rnb.set_preprocessed"))

  # Remove DISEASE samples
  disease_groups <- c("Multiple Myeloma", "Acute Lymphocytic Leukemia")
  disease_idx <- which(as.character(rnb_set@pheno$DISEASE) %in% disease_groups)
  if (length(disease_idx) > 0) {
    rnb_set <- remove.samples(rnb_set, disease_idx)
  }

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

  # Drop samples using the DISEASE column and those outside the 5 target groups
  disease_idx <- which(as.character(rnb_set@pheno$DISEASE) != "None")
  group_drop_idx <- which(is.na(cellType5Group))
  
  drop_idx <- unique(c(disease_idx, group_drop_idx))
  
  if (length(drop_idx) > 0) {
    logger.info(paste0(
      length(drop_idx), " samples excluded (disease samples or outside the 5 target groups)."
    ))
    rnb_set <- RnBeads::remove.samples(rnb_set, drop_idx)
  }

  rnbset_merged <- mergeSamples(rnb_set, "cellType5Group")
  msites <- rnb.RnBSet.to.GRangesList(rnbset_merged)
  saveRDS(msites, paste0(debug, "methylation_sites_merged_GRangesList_celltype5G_220726.rds"))
} else {
  msites <- readRDS(paste0(debug, "methylation_sites_merged_GRangesList_celltype5G_220726.rds"))
}

# -------------------------------------------------------------------
# Prepare Motifs and TFs
# -------------------------------------------------------------------
motifSet <- "JASPAR2020"
tf_bindsites <- readRDS("/scratch/icbb/igunduz/methylTFR_manuscript/blueprint/debug/methylTFRAnnotationHg38/inst/extdata/JASPAR2020_tf_bindsites.rds")
motifSet <- "JASPAR2020_distal"
gcfreqs <- readRDS("/scratch/icbb/igunduz/methylTFR_manuscript/blueprint/debug/methylTFRAnnotationHg38/inst/extdata/JASPAR2020_DISTAL_motif_gcfreq.rds")
gc_dist <- readRDS("/scratch/icbb/igunduz/methylTFR_manuscript/blueprint/debug/methylTFRAnnotationHg38/inst/extdata/genomewide_GC_hg38.rds")

tfs <- c(
  "TFAP2B", "PAX5", "PAX9", "PAX1", "NHLH2", "ASCL1", "NHLH1", "BHLHE22", "FERD3L", "PAX6",
  "EMX1", "PAX4", "EN1", "LHX1", "ELF4", "ELF2", "ETV5", "ETV6", "ELF5", "SPIC",
  "SPIB", "SPI1", "EHF", "ELF3", "IKZF1", "TFAP2C", "TFAP2B", "TFAP2A", "TFAP2E", "TFAP2C",
  "TGIF1", "CREB3L4", "PBX3", "TAL1::TCF3", "MYOG", "ATOH1", "MYF5", "BHLHA15", "ZBTB18", "EBF3",
  "EBF1", "TFAP4", "NEUROD1", "VSX2", "EGR4", "DPRX", "SOX8", "CUX1", "CUX2", "TEAD3", "CEBPB",
  "FOSL1::JUND","CEBPB"
)
tf_bindsites <- tf_bindsites[names(tf_bindsites) %in% tfs]

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
  "Monocytes"    = "#C2703D", 
  "Granulocytes" = "#E8A33D", 
  "Bcell_mem"    = "#7B1E3D", 
  "Bcell_naive"  = "#C46B8C", 
  "Tcell"        = "#4FC3D9"  
)

plot_and_save_difference <- function(samples, save_dir, base_colors, dev_mat, dev_grps) {
  for (motif in names(tf_bindsites)) {
    if (!file.exists(file.path(save_dir, paste0("TF_footprint_diff_", motif, ".pdf")))) {
      logger.start(paste("Processing motif", motif, "for all 5 groups"))

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
        
        # Calculate observed minus expected methylation
        difference_data <- plot_data$plotDF[, .(avg_methyl = avg_methyl[type == "Observed"] - avg_methyl[type == "Expected"]), by = x]
        
        # Format the legend label to include the mean deviation score
        dev_val <- round(mean_devs[[cell_type]], 2)
        dev_str <- ifelse(is.na(dev_val), "N/A", dev_val)
        label_str <- paste0(cell_type, " (Mean Dev: ", dev_str, ")")
        
        difference_data[, type := label_str]
        difference_data[, original_group := cell_type]

        # Normalize by subtracting the flanking-region baseline
        flankNorm <- 30
        flank <- max(abs(difference_data$x), na.rm = TRUE)
        idx <- abs(difference_data$x) >= flank - flankNorm
        norm_baseline <- mean(difference_data$avg_methyl[idx], na.rm = TRUE)
        difference_data[, avg_methyl := avg_methyl - norm_baseline]

        return(difference_data)
      }))
      logger.completed()

      # Map dynamic labels to base colors
      dynamic_colors <- setNames(base_colors[combined_data$original_group], combined_data$type)
      dynamic_colors <- dynamic_colors[!duplicated(names(dynamic_colors))]
      
      # Ensure factor levels keep consistent legend order
      combined_data[, type := factor(type, levels = names(dynamic_colors))]

      # Plot the combined difference data for all samples
      logger.info("Plotting combined difference data for all samples")
      p_combined <- ggplot(combined_data, aes(x = x, y = avg_methyl, color = type)) +
        geom_line() + 
        xlab("Distance from motif center") +
        ylab("Methylation difference (Observed - Expected)") +
        theme_classic() +
        ggtitle(paste("TF footprint difference for", motif)) +
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
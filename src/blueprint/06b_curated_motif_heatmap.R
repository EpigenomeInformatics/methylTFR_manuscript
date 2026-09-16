#!/usr/bin/env Rscript

#####################################################################
# 06b_curated_motif_heatmap.R
# created on 16-09-2026 by Irem B Gunduz
# Heatmap of the Z-scores for the curated AP-1, ETS, SPI and C/EBP
# motif panel (Figure C) across the main cell type groups, with the
# blueprint cell type colours and the cell labels along the bottom
#####################################################################

suppressPackageStartupMessages({
  library(grid)
  library(circlize)
  library(ComplexHeatmap)
  library(logger)
  library(methylTFR)
})
set.seed(42)

#####################################################################
# Settings
#####################################################################

motifSets <- c("jaspar2020", "jaspar2020_distal")

# Curated motif panel of Figure C, in the AP-1, ETS, SPI, C/EBP order
curated.motifs <- c(
  "JDP2", "FOS::JUNB", "FOSL2::JUNB", "FOS::JUN", "FOSL2::JUN", "FOS::JUND",
  "FOS", "JUN(var.2)", "FOSL1", "JUND", "FOSL1::JUNB", "FOSL1::JUND",
  "JUNB", "BATF3", "BATF::JUN", "BATF", "NFE2L1",
  "FLI1", "ERG", "ETS1", "ERF", "ETV2", "ETV1", "ELF1", "GABPA", "ETV4",
  "IKZF1", "ELF3", "EHF", "ELK4", "ELF4", "ETV5", "ETV6", "ELF5",
  "SPI1", "SPIB", "SPIC",
  "TEF", "DBP", "CEBPB", "CEBPE", "CEBPG", "GMEB2", "HLF", "NFIL3",
  "ATF7", "CEBPA", "CEBPD", "ATF4", "CEBPG(var.2)"
)

keep.cell.types <- c(
  "Granulocytes", "Monocytes", "Macrophages", "DC",
  "T-cells", "Thymocyte", "B-cells", "Plasma"
)

recompute.z <- TRUE

analysis.dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/blueprint/"
rnb.tag <- "RnBeads_230826"
dev.tag <- "mTFR_devs_230826"
sannot.file <- file.path(analysis.dir, rnb.tag, "reports", "data_import_data", "annotation.csv")

github.dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR_manuscript"
plot.dir <- file.path(github.dir, "figures", "blueprint", "curated_heatmaps_230826")
if (!dir.exists(plot.dir)) dir.create(plot.dir, recursive = TRUE)

group_remap <- c(
  "Bcell" = "B-cells",
  "DC" = "DC",
  "eryt" = "Erythrocytes",
  "gran" = "Granulocytes",
  "megK" = "Megakaryocytes",
  "Mf" = "Macrophages",
  "mono" = "Monocytes",
  "NK" = "NK-cells",
  "osteoclast" = "Osteoclast",
  "other" = "Other",
  "plasma" = "Plasma",
  "progenitor" = "Progenitors",
  "Tcell" = "T-cells",
  "thymocyte" = "Thymocyte"
)

cell_type_colors <- c(
  "B-cells" = "#C2377C",
  "DC" = "#8C6D3F",
  "Erythrocytes" = "#7E4B2A",
  "Granulocytes" = "#E8A33D",
  "Macrophages" = "#B5A38A",
  "Megakaryocytes" = "#6A3D9A",
  "Monocytes" = "#C2703D",
  "NK-cells" = "#2CA02C",
  "Osteoclast" = "#8C8C8C",
  "Other" = "#BCBD22",
  "Plasma" = "#7B1E3D",
  "Progenitors" = "#17BECF",
  "T-cells" = "#4FC3D9",
  "Thymocyte" = "#5B9BD5"
)

zscore.limit <- 2
heatmap_colors <- colorRamp2(
  c(-zscore.limit, -1, 0, 1, zscore.limit),
  c("#2B4B9B", "#4EC3E0", "#FFF200", "#F5A623", "#B21212")
)

#####################################################################
# Helper functions
#####################################################################

row_zscore <- function(mat) {
  centred <- mat - rowMeans(mat, na.rm = TRUE)
  sds <- apply(mat, 1, sd, na.rm = TRUE)
  out <- centred / sds
  out[!is.finite(out)] <- 0
  out
}

# Curated panel heatmap. Samples in columns, split by cell type and
# clustered within each block, motifs in rows and clustered across the
# whole matrix. The blueprint colours mark each block, and each block
# carries a colour strip and an angled cell type label along the bottom.
plot_curated_heatmap <- function(mat, cell_types, title, file) {
  if (nrow(mat) < 2) {
    log_warn(basename(file), ": only ", nrow(mat), " motif(s), skipping the heatmap")
    return(invisible(NULL))
  }
  block_order <- intersect(keep.cell.types, unique(cell_types))
  column_split_factor <- factor(cell_types, levels = block_order)

  top_ha <- HeatmapAnnotation(
    `Cell type` = cell_types,
    col = list(`Cell type` = cell_type_colors[block_order]),
    annotation_name_gp = gpar(fontsize = 8),
    show_legend = TRUE
  )
  bottom_ha <- HeatmapAnnotation(
    strip = anno_block(
      gp = gpar(fill = cell_type_colors[block_order], col = NA),
      height = unit(3, "mm")
    ),
    show_annotation_name = FALSE
  )

  ht <- Heatmap(
    mat,
    name = "methylTFR\nZ-scores",
    col = heatmap_colors,
    top_annotation = top_ha,
    bottom_annotation = bottom_ha,
    column_title = block_order,
    column_title_side = "bottom",
    column_title_rot = 45,
    column_title_gp = gpar(fontsize = 8),
    column_split = column_split_factor,
    cluster_column_slices = FALSE,
    column_gap = unit(0.8, "mm"),
    show_row_names = TRUE,
    show_column_names = FALSE,
    row_names_side = "right",
    row_names_gp = gpar(fontsize = 7),
    row_dend_side = "left",
    clustering_distance_rows = "euclidean",
    clustering_distance_columns = "euclidean",
    clustering_method_rows = "complete",
    clustering_method_columns = "complete",
    heatmap_legend_param = list(
      direction = "horizontal",
      at = seq(-zscore.limit, zscore.limit, by = 1),
      legend_width = unit(4, "cm"),
      title_position = "lefttop"
    )
  )
  pdf(file, width = 9, height = 11)
  draw(ht,
    column_title = title,
    column_title_gp = gpar(fontsize = 11),
    heatmap_legend_side = "bottom",
    annotation_legend_side = "bottom",
    merge_legend = TRUE
  )
  dev.off()
  log_info("Wrote ", file)
}

#####################################################################
# Sample annotation
#####################################################################

sannot <- read.csv(sannot.file, stringsAsFactors = FALSE)
sannot$bedFile <- as.character(sannot$bedFile)

#####################################################################
# Curated heatmap per motif set
#####################################################################

for (motifSet in motifSets) {
  dev.file <- file.path(analysis.dir, dev.tag, paste0(motifSet, "_deviations.RDS"))
  if (!file.exists(dev.file)) {
    log_warn(motifSet, ": ", dev.file, " not found, skipping this motif set")
    next
  }

  log_info("Loading deviations from ", dev.file)
  dev_obj <- readRDS(dev.file)

  matched <- sannot[match(colnames(dev_obj), sannot$bedFile), ]
  mapped <- unname(group_remap[matched$cellTypeGroup])
  keep <- which(
    !is.na(matched$DISEASE) &
      matched$DISEASE == "None" &
      matched$cellTypeGroup %in% names(group_remap) &
      mapped %in% keep.cell.types
  )
  if (length(keep) == 0) {
    log_warn(motifSet, ": no samples left after filtering, skipping")
    next
  }

  dev_obj <- dev_obj[, keep]
  cell_types <- unname(group_remap[matched$cellTypeGroup[keep]])
  absent_groups <- setdiff(keep.cell.types, unique(cell_types))
  if (length(absent_groups) > 0) {
    log_warn(motifSet, ": no sample for ", paste(absent_groups, collapse = ", "))
  }
  if (length(unique(cell_types)) < 2) {
    log_warn(motifSet, ": only one cell type group present, skipping")
    next
  }

  zscores <- if (recompute.z) row_zscore(deviations(dev_obj)) else deviationZScores(dev_obj)

  present.motifs <- intersect(curated.motifs, rownames(zscores))
  missing.motifs <- setdiff(curated.motifs, rownames(zscores))
  if (length(missing.motifs) > 0) {
    log_warn(motifSet, ": ", length(missing.motifs), " curated motif(s) not found: ",
      paste(missing.motifs, collapse = ", "))
  }
  if (length(present.motifs) < 2) {
    log_warn(motifSet, ": fewer than 2 curated motifs present, no heatmap")
    next
  }

  mat <- zscores[present.motifs, , drop = FALSE]
  finite_rows <- is.finite(rowSums(mat))
  if (any(!finite_rows)) {
    log_warn(motifSet, ": dropping ", sum(!finite_rows), " motifs with non finite Z-scores")
    mat <- mat[finite_rows, , drop = FALSE]
  }

  plot_curated_heatmap(
    mat, cell_types,
    paste0("Curated TF panel (", motifSet, ")"),
    file.path(plot.dir, paste0("heatmap_curated_motifs_", motifSet, ".pdf"))
  )
}

log_success("Finished the curated motif heatmaps")

#!/usr/bin/env Rscript

#####################################################################
# mfoot_echo_cell_type_jaspar.R
# created on 25-10-25 by Irem Gunduz
# Plot footprints for the JASPAR motifs for mapped cell-types
#####################################################################

# load libraries
suppressPackageStartupMessages({
  library(ArchR)
  library(dplyr)
  library(ggplot2)
  library(muLogR)
  library(GenomicRanges)
  library(BSgenome.Hsapiens.UCSC.hg38)
})

cell_type_colors <- c(
  "B-cell" = "#AE017E",
  "Monocyte" = "#CC4C02",
  "NK-cell" = "#A65628",
  "Th-Mem" = "#41B6C4",
  "Tc-Mem" = "#4292C6",
  "Tc-Naive" = "#888FB5",
  "Th-Naive" = "#C7E9B4"
)
set.seed(12) # set seed

# Figures go to the analysis directory rather than the Plots folder of the
# ArchR project, which is shared and not the place for manuscript panels
plot.dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/echo/mFoot_atac_230826"
if (!dir.exists(plot.dir)) dir.create(plot.dir, recursive = TRUE)

# chromVAR activity shown in the legend, one column per pseudobulk sample
# named <cellType>_<sample>. This is the adjusted Z-score matrix, so the
# label says so rather than calling it a deviation.
cvar.file <- "/icbb/projects/nitschre/methylTFR/scripts/other/cvar_zscores_adj_psuedobulk.R"
dev.label <- "Mean chromVAR Z"

outputDir <- "/icbb/projects/igunduz/archr_projects/icbb/projects/igunduz/archr_project_011023"
project <- ArchR::loadArchRProject(outputDir, showLogo = FALSE)
addArchRThreads(threads = 30) # set the cores

# Get the overlapped samples between scATAC-seq and scWGBS
annot_acc <- readRDS("/icbb/projects/igunduz/DARPA_analysis/artemis_031023/cell_matching/annot_acc_cellType.rds")
cell_annot_atac <- read.delim("/icbb/projects/igunduz/DARPA_analysis/artemis_031023/sample_annot_atac.tsv")
rownames(cell_annot_atac) <- cell_annot_atac$cellId_archr
annot_acc <- readRDS("/icbb/projects/igunduz/DARPA_analysis/artemis_031023/cell_matching/annot_acc_cellType.rds")
annot_acc$archr_id <- rownames(annot_acc)
new_annot <- merge(cell_annot_atac, annot_acc, by = "row.names")
new_annot$atac_pb_id <- paste0(new_annot$class, "_", new_annot$sample_sampleId_cminid)

# Subset the cells
idxSample <- BiocGenerics::which(project$cellNames %in% new_annot$archr_id)
cellsSample <- project$cellNames[idxSample]
project <- project[cellsSample, ]


# Add mapped cell-type information
project <- addCellColData(
  ArchRProj = project, data = new_annot$atac_pb_id,
  name = "atac_pb_id", cells = new_annot$archr_id
)
project <- addCellColData(
  ArchRProj = project, data = new_annot$mapped_cell_class,
  name = "cell_type_meth", cells = new_annot$archr_id, force = TRUE
)

# Remove the Other-cell
idxSample <- BiocGenerics::which(project$cell_type_meth %in% "Other-cell")
cellsSample <- project$cellNames[-idxSample]
project <- project[cellsSample, ]

# get motif positions
motifPositions <- getPositions(project, name = "jaspar2020")

# Motifs from MOFA2 analysis
markerMotifs <- c(
  "OLIG2", "NEUROG1", "HES1", "HES6", "HEY2", "HEY1", "KLF17", "KLF16",
  "TFAP2A", "ZNF652", "TCF21(var.2)", "TFAP4", "MSC", "TBXT", "OSR2",
  "ETV5", "ELF4", "ELF5", "SPIC", "SPIB", "SPI1", "EHF", "PRDM4", "GFI1",
  "MXI1", "ZBTB14", "E2F2", "FOXN3", "ZFP57", "ESR1", "TCF4", "ASCL1(var.2)",
  "EBF3", "POU2F3", "POU2F2", "POU5F1", "POU3F4", "LHX6", "EMX2", "TBX21",
  "TBX2", "TBR1", "EOMES", "TBX20", "TBX1", "MGA", "HIC2"
)

# add group coverages
project <- addGroupCoverages(ArchRProj = project, groupBy = "cell_type_meth")

# get footprints
seFoot <- getFootprints(
  ArchRProj = project,
  positions = motifPositions[markerMotifs],
  groupBy = "cell_type_meth"
)

#####################################################################
# chromVAR activity per cell type, for the legend
#####################################################################

# NULL when the matrix is unreachable, in which case the legend falls back
# to the plain cell type names
load_cell_type_activity <- function(file, groups) {
  if (!file.exists(file)) {
    message("No chromVAR matrix at ", file, ", the legend will omit the scores")
    return(NULL)
  }
  mat <- as.matrix(readRDS(gzfile(file)))
  cell_type <- sub("_.*", "", colnames(mat))

  out <- list()
  for (g in groups) {
    idx <- which(cell_type == g)
    if (length(idx) == 0) {
      message(g, ": not in the chromVAR matrix, the legend will show N/A")
      next
    }
    out[[g]] <- mat[, idx, drop = FALSE]
  }
  if (length(out) == 0) NULL else out
}

# Mean activity of a motif within each cell type. ArchR appends an index
# to motif names, so a trailing _<n> is stripped when the exact name is
# not in the matrix. Matching is fixed rather than regular, because motif
# names such as TCF21(var.2) carry punctuation.
motif_mean_activity <- function(motif, dev_list, groups) {
  vapply(groups, function(g) {
    mat <- dev_list[[g]]
    if (is.null(mat)) {
      return(NA_real_)
    }
    idx <- which(rownames(mat) == motif)
    if (length(idx) == 0) {
      idx <- which(rownames(mat) == sub("_[0-9]+$", "", motif))
    }
    if (length(idx) == 0) {
      idx <- grep(motif, rownames(mat), fixed = TRUE)
    }
    if (length(idx) == 0) {
      return(NA_real_)
    }
    vals <- mat[idx[1], ]
    if (all(is.na(vals))) NA_real_ else mean(vals, na.rm = TRUE)
  }, numeric(1))
}

make_label <- function(group, dev_val) {
  if (is.na(dev_val)) {
    return(group)
  }
  paste0(group, ": ", format(round(dev_val, 2), nsmall = 2))
}

activity <- load_cell_type_activity(cvar.file, names(cell_type_colors))

#####################################################################
# Plot, one PDF per motif, outside the ArchR project
#####################################################################

# plot = FALSE returns the panels instead of writing them into the Plots
# folder of the project, so they can be relabelled and saved elsewhere
foot_plots <- plotFootprints(
  seFoot = seFoot,
  ArchRProj = project,
  # Subtraction, to match the methylTFR footprints of 03_mfoot_meth.R.
  # ArchR divides by the Tn5 bias track by default, which is not on the
  # same scale as the subtracted methylTFR curves.
  normMethod = "subtract",
  pal = cell_type_colors,
  plotName = "Integrative_Footprints_JASPAR_by_cellType",
  addDOC = FALSE,
  plot = FALSE,
  force = TRUE
)

if (is.null(names(foot_plots))) names(foot_plots) <- names(motifPositions[markerMotifs])

# The scores differ per motif, so the legend is relabelled per panel.
# Only a ggplot can carry a new scale; anything else is saved untouched.
add_activity_legend <- function(p, motif) {
  if (is.null(activity) || !inherits(p, "ggplot")) {
    return(p)
  }
  means <- motif_mean_activity(motif, activity, names(cell_type_colors))
  labels <- vapply(
    names(cell_type_colors), function(g) make_label(g, means[[g]]), character(1)
  )
  suppressMessages(
    p +
      scale_color_manual(values = cell_type_colors, labels = labels, name = NULL) +
      scale_fill_manual(values = cell_type_colors, labels = labels, name = NULL)
  )
}

for (motif in names(foot_plots)) {
  p <- add_activity_legend(foot_plots[[motif]], motif)
  out.file <- file.path(plot.dir, paste0("ATAC_footprint_", make.names(motif), ".pdf"))

  pdf(out.file, width = 6, height = 8)
  if (inherits(p, "ggplot")) print(p) else grid::grid.draw(p)
  dev.off()
}
message("Wrote ", length(foot_plots), " footprint panels to ", plot.dir)

# The scores behind the legends, so a panel can be checked against a number
activity_tab <- do.call(rbind, lapply(names(foot_plots), function(motif) {
  means <- if (is.null(activity)) {
    setNames(rep(NA_real_, length(cell_type_colors)), names(cell_type_colors))
  } else {
    motif_mean_activity(motif, activity, names(cell_type_colors))
  }
  data.frame(
    motif = motif, cell_type = names(means), mean_chromvar_z = unname(means),
    stringsAsFactors = FALSE
  )
}))
write.csv(activity_tab, file.path(plot.dir, "atac_footprint_chromvar_means.csv"), row.names = FALSE)

#####################################################################

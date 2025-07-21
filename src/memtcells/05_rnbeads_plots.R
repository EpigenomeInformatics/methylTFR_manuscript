#!/usr/bin/env Rscript

#####################################################################
# 05_rnbeads_plots.R
# Created by 03-06-2025 by IBG
# Plot MA and lola plots from rnbeads
#####################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(RnBeads)
  library(grid)
})
set.seed(12)


source("/icbb/projects/igunduz/irem_github/methylTFR_manuscript/utils/utils_rnbeads.R")
analysis.dir <- "/icbb/projects/igunduz/methylTFR_manuscript/results/memoryTcells_060525/"
plot_path <- "/icbb/projects/igunduz/irem_github/methylTFR_manuscript/Figures/distal_nobatch_090725/"
if (!dir.exists(plot_path)) dir.create(plot_path, recursive = TRUE)

# Load the differential methylation data
diffMeth <- load.rnb.diffmeth(paste0(analysis.dir, "/reports_080725/differential_methylation_data/differential_rnbDiffMeth/"))
plots <- rnbeadsDensityScatter(diffMeth, "distal")

# Save all plots to files
for (name in names(plots)) {
    ggsave(filename = paste0(plot_path,name, "_scatterplot.pdf"), plot = plots[[name]])
}

plot_path <- "/icbb/projects/igunduz/irem_github/methylTFR_manuscript/Figures/distal_nobatch_090725_rankcut/"
if (!dir.exists(plot_path)) dir.create(plot_path, recursive = TRUE)
rankcut_plots <- rnbeadsDensityScatterRankCut(diffMeth, "distal")
# Save all plots to files
for (name in names(rankcut_plots)) {
    ggsave(filename = paste0(plot_path,name, "_scatterplot.pdf"), plot = rankcut_plots[[name]])
}

plot_path <- "/icbb/projects/igunduz/irem_github/methylTFR_manuscript/Figures/tiling1kb_nobatch_090725_rankcut/"
if (!dir.exists(plot_path)) dir.create(plot_path, recursive = TRUE)
rankcut_plots <- rnbeadsDensityScatterRankCut(diffMeth, "tiling1kb")
# Save all plots to files
for (name in names(rankcut_plots)) {
    ggsave(filename = paste0(plot_path,name, "_scatterplot.pdf"), plot = rankcut_plots[[name]])
}

#####################################################################
# LOLA Analysis
#######################################################################

logger.start("Running LOLA")
lolaDb_path <- "/icbb/projects/share/annotation/lolaDB/hg38/"

logger.info("Loading RnBeads objects")
diffMeth <- load.rnb.diffmeth(paste0(analysis.dir, "/reports_030625/differential_methylation_data/differential_rnbDiffMeth/"))
rnb_set <- load.rnb.set(paste0(analysis.dir, "/reports_030625/data_import_data/rnb.set_preprocessed/"))

# Run LOLA
res <- performLolaEnrichment.diffMeth(rnb_set, diffMeth, lolaDb_path)
logger.info("Saving results")
saveRDS(res, paste0(analysis.dir, "/reports_030625/differential_methylation_data/differential_rnbDiffMeth/lola_results.rds"))
logger.completed()

#######################################################################
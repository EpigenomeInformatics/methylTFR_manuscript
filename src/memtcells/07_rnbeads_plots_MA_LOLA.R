#!/usr/bin/env Rscript

#####################################################################
# 05_rnbeads_plots.R
# Created by 03-06-2025 by IBG
# Plot MA plots from rnbeads and Lola plots
#####################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(RnBeads)
  library(grid)
  library(LOLA)
})
set.seed(12)
source("/icbb/projects/nitschre/methylTFR/methylTFR_manuscript/src/utils.R")

plot_path <- "/scratch/icbb/regina/methylTFR/figures/memoryTcells/"
analysis.dir <- "/scratch/icbb/regina/data/results/memoryTcells/RnBeads/"

# Load the differential methylation data
diffMeth <- load.rnb.diffmeth(paste0(analysis.dir, "reports/differential_methylation_data/differential_rnbDiffMeth/"))

# Density scattere with rank colored (distal)
color_map <- c(
  TCM = "#84B067",
  TN = "#C8E0B4",
  TEM = "#0fb26bff",
  TEMRA = "#898FB5"
)
rankcut_plots <- rnbeadsDensityScatterRankCut(diffMeth, "distal", color_mapping = color_map)

# Save all plots to files
for (name in names(rankcut_plots)) {
    ggsave(filename = paste0(plot_path,name, "distal_rank_scatterplot.pdf"), plot = rankcut_plots[[name]])
}

#####################################################################
# LOLA Analysis
#######################################################################

logger.start("Running LOLA")
lolaDb_path <- "/icbb/projects/share/annotations/lolaDB/hg38/"

# Create tiling regions for 1kb
tiling1kb <- muRtools::getTilingRegions("hg38", width=1000L, onlyMainChrs=TRUE)%>%
  data.table::as.data.table() %>%
  dplyr::select(seqnames, start, end) %>%
  as.data.frame()
colnames(tiling1kb) <- c("Chromosome", "Start", "End")
rnb.set.annotation(type = "tiling1kb", regions = tiling1kb, assembly = "hg38")

distal <- readRDS("/icbb/projects/share/annotations/methylTFRAnnotationHg38/inst/extdata/distal_regions.RDS")
distal <- as.data.frame(distal) %>%
  dplyr::select(seqnames, start, end)
colnames(distal) <- c("Chromosome", "Start", "End")
rnb.set.annotation(type = "distal", regions = distal, assembly = "hg38")

logger.info("Loading RnBeads objects")
rnb_set <- load.rnb.set(paste0(analysis.dir, "/reports/data_import_data/rnb.set_preprocessed/"))

# Run LOLA
res <- performLolaEnrichment.diffMeth(rnb_set, diffMeth, lolaDb_path)
logger.info("Saving results")
saveRDS(res, paste0(analysis.dir, "differential_methylation_data/differential_rnbDiffMeth/lola_results.rds"))
logger.completed()

## Plot LOLA results
comparisons <- names(res$region)[1:3]

# Volcano Plot
source("/scratch/icbb/regina/methylTFR/src/lola_utils.R")
lola_obj <- paste0(analysis.dir, "differential_methylation_data/differential_rnbDiffMeth/lola_results.rds")

# Load the LOLA database
lolaDb <- loadRegionDB(lolaDb_path)

# Plot
for(comparison in comparisons){
p <- lolaVolcanoPlotC19(cell = NULL,
        lolaDb = lolaDb,
        outputDir = outputDir,
        reg = comparison,
        database = "TF_motifs"
)
sample <- sub(" \\(.*", "", comparison)
ggsave(paste0(plot_path, "LolaVolcanoPlot_distal", "_", sample, ".pdf"), plot=p$plot)
        
}


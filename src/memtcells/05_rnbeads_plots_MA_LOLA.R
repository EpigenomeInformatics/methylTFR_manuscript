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
})
set.seed(12)
source("/icbb/projects/nitschre/methylTFR/scripts/other/muScatter.R")
source("/icbb/projects/nitschre/methylTFR/scripts/blueprint/utils_rnbeads.R")

plot_path <- "/icbb/projects/nitschre/methylTFR/figures/memoryTcells/MAplots/rankcut/"
analysis.dir <- "/icbb/projects/nitschre/methylTFR/results/memoryTcells/"

if (!dir.exists(plot_path)) dir.create(plot_path, recursive = TRUE)

# Load the differential methylation data
diffMeth <- load.rnb.diffmeth(paste0(analysis.dir, "/reports_lola/differential_methylation_data/differential_rnbDiffMeth/"))

# Density scatter without rank colored
plots <- rnbeadsDensityScatter(diffMeth, "distal")

# Save all plots to files
for (name in names(plots)) {
    ggsave(filename = paste0(plot_path,name, "_scatterplot.pdf"), plot = plots[[name]])
}

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
    ggsave(filename = paste0(plot_path,name, "distal_scatterplot.pdf"), plot = rankcut_plots[[name]])
}

# Density scattere with rank colored (tiling1kb)
rankcut_plots <- rnbeadsDensityScatterRankCut(diffMeth, "tiling1kb", color_mapping = color_map)

# Save all plots to files
for (name in names(rankcut_plots)) {
    ggsave(filename = paste0(plot_path,name, "tiling1kb_scatterplot.pdf"), plot = rankcut_plots[[name]])
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
saveRDS(res, "/icbb/projects/nitschre/methylTFR/results/memoryTcells/reports_lola/differential_methylation_data/differential_rnbDiffMeth/lola_results_tiling.rds")
logger.completed()

# Plot LOLA results
comparisons <- names(res$region)[1:3]
for(comparison in comparisons){
lolaRes <- res$region[[comparison]][["tiling1kb"]]
lolaRes <- lolaRes[lolaRes$collection == "TF_motifs",]

bp <- lolaBarPlot.hyp(res$lolaDb, lolaRes, scoreCol="oddsRatio", orderCol="maxRnk", pvalCut=0.05,groupByCollection=FALSE)
ggsave(file.path(plot_path, paste0(comparison, "_1kb_lola_hyper.pdf")), bp, width = 15, height = 15, units = "cm")
}

# Remove MA before name 
for(comparison in comparisons){
res$region[[comparison]][["distal"]]$description <- 
    sub(".*_", "", res$region[[comparison]][["distal"]]$description)
}
saveRDS(res, "/icbb/projects/nitschre/methylTFR/results/memoryTcells/reports_lola/differential_methylation_data/differential_rnbDiffMeth/lola_results_tiling.rds")

# Lola Bar plot
userSets <- c("rankCut_1000_hyper", "rankCut_1000_hypo")
for(comparison in comparisons){
  for (set in userSets){
    lolaRes <- res$region[[comparison]][["tiling1kb"]]
    lolaRes_filtered <- lolaRes[lolaRes$userSet == set & lolaRes$collection == "TF_motifs",]

    # Plot
    sample <- sub(" \\(.*", "", comparison)
    suffix <- if (grepl("hyper", set)) "hyper" else "hypo"

    fn <- file.path(plot_dir, paste0("LOLA1kb", "_", sample, "_", suffix, ".pdf"))
    p <- lolaBarPlot(
      res$lolaDb,
      lolaRes_filtered,
      scoreCol = "oddsRatio",    
      orderCol = "oddsRatio",         
      pvalCut = 0.05,
      maxTerms = 50
    )
    ggsave(fn, p)
  }
}
# Volcano Plot
source("/icbb/projects/nitschre/methylTFR/scripts/blueprint/lola.R")
outputDir <- file.path("/icbb/projects/nitschre/methylTFR/results/memoryTcells/reports_lola/differential_methylation_data/differential_rnbDiffMeth/lola_results_tiling.rds")

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
ggsave(paste0("/icbb/projects/nitschre/methylTFR/figures/memoryTcells/LolaVolcanoPlot_tiling", "_", sample, ".pdf"), plot=p$plot)
        
}


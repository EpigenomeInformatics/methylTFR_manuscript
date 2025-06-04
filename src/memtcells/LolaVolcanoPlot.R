library(dplyr)
library(LOLA)
source("/icbb/projects/nitschre/methylTFR/scripts/memoryTcells/lola.R")

outputDir <- "/icbb/projects/nitschre/methylTFR/results/memoryTcells/reports/differential_methylation_data/differential_rnbDiffMeth/lola_results.rds"

# Load the LOLA database
lolaDb_path <- "/icbb/projects/share/annotations/lolaDB/hg38/"
lolaDb <- loadRegionDB(lolaDb_path)

# Read Lola results
res <- readRDS(outputDir)

# Plot
p <- lolaVolcanoPlotC19(cell = NULL,
        lolaDb = lolaDb,
        outputDir = outputDir,
        reg = "TCM vs. non.TCM (based on cellType)",
        database = "TF_motif_clusters"
)

ggsave("/icbb/projects/nitschre/methylTFR/figures/LolaVolcanoPlot.pdf", plot=p$plot)

p <- lolaVolcanoPlotC19(cell = NULL,
        lolaDb = lolaDb,
        outputDir = outputDir,
        reg = "TEM vs. non.TEM (based on cellType)",
        database = "TF_motif_clusters"
)

ggsave("/icbb/projects/nitschre/methylTFR/figures/LolaVolcanoPlotTEMvsNonTEM.pdf", plot=p$plot)
library(dplyr)
library(LOLA)
source("/icbb/projects/nitschre/methylTFR/scripts/memoryTcells/lola.R")

outputDir <- "/icbb/projects/nitschre/methylTFR/results/memoryTcells/reports_lola/differential_methylation_data/differential_rnbDiffMeth/lola_results.rds"

# Load the LOLA database
lolaDb_path <- "/icbb/projects/share/annotations/lolaDB/hg38/"
lolaDb <- loadRegionDB(lolaDb_path)

# Read Lola results
res <- readRDS(outputDir)

# Plot
comparisons <- names(res$region)[1:3]

for(comparison in comparisons){
p <- lolaVolcanoPlotC19(cell = NULL,
        lolaDb = lolaDb,
        outputDir = outputDir,
        reg = comparison,
        database = "TF_motif_clusters"
)
sample <- sub(" \\(.*", "", comparison)
ggsave(paste0("/icbb/projects/nitschre/methylTFR/figures/memoryTcells/LolaVolcanoPlot", "_", sample), plot=p$plot)
        
}

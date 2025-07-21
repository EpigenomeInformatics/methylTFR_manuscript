#####################################################################
# LOLA Analysis
#######################################################################
suppressPackageStartupMessages({
  library(dplyr)
  library(RnBeads)
  library(grid)
})
set.seed(12)

# Directory where the output should be written to
analysis.dir <- "/icbb/projects/igunduz/methylTFR_manuscript/results"
if (!dir.exists(analysis.dir)) dir.create(analysis.dir)
analysis.dir <- file.path(analysis.dir, "memoryTcells_060525")
if (!dir.exists(analysis.dir)) dir.create(analysis.dir)

logger.start("Running LOLA")
lolaDb_path <- "/icbb/projects/share/annotations/lolaDB/hg38/"

logger.info("Loading RnBeads objects")
diffMeth <- load.rnb.diffmeth(paste0(analysis.dir, "/reports_080725/differential_methylation_data/differential_rnbDiffMeth/"))
rnb_set <- load.rnb.set(paste0(analysis.dir, "/reports/data_import_data/rnb.set_preprocessed/"))

# Run LOLA
res <- performLolaEnrichment.diffMeth(rnb_set, diffMeth, lolaDb_path)
logger.info("Saving results")
saveRDS(res, paste0(analysis.dir, "/reports_080725/differential_methylation_data/differential_rnbDiffMeth/lola_results.rds"))
logger.completed()

#######################################################################
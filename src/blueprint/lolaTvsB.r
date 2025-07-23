

#####################################################################
# LOLA Analysis
#######################################################################

logger.start("Running LOLA for TvsBcells")
analysis.dir <- "/icbb/projects/nitschre/methylTFR/results/blueprint_subset_tiling/"
#cell <- cells[2] # Monocyte

# Load the LOLA database
lolaDb_path <- "/icbb/projects/share/annotations/lolaDB/hg38/"

logger.info("Loading RnBeads objects")
diffMeth <- load.rnb.diffmeth(paste0(analysis.dir, "/reports_lola/differential_methylation_data/differential_rnbDiffMeth/"))
rnb_set <- load.rnb.set(paste0(analysis.dir, "/reports/data_import_data/rnb.set_preprocessed/"))

# Run LOLA
res <- performLolaEnrichment.diffMeth(rnb_set, diffMeth, lolaDb_path)
logger.info("Saving results")
saveRDS(res, paste0(analysis.dir, "/reports/differential_methylation_data/differential_rnbDiffMeth/lola_results.rds"))
logger.completed()




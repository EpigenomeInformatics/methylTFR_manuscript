

#####################################################################
# LOLA Analysis
#######################################################################

logger.start("Running LOLA for memorTcells")
analysis.dir <- "/icbb/projects/nitschre/methylTFR/results/memoryTcells/"
#cell <- cells[2] # Monocyte

# Load the LOLA database
lolaDb_path <- "/icbb/projects/share/annotations/lolaDB/hg38/"

logger.info("Loading RnBeads objects")
diffMeth <- load.rnb.diffmeth(paste0(analysis.dir, "/reports/differential_methylation_data/differential_rnbDiffMeth/"))
rnb_set <- load.rnb.set(paste0(analysis.dir, "/reports/data_import_data/rnb.set_preprocessed/"))

# Run LOLA
res <- performLolaEnrichment.diffMeth(rnb_set, diffMeth, lolaDb_path)
logger.info("Saving results")
saveRDS(res, paste0(analysis.dir, "/reports/differential_methylation_data/differential_rnbDiffMeth/lola_results.rds"))
logger.completed()

# select the 50 most hypermethylated tiling regions in TCM compared to TN
lolaRes <- res$region[["TCM vs. non.TCM (based on cellType)"]][["tiling"]]
lolaRes_filtered <- lolaRes[lolaRes$userSet == "rankCut_100_hyper" & lolaRes$collection == "TF_motifs", ]

# plot
plot_dir <- "/icbb/projects/nitschre/methylTFR/figures"
fn <- file.path(plot_dir, "LOLAEnrichment_TCM.pdf")
pdf(fn)
lolaBarPlot(
  res$lolaDb,
  lolaRes_filtered,
  scoreCol = "oddsRatio",    
  orderCol = "oddsRatio",         
  pvalCut = 0.05,                 
  maxTerms = 50
)
dev.off()


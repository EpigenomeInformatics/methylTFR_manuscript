

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

# select the 50 most hyper and hypomethylated promotor regions in TCM compared to TN
comparisons <- names(res$region)[1:3]
userSets <- c("rankCut_1000_hyper", "rankCut_1000_hypo")

plot_dir <- "/icbb/projects/nitschre/methylTFR/figures"

for(comparison in comparisons){
  for (set in userSets){
    lolaRes <- res$region[[comparison]][["promoters"]]
    lolaRes_filtered <- lolaRes[lolaRes$userSet == set & lolaRes$collection == "TF_motif_clusters",]

    # Plot
    sample <- sub(" vs.*", "", comparison)
    suffix <- if (grepl("hyper", set)) "hyper" else "hypo"

    fn <- file.path(plot_dir, paste0("LOLA", "_", sample, "_", suffix, ".pdf"))
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



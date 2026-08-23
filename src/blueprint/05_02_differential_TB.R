#!/usr/bin/env Rscript
#####################################################################
# differentials_heatmap.R
# Created by RN on 13-11-2025
# Script to calculate differntials between T and B cells from Blueprint data
#####################################################################

suppressPackageStartupMessages({
  library(methylTFR)
})

set.seed(13)

# Set up
results_dir <- "/scratch/icbb/regina/data/blueprint/data"

# Loading deviations scores and remove problematic sample
dev_obj <- readRDS("/icbb/projects/nitschre/methylTFR/r_objects/jaspar2020_distal_deviations_TvsB.RDS")
dev_obj <- dev_obj[,!colnames(dev_obj) %in% "Bcell_naive_VB_NBC_NC11_83.bed"]
deviations <- deviations(dev_obj)

# Aggregating cell names 
colnames(deviations) <- case_when(
  grepl("Bcell", colnames(deviations), ignore.case = TRUE) ~ "Bcell",
  grepl("TC", colnames(deviations), ignore.case = TRUE) ~ "Tcell"
)

# Set groups to compare
groups <- colnames(deviations)

# Perform differential analysis
if(file.exists(paste0(results_dir, "diff_BvsT.RDS"))) {  
  diff <- readRDS(paste0(results_dir, "diff_BvsT.RDS"))
} else {
  diff <- differential_deviation_test(deviations, groups=groups, 
                                     alternative="two.sided", parametric = TRUE)
  saveRDS(diff, file=paste0(results_dir, "diff_BvsT.RDS"))
}

set.seed(42)
logger::log_info("Loading libraries...")

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(methylTFR)
  library(muLogR)
  library(gplots)
  library(muRtools)
  library(ComplexHeatmap)
  library(factoextra)
  library(ggfortify)
})

plot_dir <- "/icbb/projects/nitschre/methylTFR/figures"
if(!dir.exists(plot_dir)) dir.create(plot_dir)
deviations <- readRDS("/icbb/projects/nitschre/methylTFR/r_objects/blueprintjaspar2020_distal_deviations_subset.RDS")
motifset <- "jaspar2020_distal"
deviations <- deviations(deviations)
sannot <- fread("/icbb/projects/skumar/blueprint/bed/samples.tsv")
sannot <- fread("/icbb/projects/nitschre/methylTFR/samples_subset.tsv")

# Get cell types
tdf <- as.data.frame(t(deviations))
pca <- prcomp(tdf, scale. = T)
tdf$cell_type <- sannot$CELL_TYPE

tdf$cell_type <- case_when(
  grepl("naive B", tdf$cell_type, ignore.case = TRUE) ~ "naive B",
  grepl("precursor B cell", tdf$cell_type, ignore.case = TRUE) ~ "precursor B",
  grepl("precursor lymphocyte of B lineage", tdf$cell_type, ignore.case = TRUE) ~ "precursor B",  
  grepl("memory B cell", tdf$cell_type, ignore.case = TRUE) ~ "memory B",
  grepl("germinal center B cell", tdf$cell_type, ignore.case = TRUE) ~ "memory B",  
  grepl("CD8", tdf$cell_type, ignore.case = TRUE) ~ "CD8 T",
  grepl("CD4", tdf$cell_type, ignore.case = TRUE) ~ "CD4 T",
)


# File path
fn_pca_ind <- file.path(plot_dir, paste0("pca_blueprintSubset_", motifset, ".pdf"))

# Save as plot
pdf(fn_pca_ind)
autoplot(pca, data = tdf,
         colour = 'cell_type', 
         size = 5,
         main = "PCA - Bcell vs Tcell") +
         theme_classic() +
         theme(legend.position = "bottom",
         legend.text = element_text(size=10))
dev.off()

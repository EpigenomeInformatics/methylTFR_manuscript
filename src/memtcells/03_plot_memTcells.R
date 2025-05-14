setwd("/icbb/projects/igunduz/methylTFR_manuscript/")
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

plot_dir <- "/icbb/projects/igunduz/methylTFR_manuscript/Figures"
if (!dir.exists(plot_dir)) dir.create(plot_dir)
deviations <- readRDS("//icbb/projects/igunduz/methylTFR_manuscript/results/memoryTcells/mtfr_final_221223/jaspar2020_distal_deviations.RDS")
motifset <- "jaspar2020_distal"
deviations <- deviations(deviations)
sannot <- fread("/icbb/projects/skumar/memoryTcells/bed/samples.tsv")

# Get cell types
tdf <- as.data.frame(t(deviations))
pca <- prcomp(tdf, scale. = T)
tdf$cell_type <- sannot$cellType


# Save PCA individual plot
fn_pca_ind <- file.path(plot_dir, paste0("pca_individuals_", motifset, ".pdf"))
pdf(fn_pca_ind)
fviz_pca_ind(res.pca, repel = TRUE, title = "PCA on all memTcells samples")
dev.off()


# Re-group T-cells
# tdf$cell_type  <- ifelse(tdf$cell_type %in% c("TCD4", "TCD8"),"Tcell", tdf$cell_type)


# skip cell_type column
pca <- prcomp(tdf[, -633])
fn_pca_ind <- file.path(plot_dir, paste0("bcell_vs_tcell_memTcells_", motifset, ".pdf"))
pdf(fn_pca_ind)
autoplot(pca,
  data = tdf,
  colour = "cell_type",
  size = 5,
  main = "PCA - Bcell vs Tcell"
) +
  theme_classic() +
  theme(legend.position = "bottom")
dev.off()


# merging TCD4 and TCD8 into single group
groups <- sannot$cellType
match <- which(groups %in% c("TCD4", "TCD8"))
groups[match] <- "Tcell"
groups <- as.factor(groups)

tdf <- as.data.frame(deviations)
diff <- differential_deviation_test(tdf, groups = groups, alternative = "two.sided", parametric = FALSE)

dim(diff[diff$p_value_adjusted < 0.05, ])
head(diff[diff$p_value_adjusted < 0.05, ])

write.table(diff, file = paste0(plot_dir, "/", motifset, "_diff_devs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

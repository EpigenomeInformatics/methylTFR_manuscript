#!/usr/bin/env Rscript

#####################################################################
# 04_plot_blueprint.R
# created on 2024-01-10 by Irem Gunduz
# Plot PCA and Heatmaps for Blueprint MethylTFR results
#####################################################################

set.seed(42)
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
  library(RnBeads)
})
source("/icbb/projects/igunduz/methylTFR_manuscript/src/utils.R")

plot_dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR_manuscript/figures/blueprint"
if (!dir.exists(plot_dir)) dir.create(plot_dir)

# Full blueprint methylTFR results
full_mtfr <- readRDS("/scratch/icbb/mtfr_manuscript/BLUEPRINT/mtfr_final_251223/jaspar2020_distal_deviations.RDS")
bdevs <- deviationZScores(full_mtfr)
motifset <- "jaspar2020_distal"
sannot <- read.delim("/scratch/icbb/mtfr_manuscript//blueprint_data/samples_subset.tsv")[1:5]

tdf <- as.data.frame(t(bdevs))
rownames(sannot) <- sannot$bedFile
tdf$cell_type <- sannot[colnames(bdevs), "cellTypeGroup"]
cell_types <- unique(tdf$cell_type)
colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", "#1b9e77", "#d95f02", "#7570b3")
named_colors <- setNames(colors[1:length(cell_types)], cell_types)

# skip cell_type column
k <- ifelse(motifset == "altius", 287, 633)
pca <- prcomp(tdf[, -k])
fn_pca_ind <- file.path(plot_dir, paste0("pca_blueprint_", motifset, ".pdf"))
pdf(fn_pca_ind)
autoplot(pca,
  data = tdf,
  colour = "cell_type",
  size = 5,
  main = "PCA - Blueprint"
) +
  theme_classic() +
  theme(legend.position = "right") +
  scale_color_manual(values = colors)
dev.off()

# Plot pie chart for cell types
# Calculate the percentage of each cell type
percentages <- round(sample_counts / sum(sample_counts) * 100, 1)

# Create labels
labels <- paste(names(sample_counts), "\n", percentages, "%", sep = "")

# Create the pie chart
pdf(file.path(plot_dir, "blueprint_pie.pdf"), width = 15, height = 15, onefile = FALSE)
pie(sample_counts, main = "Number of samples per cell type", col = colors, labels = labels)
dev.off()

# RnBeads plot for BLUEPRINT
analysis.dir <- "/icbb/projects/igunduz/methylTFR_manuscript/results/BLUEPRINT/reports/differential_methylation_data/differential_rnbDiffMeth/"
diffMeth <- load.rnb.diffmeth(analysis.dir)
region <- "cpgislands"
pp <- rnbeadsDensityScatter(diffMeth, "cpgislands")
ggsave(file.path(plot_dir, "blueprint_rnb_diffmeth.pdf"), pp, width = 15, height = 15, units = "cm")

# logger.info("Loading LOLA database")
lolaDb_path <- "/icbb/projects/share/annotations/lolaDB/hg38/"
rnb_set <- load.rnb.set("/icbb/projects/igunduz/methylTFR_manuscript/results/BLUEPRINT/reports/data_import_data/rnb.set_preprocessed")

# Run LOLA
res <- performLolaEnrichment.diffMeth(rnb_set, diffMeth, lolaDb_path)
logger.info("Saving results")
saveRDS(res, paste0(analysis.dir, "lola_results.rds"))
# res <- readRDS(paste0(analysis.dir, "lola_results.rds"))

# Plot LOLA results
lolaRes <- res$region[["Bcell vs. Tcell (based on cellTypeGroup)"]][["cpgislands"]]
lolaRes <- lolaRes[lolaRes$collection == "TF_motif_clusters", ]

bp <- lolaBarPlot.hyp(res$lolaDb, lolaRes, scoreCol = "oddsRatio", orderCol = "maxRnk", pvalCut = 0.05, groupByCollection = FALSE)
ggsave(file.path(plot_dir, "blueprint_lola_hyper.pdf"), bp, width = 15, height = 15, units = "cm")


#####################################################################
mtfr <- readRDS("/icbb/projects/igunduz/methylTFR_manuscript/results/BLUEPRINT/mtfr_final_221223/altius_deviations.RDS")
motifset <- "altius"
deviations <- deviations(mtfr)

fn <- file.path(plot_dir, paste0("deviation_score_all_", motifset, ".pdf"))
pdf(fn, width = 15, height = 15, onefile = FALSE)
Heatmap(as.matrix(deviations),
  name = "deviation_score",
  column_title = "samples", row_title = "motifs",
  cluster_rows = FALSE, show_row_names = TRUE
)
dev.off()


samples <- colnames(deviations)
get_groupname <- function(x) {
  return(unlist(strsplit(x, split = "_"))[1])
}
groups <- unlist(lapply(FUN = get_groupname, X = samples))
tdf <- as.data.frame(t(deviations))
tdf$cell_type <- groups
res.pca <- prcomp(t(deviations))


# Save PCA individual plot
fn_pca_ind <- file.path(plot_dir, paste0("pca_individuals_", motifset, ".pdf"))
pdf(fn_pca_ind)
fviz_pca_ind(res.pca, repel = TRUE, title = "PCA samples")
dev.off()


# skip cell_type column
k <- ifelse(motifset == "altius", 287, 633)
pca <- prcomp(tdf[, -k])
fn_pca_ind <- file.path(plot_dir, paste0("bcell_vs_tcell_blueprint_", motifset, ".pdf"))
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
match <- which(groups %in% c("TCD4", "TCD8"))
groups[match] <- "Tcell"
groups <- as.factor(groups)

tdf <- as.data.frame(deviations)
diff <- differential_deviation_test(tdf, groups = groups, alternative = "two.sided", parametric = FALSE)

dim(diff[diff$p_value_adjusted < 0.05, ])
head(diff[diff$p_value_adjusted < 0.05, ])
write.table(diff, file = paste0(plot_dir, "/", motifset, "_diff_devs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Plot heatmap for differential deviations
deviations <- deviationZScores(mtfr)
deviations <- deviations[diff[diff$p_value_adjusted < 0.01, ]$motif, ]

# Add group information
ann <- data.frame(groups = groups)
rownames(ann) <- colnames(deviations)
group_colors <- c("Bcell" = "#08519c", "Tcell" = "#a50f15")
ann_heatmap <- HeatmapAnnotation(df = ann, col = list(groups = group_colors))
dend_cols <- cluster_within_group(deviations, ann$groups)
colors <- muRtools::colpal.cont(nrow(ann), "cptcity.arendal_temperature")

# Plot heatmap
fn <- file.path(plot_dir, paste0("deviation_score_diff_", motifset, ".pdf"))
pdf(fn, width = 15, height = 15, onefile = FALSE)
Heatmap(as.matrix(deviations),
  name = "deviation_score",
  cluster_rows = TRUE, show_row_names = TRUE,
  col = colors,
  cluster_columns = dend_cols, show_column_names = FALSE,
  top_annotation = ann_heatmap
)
dev.off()

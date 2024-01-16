setwd("/icbb/projects/igunduz/methylTFR_manuscript/")
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
})
source("/icbb/projects/igunduz/methylTFR_manuscript/src/utils.R")

plot_dir <- "/icbb/projects/igunduz/methylTFR_manuscript/Figures"
if(!dir.exists(plot_dir)) dir.create(plot_dir)
#Full blueprint methylTFR results
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
autoplot(pca, data = tdf,
         colour = 'cell_type',  
         size = 5,
         main = "PCA - Blueprint") +
         theme_classic() +
         theme(legend.position = "right")+
         scale_color_manual(values=colors)
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























#####################################################################
mtfr <- readRDS("/icbb/projects/igunduz/methylTFR_manuscript/results/BLUEPRINT/mtfr_final_221223/cisbpv2_deviations.RDS")
motifset <- "cisbpv2"
deviations <- deviations(mtfr)

fn <- file.path(plot_dir, paste0("proximal_devs_",motifset, ".pdf"))
pdf(fn, width = 15, height = 15, onefile = FALSE)
heatmap.2(as.matrix(deviations), col = bluered(100), trace = c("none"),
          dendrogram = c("column"), cexCol = 0.5, density.info = c("none"),
          keysize = 1,
          main = "Proximal motifs Deviations")
dev.off()

fn <- file.path(plot_dir, paste0("deviation_score_all_",motifset, ".pdf"))
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
pca <- prcomp(tdf, scale. = T)
tdf$cell_type <- groups
res.pca <- prcomp(t(deviations))


# Save eigenvalues plot
fn_eig <- file.path(plot_dir, paste0("eigenvalues_plot_", motifset, ".pdf"))
pdf(fn_eig)
fviz_eig(res.pca)
dev.off()

# Save PCA individual plot
fn_pca_ind <- file.path(plot_dir, paste0("pca_individuals_", motifset, ".pdf"))
pdf(fn_pca_ind)
fviz_pca_ind(res.pca, repel = TRUE, title = "PCA on all blueprint samples")
dev.off()


#Re-group T-cells
#tdf$cell_type  <- ifelse(tdf$cell_type %in% c("TCD4", "TCD8"),"Tcell", tdf$cell_type)


# skip cell_type column 
k <- ifelse(motifset == "altius", 287, 633)
pca <- prcomp(tdf[, -k])
fn_pca_ind <- file.path(plot_dir, paste0("bcell_vs_tcell_blueprint_", motifset, ".pdf"))
pdf(fn_pca_ind)
autoplot(pca, data = tdf,
         colour = 'cell_type',  
         size = 5,
         main = "PCA - Bcell vs Tcell") +
         theme_classic() +
         theme(legend.position = "bottom")
dev.off()


# merging TCD4 and TCD8 into single group
match <- which(groups %in% c("TCD4", "TCD8"))
groups[match] <- "Tcell"
groups <- as.factor(groups)

tdf <- as.data.frame(deviations)
diff <- differential_deviation_test(tdf, groups = groups,alternative = "two.sided",parametric =FALSE)

dim(diff[diff$p_value_adjusted < 0.05, ])
head(diff[diff$p_value_adjusted < 0.05, ])

write.table(diff, file = paste0(plot_dir,"/",motifset,"_diff_devs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

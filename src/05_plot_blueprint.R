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
if(!dir.exists(plot_dir)) dir.create(plot_dir)
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

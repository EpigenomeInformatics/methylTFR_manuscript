# expr: genes x samples matrix of log2(FPKM+1) (or TPM+1)
library(dorothea)
library(viper)
library(dplyr)
library(ComplexHeatmap)

# Expr matrix
source("/icbb/projects/nitschre/methylTFR/scripts/memoryTcells/06_rna_expr_matrix.R")

# Change to gene short names
expr_matrix_log <- merge(expr_matrix_log, tracking2gene, by="tracking_id")

# Duplicate tracking id replaced with mean gene exppression
expr <- expr_matrix_log %>%
  group_by(gene_short_name) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))

# Gene ids as rownames
expr <- as.data.frame(expr)
rownames(expr) <- expr$gene_short_name
expr <- expr[,!colnames(expr) %in% c("gene_short_name", "tracking_id")]

# 1) Get human regulons (A–C for higher confidence)
data(dorothea_hs, package = "dorothea")
reg <- dorothea_hs %>% filter(confidence %in% c("A","B","C"))

# 2) Build regulon object for VIPER
regulon <- dorothea::df2regulon(reg)

# 3) Single-sample TF activity (aREA)
# expr must be a matrix with gene symbols as rownames
tf_activity <- viper(as.matrix(expr), regulon, method = "scale", minsize = 5, eset.filter = FALSE)

# tf_activity: TF x samples matrix of normalized enrichment scores

### Heatmap containing differential Tfs from mtfr
mtfr_tfs <- readRDS("/icbb/projects/nitschre/methylTFR/r_objects/zscores_top50_memTcells.RDS")

tf_activity_filtered <- tf_activity[rownames(tf_activity) %in% rownames(mtfr_tfs),]

# Change sample names
colnames(tf_activity_filtered) <- c("Hf03_CM", "Hf03_EM", "Hf03_TN", "Hf04_CM", "Hf04_TN")

ha <- HeatmapAnnotation(
  celltypes=c("CM","EM", "TN", "CM", "TN"),
  col = list(celltypes = c("TN" = "#C8E0B4", "CM"="#4492C6", "EM"="#43B6C4"))
)

path <- "/icbb/projects/nitschre/methylTFR/figures/memoryTcells/heatmap_gex_tf_activity.pdf"
col_fun <- colorRamp2(seq(min(tf_activity_filtered), max(tf_activity_filtered), length.out = 100),
                      viridis(100))
pdf(path)
Heatmap(
  tf_activity_filtered, column_names_gp = gpar(fontsize = 9),
  top_annotation = ha,
  column_title = "Genexpression based TF activity",
  show_row_names = TRUE,
  column_names_side = "top",
  column_order = c("Hf03_CM","Hf04_CM", "Hf03_EM", "Hf03_TN", "Hf04_TN" ),
  col = col_fun
  )
dev.off()


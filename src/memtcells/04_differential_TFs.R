#!/usr/bin/env Rscript

#####################################################################
# 04_differential_TFs.R
# Created on 27-10-25 by IBG
# Identify differential TFs from methylTFR analysis and
# TF activity analysis based on gene expression data with dorothea/viper
#####################################################################

# load libraries
suppressPackageStartupMessages({
  library(dorothea)
  library(viper)
  library(dplyr)
  library(ComplexHeatmap)
  library(org.Hs.eg.db)
  library(circlize)
  library(viridis)
  library(methylTFR)
  library(ComplexHeatmap)
})
set.seed(13)

plot_dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR_manuscript/figures/memtcells"
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}
mem_tcell <- "/scratch/icbb/igunduz/methylTFR_manuscript/memTcell"
if (!dir.exists(mem_tcell)) {
  dir.create(mem_tcell, recursive = TRUE)
}

#########################################################################
# Get differential TFs from methylTFR analysis
#########################################################################

# Loading deviations scores
deviations_raw <- readRDS("/icbb/projects/nitschre/methylTFR/r_objects/jaspar2020_distal_deviations.RDS")
deviations <- deviations(deviations_raw)

### EM vs TN
# Subsetting deviation matrix to only naive and EM cells
deviations_em <- deviations[, grepl("TN|EM", colnames(deviations))]

# Set groups to compare
groups <- c("EM", "TN", "EM", "TN")

# differntial analysis
diff_em <- differential_deviation_test(deviations_em,
  groups = groups,
  alternative = "two.sided", parametric = TRUE
)
saveRDS(diff_em, file = file.path(mem_tcell, "diff_em_jaspar2020.RDS"))

# Filter for pval < .05
diff_em_filtered <- diff_em[diff_em$p_value_adjusted < 0.05, ]

### CM vs TN
# Subsetting deviations matrix to only naive and CM cells
deviations_cm <- deviations[, grepl("TN|CM", colnames(deviations))]

# groups
groups <- c("CM", "TN", "CM", "TN")

# differential analysis
diff_cm <- differential_deviation_test(deviations_cm,
  groups = groups,
  alternative = "two.sided", parametric = TRUE
)
saveRDS(diff_cm, file = file.path(mem_tcell, "diff_cm_jaspar2020.RDS"))

# Filter for pval < .05
diff_cm_filtered <- diff_cm[diff_cm$p_value_adjusted < 0.05, ]

# Get Z-scores
zscores <- deviationZScores(deviations_raw)

# Filter for diff TFs
zscores_filtered <- zscores[rownames(zscores) %in% c(rownames(diff_em_filtered), rownames(diff_cm_filtered)), ]
saveRDS(zscores_filtered, file = file.path(mem_tcell, "zscores_diffmotifs_jaspar2020.RDS"))

#########################################################################
# VIPER analysis for TF activity from gene expression data
#########################################################################

# Get the RNA-seq data files
rna_files <- list.files("/icbb/projects/share/datasets/memoryTcells/rna_raw_081025",
  pattern = "IHECrsem.20190606.hs38.genes.results$", full.names = TRUE
)

# Get the methylation data files
meth_files <- list.files("/icbb/projects/share/datasets/memoryTcells/",
  pattern = "bed$", full.names = TRUE
)

# Extract sample IDs from filenames
meth_sample_ids <- gsub(".*?/(\\d+_Hf\\d+_Bl\\w+_Ct)_.*", "\\1", meth_files)
rna_sample_ids <- sub("\\..*", "", basename(rna_files))

# Subset RNA files where sample ID matches methylation data
matched_rna_files <- rna_files[rna_sample_ids %in% meth_sample_ids]

# Load by tracking_id instead of gene_short_name
rna_list <- lapply(matched_rna_files, function(file) {
  df <- read.delim(file, stringsAsFactors = FALSE)
  sample_id <- tools::file_path_sans_ext(basename(file))
  df <- df[, c("gene_id", "expected_count")]
  colnames(df)[2] <- sample_id
  return(df)
})

# Merge by tracking_id (unique)
expr_matrix <- Reduce(function(x, y) full_join(x, y, by = "gene_id"), rna_list)

# Map Ensembl IDs to gene symbols
expr_matrix$gene_id_clean <- sub("\\..*$", "", expr_matrix$gene_id)
mapping <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = expr_matrix$gene_id_clean,
  keytype = "ENSEMBL",
  columns = c("SYMBOL")
)
colnames(mapping) <- c("gene_id_clean", "gene_short_name")

# Remove NAs from mapping
mapping <- mapping[!is.na(mapping$gene_short_name), ]

# Merge back with expression table:
expr_matrix <- merge(expr_matrix, mapping, "gene_id_clean")

# Gene ids as rownames
expr <- as.data.frame(expr_matrix)
rownames(expr) <- make.unique(expr$gene_short_name)
expr <- expr[, !colnames(expr) %in% c("gene_id_clean", "gene_id", "gene_short_name", "tracking_id")]

# 1) Get human regulons (A–C for higher confidence)
data(dorothea_hs, package = "dorothea")
reg <- dorothea_hs %>% filter(confidence %in% c("A", "B", "C"))

# 2) Build regulon object for VIPER
regulon <- dorothea::df2regulon(reg)

# 3) Single-sample TF activity (aREA)
# expr must be a matrix with gene symbols as rownames
tf_activity <- viper(as.matrix(expr), regulon, method = "scale", minsize = 5, eset.filter = FALSE)
saveRDS(tf_activity, file = file.path(mem_tcell, "viper_tf_activity_memTcells.RDS"))

### Heatmap containing differential Tfs from mtfr
mtfr_tfs <- readRDS(file.path(mem_tcell, "zscores_diffmotifs_jaspar2020.RDS"))
tf_activity_filtered <- tf_activity[rownames(tf_activity) %in% rownames(mtfr_tfs), ]

# Change sample names
colnames(tf_activity_filtered) <- c("Hf03_CM", "Hf03_EM", "Hf03_TN", "Hf04_CM", "Hf04_EM", "Hf04_TN")

# Rowwise zscores of gex data
tf_activity_filtered <- methylTFR:::computeRowZScore(tf_activity_filtered)

ha <- HeatmapAnnotation(
  celltypes = c("CM", "EM", "TN", "CM", "EM", "TN"),
  col = list(celltypes = c("TN" = "#C8E0B4", "CM" = "#4492C6", "EM" = "#43B6C4"))
)

path <- file.path(plot_dir, "viper_tf_activity_memTcells_mtfr_tfs.pdf")
col_fun <- colorRamp2(
  seq(-2, 2, length.out = 100),
  viridis(100)
)

# Define sample types for column splitting
sample_types <- sub("^.*_", "", colnames(tf_activity_filtered))
sample_types <- factor(sample_types, levels = c("CM", "EM", "TN"))

hm <- Heatmap(
  tf_activity_filtered,
  column_names_gp = gpar(fontsize = 9),
  top_annotation = ha,
  column_title = "VIPER TF Expression",
  show_row_names = TRUE,
  column_split = sample_types,
  cluster_columns = TRUE,
  cluster_column_slices = TRUE,
  column_order = NULL,
  col = col_fun
)
pdf(path, width = 10, height = 10)
ht <- draw(hm)
dev.off()

# Get the row order from heatmap
row_order_vector <- rownames(tf_activity_filtered)[row_order(ht)]

# Save the TF activity matrix with ordered rows
tf_activity_final <- tf_activity_filtered[row_order_vector, ]
saveRDS(tf_activity_final, file = file.path(mem_tcell, "viper_tf_activity_memTcells_mtfr_tfs.RDS"))

#########################################################################
# methylTFR heatmap for the differential TFs
#########################################################################

# load zscores for differential motifs
zscores_diffmotifs <- readRDS(file.path(mem_tcell, "zscores_diffmotifs_jaspar2020.RDS"))
zscores_diffmotifs <- zscores_diffmotifs[row_order_vector, ]
colnames(zscores_diffmotifs) <- c(
  "Hf03_CM", "Hf03_EM", "Hf03_TN", "Hf04_CM", "Hf04_EM",
  "Hf04_TN", "Hf04_TEMRA"
)
# Remove TEMRA
zscores_diffmotifs <- zscores_diffmotifs[, !colnames(zscores_diffmotifs) %in% c("Hf04_TEMRA")]

# Heatmap annotation
ha <- HeatmapAnnotation(
  celltypes = c("CM", "EM", "TN", "CM", "EM", "TN"),
  col = list(celltypes = c("TN" = "#ff0000", "CM" = "#008cff", "EM" = "#0037ff"))
)

# Define sample types for column splitting
sample_types <- sub("^.*_", "", colnames(zscores_diffmotifs))
sample_types <- factor(sample_types, levels = c("CM", "EM", "TN"))

path <- file.path(plot_dir, "mtfr_zscores_diffmotifs_memTcells.pdf")
hm <- Heatmap(
  zscores_diffmotifs,
  column_names_gp = gpar(fontsize = 9),
  top_annotation = ha,
  column_title = "methylTFR Z-scores",
  show_row_names = TRUE,
  column_split = sample_types,
  cluster_columns = TRUE,
  cluster_rows = FALSE,
  cluster_column_slices = TRUE,
  column_order = NULL
)
pdf(path, width = 10, height = 10)
ht <- draw(hm)
dev.off()

#########################################################################

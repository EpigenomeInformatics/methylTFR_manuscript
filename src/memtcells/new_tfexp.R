suppressPackageStartupMessages({
  library(dplyr)
  library(ComplexHeatmap)
  library(org.Hs.eg.db)
  library(circlize)
  library(ggplot2)
  library(ggrepel)
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
rna_files <- list.files("/icbb/projects/share/datasets/memoryTcells/rna_raw_081025",
  pattern = "IHECrsem.20190606.hs38.genes.results$", full.names = TRUE
)
meth_files <- list.files("/icbb/projects/share/datasets/memoryTcells/",
  pattern = "bed$", full.names = TRUE
)
meth_sample_ids <- gsub(".*?/(\\d+_Hf\\d+_Bl\\w+_Ct)_.*", "\\1", meth_files)
rna_sample_ids <- sub("\\..*", "", basename(rna_files))
matched_rna_files <- rna_files[rna_sample_ids %in% meth_sample_ids]
rna_list <- lapply(matched_rna_files, function(file) {
  df <- read.delim(file, stringsAsFactors = FALSE)
  sample_id <- tools::file_path_sans_ext(basename(file))
  df <- df[, c("gene_id", "expected_count")]
  colnames(df)[2] <- sample_id
  return(df)
})
expr_matrix <- Reduce(function(x, y) full_join(x, y, by = "gene_id"), rna_list)
expr_matrix$gene_id_clean <- sub("\\..*$", "", expr_matrix$gene_id)
mapping <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = expr_matrix$gene_id_clean,
  keytype = "ENSEMBL",
  columns = c("SYMBOL")
)
colnames(mapping) <- c("gene_id_clean", "gene_short_name")
mapping <- mapping[!is.na(mapping$gene_short_name), ]
expr_matrix <- merge(expr_matrix, mapping, "gene_id_clean")
expr <- as.data.frame(expr_matrix)
rownames(expr) <- make.unique(expr$gene_short_name)
expr <- expr[, !colnames(expr) %in% c("gene_id_clean", "gene_id", "gene_short_name", "tracking_id")]
mtfr_tfs <- readRDS(file.path(mem_tcell, "zscores_diffmotifs_jaspar2020.RDS"))
tf_expression_filtered <- expr[rownames(expr) %in% rownames(mtfr_tfs), ]
original_cols <- colnames(tf_expression_filtered)
clean_cols <- sub("\\..*", "", original_cols)
clean_cols <- sub("^\\d+_", "", clean_cols)
clean_cols <- sub("_BlEM", "_EM", clean_cols)
clean_cols <- sub("_BlCM", "_CM", clean_cols)
clean_cols <- sub("_BlTN", "_TN", clean_cols)
clean_cols <- sub("_Ct$", "", clean_cols)
cat("DEBUG: Column Name Mapping\n")
print(data.frame(Original = original_cols, Cleaned = clean_cols))
colnames(tf_expression_filtered) <- clean_cols
tf_expression_zscore <- methylTFR:::computeRowZScore(as.matrix(tf_expression_filtered))
sample_types <- sub("^.*_", "", colnames(tf_expression_zscore))
sample_types_factor <- factor(sample_types, levels = c("CM", "EM", "TN"))
ha <- HeatmapAnnotation(
  celltypes = sample_types,
  col = list(celltypes = c("TN" = "#C8E0B4", "CM" = "#4492C6", "EM" = "#43B6C4"))
)
path <- file.path(plot_dir, "gex_tf_expression_memTcells_mtfr_tfs.pdf")
col_fun <- colorRamp2(
  seq(-2, 2, length.out = 70),
  viridis(70)
)
hm <- Heatmap(
  tf_expression_zscore,
  column_names_gp = gpar(fontsize = 9),
  top_annotation = ha,
  column_title = "TF Gene Expression (Z-Score)",
  show_row_names = TRUE,
  column_split = sample_types_factor,
  cluster_columns = TRUE,
  cluster_column_slices = TRUE,
  column_order = NULL,
  col = col_fun
)
pdf(path, width = 10, height = 15)
ht <- draw(hm)
dev.off()
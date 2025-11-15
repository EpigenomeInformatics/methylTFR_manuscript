#!/usr/bin/env Rscript

#####################################################################
# 03_differential_TFs.R
# Created on 27-10-25 by IBG
# Identify differential TFs from methylTFR analysis and
# TF activity analysis based on gene expression data
#####################################################################

# load libraries
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
table_dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR_manuscript/tables"
if (!dir.exists(table_dir)) {
  dir.create(table_dir, recursive = TRUE)
}
cut_padj <- 0.05
cut_mean_diff <- 0.0

#########################################################################
# Get differential TFs from methylTFR analysis
#########################################################################

# Loading deviations scores
deviations_raw <- readRDS("/scratch/icbb/igunduz/methylTFR_manuscript/memoryTcells/mTFR_devs_131125/jaspar2020_distal_deviations.RDS")
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
write.csv(diff_em, file = file.path(table_dir, "diff_em_jaspar2020.csv"), row.names = TRUE)

# Filter for top 50 TFs based on p_value_adjusted
diff_em_ordered <- diff_em[order(diff_em[, "p_value_adjusted"]), ]
top_em_tfs <- rownames(diff_em_ordered)[1:min(50, nrow(diff_em_ordered))]


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
write.csv(diff_cm, file = file.path(table_dir, "diff_cm_jaspar2020.csv"), row.names = TRUE)

# Filter for top 50 TFs based on p_value_adjusted
diff_cm_ordered <- diff_cm[order(diff_cm[, "p_value_adjusted"]), ]
top_cm_tfs <- rownames(diff_cm_ordered)[1:min(50, nrow(diff_cm_ordered))]

# Get Z-scores
zscores <- deviationZScores(deviations_raw)

# Filter for diff TFs (top 50 from each comparison)
all_top_tfs <- unique(c(top_em_tfs, top_cm_tfs))
zscores_filtered <- zscores[rownames(zscores) %in% all_top_tfs, ]
saveRDS(zscores_filtered, file = file.path(mem_tcell, "zscores_diffmotifs_jaspar2020.RDS"))

#########################################################################
# GEX analysis for TF activity from gene expression data
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

### Heatmap containing differential Tfs from mtfr
mtfr_tfs <- readRDS(file.path(mem_tcell, "zscores_diffmotifs_jaspar2020.RDS"))
tf_expression_filtered <- expr[rownames(expr) %in% rownames(mtfr_tfs), ]

# Clean sample names
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

# Rowwise zscores of gex data
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

hm_gex <- Heatmap(
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
ht_gex <- draw(hm_gex)
dev.off()

# Get the row order from heatmap
row_order_vector <- rownames(tf_expression_zscore)[row_order(ht_gex)]


#########################################################################
# methylTFR heatmap for the differential TFs
#########################################################################

# load zscores for differential motifs
zscores_diffmotifs <- readRDS(file.path(mem_tcell, "zscores_diffmotifs_jaspar2020.RDS"))

# Order rows based on the GEX heatmap clustering
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
  col = list(celltypes = c("TN" = "#C8E0B4", "CM" = "#4492C6", "EM" = "#43B6C4"))
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
  cluster_rows = FALSE, # Use the order from GEX heatmap
  cluster_column_slices = TRUE,
  column_order = NULL
)
pdf(path, width = 10, height = 15)
ht <- draw(hm)
dev.off()


#########################################################################
# Mean difference plots for EM vs TN and CM vs TN
#########################################################################
cut_mean_diff <- 0.1
cut_padj <- 0.05
deviations <- deviationZScores(deviations_raw)
deviations_em <- deviations[, grepl("TN|EM", colnames(deviations))]
deviations_cm <- deviations[, grepl("TN|CM", colnames(deviations))]

zscores_em <- deviations_em
zscores_cm <- deviations_cm
# Apply Z-scores seperately : This gave similar results
# zscores_em <- methylTFR:::computeRowZScore(as.matrix(deviations_em))
# zscores_cm <- methylTFR:::computeRowZScore(as.matrix(deviations_cm))

# Compute group mean
group_means_em <- data.frame(
  EM = rowMeans(zscores_em[, grepl("TN", colnames(zscores_em))]) - rowMeans(zscores_em[, grepl("EM", colnames(zscores_em))])
)
group_means_em$motifs <- rownames(group_means_em)
group_means_cm <- data.frame(
  CM = rowMeans(zscores_em[, grepl("TN", colnames(zscores_em))]) - rowMeans(zscores_cm[, grepl("CM", colnames(zscores_cm))])
)
group_means_cm$motifs <- rownames(group_means_cm)

# Add differential information to group_means
diff_em <- readRDS(file.path(mem_tcell, "diff_em_jaspar2020.RDS"))
diff_cm <- readRDS(file.path(mem_tcell, "diff_cm_jaspar2020.RDS"))

group_means_em <- group_means_em %>%
  mutate(TF = rownames(diff_em)) %>%  # 
  left_join(diff_em, by = "motifs")  
group_means_cm <- group_means_cm %>%
  mutate(TF = rownames(diff_cm)) %>%
  left_join(diff_cm, by = "motifs")

# combine the results
comb_df <- merge(group_means_em, group_means_cm, by = "motifs", suffixes = c("_em", "_cm"))

# Add differential identifiers with three groups
comb_df$isDiff <- dplyr::case_when(
  abs(comb_df$mean_difference_em) > cut_mean_diff & comb_df$p_value_adjusted_em < cut_padj &
    abs(comb_df$mean_difference_cm) > cut_mean_diff & comb_df$p_value_adjusted_cm < cut_padj ~ "Differential Both",
  abs(comb_df$mean_difference_em) > cut_mean_diff & comb_df$p_value_adjusted_em < cut_padj ~ "Differential EM",
  abs(comb_df$mean_difference_cm) > cut_mean_diff & comb_df$p_value_adjusted_cm < cut_padj ~ "Differential CM",
  TRUE ~ "Not Differential"
)

cor_val <- cor(comb_df$CM, comb_df$EM, method = "pearson", use = "complete.obs")
plot_path <- file.path(plot_dir, "mean_difference_em_vs_cm_memTcells.pdf")

# Reset label column
comb_df$label <- NA

# Top 20 Differential EM
em_rows <- which(comb_df$isDiff == "Differential EM")
top_em_idx <- order(abs(comb_df$mean_difference_em[em_rows]), decreasing = TRUE)[1:min(20, length(em_rows))]
comb_df$label[em_rows[top_em_idx]] <- comb_df$TF_em[em_rows[top_em_idx]]

# Top 20 Differential CM
cm_rows <- which(comb_df$isDiff == "Differential CM")
top_cm_idx <- order(abs(comb_df$mean_difference_cm[cm_rows]), decreasing = TRUE)[1:min(20, length(cm_rows))]
comb_df$label[cm_rows[top_cm_idx]] <- comb_df$TF_cm[cm_rows[top_cm_idx]]
# Subset only rows to be labeled
label_df <- comb_df[!is.na(comb_df$label), ]
p <- ggplot(comb_df, aes(x = CM, y = EM, color = isDiff)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_manual(
    values = c(
      "Differential CM" = "#377EB8",
      "Differential EM" = "#E41A1C",
      "Differential Both" = "#4DAF4A",
      "Not Differential" = "grey70"
    )
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.3) +
  labs(
    x = "Mean difference (CM)",
    y = "Mean difference (EM)",
    color = "Differential status"
  ) +
  theme_classic(base_size = 13) +
  theme(legend.position = "top") +
  annotate("text", x = Inf, y = Inf, label = paste0("r = ", round(cor_val, 2)),
           hjust = 1.1, vjust = 1.5, size = 4.5, fontface = "bold") +
  geom_text_repel(
    data = label_df,
    aes(label = label, color = isDiff),
    size = 5,
    box.padding = 0.3,
    point.padding = 0.3,
    max.overlaps = Inf,
    segment.color = NA  # removes the lines
  )

ggsave(plot_path, p, width = 10, height = 10)


#########################################################################
# L2FC from LOLA vs mean difference from methylTFR
#########################################################################
# Load LOLA results
res_lola <- readRDS("/scratch/icbb/igunduz/methylTFR_manuscript/memoryTcells/reports/differential_methylation_data/differential_rnbDiffMeth/TF_motifs_lola.rds")

# Get mTFR difference for TN/TEM
group_means_em$zdiff <- group_means_em$EM
group_means_cm$zdiff <- group_means_cm$CM
em_diff <- group_means_em[, c("motifs", "zdiff","p_value_adjusted")]
cm_diff <- group_means_cm[, c("motifs", "zdiff","p_value_adjusted")]

# Get the LOLA enrichment results for TN/EM and TN/CM
em_lola <- res_lola$region[[3]]$tiling
cm_lola <- res_lola$region[[2]]$tiling

# Subset the userSet based on "rankCut_100_hyper" and, "rankCut_100_hypo"
em_lola <- em_lola[em_lola$userSet %in% c("rankCut_1000_hyper", "rankCut_1000_hypo"), ]
cm_lola <- cm_lola[cm_lola$userSet %in% c("rankCut_1000_hyper", "rankCut_1000_hypo"), ]

# Add condition column, if hyper then its "TEM" else "TN"
em_lola$condition <- ifelse(em_lola$userSet == "rankCut_1000_hyper", "EM", "TN")
cm_lola$condition <- ifelse(cm_lola$userSet == "rankCut_1000_hyper", "CM", "TN")

# Calculate log2OR
em_lola$log2OR <- log2(em_lola$oddsRatio)
cm_lola$log2OR <- log2(cm_lola$oddsRatio)
#em_lola$qValue <- -log10(em_lola$qValue)
#cm_lola$qValue <- -log10(cm_lola$qValue)

# Multiply log2OR with - if condition is TN
em_lola$log2OR <- ifelse(em_lola$condition == "TN", -em_lola$log2OR, em_lola$log2OR)
cm_lola$log2OR <- ifelse(cm_lola$condition == "TN", -cm_lola$log2OR, cm_lola$log2OR)

# Re-arrange motifs
em_lola$motifs <- gsub(".*_", "", em_lola$description)
cm_lola$motifs <- gsub(".*_", "", cm_lola$description)

em_lola <- em_lola[, c("motifs", "log2OR", "qValue")]
cm_lola <- cm_lola[, c("motifs", "log2OR", "qValue")]

# Merge with mTFR mean difference
em_merged <- merge(em_diff, em_lola, by = "motifs")
em_merged_unique <- em_merged[ave(em_merged$qValue, em_merged$motifs, FUN = rank) == 1, ]
cm_merged <- merge(cm_diff, cm_lola, by = "motifs")
cm_merged_unique <- cm_merged[ave(cm_merged$qValue, cm_merged$motifs, FUN = rank) == 1, ]

# Plotting function
plotlog2OR <- function(df){
  x_max <- max(abs(df$log2OR), na.rm = TRUE)
  y_max <- max(abs(df$zdiff), na.rm = TRUE)
  
  df <- df %>%
    mutate(
      isDiff = case_when(
        p_value_adjusted < 0.05 & qValue < 0.05 ~ "Differential in both",
        p_value_adjusted < 0.05 ~ "mTFR differential",
        qValue < 0.05 ~ "LOLA differential",
        TRUE ~ "Not differential"
      )
    )
  
  df$isDiff <- factor(df$isDiff, levels = c(
    "Differential in both", 
    "mTFR differential", 
    "LOLA differential", 
    "Not differential"
  ))
  
  # Create label column: NA if "Not differential"
  df$label <- ifelse(df$isDiff != "Not differential", df$motifs, NA)
  
  p <- ggplot(df, aes(x = log2OR, y = zdiff, color = isDiff)) +
    geom_point(alpha = 0.7, size = 3) +
    
    labs(
      title = paste0("Motif Enrichment vs. Activity Difference"),
      x = expression(Log[2]~"Odds Ratio (LOLA Enrichment)"),
      y = expression("Z-Score Difference (mTFR Activity)"),
      color = "Differential Status"
    ) +
    scale_x_continuous(limits = c(-x_max, x_max)) +
    scale_y_continuous(limits = c(-y_max, y_max)) +
    theme_classic(base_size = 14) +
    scale_color_manual(values = c(
      "Differential in both" = "red", 
      "mTFR differential" = "darkgreen", 
      "LOLA differential" = "dodgerblue", 
      "Not differential" = "gray50"
    )) +
    # Use the new 'label' column for repel
    geom_text_repel(aes(label = label), size = 3.5, box.padding = 0.5, max.overlaps = 15) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "black") +
    geom_vline(xintercept = 0, linetype = "dotted", color = "black")
  
  return(p)
}

# Plot for EM vs TN
p_em <- plotlog2OR(em_merged_unique)
ggsave(file.path(plot_dir, "lola_log2OR_vs_mtfr_meanDiff_em_vs_tn.pdf"), p_em, width = 8, height = 6)

# Plot for CM vs TN
p_cm <- plotlog2OR(cm_merged_unique)
ggsave(file.path(plot_dir, "lola_log2OR_vs_mtfr_meanDiff_cm_vs_tn.pdf"), p_cm, width = 8, height = 6)

#########################################################################
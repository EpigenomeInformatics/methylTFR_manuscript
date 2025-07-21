
#!/usr/bin/env Rscript

#####################################################################
# 06.rna_process.R
# created on 08-07-25 by Irem Gunduz
# Load RNA-seq data from memTcells and correalate with TF devs
#####################################################################
set.seed(12)
suppressPackageStartupMessages({
  library(dplyr)
  library(dorothea)
  library(decoupleR)
  library(methylTFR)
})

# Get the RNA-seq data files
rna_files <- list.files("/icbb/projects/share/datasets/memoryTcells/rna-seq", 
                    pattern = "genes.fpkm_tracking$", full.names = TRUE)

# Get the methylation data files
meth_files <- list.files("/icbb/projects/share/datasets/memoryTcells/", 
                         pattern = "bed$", full.names = TRUE)

# Extract sample IDs from filenames
meth_sample_ids <- gsub(".*?/(\\d+_Hf\\d+_Bl\\w+_Ct)_.*", "\\1", meth_files)
rna_sample_ids <- gsub("_mRNA$", "", rna_sample_ids)

# Subset RNA files where sample ID matches methylation data
matched_rna_files <- rna_files[rna_sample_ids %in% meth_sample_ids]

# Load by tracking_id instead of gene_short_name
fpkm_list <- lapply(matched_rna_files, function(file) {
  df <- read.delim(file, stringsAsFactors = FALSE)
  sample_id <- tools::file_path_sans_ext(basename(file))
  df <- df[, c("tracking_id", "FPKM")]
  colnames(df)[2] <- sample_id
  return(df)
})

# Merge by tracking_id (unique)
expr_matrix <- Reduce(function(x, y) full_join(x, y, by = "tracking_id"), fpkm_list)

# Log2-transform for correlation
expr_matrix_log <- expr_matrix
#expr_matrix_log[,-1] <- log2(expr_matrix_log[,-1] + 1)

# Link tracking_id to gene names
lookup <- read.delim(matched_rna_files[1], stringsAsFactors = FALSE)
tracking2gene <- lookup[, c("tracking_id", "gene_short_name")]
tracking2gene$gene_short_name <- toupper(tracking2gene$gene_short_name)

# Get the dorothea regulons
dorothea_hs <- decoupleR::get_dorothea(levels = c('A', 'B', 'C', 'D'))

# Keep only high-confidence interactions
dorothea_regulon <- subset(dorothea_hs, confidence %in% c("A", "B", 'C', 'D'))
dorothea_regulon$target <- toupper(dorothea_regulon$target)
dorothea_regulon$tf <- toupper(dorothea_regulon$source)

# Merge to get tracking_ids for DoRothEA targets
dorothea_with_tracking <- left_join(dorothea_regulon, tracking2gene, 
                                    by = c("target" = "gene_short_name"))

# Remove entries with no match
dorothea_with_tracking <- dorothea_with_tracking[!is.na(dorothea_with_tracking$tracking_id), ]

# Load memTcell deviation scores
deviation_scores <- readRDS("/icbb/projects/nitschre/methylTFR/r_objects/jaspar2020_distal_deviations.RDS")
deviation_scores <- deviationZscores(deviation_scores)

# From RNA-seq file names
rna_sample_ids <- gsub("_mRNA_.*", "", colnames(expr_matrix_log)[-1])
colnames(expr_matrix_log)[-1] <- rna_sample_ids

# From deviation score columns
dev_sample_ids <- gsub("_WGBS_.*", "", colnames(deviation_scores))
colnames(deviation_scores) <- dev_sample_ids

# Create lookup table for mapping
rna_map <- data.frame(
  full_rna_name = colnames(expr_matrix_log)[-1],
  core_id = rna_sample_ids,
  stringsAsFactors = FALSE
)

dev_map <- data.frame(
  full_dev_name = colnames(deviation_scores),
  core_id = dev_sample_ids,
  stringsAsFactors = FALSE
)

# Inner join to find common samples
matched_samples <- inner_join(rna_map, dev_map, by = "core_id")

# Reorder expression and deviation matrices to match
expr_matrix_log_matched <- expr_matrix_log[, c("tracking_id", matched_samples$full_rna_name)]
deviation_scores_matched <- deviation_scores[, matched_samples$full_dev_name]
deviation_scores_matched <- methylTFR:::computeRowZScore(deviation_scores_matched)


##################################################################################################
# Correlate expression with deviation scores
##################################################################################################

results <- list()

for (tf in rownames(deviation_scores_matched)) {
  # Get all known tracking_ids for this TF
  tf_targets <- dorothea_with_tracking %>%
    filter(tf == !!tf) %>%
    pull(tracking_id) %>%
    unique()

  # Filter for those that exist in the expression matrix
  valid_targets <- intersect(tf_targets, expr_matrix_log_matched$tracking_id)
  if (length(valid_targets) < 3) next  # skip TFs with too few targets

  tf_dev <- as.numeric(deviation_scores_matched[tf, ])

  tf_exprs <- expr_matrix_log_matched %>%
    filter(tracking_id %in% valid_targets)

  tf_result <- data.frame(
    tf = tf,
    tracking_id = tf_exprs$tracking_id,
    cor = NA,
    pval = NA
  )

  for (i in seq_len(nrow(tf_exprs))) {
    expr <- as.numeric(tf_exprs[i, -1])
    test <- suppressWarnings(cor.test(expr, tf_dev, method = "spearman"))
    tf_result$cor[i] <- test$estimate
    tf_result$pval[i] <- test$p.value
  }

  # Ensure tf_result is a data frame before adding to results
  results[[tf]] <- as.data.frame(tf_result)
}

# Combine all data frames into one
correlation_df <- bind_rows(results)

# Adjust p-values
correlation_df$p_adj <- p.adjust(correlation_df$pval, method = "fdr")

top_hits <- correlation_df %>%
  filter(p_adj < 0.05) %>%
  arrange(desc(abs(cor))) %>%
  head(20)

print(top_hits)

top_hits <- correlation_df %>%
  filter(pval < 0.05) %>%
  arrange(desc(abs(cor))) %>%
  head(20)

print(top_hits)


##################################################################################################

# Create group map: match your actual sample names here
sample_groups <- c(
  "51_Hf03_BlCM_Ct_mRNA_M_1.LXPv1.20150708_genes.fpkm_tracking" = "TCM",
  "51_Hf03_BlEM_Ct_mRNA_M_1.LXPv1.20150708_genes.fpkm_tracking" = "TEM",
  "51_Hf03_BlTN_Ct_mRNA_M_1.LXPv1.20150708_genes.fpkm_tracking" = "Tnaive",
  "51_Hf04_BlCM_Ct_mRNA_M_1.LXPv1.20150708_genes.fpkm_tracking" = "TCM",
  "51_Hf04_BlTN_Ct_mRNA_M_1.LXPv1.20150708_genes.fpkm_tracking" = "Tnaive"
)

expr_matrix_log_matched$tracking_id <- make.unique(expr_matrix_log_matched$tracking_id)

expr_mat <- expr_matrix_log_matched %>%
  column_to_rownames("tracking_id") %>%
  as.matrix()

get_tf_target_expression <- function(tf) {
  targets <- dorothea_with_tracking %>%
    filter(tf == !!tf) %>%
    pull(tracking_id) %>%
    unique()

  targets <- intersect(targets, rownames(expr_mat))
  if (length(targets) < 3) return(NULL)

  expr_subset <- expr_mat[targets, , drop = FALSE]
  means <- colMeans(expr_subset, na.rm = TRUE)

  # Build and return a data.frame
  df <- data.frame(
    tf = tf,
    sample = names(means),
    mean_target_expr = means,
    group = sample_groups[names(means)],
    stringsAsFactors = FALSE
  )

  return(df)
}

expr_profiles_list <- lapply(rownames(deviation_scores_matched), get_tf_target_expression)

# Filter out NULLs just in case
expr_profiles_list <- Filter(Negate(is.null), expr_profiles_list)

# Combine safely
expr_profiles <- do.call(rbind, expr_profiles_list)

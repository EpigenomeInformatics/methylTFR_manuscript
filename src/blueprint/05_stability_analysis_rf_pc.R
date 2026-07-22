#!/usr/bin/env Rscript

#####################################################################
# 05_stability_analysis_rf_pc.R
# Written by IBG on 02-11-2025
# Assess classification accuracy using Random Forests on PCA-reduced data
#
# Modified on 15-11-2025:
# - Select number of PCs based on cumulative variance explained
#   instead of a fixed number (n_pc = 20).
# - Added a cap (max_pc) to limit the total number of PCs selected.
#
# Updated: Fixed disease filtering to use the DISEASE column and
# corrected variable names during RnBeads sample removal.
#####################################################################

suppressPackageStartupMessages({
  library(methylTFR)
  library(RnBeads)
  library(mclust)
  library(dplyr)
  library(ComplexHeatmap)
  library(Cairo)
  library(randomForest)
})

set.seed(13)

source("/icbb/projects/nitschre/methylTFR/scripts/other/helpers.R")
fig_dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR_manuscript/figures/blueprint/"
table_dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR_manuscript/tables/"
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

sannot <- read.csv("/scratch/icbb/igunduz/methylTFR_manuscript/blueprint/RnBeads_291025/reports/data_import_data/annotation.csv", stringsAsFactors = FALSE)
sannot$bedFile <- as.character(sannot$bedFile)

var_threshold <- 0.9 # Use PCs explaining 90% of variance
max_pc <- 20 # Set a maximum cap for number of PCs

# Remap to cleaner group names
group_remap <- c(
  "Bcell" = "B-cells",
  "DC" = "DC",
  "eryt" = "Erythrocytes",
  "gran" = "Granulocytes",
  "megK" = "Megakaryocytes",
  "Mf" = "Mf",
  "mono" = "Monocytes",
  "NK" = "NK",
  "osteoclast" = "Osteoclast",
  "other" = "Other",
  "plasma" = "Plasma",
  "progenitor" = "Progenitors",
  "Tcell" = "T-cells",
  "thymocyte" = "Thymocyte"
)

#####################################################################
# Get the data
#####################################################################

# 1. Get mTFR deviation matrix and filter
mtfr_obj <- readRDS("/scratch/icbb/igunduz/methylTFR_manuscript/blueprint/mTFR_devs_121125/JASPAR2020_distal_deviations.RDS")
deviations_mat <- deviations(mtfr_obj)

# Match annotations to deviation matrix columns
matched_sannot_dev <- sannot[match(colnames(deviations_mat), sannot$bedFile), ]

# Remove deviations for samples with disease
keep_dev <- matched_sannot_dev$DISEASE == "None"
mtfr <- deviations_mat[, keep_dev, drop = FALSE]

# 2. Load RnBeads objects and filter
rnbeads <- load.rnb.set("/scratch/icbb/igunduz/methylTFR_manuscript/blueprint/RnBeads_291025/reports/data_import_data/rnb.set_preprocessed")

# Remove DISEASE samples (fixed variable name from rnb_set to rnbeads)
disease_idx <- which(as.character(rnbeads@pheno$DISEASE) != "None")
if (length(disease_idx) > 0) {
  rnbeads <- remove.samples(rnbeads, disease_idx)
}

# Extract methylation matrices from the filtered RnBeads object
distal <- meth(rnbeads, type = "distal")
tiling <- meth(rnbeads, type = "tiling1kb")

# 3. Extract and remap cell types for the remaining samples
sample_names <- colnames(mtfr)
matched_sannot_final <- sannot[match(sample_names, sannot$bedFile), ]

# Remap cell type group names using the group_remap dictionary
cell_types <- unname(group_remap[matched_sannot_final$cellTypeGroup])

#####################################################################
# Functions for cross-validation with RF
#####################################################################

cv_rf_accuracy_safe <- function(mat, labels, k = 5, ntree = 500) {
  stopifnot(ncol(mat) == length(labels))
  labels <- factor(labels)
  tbl <- table(labels)

  # keep only classes with >=2 samples (cannot safely CV otherwise)
  keep <- names(tbl)[tbl >= 2]
  removed <- setdiff(names(tbl), keep)
  if (length(removed) > 0) {
    message("Removing classes with <2 samples (cannot reliably CV them): ", paste(removed, collapse = ", "))
  }

  keep_idx <- which(labels %in% keep)
  if (length(keep_idx) < 2) stop("Not enough samples after removing rare classes.")

  mat2 <- mat[, keep_idx, drop = FALSE]
  labels2 <- droplevels(labels[keep_idx])

  # choose k no larger than min class size
  min_count <- min(table(labels2))
  k_use <- min(k, min_count)
  if (k_use < 2) k_use <- 2

  # stratified folds by assigning fold ids per class
  n <- ncol(mat2)
  folds <- integer(n)
  classes <- levels(labels2)
  for (cl in classes) {
    idx <- which(labels2 == cl)
    # shuffle then assign folds evenly
    idx <- sample(idx)
    folds[idx] <- rep(1:k_use, length.out = length(idx))
  }

  acc <- numeric(k_use)
  x <- t(mat2)
  for (i in 1:k_use) {
    test_idx <- which(folds == i)
    train_idx <- setdiff(seq_len(n), test_idx)
    rf <- randomForest(x[train_idx, , drop = FALSE], labels2[train_idx], ntree = ntree)
    pred <- predict(rf, x[test_idx, , drop = FALSE])
    acc[i] <- mean(pred == labels2[test_idx])
  }
  mean(acc)
}

cv_rf_stratified_splits <- function(mat, labels, repeats = 10, test_frac = 0.2, ntree = 500, seed = 13) {
  set.seed(seed) # top-level seed
  labels <- factor(labels)
  n <- ncol(mat)
  classes <- levels(labels)
  tbl <- table(labels)

  singletons <- names(tbl)[tbl == 1]
  if (length(singletons) > 0) message("Singleton classes kept only in train: ", paste(singletons, collapse = ", "))

  accs <- numeric(repeats)
  x <- t(mat)

  for (r in 1:repeats) {
    set.seed(seed + r) # ensures deterministic repeat
    test_idx <- integer(0)
    for (cl in classes) {
      idx <- which(labels == cl)
      if (length(idx) == 1) next
      nt <- max(1, floor(length(idx) * test_frac))
      test_idx <- c(test_idx, sample(idx, nt))
    }
    train_idx <- setdiff(seq_len(n), test_idx)
    rf <- randomForest(x[train_idx, , drop = FALSE], labels[train_idx], ntree = ntree)
    pred <- predict(rf, x[test_idx, , drop = FALSE])
    accs[r] <- mean(pred == labels[test_idx])
  }
  mean(accs, na.rm = TRUE)
}

# Helper function to select PCs based on cumulative variance, with a max cap
select_pcs_by_variance <- function(pca_obj, threshold = 0.90, max_pcs = 30) {
  variances <- pca_obj$sdev^2
  cum_var_prop <- cumsum(variances) / sum(variances)

  # Find the first component that meets or exceeds the threshold
  n_pcs_variance_candidates <- which(cum_var_prop >= threshold)

  n_pcs_variance <- if (length(n_pcs_variance_candidates) > 0) {
    n_pcs_variance_candidates[1] # Take the first one
  } else {
    length(pca_obj$sdev) # Use all PCs if threshold is never met
  }

  # Apply the cap
  n_pcs <- min(n_pcs_variance, max_pcs)

  # Also cap by the total number of available PCs
  n_pcs <- min(n_pcs, length(pca_obj$sdev))

  return(n_pcs)
}

#####################################################################
# Execute PCA and Random Forest
#####################################################################

# mtfr PCA
pca_mtfr <- prcomp(t(mtfr), center = FALSE, scale. = FALSE)
n_pc_mtfr <- select_pcs_by_variance(pca_mtfr, var_threshold, max_pcs = max_pc)
message(paste("mTFR: Using", n_pc_mtfr, "PCs (variance threshold", var_threshold, ", max", max_pc, ")"))
pcs_mtfr <- pca_mtfr$x[, 1:n_pc_mtfr]

# distal PCA
pca_distal <- prcomp(t(distal), center = FALSE, scale. = FALSE)
n_pc_distal <- select_pcs_by_variance(pca_distal, var_threshold, max_pcs = max_pc)
message(paste("Distal: Using", n_pc_distal, "PCs (variance threshold", var_threshold, ", max", max_pc, ")"))
pcs_distal <- pca_distal$x[, 1:n_pc_distal]

# tiling PCA
pca_tiling <- prcomp(t(tiling), center = FALSE, scale. = FALSE)
n_pc_tiling <- select_pcs_by_variance(pca_tiling, var_threshold, max_pcs = max_pc)
message(paste("Tiling: Using", n_pc_tiling, "PCs (variance threshold", var_threshold, ", max", max_pc, ")"))
pcs_tiling <- pca_tiling$x[, 1:n_pc_tiling]

# RF with stratified splits
acc_mtfr_pc <- cv_rf_stratified_splits(t(pcs_mtfr), cell_types, repeats = 100, test_frac = 0.2, seed = 42)
acc_distal_pc <- cv_rf_stratified_splits(t(pcs_distal), cell_types, repeats = 100, test_frac = 0.2, seed = 42)
acc_tiling_pc <- cv_rf_stratified_splits(t(pcs_tiling), cell_types, repeats = 100, test_frac = 0.2, seed = 42)

results <- data.frame(
  Representation = c("mTFR", "Distal", "Tiling1kb"),
  Accuracy = c(acc_mtfr_pc, acc_distal_pc, acc_tiling_pc),
  NumPCs = c(n_pc_mtfr, n_pc_distal, n_pc_tiling)
)
write.csv(results, file = paste0(table_dir, "rf_classification_accuracy_pc_stratified_splits_var90.csv"), row.names = FALSE)

# RF with 5-fold CV
set.seed(42)
acc_mtfr_pc_cv <- cv_rf_accuracy_safe(t(pcs_mtfr), cell_types, k = 5, ntree = 500)
acc_distal_pc_cv <- cv_rf_accuracy_safe(t(pcs_distal), cell_types, k = 5, ntree = 500)
acc_tiling_pc_cv <- cv_rf_accuracy_safe(t(pcs_tiling), cell_types, k = 5, ntree = 500)

res <- data.frame(
  Representation = c("mTFR", "Distal", "Tiling1kb"),
  Accuracy = c(acc_mtfr_pc_cv, acc_distal_pc_cv, acc_tiling_pc_cv),
  NumPCs = c(n_pc_mtfr, n_pc_distal, n_pc_tiling)
)
write.csv(res, file = paste0(table_dir, "rf_classification_accuracy_pc_5foldcv_var90.csv"), row.names = FALSE)
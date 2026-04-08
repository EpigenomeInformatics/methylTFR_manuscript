
#!/usr/bin/env Rscript

#####################################################################
# 05_stability_analysis_rf_pc.R
# Written by IBG on 02-11-2025
# Assess classification accuracy using Random Forests on PCA-reduced data
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
if(!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)
sannot <- read.csv("/scratch/icbb/igunduz/methylTFR_manuscript/blueprint/RnBeads_291025/reports/data_import_data/annotation.csv", stringsAsFactors = FALSE)
n_pc <- 20

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

# Get mtfr, distal and 1kbtiling matrix
mtfr <- readRDS("/scratch/icbb/igunduz/methylTFR_manuscript/blueprint/mTFR_devs_071125/jaspar2020_distal_deviations.RDS")
mtfr <- deviationZScores(mtfr)

# Load RnBeads objects
rnbeads <- load.rnb.set("/scratch/icbb/igunduz/methylTFR_manuscript/blueprint/RnBeads_291025/reports/data_import_data/rnb.set_preprocessed")

# Extract methylation matrices
distal <- meth(rnbeads, type = "distal")
tiling <- meth(rnbeads, type = "tiling1kb")

# Remap cell type group names
sannot$cellTypeGroup <- group_remap[sannot$cellTypeGroup]
sannot$bedFile <- as.character(sannot$bedFile)
sample_names <- colnames(mtfr)

# Match and extract cellTypeGroup
cell_types <- sannot$cellTypeGroup[match(sample_names, sannot$bedFile)]

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
  if(length(removed) > 0){
    message("Removing classes with <2 samples (cannot reliably CV them): ", paste(removed, collapse = ", "))
  }
  
  keep_idx <- which(labels %in% keep)
  if(length(keep_idx) < 2) stop("Not enough samples after removing rare classes.")
  
  mat2 <- mat[, keep_idx, drop = FALSE]
  labels2 <- droplevels(labels[keep_idx])
  
  # choose k no larger than min class size
  min_count <- min(table(labels2))
  k_use <- min(k, min_count)
  if(k_use < 2) k_use <- 2
  
  # stratified folds by assigning fold ids per class
  n <- ncol(mat2)
  folds <- integer(n)
  classes <- levels(labels2)
  for(cl in classes){
    idx <- which(labels2 == cl)
    # shuffle then assign folds evenly
    idx <- sample(idx)
    folds[idx] <- rep(1:k_use, length.out = length(idx))
  }
  
  acc <- numeric(k_use)
  x <- t(mat2)
  for(i in 1:k_use){
    test_idx <- which(folds == i)
    train_idx <- setdiff(seq_len(n), test_idx)
    rf <- randomForest(x[train_idx, , drop = FALSE], labels2[train_idx], ntree = ntree)
    pred <- predict(rf, x[test_idx, , drop = FALSE])
    acc[i] <- mean(pred == labels2[test_idx])
  }
  mean(acc)
}

cv_rf_stratified_splits <- function(mat, labels, repeats = 10, test_frac = 0.2, ntree = 500, seed = 13){
  set.seed(seed)  # top-level seed
  labels <- factor(labels)
  n <- ncol(mat)
  classes <- levels(labels)
  tbl <- table(labels)
  
  singletons <- names(tbl)[tbl == 1]
  if(length(singletons) > 0) message("Singleton classes kept only in train: ", paste(singletons, collapse = ", "))
  
  accs <- numeric(repeats)
  x <- t(mat)
  
  for(r in 1:repeats){
    set.seed(seed + r)  # ensures deterministic repeat
    test_idx <- integer(0)
    for(cl in classes){
      idx <- which(labels == cl)
      if(length(idx) == 1) next
      nt <- max(1, floor(length(idx) * test_frac))
      test_idx <- c(test_idx, sample(idx, nt))
    }
    train_idx <- setdiff(seq_len(n), test_idx)
    rf <- randomForest(x[train_idx, , drop=FALSE], labels[train_idx], ntree = ntree)
    pred <- predict(rf, x[test_idx, , drop=FALSE])
    accs[r] <- mean(pred == labels[test_idx])
  }
  mean(accs, na.rm = TRUE)
}

#####################################################################

# mtfr PCA (632 features)
pca_mtfr <- prcomp(t(mtfr), center = FALSE, scale. = FALSE)
pcs_mtfr <- pca_mtfr$x[, 1:min(n_pc, ncol(pca_mtfr$x))]

# distal PCA 
pca_distal <- prcomp(t(distal), center = FALSE, scale. = FALSE)
pcs_distal <- pca_distal$x[, 1:min(n_pc, ncol(pca_distal$x))]

# tiling PCA
pca_tiling <- prcomp(t(tiling), center = FALSE, scale. = FALSE)
pcs_tiling <- pca_tiling$x[, 1:min(n_pc, ncol(pca_tiling$x))]

acc_mtfr_pc   <- cv_rf_stratified_splits(t(pcs_mtfr), cell_types, repeats = 100, test_frac = 0.2,seed=42)
acc_distal_pc <- cv_rf_stratified_splits(t(pcs_distal), cell_types, repeats = 100, test_frac = 0.2,seed=42)
acc_tiling_pc <- cv_rf_stratified_splits(t(pcs_tiling), cell_types, repeats = 100, test_frac = 0.2,seed=42)

results <- data.frame(
  Representation = c("mTFR", "Distal", "Tiling1kb"),
  Accuracy = c(acc_mtfr_pc, acc_distal_pc, acc_tiling_pc)
)
#  Representation  Accuracy
#1           mTFR 0.8703333
#2         Distal 0.9076667
#3      Tiling1kb 0.8453333
write.csv(results, file = paste0(table_dir, "rf_classification_accuracy_pc_stratified_splits.csv"), row.names = FALSE)

set.seed(42)
acc_mtfr_pc <- cv_rf_accuracy_safe(t(pcs_mtfr), cell_types, k = 5, ntree = 500)
acc_distal_pc <- cv_rf_accuracy_safe(t(pcs_distal), cell_types, k = 5, ntree = 500)
acc_tiling_pc <- cv_rf_accuracy_safe(t(pcs_tiling), cell_types, k = 5, ntree = 500)

res <- data.frame(
  Representation = c("mTFR", "Distal", "Tiling1kb"),
  Accuracy = c(acc_mtfr_pc, acc_distal_pc, acc_tiling_pc)
)

#  Representation  Accuracy
#1           mTFR 0.8726874
#2         Distal 0.9108244
#3      Tiling1kb 0.8216488
 
write.csv(res, file = paste0(table_dir, "rf_classification_accuracy_pc_5foldcv.csv"), row.names = FALSE)
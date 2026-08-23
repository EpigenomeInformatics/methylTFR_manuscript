#!/usr/bin/env Rscript

#####################################################################
# 03_pca_stability_blueprint.R
# created on 21-07-25 by Irem B Gunduz
# Updated by IBG on 23-08-2026
# Merged 04_plot_blueprint_pca.R and 05_stability_analysis_rf_pc.R
# PCA of the mTFR deviations (genome wide and distal restricted),
# tiling1kb and distal methylation, followed by a Random Forest
# stability analysis on the same PCs
#####################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggfortify)
  library(logger)
  library(randomForest)
  library(RnBeads)
  library(methylTFR)
})
set.seed(42)

#####################################################################
# Settings
#####################################################################

# Motif sets whose deviations enter the PCA and the stability analysis
motifSets <- c("jaspar2020", "jaspar2020_distal")

# PCA is run on the raw matrices, as in the original scripts
pca.center <- FALSE
pca.scale <- FALSE

# Number of PCs kept for the Random Forest
var.threshold <- 0.9 # Use PCs explaining 90% of the variance
max.pc <- 20 # Cap for the number of PCs

# Random Forest settings
rf.ntree <- 500
rf.repeats <- 100
rf.test.frac <- 0.2
rf.folds <- 5
rf.seed <- 42

# Directories
analysis.dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/blueprint/"
rnb.tag <- "RnBeads_230826"
dev.tag <- "mTFR_devs_230826"
rnb.set.path <- file.path(analysis.dir, rnb.tag, "reports", "data_import_data", "rnb.set_preprocessed")
sannot.file <- file.path(analysis.dir, rnb.tag, "reports", "data_import_data", "annotation.csv")
dev.files <- setNames(
  file.path(analysis.dir, dev.tag, paste0(motifSets, "_deviations.RDS")),
  paste0("mTFR_", motifSets)
)

github.dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR_manuscript/"
fig.dir <- file.path(github.dir, "figures", "blueprint")
table.dir <- file.path(github.dir, "tables")
if (!dir.exists(fig.dir)) dir.create(fig.dir, recursive = TRUE)
if (!dir.exists(table.dir)) dir.create(table.dir, recursive = TRUE)

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

cell_type_colors <- c(
  "B-cells" = "#1f77b4",
  "DC" = "#ff7f0e",
  "Erythrocytes" = "#2ca02c",
  "Granulocytes" = "#d62728",
  "Megakaryocytes" = "#9467bd",
  "Mf" = "#8c564b",
  "Monocytes" = "#e377c2",
  "NK" = "#7f7f7f",
  "Osteoclast" = "#1b9e77",
  "Other" = "#bcbd22",
  "Plasma" = "#17becf",
  "Progenitors" = "#2ca4a2",
  "T-cells" = "#393b79",
  "Thymocyte" = "#6a5acd"
)

#####################################################################
# Helper functions
#####################################################################

# Drop features with non finite values, prcomp cannot handle them
clean_matrix <- function(mat, label) {
  mat <- as.matrix(mat)
  keep <- is.finite(rowSums(mat))
  if (any(!keep)) {
    log_info(label, ": dropping ", sum(!keep), " features with non finite values")
  }
  mat[keep, , drop = FALSE]
}

# PCA on samples, matrices are features x samples
run_pca <- function(mat, label) {
  log_info(label, ": running PCA on ", nrow(mat), " features x ", ncol(mat), " samples")
  prcomp(t(mat), center = pca.center, scale. = pca.scale)
}

plot_pca <- function(pca_obj, groups, title, file) {
  # Only the grouping column is passed, ggfortify binds the scores itself
  plot_data <- data.frame(groups = groups, stringsAsFactors = FALSE)
  pdf(file, width = 10, height = 10)
  print(
    autoplot(pca_obj,
      data = plot_data,
      colour = "groups",
      main = title,
      size = 5
    ) +
      theme_classic() +
      scale_color_manual(values = cell_type_colors) +
      theme(legend.position = "bottom")
  )
  dev.off()
  log_info("Wrote ", file)
}

# Select PCs based on cumulative variance, with a cap
select_pcs_by_variance <- function(pca_obj, threshold = 0.9, max_pcs = 20) {
  variances <- pca_obj$sdev^2
  cum_var_prop <- cumsum(variances) / sum(variances)

  # First component that meets or exceeds the threshold
  candidates <- which(cum_var_prop >= threshold)
  n_pcs <- if (length(candidates) > 0) candidates[1] else length(pca_obj$sdev)

  # Cap by max_pcs and by the number of available PCs
  as.integer(min(n_pcs, max_pcs, ncol(pca_obj$x)))
}

# Return the first n PCs as a features x samples matrix
get_pc_matrix <- function(pca_obj, n_pcs) {
  t(pca_obj$x[, seq_len(n_pcs), drop = FALSE])
}

# Repeated stratified hold out, singleton classes stay in the training set
cv_rf_stratified_splits <- function(mat, labels, repeats = 10, test_frac = 0.2, ntree = 500, seed = 13) {
  stopifnot(ncol(mat) == length(labels))
  labels <- factor(labels)
  n <- ncol(mat)
  classes <- levels(labels)
  tbl <- table(labels)

  singletons <- names(tbl)[tbl == 1]
  if (length(singletons) > 0) {
    log_info("Singleton classes kept only in train: ", paste(singletons, collapse = ", "))
  }

  accs <- numeric(repeats)
  x <- t(mat)

  for (r in seq_len(repeats)) {
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
  c(mean = mean(accs, na.rm = TRUE), sd = sd(accs, na.rm = TRUE))
}

# Stratified k fold cross validation, classes with a single sample are dropped
cv_rf_accuracy_safe <- function(mat, labels, k = 5, ntree = 500, seed = 13) {
  stopifnot(ncol(mat) == length(labels))
  set.seed(seed)
  labels <- factor(labels)
  tbl <- table(labels)

  # Keep only classes with at least 2 samples, they cannot be split otherwise
  keep <- names(tbl)[tbl >= 2]
  removed <- setdiff(names(tbl), keep)
  if (length(removed) > 0) {
    log_info("Removing classes with <2 samples: ", paste(removed, collapse = ", "))
  }

  keep_idx <- which(labels %in% keep)
  if (length(keep_idx) < 2) stop("Not enough samples after removing rare classes.")

  mat2 <- mat[, keep_idx, drop = FALSE]
  labels2 <- droplevels(labels[keep_idx])

  # k cannot be larger than the smallest class
  k_use <- max(2, min(k, min(table(labels2))))

  # Stratified folds, fold ids are assigned per class
  n <- ncol(mat2)
  folds <- integer(n)
  for (cl in levels(labels2)) {
    idx <- sample(which(labels2 == cl))
    folds[idx] <- rep(seq_len(k_use), length.out = length(idx))
  }

  acc <- numeric(k_use)
  x <- t(mat2)
  for (i in seq_len(k_use)) {
    test_idx <- which(folds == i)
    train_idx <- setdiff(seq_len(n), test_idx)
    rf <- randomForest(x[train_idx, , drop = FALSE], labels2[train_idx], ntree = ntree)
    pred <- predict(rf, x[test_idx, , drop = FALSE])
    acc[i] <- mean(pred == labels2[test_idx])
  }
  c(mean = mean(acc), sd = sd(acc))
}

#####################################################################
# Get the data
#####################################################################

# Sample annotation
sannot <- read.csv(sannot.file, stringsAsFactors = FALSE)
sannot$bedFile <- as.character(sannot$bedFile)

# mTFR deviations, one matrix per motif set
dev_mats <- list()
for (nm in names(dev.files)) {
  if (!file.exists(dev.files[[nm]])) {
    log_warn(nm, ": ", dev.files[[nm]], " not found, skipping this motif set")
    next
  }
  log_info("Loading deviations from ", dev.files[[nm]])
  dev_mats[[nm]] <- deviations(readRDS(dev.files[[nm]]))
}
if (length(dev_mats) == 0) {
  stop("None of the deviation files were found, run 02 first")
}

# RnBeads object, disease samples are removed here so that all
# representations are built from the same set of samples
log_info("Loading RnBeads object from ", rnb.set.path)
rnbeads <- load.rnb.set(rnb.set.path)
disease_idx <- which(as.character(pheno(rnbeads)$DISEASE) != "None")
if (length(disease_idx) > 0) {
  log_info("Removing ", length(disease_idx), " disease samples from the RnBeads object")
  rnbeads <- remove.samples(rnbeads, disease_idx)
}

distal <- meth(rnbeads, type = "distal")
tiling <- meth(rnbeads, type = "tiling1kb")

#####################################################################
# Align the samples across all representations
#####################################################################

# Every matrix that goes into the PCA and the Random Forest
mats <- c(dev_mats, list(Distal = distal, Tiling1kb = tiling))

# Healthy samples that carry a cell type we can map, NAs are dropped here
ann_ok <- !is.na(sannot$DISEASE) & sannot$DISEASE == "None" &
  sannot$cellTypeGroup %in% names(group_remap)
healthy <- sannot$bedFile[which(ann_ok)]
common.samples <- Reduce(intersect, c(list(healthy), lapply(mats, colnames)))

if (length(common.samples) == 0) {
  stop("No overlapping samples between the deviations, the RnBeads object and the annotation")
}
log_info("Using ", length(common.samples), " samples shared by all representations")

# Same samples, same order, in every representation
mats <- mapply(function(m, nm) clean_matrix(m[, common.samples, drop = FALSE], nm),
  mats, names(mats),
  SIMPLIFY = FALSE
)

# Cell types in the same order as the matrix columns
cell_types <- unname(group_remap[sannot$cellTypeGroup[match(common.samples, sannot$bedFile)]])

#####################################################################
# PCA and plots
#####################################################################

pca_list <- mapply(run_pca, mats, names(mats), SIMPLIFY = FALSE)

for (nm in names(pca_list)) {
  plot_pca(pca_list[[nm]], cell_types,
    paste("PCA", nm),
    file.path(fig.dir, paste0("PCA_", nm, "_blueprint.pdf"))
  )
}

#####################################################################
# Stability analysis with Random Forests on the PCs
#####################################################################

n_pcs <- vapply(names(pca_list), function(nm) {
  n <- select_pcs_by_variance(pca_list[[nm]], var.threshold, max_pcs = max.pc)
  log_info(nm, ": using ", n, " PCs (variance threshold ", var.threshold, ", max ", max.pc, ")")
  n
}, integer(1))

pcs <- lapply(names(pca_list), function(nm) get_pc_matrix(pca_list[[nm]], n_pcs[[nm]]))
names(pcs) <- names(pca_list)

# Repeated stratified hold out splits
log_info("Random Forest with ", rf.repeats, " stratified splits")
acc_splits <- vapply(pcs, function(m) {
  cv_rf_stratified_splits(m, cell_types,
    repeats = rf.repeats, test_frac = rf.test.frac,
    ntree = rf.ntree, seed = rf.seed
  )
}, numeric(2))

results_splits <- data.frame(
  Representation = names(pcs),
  Accuracy = unname(acc_splits["mean", ]),
  SD = unname(acc_splits["sd", ]),
  NumPCs = as.integer(n_pcs),
  row.names = NULL
)
write.csv(results_splits,
  file = file.path(table.dir, "rf_classification_accuracy_pc_stratified_splits_var90.csv"),
  row.names = FALSE
)

# k fold cross validation
log_info("Random Forest with ", rf.folds, " fold cross validation")
acc_cv <- vapply(pcs, function(m) {
  cv_rf_accuracy_safe(m, cell_types,
    k = rf.folds, ntree = rf.ntree, seed = rf.seed
  )
}, numeric(2))

results_cv <- data.frame(
  Representation = names(pcs),
  Accuracy = unname(acc_cv["mean", ]),
  SD = unname(acc_cv["sd", ]),
  NumPCs = as.integer(n_pcs),
  row.names = NULL
)
write.csv(results_cv,
  file = file.path(table.dir, "rf_classification_accuracy_pc_5foldcv_var90.csv"),
  row.names = FALSE
)

log_success("Finished PCA and stability analysis")
print(results_splits)
print(results_cv)

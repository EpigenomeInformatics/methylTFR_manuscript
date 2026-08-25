#!/usr/bin/env Rscript

#####################################################################
# 03_pca_stability_blueprint.R
# created on 21-07-25 by Irem B Gunduz
# Updated by IBG on 23-08-2026
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

# Cell type groups excluded from the PCA and the stability analysis.
# "Other" collects samples with no assigned cell type, matching 08_bp_mofa.R
drop.cell.types <- "Other"

# Uncorrected control. The deviations assay holds observed minus expected
# methylation, so adding the expected assay back recovers the observed
# methylation before the GC correction. Kept out of the stability analysis
# by default, so the accuracy table stays comparable across runs.
add.uncorrected <- TRUE
uncorrected.motifSet <- "jaspar2020"
uncorrected.in.stability <- FALSE

# PCA is run on the raw matrices, without centring or scaling
pca.center <- FALSE
pca.scale <- FALSE

# Number of PCs kept for the Random Forest.
# On uncentred matrices the first component is the mean methylation profile
# and absorbs almost all of the total sum of squares, so a cumulative
# variance rule always collapses to a single component that carries no cell
# type information. A fixed number of components is used instead.
pc.selection <- "fixed" # "fixed" or "variance"
n.pc <- 20 # Components kept when pc.selection is "fixed"
var.threshold <- 0.9 # Cumulative variance when pc.selection is "variance"
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

# Accuracy tables, named after the PC selection they were produced with
pc.tag <- if (pc.selection == "fixed") paste0(n.pc, "pc") else paste0("var", 100 * var.threshold)
pc.tag_splits <- paste0("rf_classification_accuracy_pc_stratified_splits_", pc.tag, ".csv")
pc.tag_cv <- paste0("rf_classification_accuracy_pc_5foldcv_", pc.tag, ".csv")

# Display names and palette, identical to 06_differential_heatmap.R so that
# every Blueprint panel uses the same colour per cell type
group_remap <- c(
  "Bcell" = "B-cells",
  "DC" = "Dendritic cells",
  "eryt" = "Erythrocytes",
  "gran" = "Granulocytes",
  "megK" = "Megakaryocytes",
  "Mf" = "Macrophages",
  "mono" = "Monocytes",
  "NK" = "NK-cells",
  "osteoclast" = "Osteoclast",
  "other" = "Other",
  "plasma" = "Plasma",
  "progenitor" = "Progenitors",
  "Tcell" = "T-cells",
  "thymocyte" = "Thymocyte"
)

cell_type_colors <- c(
  "B-cells" = "#C2377C",
  "Dendritic cells" = "#8C6D3F",
  "Erythrocytes" = "#7E4B2A",
  "Granulocytes" = "#E8A33D",
  "Macrophages" = "#B5A38A",
  "Megakaryocytes" = "#6A3D9A",
  "Monocytes" = "#C2703D",
  "NK-cells" = "#2CA02C",
  "Osteoclast" = "#8C8C8C",
  "Other" = "#BCBD22",
  "Plasma" = "#7B1E3D",
  "Progenitors" = "#17BECF",
  "T-cells" = "#4FC3D9",
  "Thymocyte" = "#5B9BD5"
)

# Fixed legend order, so the key reads the same in every panel
cell_type_levels <- names(cell_type_colors)

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
  present <- intersect(cell_type_levels, unique(groups))
  plot_data <- data.frame(groups = factor(groups, levels = present))
  pdf(file, width = 10, height = 10)
  print(
    autoplot(pca_obj,
      data = plot_data,
      colour = "groups",
      main = title,
      size = 5
    ) +
      theme_classic() +
      scale_color_manual(values = cell_type_colors[present], name = "Cell type") +
      theme(legend.position = "bottom")
  )
  dev.off()
  log_info("Wrote ", file)
}

# Number of PCs handed to the Random Forest, capped by the number available
select_pcs <- function(pca_obj, label) {
  available <- ncol(pca_obj$x)

  if (pc.selection == "fixed") {
    n_pcs <- min(n.pc, available)
  } else {
    variances <- pca_obj$sdev^2
    cum_var_prop <- cumsum(variances) / sum(variances)
    candidates <- which(cum_var_prop >= var.threshold)
    n_pcs <- if (length(candidates) > 0) candidates[1] else length(pca_obj$sdev)
    n_pcs <- min(n_pcs, max.pc, available)
    if (n_pcs <= 2) {
      log_warn(
        label, ": the variance rule selected ", n_pcs, " PC(s). ",
        "PC1 explains ", round(100 * cum_var_prop[1], 2),
        "% of the total sum of squares, which is expected on uncentred ",
        "matrices. Set pc.selection to \"fixed\" for the stability analysis."
      )
    }
  }

  n_pcs <- as.integer(n_pcs)
  log_info(label, ": using ", n_pcs, " PCs (selection: ", pc.selection, ")")
  n_pcs
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
  dev_obj <- readRDS(dev.files[[nm]])
  dev_mats[[nm]] <- deviations(dev_obj)

  # Observed methylation, before the GC expectation is subtracted
  if (add.uncorrected && nm == paste0("mTFR_", uncorrected.motifSet)) {
    if (!"expected" %in% SummarizedExperiment::assayNames(dev_obj)) {
      log_warn(nm, ": no expected assay, the uncorrected representation is skipped")
    } else {
      dev_mats[[paste0(nm, "_uncorrected")]] <- deviations(dev_obj) +
        SummarizedExperiment::assay(dev_obj, "expected")
      log_info(nm, ": added the uncorrected representation")
    }
  }
  rm(dev_obj)
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

# Healthy samples carrying a mappable cell type, NAs are dropped here
ann_ok <- !is.na(sannot$DISEASE) & sannot$DISEASE == "None" &
  sannot$cellTypeGroup %in% names(group_remap)

# Groups that are not a defined cell type are excluded from both the PCA and
# the classifier, so the accuracy is over labelled cell types only
if (length(drop.cell.types) > 0) {
  mapped <- unname(group_remap[sannot$cellTypeGroup])
  drop_idx <- which(ann_ok & mapped %in% drop.cell.types)
  if (length(drop_idx) > 0) {
    log_info(
      "Dropping ", length(drop_idx), " samples from ",
      paste(sort(unique(mapped[drop_idx])), collapse = ", ")
    )
    ann_ok[drop_idx] <- FALSE
  }
}

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

# The uncorrected control is plotted but, by default, left out of the classifier
stability.reps <- names(pca_list)
if (!uncorrected.in.stability) {
  stability.reps <- grep("_uncorrected$", stability.reps, value = TRUE, invert = TRUE)
}

n_pcs <- vapply(stability.reps, function(nm) select_pcs(pca_list[[nm]], nm), integer(1))

pcs <- lapply(stability.reps, function(nm) get_pc_matrix(pca_list[[nm]], n_pcs[[nm]]))
names(pcs) <- stability.reps

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
  file = file.path(table.dir, pc.tag_splits),
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
  file = file.path(table.dir, pc.tag_cv),
  row.names = FALSE
)

log_success("Finished PCA and stability analysis")
print(results_splits)
print(results_cv)

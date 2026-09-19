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
  library(patchwork)
  library(RnBeads)
  library(methylTFR)
})
set.seed(42)

#####################################################################
# Settings
#####################################################################

# Motif sets whose deviations enter the PCA and the stability analysis
motifSets <- c("jaspar2020", "jaspar2020_distal")

# Cell type groups excluded from the PCA and the stability analysis
drop.cell.types <- c("Other", "Thymocyte")

# Uncorrected control, kept out of the stability analysis by default
add.uncorrected <- TRUE
uncorrected.motifSet <- "jaspar2020"
uncorrected.in.stability <- FALSE

# PCA is run on the raw matrices, without centring or scaling
pca.center <- FALSE
pca.scale <- FALSE

# Number of PCs kept for the Random Forest, fixed rather than variance based
pc.selection <- "fixed" # "fixed" or "variance"
n.pc <- 10 # Components kept when pc.selection is "fixed"
var.threshold <- 0.9 # Cumulative variance when pc.selection is "variance"
max.pc <- 10 # Cap for the number of PCs

# Random Forest settings
rf.ntree <- 500
rf.repeats <- 100
rf.test.frac <- 0.2
rf.folds <- 5
rf.seed <- 42

# Cache the forests, set rf.recompute to force a refit
rf.recompute <- FALSE
# Representation whose motif matrix is used for the feature importance
importance.rep <- "mTFR_jaspar2020_distal"
importance.top <- 20
# Cell types left out of the confusion matrices
confusion.drop <- character(0)

# Directories
analysis.dir <- "/icbb_triton/scratch/igunduz/methylTFR_manuscript/blueprint/"
rnb.tag <- "RnBeads_230826"
dev.tag <- "mTFR_devs_230826"
rnb.set.path <- file.path(analysis.dir, rnb.tag, "reports", "data_import_data", "rnb.set_preprocessed")
sannot.file <- file.path(analysis.dir, rnb.tag, "reports", "data_import_data", "annotation.csv")
dev.files <- setNames(
  file.path(analysis.dir, dev.tag, paste0(motifSets, "_deviations.RDS")),
  paste0("mTFR_", motifSets)
)

github.dir <- "/icbb_triton/scratch/igunduz/methylTFR_manuscript/github/methylTFR_manuscript/"
fig.dir <- file.path(github.dir, "figures", "blueprint")
table.dir <- file.path(github.dir, "tables")
if (!dir.exists(fig.dir)) dir.create(fig.dir, recursive = TRUE)
if (!dir.exists(table.dir)) dir.create(table.dir, recursive = TRUE)

# Accuracy tables, named after the PC selection they were produced with
pc.tag <- if (pc.selection == "fixed") paste0(n.pc, "pc") else paste0("var", 100 * var.threshold)
pc.tag_splits <- paste0("rf_classification_accuracy_pc_stratified_splits_", pc.tag, ".csv")
pc.tag_cv <- paste0("rf_classification_accuracy_pc_5foldcv_", pc.tag, ".csv")

cache.dir <- file.path(analysis.dir, "debug")
if (!dir.exists(cache.dir)) dir.create(cache.dir, recursive = TRUE)
rf.cache.file <- file.path(cache.dir, paste0(
  "rf_stability_", pc.tag, "_", rf.repeats, "splits_", rf.folds, "fold.rds"
))

# Display names and palette, identical to 06_differential_heatmap.R
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

# Confusion matrix from pooled predictions, cells shaded by the row fraction
# Variance explained by the first PCs, on a centred PCA so PC1 is not just
# the mean methylation level that dominates the uncentred classifier PCA
build_scree <- function(mat, title, n.pc = 10) {
  pca <- prcomp(t(clean_matrix(mat, "scree")), center = TRUE, scale. = FALSE)
  pct <- 100 * pca$sdev^2 / sum(pca$sdev^2)
  k <- min(n.pc, length(pct))
  df <- data.frame(pc = factor(seq_len(k)), variance = pct[seq_len(k)])
  ggplot(df, aes(x = pc, y = variance)) +
    geom_col(fill = "#7B1E3D", width = 0.7) +
    labs(title = title, x = "PC", y = "Variance explained (%)") +
    theme_classic(base_size = 9) +
    theme(plot.margin = ggplot2::margin(3, 3, 3, 3))
}

build_confusion <- function(preds, title,
                            drop = if (exists("confusion.drop")) confusion.drop else character(0)) {
  present <- intersect(cell_type_levels, unique(as.character(preds$true)))
  present <- setdiff(present, drop)
  cm <- as.data.frame(table(
    True = factor(preds$true, levels = present),
    Predicted = factor(preds$pred, levels = present)
  ))
  totals <- tapply(cm$Freq, cm$True, sum)
  cm$frac <- cm$Freq / totals[as.character(cm$True)]
  cm$frac[!is.finite(cm$frac)] <- 0

  ggplot(cm, aes(x = Predicted, y = True, fill = frac)) +
    geom_tile(colour = "grey85") +
    geom_text(aes(label = sprintf("%.2f", frac), colour = frac > 0.5), size = 2.4) +
    scale_fill_gradient(low = "white", high = "#2B4B9B", limits = c(0, 1),
      name = "Row fraction") +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey20"), guide = "none") +
    scale_x_discrete(limits = present, expand = c(0, 0)) +
    scale_y_discrete(limits = rev(present), expand = c(0, 0)) +
    labs(title = title, x = "Predicted", y = "True") +
    theme_classic(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank(), plot.margin = ggplot2::margin(3, 3, 3, 3))
}

# Per cell type permutation importance from a forest on the motif matrix
rf_importance <- function(mat, labels, ntree = rf.ntree, seed = rf.seed) {
  x <- t(clean_matrix(mat, "importance"))
  labels <- factor(labels)
  set.seed(seed)
  rf <- randomForest(x, labels, ntree = ntree, importance = TRUE)
  imp <- importance(rf)
  class_cols <- intersect(levels(labels), colnames(imp))
  imp <- imp[, class_cols, drop = FALSE]
  data.frame(
    motif = rep(rownames(imp), times = ncol(imp)),
    celltype = rep(colnames(imp), each = nrow(imp)),
    importance = as.numeric(imp),
    row.names = NULL
  )
}

# Top predictors per cell type, one facet each
build_importance <- function(imp_df, title, top = importance.top) {
  present <- intersect(cell_type_levels, unique(imp_df$celltype))
  top_df <- imp_df %>%
    group_by(celltype) %>%
    slice_max(importance, n = top, with_ties = FALSE) %>%
    ungroup() %>%
    as.data.frame()
  top_df$celltype <- factor(top_df$celltype, levels = present)
  top_df$key <- reorder(interaction(top_df$motif, top_df$celltype, sep = "@@"), top_df$importance)
  ggplot(top_df, aes(x = importance, y = key, fill = celltype)) +
    geom_col(width = 0.7) +
    facet_wrap(~celltype, scales = "free_y") +
    scale_fill_manual(values = cell_type_colors, guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
    scale_y_discrete(labels = function(x) sub("@@.*$", "", x)) +
    labs(title = title, x = "Mean decrease in accuracy", y = NULL) +
    theme_classic(base_size = 9) +
    theme(plot.margin = ggplot2::margin(3, 3, 3, 3))
}

# PCA scatter coloured by cell type, on a centred PCA so the axis variance
# matches the scree panels
build_pca_panel <- function(mat, groups, title) {
  present <- intersect(cell_type_levels, unique(groups))
  df <- data.frame(groups = factor(groups, levels = present))
  pca <- prcomp(t(clean_matrix(mat, "pca")), center = TRUE, scale. = FALSE)
  autoplot(pca, data = df, colour = "groups", size = 1.4) +
    scale_color_manual(values = cell_type_colors[present], name = "Cell type") +
    labs(title = title) +
    theme_classic(base_size = 9) +
    theme(plot.margin = ggplot2::margin(3, 3, 3, 3))
}

# Random Forest classification accuracy across representations
build_accuracy <- function(res, title) {
  res <- res[order(res$Accuracy), ]
  res$Representation <- factor(res$Representation, levels = res$Representation)
  ggplot(res, aes(x = Accuracy, y = Representation)) +
    geom_col(fill = "#4FC3D9", width = 0.7) +
    geom_errorbarh(aes(xmin = Accuracy - SD, xmax = Accuracy + SD), height = 0.25) +
    geom_text(aes(label = sprintf("%.2f", Accuracy)), hjust = -0.25, size = 3) +
    coord_cartesian(xlim = c(0, 1.08)) +
    labs(title = title, x = "Classification accuracy", y = NULL) +
    theme_classic(base_size = 9) +
    theme(plot.margin = ggplot2::margin(3, 3, 3, 3))
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
  preds_list <- vector("list", repeats)
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
    preds_list[[r]] <- data.frame(true = labels[test_idx], pred = pred)
  }
  list(
    mean = mean(accs, na.rm = TRUE), sd = sd(accs, na.rm = TRUE),
    preds = do.call(rbind, preds_list)
  )
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
  preds_list <- vector("list", k_use)
  x <- t(mat2)
  for (i in seq_len(k_use)) {
    test_idx <- which(folds == i)
    train_idx <- setdiff(seq_len(n), test_idx)
    rf <- randomForest(x[train_idx, , drop = FALSE], labels2[train_idx], ntree = ntree)
    pred <- predict(rf, x[test_idx, , drop = FALSE])
    acc[i] <- mean(pred == labels2[test_idx])
    preds_list[[i]] <- data.frame(true = labels2[test_idx], pred = pred)
  }
  list(mean = mean(acc), sd = sd(acc), preds = do.call(rbind, preds_list))
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

# RnBeads object, disease samples removed so all representations share samples
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

# Non cell type groups excluded from the PCA and the classifier
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

# Cache the forests; recompute if missing, forced, or the importance format is stale
rf_cache <- if (!rf.recompute && file.exists(rf.cache.file)) readRDS(rf.cache.file) else NULL
cache_ok <- !is.null(rf_cache) &&
  (is.null(rf_cache$importance) || "celltype" %in% names(rf_cache$importance))
if (cache_ok) {
  log_info("Loading cached Random Forest results from ", rf.cache.file)
} else {
  log_info("Random Forest with ", rf.repeats, " stratified splits")
  split_res <- lapply(pcs, function(m) {
    cv_rf_stratified_splits(m, cell_types,
      repeats = rf.repeats, test_frac = rf.test.frac, ntree = rf.ntree, seed = rf.seed
    )
  })
  log_info("Random Forest with ", rf.folds, " fold cross validation")
  cv_res <- lapply(pcs, function(m) {
    cv_rf_accuracy_safe(m, cell_types, k = rf.folds, ntree = rf.ntree, seed = rf.seed)
  })
  importance_df <- if (importance.rep %in% names(mats)) {
    log_info("Random Forest importance on ", importance.rep)
    rf_importance(mats[[importance.rep]], cell_types)
  } else {
    log_warn(importance.rep, " is not among the representations, skipping importance")
    NULL
  }
  rf_cache <- list(split_res = split_res, cv_res = cv_res, importance = importance_df)
  saveRDS(rf_cache, rf.cache.file)
  log_info("Cached Random Forest results to ", rf.cache.file)
}
split_res <- rf_cache$split_res
cv_res <- rf_cache$cv_res
importance_df <- rf_cache$importance

results_splits <- data.frame(
  Representation = names(split_res),
  Accuracy = vapply(split_res, function(x) x$mean, numeric(1)),
  SD = vapply(split_res, function(x) x$sd, numeric(1)),
  NumPCs = as.integer(n_pcs[names(split_res)]),
  row.names = NULL
)
write.csv(results_splits, file.path(table.dir, pc.tag_splits), row.names = FALSE)

results_cv <- data.frame(
  Representation = names(cv_res),
  Accuracy = vapply(cv_res, function(x) x$mean, numeric(1)),
  SD = vapply(cv_res, function(x) x$sd, numeric(1)),
  NumPCs = as.integer(n_pcs[names(cv_res)]),
  row.names = NULL
)
write.csv(results_cv, file.path(table.dir, pc.tag_cv), row.names = FALSE)

# Supplementary 1: A PCA and B scree share the top row, C confusion, D importance
# Panel A: only the uncorrected and bias corrected JASPAR2020 PCAs
pca.panel.reps <- c(
  "mTFR_jaspar2020_uncorrected" = "Uncorrected (JASPAR2020)",
  "mTFR_jaspar2020" = "Bias corrected (JASPAR2020)"
)
pca.panel.reps <- pca.panel.reps[names(pca.panel.reps) %in% names(pca_list)]
pca_plots <- lapply(names(pca.panel.reps), function(nm)
  build_pca_panel(mats[[nm]], cell_types, pca.panel.reps[[nm]]))
pca_plots[[1]] <- pca_plots[[1]] + labs(tag = "A")
pca_grid <- wrap_plots(pca_plots, ncol = 2, guides = "collect")

# Panel B: scree for every classified representation loaded this run
show.reps <- intersect(names(cv_res), names(pca_list))
scree_plots <- lapply(show.reps, function(nm) build_scree(mats[[nm]], nm))
scree_plots[[1]] <- scree_plots[[1]] + labs(tag = "B")
scree_grid <- wrap_plots(scree_plots, ncol = 2)

conf_plots <- lapply(names(cv_res), function(nm) build_confusion(cv_res[[nm]]$preds, nm))
conf_plots[[1]] <- conf_plots[[1]] + labs(tag = "C")
conf_grid <- wrap_plots(conf_plots, ncol = 2, guides = "collect")

row_ab <- wrap_plots(pca_grid, scree_grid, ncol = 2)

panels <- list(row_ab, conf_grid)
heights <- c(1, 1.1)
if (!is.null(importance_df)) {
  write.csv(importance_df,
    file.path(table.dir, paste0("rf_importance_", importance.rep, ".csv")),
    row.names = FALSE
  )
  imp_panel <- build_importance(
    importance_df, paste0("Top ", importance.top, " methylTFR predictors per cell type")
  ) + labs(tag = "D")
  panels <- c(panels, list(imp_panel))
  heights <- c(1, 1.1, 1.5)
}

combined <- wrap_plots(panels, ncol = 1, heights = heights) &
  theme(plot.margin = ggplot2::margin(3, 3, 3, 3),
    plot.tag = element_text(face = "bold", size = 14))

rf.figure.file <- file.path(fig.dir, paste0("blueprint_supplementary1_", pc.tag, ".pdf"))
ggsave(rf.figure.file, combined, width = 14, height = 20, bg = "white")
log_info("Wrote ", rf.figure.file)

# RF classification accuracy, saved on its own for the main figure
p_acc <- build_accuracy(results_splits, "RF classification accuracy")
acc.figure.file <- file.path(fig.dir, paste0("rf_accuracy_", pc.tag, ".pdf"))
ggsave(acc.figure.file, p_acc, width = 6, height = 4, bg = "white")
log_info("Wrote ", acc.figure.file)

log_success("Finished PCA and stability analysis")


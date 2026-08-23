#!/usr/bin/env Rscript

#####################################################################
# 06_differential_TFs.R
# created on 27-10-2025 by Irem B Gunduz
# Updated by IBG on 23-08-2026
# Differential TFs from the methylTFR deviations of the CD4 memory
# T cell subtypes, the matching TF gene expression heatmap, the mean
# difference comparison of the two contrasts, and the comparison
# against the LOLA motif enrichment of 01
#####################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(grid)
  library(circlize)
  library(viridis)
  library(ComplexHeatmap)
  library(org.Hs.eg.db)
  library(logger)
  library(SummarizedExperiment)
  library(methylTFR)
})
set.seed(13)

#####################################################################
# Settings
#####################################################################

motifSet <- "jaspar2020_distal"

# Contrasts, each memory subtype against naive. The cell type values are
# the ones carried by the cellType column of the sample annotation.
comparisons <- list(
  EM = list(test = "TEM", ref = "TN"),
  CM = list(test = "TCM", ref = "TN")
)

# TEMRA is unreplicated and excluded everywhere
drop.cell.types <- "TEMRA"

top.n <- 50
cut.padj <- 0.05
cut.mean.diff <- 0.1

cell_type_colors <- c(
  "TN" = "#C8E0B4",
  "TCM" = "#4492C6",
  "TEM" = "#43B6C4",
  "TEMRA" = "#898FB5"
)

# Which element of the LOLA result holds each contrast. These are
# positional indices into res_lola$region, they are checked below.
lola.index <- c(EM = 3, CM = 2)
lola.region <- "tiling"
lola.userSets <- c("rankCut_1000_hyper", "rankCut_1000_hypo")

# Directories
analysis.dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/memoryTcells"
dev.tag <- "mTFR_devs_230826"
dev.file <- file.path(analysis.dir, dev.tag, paste0(motifSet, "_deviations.RDS"))
lola.file <- file.path(
  analysis.dir, "reports", "differential_methylation_data",
  "differential_rnbDiffMeth", "TF_motifs_lola.rds"
)

rna.dir <- "/icbb/projects/share/datasets/memoryTcells/rna_raw_081025"
meth.dir <- "/icbb/projects/share/datasets/memoryTcells"

plot.dir <- file.path(analysis.dir, "diff_TFs_230826")
if (!dir.exists(plot.dir)) dir.create(plot.dir, recursive = TRUE)

table.dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR_manuscript/tables"
if (!dir.exists(table.dir)) dir.create(table.dir, recursive = TRUE)

#####################################################################
# Helper functions
#####################################################################

row_zscore <- function(mat) {
  mat <- as.matrix(mat)
  centred <- mat - rowMeans(mat, na.rm = TRUE)
  sds <- apply(mat, 1, sd, na.rm = TRUE)
  out <- centred / sds
  out[!is.finite(out)] <- 0
  out
}

# Donor identifier, Hf03 or Hf04, parsed from a sample name
sample_donors <- function(ids) {
  donors <- rep(NA_character_, length(ids))
  hit <- grepl("Hf[0-9]+", ids)
  donors[hit] <- sub(".*?(Hf[0-9]+).*", "\\1", ids[hit])
  donors
}

# RNA file names carry the subtype as BlCM, BlEM, BlTN, BlTR. They are
# rewritten to the cellType vocabulary so the expression columns and the
# deviation columns can be matched by donor and subtype.
clean_gex_names <- function(x) {
  x <- sub("\\..*$", "", x)
  x <- sub("^[0-9]+_", "", x)
  x <- sub("_Ct$", "", x)
  x <- sub("_BlCM$", "_TCM", x)
  x <- sub("_BlEM$", "_TEM", x)
  x <- sub("_BlTN$", "_TN", x)
  x <- sub("_BlTR$", "_TEMRA", x)
  x
}

#####################################################################
# Load the deviations
#####################################################################

if (!file.exists(dev.file)) {
  stop("Deviations not found, run 02 first: ", dev.file)
}
log_info("Loading deviations from ", dev.file)
dev_obj <- readRDS(dev.file)

cd <- as.data.frame(colData(dev_obj), stringsAsFactors = FALSE)
if (!"cellType" %in% colnames(cd)) {
  stop(
    "The deviations object carries no cellType column. Available: ",
    paste(colnames(cd), collapse = ", ")
  )
}
cell_types <- as.character(cd$cellType)

keep <- which(!is.na(cell_types) & !cell_types %in% drop.cell.types)
dev_obj <- dev_obj[, keep]
cell_types <- cell_types[keep]

# Donor_subtype labels, derived rather than typed out in column order
sample_labels <- paste(sample_donors(colnames(dev_obj)), cell_types, sep = "_")
log_info(ncol(dev_obj), " samples: ", paste(sample_labels, collapse = ", "))

deviations_mat <- deviations(dev_obj)
zscores_all <- deviationZScores(dev_obj)
colnames(zscores_all) <- sample_labels

#####################################################################
# Differential TFs, one contrast at a time
#####################################################################

results <- list()

for (nm in names(comparisons)) {
  test <- comparisons[[nm]]$test
  ref <- comparisons[[nm]]$ref

  idx <- which(cell_types %in% c(test, ref))
  groups <- cell_types[idx]
  if (length(unique(groups)) < 2) {
    stop(nm, ": need both ", test, " and ", ref, " but found ",
      paste(unique(groups), collapse = ", "))
  }
  log_info(nm, ": ", test, " (n = ", sum(groups == test), ") vs ",
    ref, " (n = ", sum(groups == ref), ")")

  # Groups come from the annotation, not from the column order
  diff <- differential_deviation_test(
    deviations_mat[, idx, drop = FALSE],
    groups = groups,
    alternative = "two.sided",
    parametric = TRUE
  )

  saveRDS(diff, file.path(table.dir, paste0("diff_", tolower(nm), "_", motifSet, ".RDS")))
  write.csv(diff,
    file.path(table.dir, paste0("diff_", tolower(nm), "_", motifSet, ".csv")),
    row.names = FALSE
  )

  # Mean Z-score difference, reference minus test, as in the original
  z <- zscores_all[, idx, drop = FALSE]
  zdiff <- rowMeans(z[, groups == ref, drop = FALSE]) -
    rowMeans(z[, groups == test, drop = FALSE])

  diff$zdiff <- zdiff[diff$motifs]
  diff <- diff[order(diff$p_value_adjusted, -diff$mean_difference), ]

  results[[nm]] <- diff
}

# Top motifs of either contrast drive both heatmaps
top_tfs <- unique(unlist(lapply(results, function(d) head(d$motifs, top.n))))
log_info(length(top_tfs), " motifs in the union of the top ", top.n, " of each contrast")

zscores_top <- zscores_all[rownames(zscores_all) %in% top_tfs, , drop = FALSE]
saveRDS(zscores_top, file.path(table.dir, paste0("zscores_diffmotifs_", motifSet, ".RDS")))

#####################################################################
# TF gene expression heatmap
#####################################################################

rna_files <- list.files(rna.dir,
  pattern = "IHECrsem.20190606.hs38.genes.results$", full.names = TRUE
)
meth_files <- list.files(meth.dir, pattern = "bed$", full.names = TRUE)

# Keep the RNA samples that also have methylation data
meth_sample_ids <- gsub(".*?/(\\d+_Hf\\d+_Bl\\w+_Ct)_.*", "\\1", meth_files)
rna_sample_ids <- sub("\\..*", "", basename(rna_files))
matched_rna_files <- rna_files[rna_sample_ids %in% meth_sample_ids]
if (length(matched_rna_files) == 0) {
  stop("No RNA files matched the methylation samples under ", rna.dir)
}
log_info(length(matched_rna_files), " RNA samples matched")

rna_list <- lapply(matched_rna_files, function(file) {
  df <- read.delim(file, stringsAsFactors = FALSE)
  df <- df[, c("gene_id", "expected_count")]
  colnames(df)[2] <- tools::file_path_sans_ext(basename(file))
  df
})
expr_matrix <- Reduce(function(x, y) full_join(x, y, by = "gene_id"), rna_list)

# Map Ensembl IDs to gene symbols
expr_matrix$gene_id_clean <- sub("\\..*$", "", expr_matrix$gene_id)
mapping <- AnnotationDbi::select(org.Hs.eg.db,
  keys = expr_matrix$gene_id_clean,
  keytype = "ENSEMBL", columns = "SYMBOL"
)
colnames(mapping) <- c("gene_id_clean", "gene_short_name")
mapping <- mapping[!is.na(mapping$gene_short_name), ]

expr <- merge(expr_matrix, mapping, by = "gene_id_clean")
rownames(expr) <- make.unique(expr$gene_short_name)
expr <- expr[, !colnames(expr) %in% c("gene_id_clean", "gene_id", "gene_short_name")]
colnames(expr) <- clean_gex_names(colnames(expr))

# Drop the excluded subtypes from the expression side as well
gex_cell_types <- sub("^.*_", "", colnames(expr))
expr <- expr[, !gex_cell_types %in% drop.cell.types, drop = FALSE]
gex_cell_types <- sub("^.*_", "", colnames(expr))

# Only the motifs that are plain gene symbols have an expression row.
# Dimer motifs such as FOSL1::JUND cannot be matched and are reported.
tf_expr <- expr[rownames(expr) %in% rownames(zscores_top), , drop = FALSE]
unmatched <- setdiff(rownames(zscores_top), rownames(tf_expr))
if (length(unmatched) > 0) {
  log_warn(length(unmatched), " motifs have no expression row and are dropped ",
    "from both heatmaps: ", paste(head(unmatched, 10), collapse = ", "),
    if (length(unmatched) > 10) ", ..." else "")
}
if (nrow(tf_expr) < 2) {
  stop("Fewer than two motifs could be matched to gene expression")
}

tf_expr_z <- row_zscore(tf_expr)
gex_types_factor <- factor(gex_cell_types, levels = intersect(
  names(cell_type_colors), unique(gex_cell_types)
))

ha_gex <- HeatmapAnnotation(
  celltypes = gex_cell_types,
  col = list(celltypes = cell_type_colors[levels(gex_types_factor)])
)
hm_gex <- Heatmap(tf_expr_z,
  name = "GEX Z-score",
  column_names_gp = gpar(fontsize = 9),
  top_annotation = ha_gex,
  column_title = "TF Gene Expression (Z-Score)",
  show_row_names = TRUE,
  column_split = gex_types_factor,
  cluster_columns = TRUE,
  cluster_column_slices = TRUE,
  col = colorRamp2(seq(-2, 2, length.out = 70), viridis(70))
)

file <- file.path(plot.dir, "gex_tf_expression_memTcells_mtfr_tfs.pdf")
pdf(file, width = 10, height = 15)
ht_gex <- draw(hm_gex)
dev.off()
log_info("Wrote ", file)

# Row order of the expression heatmap, reused for the methylTFR heatmap
row_order_vector <- rownames(tf_expr_z)[row_order(ht_gex)]

#####################################################################
# methylTFR heatmap, rows ordered by the expression clustering
#####################################################################

zscores_ordered <- zscores_top[row_order_vector, , drop = FALSE]
mtfr_cell_types <- sub("^.*_", "", colnames(zscores_ordered))
mtfr_types_factor <- factor(mtfr_cell_types, levels = intersect(
  names(cell_type_colors), unique(mtfr_cell_types)
))

ha <- HeatmapAnnotation(
  celltypes = mtfr_cell_types,
  col = list(celltypes = cell_type_colors[levels(mtfr_types_factor)])
)
hm <- Heatmap(zscores_ordered,
  name = "methylTFR\nZ-scores",
  column_names_gp = gpar(fontsize = 9),
  top_annotation = ha,
  column_title = "methylTFR Z-scores",
  show_row_names = TRUE,
  column_split = mtfr_types_factor,
  cluster_columns = TRUE,
  cluster_rows = FALSE, # order comes from the expression heatmap
  cluster_column_slices = TRUE
)

file <- file.path(plot.dir, "mtfr_zscores_diffmotifs_memTcells.pdf")
pdf(file, width = 10, height = 15)
draw(hm)
dev.off()
log_info("Wrote ", file)

#####################################################################
# Mean difference of the two contrasts
#####################################################################

comb_df <- results[["EM"]] %>%
  select(motifs, zdiff_em = zdiff,
    mean_difference_em = mean_difference,
    p_value_adjusted_em = p_value_adjusted) %>%
  inner_join(
    results[["CM"]] %>%
      select(motifs, zdiff_cm = zdiff,
        mean_difference_cm = mean_difference,
        p_value_adjusted_cm = p_value_adjusted),
    by = "motifs"
  )

comb_df$isDiff <- case_when(
  abs(comb_df$mean_difference_em) > cut.mean.diff & comb_df$p_value_adjusted_em < cut.padj &
    abs(comb_df$mean_difference_cm) > cut.mean.diff & comb_df$p_value_adjusted_cm < cut.padj ~ "Differential Both",
  abs(comb_df$mean_difference_em) > cut.mean.diff & comb_df$p_value_adjusted_em < cut.padj ~ "Differential EM",
  abs(comb_df$mean_difference_cm) > cut.mean.diff & comb_df$p_value_adjusted_cm < cut.padj ~ "Differential CM",
  TRUE ~ "Not Differential"
)

# Label the strongest motifs of each single contrast. The label is taken
# from the motifs column itself, the original carried a separate TF column
# built from the row names of a differently ordered table.
comb_df$label <- NA_character_
for (grp in c("Differential EM", "Differential CM")) {
  col <- if (grp == "Differential EM") "mean_difference_em" else "mean_difference_cm"
  rows <- which(comb_df$isDiff == grp)
  if (length(rows) == 0) next
  top_idx <- order(abs(comb_df[[col]][rows]), decreasing = TRUE)[seq_len(min(20, length(rows)))]
  comb_df$label[rows[top_idx]] <- comb_df$motifs[rows[top_idx]]
}

cor_val <- cor(comb_df$zdiff_cm, comb_df$zdiff_em, method = "pearson", use = "complete.obs")

p <- ggplot(comb_df, aes(x = zdiff_cm, y = zdiff_em, color = isDiff)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_manual(values = c(
    "Differential CM" = "#377EB8",
    "Differential EM" = "#E41A1C",
    "Differential Both" = "#4DAF4A",
    "Not Differential" = "grey70"
  )) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "black", linewidth = 0.3) +
  labs(
    x = "Mean difference (CM)",
    y = "Mean difference (EM)",
    color = "Differential status"
  ) +
  theme_classic(base_size = 13) +
  theme(legend.position = "top") +
  annotate("text",
    x = Inf, y = Inf, label = paste0("r = ", round(cor_val, 2)),
    hjust = 1.1, vjust = 1.5, size = 4.5, fontface = "bold"
  ) +
  geom_text_repel(
    data = comb_df[!is.na(comb_df$label), ],
    aes(label = label, color = isDiff),
    size = 5, box.padding = 0.3, point.padding = 0.3,
    max.overlaps = Inf, segment.color = NA
  )

file <- file.path(plot.dir, "mean_difference_em_vs_cm_memTcells.pdf")
ggsave(file, p, width = 10, height = 10)
log_info("Wrote ", file)

#####################################################################
# LOLA enrichment against the methylTFR mean difference
#####################################################################

if (!file.exists(lola.file)) {
  log_warn("LOLA results not found, skipping the last figure: ", lola.file)
} else {
  res_lola <- readRDS(lola.file)

  plotlog2OR <- function(df, title) {
    x_max <- max(abs(df$log2OR), na.rm = TRUE)
    y_max <- max(abs(df$zdiff), na.rm = TRUE)

    df <- df %>%
      mutate(isDiff = case_when(
        p_value_adjusted < cut.padj & qValue < cut.padj ~ "Differential in both",
        p_value_adjusted < cut.padj ~ "mTFR differential",
        qValue < cut.padj ~ "LOLA differential",
        TRUE ~ "Not differential"
      ))
    df$isDiff <- factor(df$isDiff, levels = c(
      "Differential in both", "mTFR differential",
      "LOLA differential", "Not differential"
    ))
    df$label <- ifelse(df$isDiff != "Not differential", df$motifs, NA)

    ggplot(df, aes(x = log2OR, y = zdiff, color = isDiff)) +
      geom_point(alpha = 0.7, size = 3) +
      labs(
        title = title,
        x = expression(Log[2] ~ "Odds Ratio (LOLA Enrichment)"),
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
      geom_text_repel(aes(label = label), size = 3.5, box.padding = 0.5, max.overlaps = 15) +
      geom_hline(yintercept = 0, linetype = "dotted", colour = "black") +
      geom_vline(xintercept = 0, linetype = "dotted", colour = "black")
  }

  for (nm in names(comparisons)) {
    i <- lola.index[[nm]]
    if (is.null(res_lola$region[[i]]) || is.null(res_lola$region[[i]][[lola.region]])) {
      log_warn(nm, ": res_lola$region[[", i, "]]$", lola.region,
        " is missing, check lola.index. Available: ",
        paste(names(res_lola$region), collapse = ", "))
      next
    }
    lola <- res_lola$region[[i]][[lola.region]]
    lola <- lola[lola$userSet %in% lola.userSets, ]
    if (nrow(lola) == 0) {
      log_warn(nm, ": no rows for userSets ", paste(lola.userSets, collapse = ", "))
      next
    }

    # Hypermethylated in the memory subtype, hypomethylated in naive
    lola$condition <- ifelse(grepl("hyper$", lola$userSet), comparisons[[nm]]$test, comparisons[[nm]]$ref)
    lola$log2OR <- log2(lola$oddsRatio)
    lola$log2OR <- ifelse(lola$condition == comparisons[[nm]]$ref, -lola$log2OR, lola$log2OR)
    lola$motifs <- gsub(".*_", "", lola$description)
    lola <- lola[, c("motifs", "log2OR", "qValue")]

    merged <- merge(
      results[[nm]][, c("motifs", "zdiff", "p_value_adjusted")],
      lola,
      by = "motifs"
    )
    # One row per motif, the most significant LOLA hit
    merged <- merged[order(merged$motifs, merged$qValue), ]
    merged <- merged[!duplicated(merged$motifs), ]
    if (nrow(merged) == 0) {
      log_warn(nm, ": no motif names overlap between LOLA and methylTFR")
      next
    }

    file <- file.path(plot.dir, paste0(
      "lola_log2OR_vs_mtfr_meanDiff_", tolower(nm), "_vs_tn.pdf"
    ))
    ggsave(file,
      plotlog2OR(merged, paste0("Motif enrichment vs activity difference, ", nm, " vs TN")),
      width = 8, height = 6
    )
    log_info("Wrote ", file)
  }
}

log_success("Finished the differential TF analysis")

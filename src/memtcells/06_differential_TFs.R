#!/usr/bin/env Rscript

#####################################################################
# 06_differential_TFs.R
# created on 27-10-2025 by Irem B Gunduz
# Updated by IBG on 23-08-2026
# Differential TFs of the CD4 memory T cell subtypes:
#   - paired heatmap of the methylTFR Z-scores and the matching VIPER
#     TF activity, annotated with the row-wise correlation
#   - LOLA volcano plots and LOLA against methylTFR
#   - density scatter and MA plot of the RnBeads differential methylation
#####################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(grid)
  library(circlize)
  library(ComplexHeatmap)
  library(org.Hs.eg.db)
  library(logger)
  library(muRtools)
  library(RnBeads)
  library(SummarizedExperiment)
  library(viper)
  library(methylTFR)
})
set.seed(13)

src.dir <- "/icbb_triton/scratch/igunduz/methylTFR_manuscript/github/methylTFR_manuscript/src"
source(file.path(src.dir, "utils.R"))
source(file.path(src.dir, "lola_utils.R"))

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

# TFs shown side by side on the paired heatmap. These are the priority rows,
# kept whenever the motif carries both a deviation and a VIPER activity. Set
# to NULL to fill the heatmap from the contrasts alone.
heatmap.tfs <- c(
  "XBP1", "ZNF449", "IRF9", "RFX3", "NFE2L1", "RBPJ", "HEY1", "RORB",
  "MEF2C", "TFEC", "BATF3", "NFE2", "HEY2", "HES7", "ZFP57", "EBF1",
  "MNT", "FOSL2", "BHLHE40", "BATF", "FOS", "JUN", "RORA", "KLF6",
  "BCL6", "FOSL1", "JUND", "JUNB", "IRF7", "MYC", "GLIS3", "IRF4"
)

# How many rows the heatmap should carry. DoRothEA covers far fewer factors
# than the expression matrix did, so the curated list alone leaves it short.
# Whatever is missing is filled with the motifs the contrasts rank highest.
# Set to NULL to plot only the curated ones.
heatmap.rows <- 40

# VIPER. The right half of the paired heatmap is TF activity inferred from
# the RNA rather than the expression of the TF gene itself, so a factor
# whose own transcript barely moves can still show a shifted regulon.
# Settings match 07b of the Blueprint pipeline.
dorothea.confidence <- c("A", "B", "C")
viper.method <- "scale"
viper.minsize <- 4
viper.eset.filter <- FALSE
viper.cores <- 1

# Significance cuts of the LOLA against methylTFR panel, one per assay
cut.padj.mtfr <- 0.05 # methylTFR significance
cut.qval.lola <- 0.05 # LOLA q-value

# Which methylTFR column that cut applies to. With two donors per group a
# Welch t-test on four values bottoms out around 1e-4, and BH over ~630
# motifs multiplies that by up to 630, so the adjusted p-value can be
# unable to reach 0.05 at all for a given contrast. The script reports the
# smallest attainable value below. Switch to "p_value" to categorise on the
# nominal p-value instead, which must then be described as nominal.
mtfr.signif.col <- "p_value_adjusted"

# Optional effect size floor applied together with the significance cut
cut.effect.mtfr <- 0

# Paired test blocking on donor. Each donor contributes TN, TCM and TEM, so
# a paired test removes the donor effect that the unpaired test leaves in
# the residual. With only two donors it has one degree of freedom, so it is
# more sensitive when the donors agree and less when they do not.
diff.paired <- FALSE

# Motifs labelled on the standalone LOLA volcano. A fixed list rather
# than the top n per direction, so the same motifs are named in the TCM
# and the TEM panel and the two can be read side by side. Set to NULL to
# fall back to the top top.label.lola per direction.
top.label.lola <- 10
lola.label.motifs <- c(
  # AP-1, enriched on the naive side
  "BATF::JUN", "FOS", "FOS::JUN", "FOS::JUNB", "FOS::JUND",
  "FOSB::JUNB", "FOSL1::JUN", "FOSL1::JUNB", "FOSL1::JUND",
  "FOSL2::JUND", "JUND", "JUN(var.2)",
  # ETS, RUNX and TCF, enriched on the memory side
  "ELK4", "ERF", "ERG", "ETS1", "ETV1", "ETV2", "ETV6",
  "EWSR1-FLI1", "FLI1", "GATA1::TAL1", "LEF1", "RUNX2", "RUNX3",
  "SPI1", "SPIB", "TCF7L1", "TCF7L2", "ZBTB7A"
)

# Motifs labelled per category on that panel. Everything differential in
# both assays is labelled by default; lower it to Inf-free numbers if the
# top band is still too crowded to read.
label.n <- c(
  "Differential in both" = 15,
  "mTFR differential" = 15,
  "LOLA differential" = 15,
  "Not differential" = 0
)

# Point label size on the LOLA against methylTFR panel
label.size <- 5.5

# The correlation printed on that panel. Over every motif the two assays
# share, most of which are differential in neither, the number mostly
# measures the noise cloud. Restricting it to the motifs at least one assay
# calls differential asks the question the panel is actually about, and
# applies the same rule to both contrasts.
cor.on.differential <- FALSE

# Differential status colours of that panel
status_colors <- c(
  "Differential in both" = "#E4674E",
  "mTFR differential" = "#2D7A4F",
  "LOLA differential" = "#6BAED6",
  "Not differential" = "#BDBDBD"
)

# Shared green scheme, also defined in utils.R as CELL_TYPE_COLORS
cell_type_colors <- CELL_TYPE_COLORS

# Diverging ramps of the paired heatmap, capped at +/- 2
zscore.limit <- 2
# same methylTFR ramp as the BLUEPRINT (08) and ECHO (04) heatmaps
mtfr_colors <- colorRamp2(
  c(-2, -1, 0, 1, 2),
  c("#2B4B9B", "#4EC3E0", "#FFF200", "#F5A623", "#B21212")
)
viper_colors <- colorRamp2(
  c(-2, -1, 0, 1, 2),
  c("#40004B", "#9970AB", "#FFFFFF", "#5AAE61", "#00441B")
)
cor_colors <- colorRamp2(
  c(-1, -0.5, 0, 0.5, 1),
  c("#8E0152", "#DE77AE", "#FFFFFF", "#7FBC41", "#276419")
)

# The comparison inside the LOLA result is resolved by name, not by
# position: res_lola$region is named after the RnBeads comparisons, e.g.
# "TEM vs. TN (based on cellType)". Set an entry here to override the
# lookup with an exact name or an index.
lola.comparison <- c(EM = NA, CM = NA)

# Region type inside each comparison. 01 registers tiling1kb and distal, so
# the "tiling" of the older runs no longer exists. The first name in this
# vector that is actually present wins.
lola.region <- c("distal", "tiling1kb", "tiling", "sites")
lola.userSets <- c("rankCut_1000_hyper", "rankCut_1000_hypo")

# Region types used for the RnBeads density scatter and MA plots
diffmeth.regions <- c("tiling1kb", "distal", "sites")

# The two RnBeads panels are slow and the supplementary script draws the MA
# plot it needs, so both are off by default. With both FALSE the RnBeads
# differential methylation object is never loaded at all.
run.density.scatter <- FALSE
run.ma.plot <- FALSE

# Differential regions on those two panels come from the RnBeads combined
# rank, not from the adjusted p-value. The density scatter lets RnBeads
# choose the cut per comparison at this alpha.
diffmeth.alpha <- 0.1

# The MA plot uses the same automatic cut, falling back to this fixed one
# if auto.select.rank.cut fails
diffmeth.auto.rank.cut <- TRUE
diffmeth.rank.cut <- 1000

# Directories
analysis.dir <- "/icbb_triton/scratch/igunduz/methylTFR_manuscript/memoryTcells"
dev.tag <- "mTFR_devs_230826"
dev.file <- file.path(analysis.dir, dev.tag, paste0(motifSet, "_deviations.RDS"))
diffmeth.dir <- file.path(
  analysis.dir, "reports", "differential_methylation_data",
  "differential_rnbDiffMeth"
)
lola.file <- file.path(diffmeth.dir, "TF_motifs_lola.rds")
viper.out <- file.path(analysis.dir, "memT_viper_activity.RDS")

rna.dir <- "/icbb/projects/share/datasets/memoryTcells/rna_raw_081025"
meth.dir <- "/icbb/projects/share/datasets/memoryTcells"

# Figures live in the repository, next to the tables they belong with.
# Only the footprints stay under analysis.dir, they are too many for git.
github.dir <- "/icbb_triton/scratch/igunduz/methylTFR_manuscript/github/methylTFR_manuscript"
plot.dir <- file.path(github.dir, "figures", "memtcells", "diff_TFs_230826")
if (!dir.exists(plot.dir)) dir.create(plot.dir, recursive = TRUE)

table.dir <- file.path(github.dir, "tables")
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

# resolve_lola_comparison(), resolve_lola_region() and
# assert_distinct_comparisons() come from lola_utils.R, so 06 and 07 resolve
# contrasts the same way.

# Paired differential test blocking on donor. differential_deviation_test
# is unpaired, which for this design leaves the donor effect in the residual
# even though every donor contributes both subtypes.
paired_deviation_test <- function(mat, groups, donors, test, ref, padjMethod = "BH") {
  d_test <- donors[groups == test]
  d_ref <- donors[groups == ref]
  common <- intersect(d_test, d_ref)
  if (length(common) < 2) {
    stop("Need at least two donors with both ", test, " and ", ref)
  }
  m_test <- mat[, which(groups == test)[match(common, d_test)], drop = FALSE]
  m_ref <- mat[, which(groups == ref)[match(common, d_ref)], drop = FALSE]
  diffs <- m_test - m_ref

  p_val <- apply(diffs, 1, function(x) {
    x <- x[is.finite(x)]
    if (length(x) < 2 || stats::sd(x) == 0) {
      return(NA_real_)
    }
    stats::t.test(x)$p.value
  })
  data.frame(
    motifs = rownames(mat),
    p_value = p_val,
    p_value_adjusted = stats::p.adjust(p_val, method = padjMethod),
    mean_difference = abs(rowMeans(diffs)),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

# Cell type annotation used on top of both halves of the paired heatmap
cell_type_annotation <- function(cell_types) {
  present <- intersect(names(cell_type_colors), unique(cell_types))
  HeatmapAnnotation(
    celltypes = cell_types,
    col = list(celltypes = cell_type_colors[present]),
    annotation_name_gp = gpar(fontsize = 8),
    show_legend = FALSE
  )
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
  diff <- if (diff.paired) {
    log_info(nm, ": paired test blocking on donor")
    paired_deviation_test(
      deviations_mat[, idx, drop = FALSE],
      groups = groups,
      donors = sample_donors(colnames(dev_obj))[idx],
      test = test, ref = ref
    )
  } else {
    differential_deviation_test(
      deviations_mat[, idx, drop = FALSE],
      groups = groups,
      alternative = "two.sided",
      parametric = TRUE
    )
  }

  saveRDS(diff, file.path(table.dir, paste0("diff_", tolower(nm), "_", motifSet, ".RDS")))
  write.csv(diff,
    file.path(table.dir, paste0("diff_", tolower(nm), "_", motifSet, ".csv")),
    row.names = FALSE
  )

  # Mean Z-score difference, test minus reference (memory positive)
  z <- zscores_all[, idx, drop = FALSE]
  zdiff <- rowMeans(z[, groups == test, drop = FALSE]) -
    rowMeans(z[, groups == ref, drop = FALSE])

  diff$zdiff <- zdiff[diff$motifs]
  diff <- diff[order(diff$p_value_adjusted, -diff$mean_difference), ]

  # The smallest adjusted p-value the design can produce. If it exceeds the
  # cut, no motif can be called differential no matter how strong it is.
  min_padj <- min(diff$p_value_adjusted, na.rm = TRUE)
  log_info(nm, ": min raw p = ", signif(min(diff$p_value, na.rm = TRUE), 3),
    ", min adjusted p = ", signif(min_padj, 3), ", ",
    sum(diff$p_value_adjusted < cut.padj.mtfr, na.rm = TRUE),
    " motifs below ", cut.padj.mtfr)
  if (min_padj > cut.padj.mtfr) {
    log_warn(
      nm, ": no motif can reach an adjusted p of ", cut.padj.mtfr,
      " with ", length(idx), " samples and ", nrow(diff), " tests. ",
      "The mTFR side of the panel will be empty. Consider ",
      "mtfr.signif.col = \"p_value\" with an effect floor, or diff.paired = TRUE."
    )
  }

  results[[nm]] <- diff
}

# Top motifs of either contrast drive the paired heatmap
top_tfs <- unique(unlist(lapply(results, function(d) head(d$motifs, top.n))))
log_info(length(top_tfs), " motifs in the union of the top ", top.n, " of each contrast")

zscores_top <- zscores_all[rownames(zscores_all) %in% top_tfs, , drop = FALSE]
saveRDS(zscores_top, file.path(table.dir, paste0("zscores_diffmotifs_", motifSet, ".RDS")))

#####################################################################
# TF gene expression matrix
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

#####################################################################
# VIPER activity from the same expression matrix
#####################################################################

rna <- as.matrix(expr)
keep <- rowSums(rna, na.rm = TRUE) > 0
if (any(!keep)) log_info("Dropping ", sum(!keep), " genes with no counts")
rna <- rna[keep, , drop = FALSE]

lib_size <- colSums(rna, na.rm = TRUE)
if (any(lib_size == 0)) {
  stop("Sample(s) with an empty library: ",
    paste(colnames(rna)[lib_size == 0], collapse = ", "))
}
rna <- log2(t(t(rna) / lib_size) * 1e6 + 1)

# viper scores a gene against its spread across samples, so a gene with no
# variance cannot contribute and makes the scaling undefined
row_sd <- apply(rna, 1, stats::sd, na.rm = TRUE)
rna <- rna[is.finite(row_sd) & row_sd > 0, , drop = FALSE]
log_info("Expression into viper: ", nrow(rna), " genes x ", ncol(rna), " samples")

utils::data("dorothea_hs", package = "dorothea", envir = environment())
dorothea_hs <- get("dorothea_hs", envir = environment())
net <- dorothea_hs[dorothea_hs$confidence %in% dorothea.confidence, ]
regulon <- dorothea::df2regulon(net)
if (length(regulon) == 0) stop("The regulon is empty")

measured <- vapply(regulon, function(r) sum(names(r$tfmode) %in% rownames(rna)), numeric(1))
log_info(
  sum(measured >= viper.minsize), " of ", length(regulon),
  " TFs have at least ", viper.minsize, " measured targets"
)
if (sum(measured >= viper.minsize) == 0) {
  stop(
    "No TF reaches minsize. The counts use '", rownames(rna)[1],
    "', the regulon targets look like '", names(regulon[[1]]$tfmode)[1], "'."
  )
}

activity <- as.matrix(viper(
  eset = rna, regulon = regulon, method = viper.method,
  minsize = viper.minsize, eset.filter = viper.eset.filter,
  cores = viper.cores, verbose = FALSE
))
log_info("VIPER activity: ", nrow(activity), " TFs x ", ncol(activity), " samples")
saveRDS(activity, viper.out)
log_info("Wrote ", viper.out)

# Only the motifs that are plain gene symbols have a VIPER row. Dimer motifs
# such as FOSL1::JUND cannot be matched and are reported. The match is made
# against every motif rather than the top ones, so a curated TF that is not
# among the top of either contrast can still be shown.
tf_expr <- activity[rownames(activity) %in% rownames(zscores_all), , drop = FALSE]
unmatched <- setdiff(rownames(zscores_top), rownames(tf_expr))
if (length(unmatched) > 0) {
  log_warn(length(unmatched), " of the top motifs have no VIPER activity: ",
    paste(head(unmatched, 10), collapse = ", "),
    if (length(unmatched) > 10) ", ..." else "")
}
if (nrow(tf_expr) < 2) {
  stop("Fewer than two motifs could be matched to a VIPER activity")
}
tf_expr_z <- row_zscore(tf_expr)

#####################################################################
# Paired heatmap, methylTFR on the left, VIPER activity on the right
#####################################################################

# Every motif that carries both a deviation and a VIPER activity is eligible
shared_tfs <- intersect(rownames(zscores_all), rownames(tf_expr_z))
log_info(length(shared_tfs), " motifs carry both a deviation and a VIPER activity")

if (!is.null(heatmap.tfs)) {
  wanted <- intersect(heatmap.tfs, shared_tfs)
  absent <- setdiff(heatmap.tfs, shared_tfs)
  if (length(absent) > 0) {
    log_warn(
      length(absent), " of the requested TFs have no deviation or no VIPER ",
      "activity and are not plotted: ", paste(absent, collapse = ", ")
    )
  }
  if (length(wanted) < 2) {
    stop("Fewer than two of heatmap.tfs could be matched on both sides")
  }

  # The curated rows first, then the motifs the contrasts themselves rank
  # highest, so the heatmap reaches its target height without the list
  # having to be extended by hand every time the regulon coverage changes
  if (!is.null(heatmap.rows) && length(wanted) < heatmap.rows) {
    # Every contrast in results is already sorted by adjusted p and then by
    # effect size, so the union of their motif orders is the ranking. An
    # aggregation over the two would return Inf wherever a motif has only
    # missing adjusted p-values, which is how this came back empty before.
    ranked_motifs <- unique(unlist(lapply(results, function(d) as.character(d$motifs))))
    eligible <- setdiff(shared_tfs, wanted)
    fill_pool <- intersect(ranked_motifs, eligible)
    if (length(fill_pool) == 0 && length(eligible) > 0) {
      log_warn(
        "None of the ", length(ranked_motifs), " ranked motifs is eligible, ",
        "filling from the remaining ", length(eligible), " unranked"
      )
      fill_pool <- eligible
    }
    log_info(length(fill_pool), " motifs available to fill the heatmap with")

    add <- head(fill_pool, heatmap.rows - length(wanted))
    if (length(add) > 0) {
      log_info(
        "Topping the heatmap up with ", length(add), " differential motif(s): ",
        paste(add, collapse = ", ")
      )
      wanted <- c(wanted, add)
    }
    if (length(wanted) < heatmap.rows) {
      log_warn(
        "Only ", length(wanted), " of the ", heatmap.rows, " requested rows ",
        "could be filled, ", length(shared_tfs),
        " motifs carry both measurements in total"
      )
    }
  }

  shared_tfs <- wanted
}
log_info("Paired heatmap on ", length(shared_tfs), " TFs")

mtfr_mat <- zscores_all[shared_tfs, , drop = FALSE]
gex_mat <- tf_expr_z[shared_tfs, , drop = FALSE]

mtfr_types <- sub("^.*_", "", colnames(mtfr_mat))
gex_types <- sub("^.*_", "", colnames(gex_mat))
mtfr_split <- factor(mtfr_types, levels = intersect(names(cell_type_colors), unique(mtfr_types)))
gex_split <- factor(gex_types, levels = intersect(names(cell_type_colors), unique(gex_types)))

# Row-wise correlation between the two assays, over the samples that both
# sides share. Donor_subtype labels make that intersection explicit rather
# than assuming the two matrices are in the same column order.
shared_samples <- intersect(colnames(mtfr_mat), colnames(gex_mat))
log_info(length(shared_samples), " samples shared by methylTFR and VIPER: ",
  paste(shared_samples, collapse = ", "))
if (length(shared_samples) < 3) {
  log_warn("Fewer than three shared samples, the row correlations are unstable")
}
row_cor <- vapply(shared_tfs, function(tf) {
  suppressWarnings(cor(
    as.numeric(mtfr_mat[tf, shared_samples]),
    as.numeric(gex_mat[tf, shared_samples]),
    method = "pearson", use = "complete.obs"
  ))
}, numeric(1))
row_cor[!is.finite(row_cor)] <- 0

ht_mtfr <- Heatmap(mtfr_mat,
  name = "methylTFR\nZ-scores",
  col = mtfr_colors,
  top_annotation = cell_type_annotation(mtfr_types),
  column_split = mtfr_split,
  cluster_columns = TRUE,
  cluster_column_slices = FALSE,
  cluster_rows = TRUE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  row_dend_side = "left",
  column_title_gp = gpar(fontsize = 9),
  heatmap_legend_param = list(at = seq(-zscore.limit, zscore.limit, by = 1))
)

ht_gex <- Heatmap(gex_mat,
  name = "VIPER\nZ-scores",
  col = viper_colors,
  top_annotation = cell_type_annotation(gex_types),
  column_split = gex_split,
  cluster_columns = TRUE,
  cluster_column_slices = FALSE,
  cluster_rows = FALSE, # row order comes from the methylTFR half
  show_row_names = FALSE,
  show_column_names = FALSE,
  column_title_gp = gpar(fontsize = 9),
  heatmap_legend_param = list(at = seq(-zscore.limit, zscore.limit, by = 1))
)

# Correlation strip and the TF labels, both to the right of the two halves
ha_cor <- rowAnnotation(
  `Row-Wise Pearson Correlation` = row_cor,
  col = list(`Row-Wise Pearson Correlation` = cor_colors),
  TF = anno_text(shared_tfs, gp = gpar(fontsize = 8)),
  annotation_name_gp = gpar(fontsize = 0),
  width = unit(30, "mm")
)

# The canvas grows with the number of rows, so a taller heatmap does not
# squeeze the labels
file <- file.path(plot.dir, "paired_heatmap_mtfr_viper_memTcells.pdf")
pdf(file, width = 13, height = max(8, 0.28 * length(shared_tfs) + 3))
draw(ht_mtfr + ht_gex + ha_cor,
  merge_legend = TRUE,
  heatmap_legend_side = "right"
)
dev.off()
log_info("Wrote ", file)

#####################################################################
# LOLA: volcano plots, and enrichment against the methylTFR difference
#####################################################################

if (!file.exists(lola.file)) {
  log_warn("LOLA results not found, skipping the LOLA figures: ", lola.file)
} else {
  res_lola <- readRDS(lola.file)

  # Resolve each contrast once, by name, and reuse it for both LOLA figures
  lola_hits <- list()
  for (nm in names(comparisons)) {
    hit <- tryCatch(
      {
        h <- resolve_lola_comparison(res_lola,
          test = comparisons[[nm]]$test, ref = comparisons[[nm]]$ref,
          override = lola.comparison[[nm]]
        )
        h$region <- resolve_lola_region(res_lola, h$index, lola.region)
        h
      },
      error = function(e) {
        log_warn(nm, ": ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(hit)) next
    log_info(nm, " -> \"", hit$name, "\", region ", hit$region,
      " (hyper means hypermethylated in ", hit$grp1, ")")
    lola_hits[[nm]] <- hit
  }
  assert_distinct_comparisons(lola_hits)

  # Volcano per contrast, from lola_utils.R
  for (nm in names(lola_hits)) {
    hit <- lola_hits[[nm]]
    ok <- tryCatch(
      {
        lolaVolcanoPlot(
          lolaRes = res_lola, outputDir = plot.dir,
          comparison = hit$index, region = hit$region,
          label = paste0(nm, "_vs_", comparisons[[nm]]$ref),
          grp1 = hit$grp1, grp2 = hit$grp2,
          n.label = top.label.lola,
          motifs = lola.label.motifs,
          userSets = lola.userSets,
          cell_colors = cell_type_colors
        )
        TRUE
      },
      error = function(e) {
        log_warn(nm, " volcano: ", conditionMessage(e),
          " [call: ", paste(deparse(conditionCall(e)), collapse = " "), "]")
        FALSE
      }
    )
    if (ok) log_info("Wrote the LOLA volcano for ", nm)
  }

  # Enrichment against the methylTFR Z-score difference.
  # grp1 is the memory subtype, grp2 the naive reference, and the header
  # arrow says which side of the x axis belongs to which.
  plotlog2OR <- function(df, grp1, grp2) {
    df <- df %>%
      dplyr::mutate(
        mtfr_sig = .data[[mtfr.signif.col]] < cut.padj.mtfr &
          abs(zdiff) >= cut.effect.mtfr,
        lola_sig = qValue < cut.qval.lola,
        isDiff = dplyr::case_when(
          mtfr_sig & lola_sig ~ "Differential in both",
          mtfr_sig ~ "mTFR differential",
          lola_sig ~ "LOLA differential",
          TRUE ~ "Not differential"
        )
      )
    df$isDiff <- factor(df$isDiff, levels = names(status_colors))
    log_info(
      "  categories: ",
      paste(names(table(df$isDiff)), table(df$isDiff), sep = " = ", collapse = ", ")
    )

    # Labels per category, ranked by the quantity that defines it: the
    # combined effect for the motifs differential in both assays, the mTFR
    # effect or the LOLA odds ratio for the single assay ones.
    df$label <- NA_character_
    for (cat in names(label.n)) {
      n_cat <- label.n[[cat]]
      if (n_cat <= 0) next
      rows <- which(df$isDiff == cat)
      if (length(rows) == 0) next
      rank_by <- switch(cat,
        "Differential in both" = sqrt(df$zdiff[rows]^2 + df$log2OR[rows]^2),
        "mTFR differential" = abs(df$zdiff[rows]),
        "LOLA differential" = abs(df$log2OR[rows]),
        rep(0, length(rows))
      )
      keep <- rows[order(rank_by, decreasing = TRUE)[seq_len(min(n_cat, length(rows)))]]
      df$label[keep] <- df$motifs[keep]
    }
    log_info("  labelling ", sum(!is.na(df$label)), " motifs")

    x_max <- max(abs(df$log2OR), na.rm = TRUE)
    y_max <- max(abs(df$zdiff), na.rm = TRUE)

    cor_set <- if (cor.on.differential) {
      df[df$isDiff != "Not differential", , drop = FALSE]
    } else {
      df
    }
    # Correlation on the signed (-1 multiplied) log2OR against zdiff
    cor_val <- if (nrow(cor_set) < 3) {
      NA_real_
    } else {
      suppressWarnings(cor(cor_set$log2OR, cor_set$zdiff,
        method = "pearson", use = "complete.obs"
      ))
    }
    cor_label <- if (is.na(cor_val)) {
      "Correlation: too few motifs"
    } else if (cor.on.differential) {
      paste0(
        "Correlation over ", nrow(cor_set), " differential motifs: ",
        format(round(cor_val, 2), nsmall = 2)
      )
    } else {
      paste0("Correlation: ", format(round(cor_val, 2), nsmall = 2))
    }
    log_info("  ", cor_label)

    # Room for the labels to escape into. The points saturate near the top
    # of the y range, so the panel is padded and the header sits above it.
    x_lim <- x_max * 1.22
    y_lim <- y_max * 1.20
    y_head <- y_max * 1.34

    # Non differential points first so the coloured ones sit on top
    df <- df[order(df$isDiff, decreasing = TRUE), ]

    ggplot(df, aes(x = log2OR, y = zdiff, color = isDiff)) +
      geom_hline(yintercept = 0, linetype = "dotted", colour = "black") +
      geom_vline(xintercept = 0, linetype = "dotted", colour = "black") +
      geom_point(shape = 16, size = 2.4, alpha = 0.9) +
      geom_text_repel(aes(label = label),
        size = label.size,
        box.padding = 0.7,
        point.padding = 0.3,
        # No leader lines, the labels sit next to their own points
        min.segment.length = Inf,
        force = 6,
        force_pull = 0.4,
        max.overlaps = Inf,
        max.iter = 30000,
        max.time = 3,
        seed = 42,
        show.legend = FALSE
      ) +
      scale_color_manual(values = status_colors, drop = FALSE) +
      # Direction header, drawn above the panel
      annotate("point", x = -x_lim, y = y_head, colour = cell_type_colors[[grp2]], size = 3.5) +
      annotate("text", x = -x_lim * 0.90, y = y_head, label = grp2, hjust = 0, size = 4.5) +
      annotate("segment",
        x = -x_lim * 0.70, xend = x_lim * 0.70, y = y_head, yend = y_head,
        arrow = arrow(ends = "both", length = unit(2, "mm")), linewidth = 0.4
      ) +
      annotate("text", x = x_lim * 0.90, y = y_head, label = grp1, hjust = 1, size = 4.5) +
      annotate("point", x = x_lim, y = y_head, colour = cell_type_colors[[grp1]], size = 3.5) +
      annotate("text",
        x = -x_lim, y = -y_lim, hjust = 0, vjust = 0, size = 4,
        label = cor_label
      ) +
      labs(
        title = paste(grp1, "vs", grp2, "differentials"),
        x = expression(Log[2] ~ "Odds Ratio (LOLA Enrichment)"),
        y = "Z-Score Difference (mTFR Activity)",
        color = "Differential Status"
      ) +
      coord_cartesian(
        xlim = c(-x_lim, x_lim), ylim = c(-y_lim, y_lim), clip = "off"
      ) +
      theme_classic(base_size = 14) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14),
        plot.margin = margin(t = 30, r = 14, b = 6, l = 14),
        legend.position = "bottom"
      ) +
      guides(colour = guide_legend(override.aes = list(size = 3.5, label = "")))
  }

  for (nm in names(lola_hits)) {
    hit <- lola_hits[[nm]]
    lola <- as.data.frame(res_lola$region[[hit$index]][[hit$region]])
    lola <- lola[lola$userSet %in% lola.userSets, ]
    if (nrow(lola) == 0) {
      log_warn(nm, ": no rows for userSets ", paste(lola.userSets, collapse = ", "))
      next
    }

    # hyper refers to group1 of the comparison name, which is not always the
    # memory subtype, so the direction is taken from the name rather than
    # assumed to be the test group
    lola$condition <- ifelse(grepl("hyper$", lola$userSet), hit$grp1, hit$grp2)
    lola$log2OR <- log2(lola$oddsRatio)
    lola$log2OR <- ifelse(lola$condition == comparisons[[nm]]$ref, -lola$log2OR, lola$log2OR)
    # Rows with no description cannot be named, gsub would make them NA
    n_before <- nrow(lola)
    lola <- lola[!is.na(lola$description), ]
    if (nrow(lola) < n_before) {
      log_warn(nm, ": dropped ", n_before - nrow(lola),
        " LOLA rows without a description")
    }
    lola$motifs <- gsub(".*_", "", lola$description)
    lola <- lola[!is.na(lola$motifs) & nzchar(lola$motifs), c("motifs", "log2OR", "qValue")]

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
      plotlog2OR(merged, grp1 = hit$grp1, grp2 = hit$grp2),
      width = 10, height = 8.5
    )
    log_info("Wrote ", file)
  }
}

#####################################################################
# RnBeads differential methylation: density scatter and MA plots
#####################################################################

if (!run.density.scatter && !run.ma.plot) {
  log_info(
    "run.density.scatter and run.ma.plot are FALSE, the RnBeads section is ",
    "skipped and the differential methylation object is not loaded"
  )
} else if (!dir.exists(diffmeth.dir)) {
  log_warn("Differential methylation results not found, skipping: ", diffmeth.dir)
} else {
  if (!run.density.scatter) {
    log_info("run.density.scatter is FALSE, the density scatters are skipped")
  }
  if (!run.ma.plot) {
    log_info("run.ma.plot is FALSE, the MA plots are skipped")
  }

  diffMeth <- load.rnb.diffmeth(diffmeth.dir)
  cmps <- get.comparisons(diffMeth)
  log_info(length(cmps), " comparisons: ", paste(names(cmps), collapse = " | "))

  for (region in diffmeth.regions) {
    # One call per region: the rank cut version returns a named list with
    # one plot per comparison, each with its own automatically chosen cut
    if (run.density.scatter) {
      dens_plots <- tryCatch(
        rnbeadsDensityScatterRankCut(diffMeth, region, alpha = diffmeth.alpha),
        error = function(e) {
          log_warn("density scatter ", region, ": ", conditionMessage(e))
          NULL
        }
      )
      for (nm in names(dens_plots)) {
        tag <- gsub("[^A-Za-z0-9]+", "_", nm)
        file <- file.path(plot.dir, paste0("density_scatter_", region, "_", tag, ".pdf"))
        ggsave(file, dens_plots[[nm]], width = 6, height = 6)
        log_info("Wrote ", file)
      }
    }

    if (!run.ma.plot) next

    for (i in seq_along(cmps)) {
      tag <- gsub("[^A-Za-z0-9]+", "_", names(cmps)[i])
      p_ma <- tryCatch(
        maPlot(diffMeth, region,
          comparison = i,
          rank.cut = diffmeth.rank.cut, auto.rank.cut = diffmeth.auto.rank.cut
        ),
        error = function(e) {
          log_warn("MA plot ", region, " / ", tag, ": ", conditionMessage(e))
          NULL
        }
      )
      if (!is.null(p_ma)) {
        file <- file.path(plot.dir, paste0("ma_plot_", region, "_", tag, ".pdf"))
        ggsave(file, p_ma, width = 7, height = 6)
        log_info("Wrote ", file)
      }
    }
  }
}

log_success("Finished the differential TF analysis")
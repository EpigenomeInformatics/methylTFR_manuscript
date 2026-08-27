#!/usr/bin/env Rscript

#####################################################################
# 07b_viper_prep.R
# created on 27-08-2026 by Irem B Gunduz
# VIPER transcription factor activity for the Blueprint samples that
# 07 paired with WGBS, so that 09 can integrate activity rather than
# raw expression
#####################################################################

suppressPackageStartupMessages({
  library(logger)
  library(viper)
})
set.seed(12)

#####################################################################
# Settings
#####################################################################

# Regulon. "dorothea" builds one from the DoRothEA collection, "file"
# reads a regulon object that was saved earlier, for instance an ARACNe
# network built on this cohort.
regulon.source <- "dorothea"
regulon.file <- NULL

# DoRothEA confidence levels. A to C is the usual choice for bulk data:
# A and B are curated and ChIP supported, C adds inferred interactions.
# D and E are noisier and are left out.
dorothea.confidence <- c("A", "B", "C")

# viper settings. A regulon with fewer targets than minsize is dropped,
# which is what keeps a TF with two measured targets out of the result.
viper.method <- "scale"
viper.minsize <- 20
viper.eset.filter <- FALSE
viper.cores <- 1

# Counts are library size corrected and log transformed before viper,
# which expects a roughly symmetric expression matrix
normalise.counts <- TRUE

# Directories
analysis.dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/blueprint"
mofa.dir <- file.path(analysis.dir, "mofa_230826")
counts.file <- file.path(mofa.dir, "bp_rawRNAcounts.RDS")
viper.out <- file.path(mofa.dir, "bp_viper_activity.RDS")
regulon.out <- file.path(mofa.dir, "bp_viper_regulon.RDS")

github.dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR_manuscript"
table.dir <- file.path(github.dir, "tables")
if (!dir.exists(table.dir)) dir.create(table.dir, recursive = TRUE)

#####################################################################
# Expression matrix from 07
#####################################################################

if (!file.exists(counts.file)) {
  stop("RNA counts not found, run 07 first: ", counts.file)
}
rna <- readRDS(counts.file)
rna <- as.matrix(rna)
log_info("Counts: ", nrow(rna), " genes x ", ncol(rna), " samples")

if (normalise.counts) {
  # Counts per million, then log. Genes that are zero everywhere carry no
  # information and would only widen the row centring below.
  keep <- rowSums(rna, na.rm = TRUE) > 0
  if (any(!keep)) log_info("Dropping ", sum(!keep), " genes with no counts")
  rna <- rna[keep, , drop = FALSE]

  lib_size <- colSums(rna, na.rm = TRUE)
  if (any(lib_size == 0)) stop("Sample(s) with an empty library: ",
    paste(colnames(rna)[lib_size == 0], collapse = ", "))
  rna <- log2(t(t(rna) / lib_size) * 1e6 + 1)
  log_info("Counts normalised to log2 CPM")
}

# viper scores a gene against the spread of that gene across samples, so
# a gene with no variance cannot contribute and makes the scaling undefined
row_sd <- apply(rna, 1, stats::sd, na.rm = TRUE)
if (any(!is.finite(row_sd) | row_sd == 0)) {
  log_info("Dropping ", sum(!is.finite(row_sd) | row_sd == 0), " genes with zero variance")
  rna <- rna[is.finite(row_sd) & row_sd > 0, , drop = FALSE]
}
log_info("Expression matrix into viper: ", nrow(rna), " genes x ", ncol(rna), " samples")

#####################################################################
# Regulon
#####################################################################

if (regulon.source == "file") {
  if (is.null(regulon.file) || !file.exists(regulon.file)) {
    stop("regulon.source is \"file\" but regulon.file is not readable")
  }
  log_info("Reading the regulon from ", regulon.file)
  regulon <- readRDS(regulon.file)
} else {
  if (!requireNamespace("dorothea", quietly = TRUE)) {
    stop(
      "The dorothea package is needed for regulon.source = \"dorothea\". ",
      "Install it with BiocManager::install(\"dorothea\"), or set ",
      "regulon.source to \"file\" and point regulon.file at your own regulon."
    )
  }
  utils::data("dorothea_hs", package = "dorothea", envir = environment())
  dorothea_hs <- get("dorothea_hs", envir = environment())

  net <- dorothea_hs[dorothea_hs$confidence %in% dorothea.confidence, ]
  log_info(
    nrow(net), " interactions across ", length(unique(net$tf)),
    " TFs at confidence ", paste(dorothea.confidence, collapse = "")
  )

  # df2regulon returns the nested tfmode / likelihood structure viper wants
  regulon <- dorothea::df2regulon(net)
}

if (length(regulon) == 0) stop("The regulon is empty")
saveRDS(regulon, regulon.out)
log_info("Regulon: ", length(regulon), " TFs, written to ", regulon.out)

#####################################################################
# VIPER
#####################################################################

# A TF is only scored where enough of its targets are measured, so the
# overlap is reported before the run rather than inferred from the result
measured <- vapply(regulon, function(r) sum(names(r$tfmode) %in% rownames(rna)), numeric(1))
log_info(
  sum(measured >= viper.minsize), " of ", length(regulon),
  " TFs have at least ", viper.minsize, " measured targets"
)
if (sum(measured >= viper.minsize) == 0) {
  stop(
    "No TF reaches minsize. The regulon and the counts are probably keyed ",
    "on different identifiers: the counts use '", rownames(rna)[1],
    "', the regulon targets look like '", names(regulon[[1]]$tfmode)[1], "'."
  )
}

log_info("Running viper, method = ", viper.method)
activity <- viper(
  eset = rna,
  regulon = regulon,
  method = viper.method,
  minsize = viper.minsize,
  eset.filter = viper.eset.filter,
  cores = viper.cores,
  verbose = FALSE
)
activity <- as.matrix(activity)
log_info("VIPER activity: ", nrow(activity), " TFs x ", ncol(activity), " samples")

if (!identical(colnames(activity), colnames(rna))) {
  stop("viper returned a different sample set than it was given")
}

saveRDS(activity, viper.out)
log_success("Wrote ", viper.out)

summary_tab <- data.frame(
  tf = rownames(activity),
  n_targets = unname(measured[rownames(activity)]),
  mean_activity = rowMeans(activity, na.rm = TRUE),
  sd_activity = apply(activity, 1, stats::sd, na.rm = TRUE),
  stringsAsFactors = FALSE
)
summary_tab <- summary_tab[order(-summary_tab$sd_activity), ]
write.csv(summary_tab, file.path(table.dir, "bp_viper_tf_summary.csv"), row.names = FALSE)
log_success("Wrote ", file.path(table.dir, "bp_viper_tf_summary.csv"))

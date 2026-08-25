#!/usr/bin/env Rscript

#####################################################################
# 00b_merge_pseudobulk_devs.R
# created on 25-08-2026 by Irem B Gunduz
# Merge the per cell type pseudobulk deviations into one
# methylTFRdeviations object carrying the sample annotation
#####################################################################

suppressPackageStartupMessages({
  library(logger)
  library(methylTFR)
  library(S4Vectors)
  library(SummarizedExperiment)
})
set.seed(42)

#####################################################################
# Settings
#####################################################################

motifSet <- "jaspar2020_distal"
genome <- "hg38"

pb.dir <- "/icbb/projects/igunduz/irem_github/exposure_atlas_manuscript/data/sample_pseudobulks"

cell.types <- c(
  "Th-Naive", "Th-Mem", "Tc-Naive", "Tc-Mem",
  "B-cell", "NK-cell", "Monocyte"
)

annot.cols <- c("CommonMinID", "condition")

# A motif has no deviation score in a sample where too few of its binding
# sites are covered. Such motifs are dropped so that every column of the
# merged matrix is defined on the same motifs.
drop.nonfinite.motifs <- TRUE

# Directories
analysis.dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/echo"
out.dir <- file.path(analysis.dir, "pseudobulk_devs_230826")
if (!dir.exists(out.dir)) dir.create(out.dir, recursive = TRUE)

github.dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR_manuscript"
table.dir <- file.path(github.dir, "tables")
if (!dir.exists(table.dir)) dir.create(table.dir, recursive = TRUE)

cache.file <- file.path(out.dir, "pseudobulk_devs_list.RDS")
merged.file <- file.path(out.dir, paste0(motifSet, "_deviations.RDS"))
annot.file <- file.path(table.dir, "echo_pseudobulk_sample_annot.tsv")

#####################################################################
# Helper functions
#####################################################################

get_groupname <- function(x) {
  sub("\\.bedGraph.*$", "", vapply(
    strsplit(x, split = "_"), `[`, character(1), 1L
  ))
}

# Row-wise Z-score. Recomputed here rather than carried over, because a
# row-wise Z-score depends on which samples are present and the merged
# object holds a different sample set than any single cell type did.
row_zscore <- function(mat) {
  out <- (mat - rowMeans(mat, na.rm = TRUE)) / apply(mat, 1, stats::sd, na.rm = TRUE)
  out[!is.finite(out)] <- 0
  out
}

#####################################################################
# Load the per cell type objects
#####################################################################

if (file.exists(cache.file)) {
  log_info("Reusing the cached objects from ", cache.file)
  objs <- readRDS(cache.file)
} else {
  log_info("No cache found, reading from ", pb.dir)
  files <- setNames(
    file.path(pb.dir, paste0(cell.types, "_deviations.RDS")),
    cell.types
  )
  missing <- names(files)[!file.exists(files)]
  if (length(missing) > 0) {
    stop("Missing pseudobulk file(s): ", paste(missing, collapse = ", "), ". Run 00a first")
  }
  objs <- lapply(files, readRDS)
}

objs <- objs[intersect(cell.types, names(objs))]
if (length(objs) == 0) stop("None of the requested cell types were loaded")

motifs <- rownames(objs[[1]])
if (!all(vapply(objs, function(x) identical(rownames(x), motifs), logical(1)))) {
  stop("The pseudobulk objects do not share a motif set")
}

#####################################################################
# Merge the deviations
#####################################################################

dev_mat <- do.call(base::cbind, lapply(objs, function(x) as.matrix(deviations(x))))
log_info("Merged ", nrow(dev_mat), " motifs x ", ncol(dev_mat), " samples")

if (anyDuplicated(colnames(dev_mat)) > 0) {
  stop(
    "Duplicated sample names after merging: ",
    paste(unique(colnames(dev_mat)[duplicated(colnames(dev_mat))]), collapse = ", ")
  )
}

#####################################################################
# Sample annotation
#####################################################################

# The annotation travels with the objects, so it is taken from their
# colData rather than rebuilt from a separate table
sample_annot <- do.call(rbind, unname(lapply(objs, function(x) {
  cd <- as.data.frame(colData(x))
  present <- intersect(annot.cols, colnames(cd))
  out <- data.frame(row.names = colnames(x), stringsAsFactors = FALSE)
  for (col in annot.cols) {
    out[[col]] <- if (col %in% present) as.character(cd[[col]]) else NA_character_
  }
  out
})))
rownames(sample_annot) <- colnames(dev_mat)

# Cell type comes from the sample name, which is what the grouping in the
# downstream scripts relies on
sample_annot$cell_type <- factor(get_groupname(colnames(dev_mat)), levels = names(objs))

if (anyNA(sample_annot$cell_type)) {
  stop(
    "Sample name(s) parse to an unlisted cell type: ",
    paste(head(colnames(dev_mat)[is.na(sample_annot$cell_type)], 5), collapse = ", ")
  )
}

for (col in annot.cols) {
  n_missing <- sum(is.na(sample_annot[[col]]))
  if (n_missing > 0) log_warn(col, ": missing for ", n_missing, " sample(s)")
}

log_info("Samples per cell type:")
print(table(sample_annot$cell_type))

#####################################################################
# Motif filter
#####################################################################

if (drop.nonfinite.motifs) {
  finite_motif <- apply(dev_mat, 1, function(r) all(is.finite(r)))
  if (any(!finite_motif)) {
    log_info(
      "Dropping ", sum(!finite_motif), " motif(s) undefined in at least one sample: ",
      paste(head(rownames(dev_mat)[!finite_motif], 10), collapse = ", ")
    )
    dev_mat <- dev_mat[finite_motif, , drop = FALSE]
  }
}
if (nrow(dev_mat) == 0) stop("No motifs left after the filter")

#####################################################################
# Assemble and save
#####################################################################

se <- SummarizedExperiment(
  assays = list(deviations = dev_mat, z = row_zscore(dev_mat)),
  colData = DataFrame(sample_annot),
  rowData = DataFrame(motifs = rownames(dev_mat))
)
merged <- methods::new("methylTFRdeviations", se)

metadata(merged) <- list(
  motifSet = motifSet,
  genome = genome,
  source = "Sample pseudobulks of the single cell immune methylome atlas",
  cellTypes = names(objs)
)

saveRDS(merged, merged.file)
log_success(
  "Wrote ", merged.file, " with ", nrow(merged), " motifs x ", ncol(merged), " samples"
)

annot_out <- data.frame(sample = rownames(sample_annot), sample_annot, stringsAsFactors = FALSE)
write.table(annot_out, annot.file, sep = "\t", quote = FALSE, row.names = FALSE)
log_success("Wrote ", annot.file)

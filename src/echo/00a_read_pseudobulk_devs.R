#!/usr/bin/env Rscript

#####################################################################
# 00a_read_pseudobulk_devs.R
# created on 25-08-2026 by Irem B Gunduz
# Read the per cell type sample pseudobulk deviation objects, check
# that they share one motif set, and cache them for 00b
#####################################################################

suppressPackageStartupMessages({
  library(logger)
  library(methylTFR)
  library(SummarizedExperiment)
})
set.seed(42)

#####################################################################
# Settings
#####################################################################

# Motif set the pseudobulks were scored against
motifSet <- "jaspar2020_distal"

# One object per cell type, columns are donor samples
pb.dir <- "/icbb/projects/igunduz/irem_github/exposure_atlas_manuscript/data/sample_pseudobulks"

# Cell types to read, in the order they should appear downstream.
# Other-cell is excluded, as in the atlas manuscript.
cell.types <- c(
  "Th-Naive", "Th-Mem", "Tc-Naive", "Tc-Mem",
  "B-cell", "NK-cell", "Monocyte"
)

# colData columns carried through to the merged annotation in 00b
annot.cols <- c("CommonMinID", "condition")

# Directories
analysis.dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/echo"
out.dir <- file.path(analysis.dir, "pseudobulk_devs_230826")
if (!dir.exists(out.dir)) dir.create(out.dir, recursive = TRUE)

github.dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR_manuscript"
table.dir <- file.path(github.dir, "tables")
if (!dir.exists(table.dir)) dir.create(table.dir, recursive = TRUE)

cache.file <- file.path(out.dir, "pseudobulk_devs_list.RDS")
inventory.file <- file.path(table.dir, "echo_pseudobulk_inventory.csv")

#####################################################################
# Helper functions
#####################################################################

# Cell type is the token before the first underscore of the sample name,
# so Tc-Mem_OP_S4_Long_D1.bedGraph.bed becomes Tc-Mem. This is the
# convention the atlas manuscript uses.
get_groupname <- function(x) {
  sub("\\.bedGraph.*$", "", vapply(
    strsplit(x, split = "_"), `[`, character(1), 1L
  ))
}

#####################################################################
# Read the objects
#####################################################################

if (!dir.exists(pb.dir)) stop("Pseudobulk directory not found: ", pb.dir)

files <- setNames(
  file.path(pb.dir, paste0(cell.types, "_deviations.RDS")),
  cell.types
)
missing <- names(files)[!file.exists(files)]
if (length(missing) > 0) {
  stop("Missing pseudobulk file(s): ", paste(missing, collapse = ", "))
}

objs <- list()
for (ct in cell.types) {
  log_info("Reading ", basename(files[[ct]]))
  obj <- readRDS(files[[ct]])

  if (!methods::is(obj, "methylTFRdeviations")) {
    stop(ct, ": expected a methylTFRdeviations object, got ", class(obj)[1])
  }

  # The cell type parsed from the sample names is what the downstream
  # grouping relies on, so a disagreement with the colData is reported
  parsed <- get_groupname(colnames(obj))
  cd <- as.data.frame(colData(obj))
  if ("cell_type" %in% colnames(cd)) {
    recorded <- as.character(cd$cell_type)
    if (!identical(parsed, recorded)) {
      log_warn(ct, ": parsed cell types disagree with the cell_type column")
    }
  }
  if (any(parsed != ct)) {
    log_warn(
      ct, ": ", sum(parsed != ct), " column(s) parse to another cell type: ",
      paste(unique(parsed[parsed != ct]), collapse = ", ")
    )
  }

  missing_annot <- setdiff(annot.cols, colnames(cd))
  if (length(missing_annot) > 0) {
    log_warn(
      ct, ": colData is missing ", paste(missing_annot, collapse = ", "),
      ". Present: ", paste(colnames(cd), collapse = ", ")
    )
  }

  objs[[ct]] <- obj
  log_info(ct, ": ", nrow(obj), " motifs x ", ncol(obj), " samples")
}

#####################################################################
# One motif set across all cell types
#####################################################################

motifs <- rownames(objs[[1]])
same_motifs <- vapply(objs, function(x) identical(rownames(x), motifs), logical(1))
if (!all(same_motifs)) {
  stop(
    "The pseudobulk objects do not share a motif set. Differing: ",
    paste(names(same_motifs)[!same_motifs], collapse = ", ")
  )
}
log_info("All cell types share ", length(motifs), " motifs")

#####################################################################
# Inventory and cache
#####################################################################

inventory <- do.call(rbind, lapply(cell.types, function(ct) {
  mat <- deviations(objs[[ct]])
  data.frame(
    cell_type = ct,
    n_samples = ncol(mat),
    n_motifs = nrow(mat),
    n_motifs_nonfinite = sum(!apply(mat, 1, function(r) all(is.finite(r)))),
    stringsAsFactors = FALSE
  )
}))
write.csv(inventory, inventory.file, row.names = FALSE)
log_success("Wrote ", inventory.file)
print(inventory)

saveRDS(objs, cache.file)
log_success(
  "Wrote ", cache.file, " with ", length(objs), " cell types and ",
  sum(inventory$n_samples), " samples"
)

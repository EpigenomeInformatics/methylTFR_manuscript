#!/usr/bin/env Rscript

# 05_scECHO.R  --  Irem B Gunduz, 23-08-2026
# Run methylTFR on the ECHO single-cell allc files, restricted to the
# methylation cells overlapping the scATAC data (overlap defined in
# integration/01_prepare_sampleannot.R -> sample_annot.tsv).
# One deviations object is written per cell type so an interrupted run resumes.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(GenomicRanges)
  library(logger)
  library(methylTFR)
  library(methylTFRAnnotationHg38)
})
set.seed(42)

#####################################################################
# Settings
#####################################################################

# Distal motif set, matched to the distal GC frequency tables
motifSet <- "jaspar2020_distal"
tfSet <- "jaspar2020" # binding sites are shared with the genome-wide set

num.cores <- 30
chunk.size <- 10
cov.threshold <- 1

drop.cell.types <- "Other-cell"
min.cells <- 10

# Inf keeps a cell type in one call (cheapest); lower to bisect a failing type
batch.size <- Inf

# gzipped allc size as a depth proxy; 0 disables the filter
min.file.mb <- 0

# Columns of the sample annotation
cell.id.col <- "Cell_UID"
cell.type.col <- "cell_type"
allc.path.col <- "allC_FilePath_fixed"

# allc files moved from the recorded root to the shared dataset root; paths are
# rehomed by keeping everything after the last allcFiles/ under allc.root
allc.root <- "/icbb/projects/share/datasets/ECHO/allcFiles"
if (!dir.exists(allc.root)) {
  stop("allc directory does not exist: ", allc.root)
}

# Distal regions, Ensembl Regulatory Build v104 filtered to hg38 distal
distal.file <- "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFRAnnotationHg38_old/inst/extdata/distal_regions.RDS"

# Directories
sannot.file <- "/icbb/projects/igunduz/DARPA_analysis/artemis_031023/sample_annot.tsv"
out.dir <- file.path("/scratch/icbb/igunduz/methylTFR_manuscript/echo", paste0("mTFR_sc_150926_", motifSet))
annot.dir <- file.path(out.dir, "annot")
for (d in c(out.dir, annot.dir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# MOFA Motifs
motifs <- c(
  "TFAP4", "OSR2", "ELF5", "SPIC", "SPIB", "SPI1", "EHF", "TFAP2A", "MXI1", "HIC2",
  "TCF4", "ASCL1(var.2)", "MSC", "NHLH1", "EBF3", "EBF1", "POU2F3", "POU2F2", "POU5F1", "POU3F4",
  "OLIG2", "NEUROG1", "ZNF652", "TCF21(var.2)", "PRDM4", "HES1", "HES6", "HEY2", "HEY1", "KLF16",
  "ZBTB14", "E2F2", "FOXN3", "ZFP57", "ESR1", "LHX6", "EMX2", "ELF4", "ETV5", "PAX5",
  "KLF17", "TBXT", "GFI1", "TBX21", "TBX2", "TBR1", "EOMES", "TBX20", "TBX1", "MGA"
)
motifs <- c("POU2F3", "SPIB","EOMES", "TBX21")

#####################################################################
# Sample annotation of the cells overlapping the ATAC data
#####################################################################

if (!file.exists(sannot.file)) {
  stop("Sample annotation does not exist: ", sannot.file)
}
sannot <- fread(sannot.file, data.table = FALSE)
log_info(nrow(sannot), " cells in ", sannot.file)

required <- c(cell.id.col, cell.type.col, allc.path.col)
missing_cols <- setdiff(required, colnames(sannot))
if (length(missing_cols) > 0) {
  stop("Missing column(s) in the annotation: ", paste(missing_cols, collapse = ", "))
}

sannot <- sannot[!sannot[[cell.type.col]] %in% drop.cell.types, ]
log_info(nrow(sannot), " cells after dropping ", paste(drop.cell.types, collapse = ", "))

# Resolve each allc path: recorded first, then rehomed under allc.root; the
# variant that resolves more files wins, counts logged to catch layout changes
recorded <- as.character(sannot[[allc.path.col]])
below <- sub("^.*/allcFiles/", "", recorded)
rehomed <- file.path(allc.root, below)

n_recorded <- sum(file.exists(recorded))
n_rehomed <- sum(file.exists(rehomed))
log_info("allc paths that resolve: ", n_recorded, " as recorded, ",
  n_rehomed, " under ", allc.root)

if (n_recorded == 0 && n_rehomed == 0) {
  stop(
    "No allc file could be located.\n  recorded: ", recorded[1],
    "\n  rehomed:  ", rehomed[1],
    "\nCheck allc.root and the ", allc.path.col, " column."
  )
}
sannot$allC_FilePathfull <- if (n_rehomed >= n_recorded) rehomed else recorded

# Keep only cells whose allc file is on disk
exists <- file.exists(sannot$allC_FilePathfull)
if (any(!exists)) {
  log_warn(sum(!exists), " of ", length(exists),
    " allc files are missing and their cells are skipped, first: ",
    sannot$allC_FilePathfull[which(!exists)[1]])
  sannot <- sannot[exists, ]
}

# Basenames become column names downstream, so they must be unique
dupes <- duplicated(basename(sannot$allC_FilePathfull))
if (any(dupes)) {
  stop(
    sum(dupes), " allc file names are not unique, for example: ",
    basename(sannot$allC_FilePathfull)[which(dupes)[1]]
  )
}

# Depth proxy: a tiny allc file has too few covered CpGs and can kill a batch
sannot$allc_mb <- file.size(sannot$allC_FilePathfull) / 1024^2
log_info("allc file size in MB, quantiles:")
print(round(quantile(sannot$allc_mb, c(0, 0.01, 0.05, 0.25, 0.5, 0.75, 1)), 2))
if (min.file.mb > 0) {
  small_files <- sannot$allc_mb < min.file.mb
  if (any(small_files)) {
    log_warn("Dropping ", sum(small_files), " cells below ", min.file.mb, " MB")
    sannot <- sannot[!small_files, ]
  }
}

cell.types <- sort(unique(as.character(sannot[[cell.type.col]])))
counts <- table(sannot[[cell.type.col]])
small <- names(counts)[counts < min.cells]
if (length(small) > 0) {
  log_warn("Cell types with fewer than ", min.cells, " cells are skipped: ",
    paste(small, collapse = ", "))
  cell.types <- setdiff(cell.types, small)
}
log_info("Running ", length(cell.types), " cell types: ", paste(cell.types, collapse = ", "))

#####################################################################
# Motif annotation
#####################################################################

# Distal set has its own GC frequencies but shares the genome-wide binding sites
gcfreqs <- getGCfreq(motifSet = motifSet)
tf_bindsites <- getTFbindsites(motifSet = tfSet)
gc_dist <- getGenomeGC()
log_info("Number of motifs in gcfreqs: ", length(gcfreqs))

enhancer <- NULL
if (motifSet == "jaspar2020_distal") {
  if (!file.exists(distal.file)) {
    stop("Distal regions file does not exist: ", distal.file)
  }
  enhancer <- readRDS(distal.file)
  log_info("Loaded ", length(enhancer), " distal regions from ", distal.file)
}

# Subset the motif sets into MOFA motifs
gcfreqs <- gcfreqs[names(gcfreqs) %in% motifs]
tf_bindsites <- tf_bindsites[names(tf_bindsites) %in% motifs]

#####################################################################
# Pre-flight: GC bin coverage of each cell under the distal restriction
#####################################################################

# methyltfr_core restricts the GC table to enhancer regions before binning, and
# computeExpectations needs every GC bin present. A sparse single cell can miss
# a bin inside the small distal slice, which fails every motif and kills the
# batch; such cells are dropped here (bulk samples always cover all bins).

gc_binning <- if (is.null(enhancer)) {
  gc_dist
} else {
  log_info("Restricting the GC table to the distal regions")
  subsetByOverlaps(gc_dist, enhancer, ignore.strand = TRUE)
}
bins.needed <- sort(unique(gc_binning$GC_bin))
log_info(
  "GC table for binning: ", length(gc_binning), " ranges, bins ",
  paste(bins.needed, collapse = ", ")
)

qc.file <- file.path(out.dir, "cell_gcbin_qc.tsv")
if (file.exists(qc.file)) {
  log_info("Reusing the cell QC table ", qc.file)
  qc <- read.delim(qc.file, stringsAsFactors = FALSE)
} else {
  log_info("Checking GC bin coverage of ", nrow(sannot), " cells on ", num.cores, " cores")
  qc_list <- parallel::mclapply(seq_len(nrow(sannot)), function(i) {
    ms <- tryCatch(
      read_methylome(sannot$allC_FilePathfull[i],
        type = "allc", cov_threshold = cov.threshold
      ),
      error = function(e) NULL
    )
    if (is.null(ms)) {
      return(data.frame(n_sites = NA_integer_, n_distal_sites = NA_integer_,
        n_bins = NA_integer_))
    }
    hits <- suppressWarnings(
      findOverlaps(ms, gc_binning, ignore.strand = TRUE)
    )
    data.frame(
      n_sites = length(ms),
      n_distal_sites = length(unique(hits@from)),
      n_bins = length(unique(gc_binning$GC_bin[hits@to]))
    )
  }, mc.cores = num.cores)

  ok <- !vapply(qc_list, inherits, logical(1), "try-error")
  qc_list[!ok] <- list(data.frame(
    n_sites = NA_integer_, n_distal_sites = NA_integer_, n_bins = NA_integer_
  ))
  qc <- do.call(rbind, qc_list)
  qc[[cell.id.col]] <- sannot[[cell.id.col]]
  qc[[cell.type.col]] <- sannot[[cell.type.col]]
  qc$allc <- basename(sannot$allC_FilePathfull)
  write.table(qc, qc.file, sep = "\t", quote = FALSE, row.names = FALSE)
  log_info("Wrote ", qc.file)
}

log_info("Sites per cell inside the distal regions, quantiles:")
print(quantile(qc$n_distal_sites, c(0, 0.01, 0.05, 0.25, 0.5, 0.75, 1), na.rm = TRUE))
log_info("GC bins covered per cell:")
print(table(qc$n_bins, useNA = "ifany"))

usable <- !is.na(qc$n_bins) & qc$n_bins == length(bins.needed)
if (!any(usable)) {
  stop(
    "No cell covers all ", length(bins.needed), " GC bins inside the distal ",
    "regions. Single cell coverage is too sparse for the distal motif set. ",
    "Either run motifSet = \"jaspar2020\" without the distal restriction, or ",
    "aggregate the cells into pseudobulks first (see 02_echo_cell_pseudobulks.R)."
  )
}
if (any(!usable)) {
  log_warn(
    "Dropping ", sum(!usable), " of ", nrow(sannot),
    " cells that do not cover all ", length(bins.needed), " GC bins"
  )
}
sannot <- sannot[usable, ]

# Cell types can fall below the threshold once the sparse cells are gone
counts <- table(sannot[[cell.type.col]])
cell.types <- intersect(cell.types, names(counts)[counts >= min.cells])
log_info("After QC: ", nrow(sannot), " cells across ", length(cell.types), " cell types")

#####################################################################
# methylTFR, one cell type at a time
#####################################################################

# Run one batch of cells and return the deviations, or NULL if the batch failed
run_batch <- function(sub_annot, tag) {
  annfile <- file.path(annot.dir, paste0(tag, "_sannot.tsv"))
  write.table(sub_annot, annfile, sep = "\t", quote = FALSE, row.names = FALSE)

  deviations <- tryCatch(
    run_methyltfr(
      annfile = annfile,
      sampleColName = "allC_FilePathfull",
      full_path = TRUE,
      filetype = "allc",
      tf_bindsites = tf_bindsites,
      gcfreqs = gcfreqs,
      gc_dist = gc_dist,
      enhancer = enhancer,
      chunkSize = chunk.size,
      threads = num.cores,
      ignoreStrand = TRUE,
      cov_threshold = cov.threshold
    ),
    error = function(e) {
      # Only log full details at single-cell level, to avoid spamming during bisection
      if (nrow(sub_annot) == 1) {
        log_error("Single cell failed: ", basename(sub_annot$allC_FilePathfull),
                  " - Reason: ", conditionMessage(e))
      }
      NULL
    }
  )
  if (is.null(deviations)) {
    return(NULL)
  }

  # Columns come back named after the allc file, relabel them with the cell id
  if (!identical(colnames(deviations), basename(sub_annot$allC_FilePathfull))) {
    log_error(tag, ": column order does not match the annotation, not relabelling")
    return(NULL)
  }
  colnames(deviations) <- as.character(sub_annot[[cell.id.col]])
  deviations
}

# Recursive bisection: isolate and exclude bad cells without dropping the batch
run_batch_with_fallback <- function(sub_annot, tag) {
  # Attempt to run the whole batch at once (fastest)
  res <- run_batch(sub_annot, tag)
  if (!is.null(res)) return(res) # Success!

  n <- nrow(sub_annot)
  # Base case: failure narrowed to a single cell
  if (n == 1) {
    log_warn("Cell ", sub_annot[[cell.id.col]], " is missing required motif sites and will be excluded.")
    return(NULL)
  }

  # Recursive step: split in half to isolate the bad cells
  log_info("Batch '", tag, "' failed. Bisecting ", n, " cells into two halves to isolate failures...")
  mid <- floor(n / 2)

  # Process left half
  res1 <- run_batch_with_fallback(sub_annot[1:mid, , drop = FALSE], paste0(tag, "_L"))
  # Process right half
  res2 <- run_batch_with_fallback(sub_annot[(mid + 1):n, , drop = FALSE], paste0(tag, "_R"))

  # Combine successful halves (ignoring any isolated cells that returned NULL)
  out <- list()
  if (!is.null(res1)) out[[length(out) + 1]] <- res1
  if (!is.null(res2)) out[[length(out) + 1]] <- res2

  if (length(out) == 0) return(NULL)
  if (length(out) == 1) return(out[[1]])
  return(do.call(cbind, out))
}

for (cell.type in cell.types) {
  out.file <- file.path(out.dir, paste0(make.names(cell.type), "_deviations.RDS"))
  if (file.exists(out.file)) {
    # Motif-aware skip: only skip if the saved file already has every requested motif
    have <- tryCatch(rownames(methylTFR::deviations(readRDS(out.file))),
      error = function(e) NULL)
    if (!is.null(have) && all(motifs %in% have)) {
      log_info("Skipping ", cell.type, ", deviations already cover all requested motifs")
      next
    }
    log_info(cell.type, ": recomputing, missing motif(s): ",
      paste(setdiff(motifs, have), collapse = ", "))
  }

  sub_annot <- sannot[sannot[[cell.type.col]] == cell.type, ]
  n <- nrow(sub_annot)
  size <- if (is.finite(batch.size)) min(batch.size, n) else n
  batches <- split(seq_len(n), ceiling(seq_len(n) / size))
  log_info("Running methylTFR for ", cell.type, " (", n, " cells in ",
    length(batches), " batch(es))")

  devs <- list()
  for (b in seq_along(batches)) {
    tag <- paste0(make.names(cell.type), if (length(batches) > 1) paste0("_batch", b) else "")

    # Use the smart fallback function instead of the standard batch runner
    res <- run_batch_with_fallback(sub_annot[batches[[b]], , drop = FALSE], tag)

    if (!is.null(res)) devs[[length(devs) + 1]] <- res
    gc()
  }

  if (length(devs) == 0) {
    log_warn(cell.type, ": every cell failed, no deviations written")
    next
  }

  deviations <- if (length(devs) == 1) devs[[1]] else do.call(cbind, devs)
  saveRDS(deviations, out.file)
  log_success("Wrote ", out.file, " with ", ncol(deviations), " of ", n, " cells")

  rm(deviations, devs)
  gc()
}

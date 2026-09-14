#!/usr/bin/env Rscript

#####################################################################
# 12_bp_tf_chip_peaks.R
# created on 09-09-2026 by Irem B Gunduz
#
# Fetches and caches the TF ChIP-seq peak sets used to test whether the
# methylTFR deviation score tracks actual TF occupancy. This script only
# downloads and normalises, nothing is analysed here: 13 does the analysis
# and 14 draws the panels.
#
# Every peak set is written to cache.dir as a plain RDS holding a GRanges
# in hg38. The script is resumable: an entry whose RDS already exists is
# skipped, so it can be stopped and restarted on a flaky mount.
#
# Sources, in the order they are preferred:
#   ReMap2022   uniformly reprocessed MACS2 peaks, hg38 native, per dataset,
#               cell type in the file name. CC BY-NC 4.0.
#   UniBind     ChIP-eat TFBS, hg38, per dataset and per JASPAR motif, so the
#               peak set already carries direct motif evidence. CC BY 4.0.
#   ENCODE      IDR thresholded peaks, GRCh38, best QC but no primary immune
#               cells for any of our factors.
#
# See notes/tf_chip_validation.md for why each dataset is in the manifest and
# what is missing. tables/chip_validation/tf_chip_availability.csv is the
# full survey including the gaps.
#####################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(logger)
  library(GenomicRanges)
  library(rtracklayer)
  library(jsonlite)
})
set.seed(42)

#####################################################################
# Settings
#
# All thresholds are fixed here, before anything is looked at. Nothing
# below this block is tuned against the result.
#####################################################################

# A dataset is dropped outright if it has fewer than this many peaks. Small
# peak sets are almost always failed ChIPs and they would dominate the
# bound-fraction axis of analysis B with noise.
min.peaks.per.dataset <- 5000L

# When a factor x cell type has two or more datasets, the peak set used is
# the set of peaks seen in at least this many of them. With one dataset the
# single peak set is used and the fact is recorded.
min.datasets.for.consensus <- 2L

# Peaks are trimmed to a fixed width around their centre so that a deeply
# sequenced dataset with broad peaks does not gain an unfair overlap
# advantage over a shallow one. Set to NA to keep peaks as called.
peak.halfwidth <- 250L

# Standard chromosomes only. Peaks on scaffolds cannot be matched to the
# distal region set anyway.
keep.seqlevels <- paste0("chr", c(1:22, "X", "Y"))

analysis.dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/blueprint/"
cache.dir <- file.path(analysis.dir, "tf_chip_peaks")
if (!dir.exists(cache.dir)) dir.create(cache.dir, recursive = TRUE)

github.dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR_manuscript/"
tab.dir <- file.path(github.dir, "tables", "chip_validation")
if (!dir.exists(tab.dir)) dir.create(tab.dir, recursive = TRUE)
manifest.file <- file.path(tab.dir, "chip_manifest.csv")
summary.out <- file.path(tab.dir, "chip_peak_summary.csv")

# Set to TRUE on a node without outbound network: the script then only
# reports which cache entries are missing and exits.
offline <- FALSE

download.timeout <- 600
options(timeout = download.timeout)

#####################################################################
# Helpers
#####################################################################

# ReMap 2022, one file per dataset. Column 4 is <experiment>.<TF>.<biotype>.
remap_url <- function(acc, tf, biotype) {
  sprintf(
    "https://remap.univ-amu.fr/storage/remap2022/hg38/MACS2/DATASET/%s.%s.%s.bed.gz",
    acc, tf, biotype
  )
}

# ENCODE, the preferred_default IDR thresholded peak file of an experiment.
encode_url <- function(encff) {
  sprintf("https://www.encodeproject.org/files/%s/@@download/%s.bed.gz", encff, encff)
}

# UniBind does not expose a stable per-dataset BED path without knowing which
# JASPAR motif the dataset was called against, so the API is asked. One
# dataset can carry several motifs (a TF with more than one JASPAR profile);
# all of them are taken and merged, which is what "sites of this factor with
# ChIP support" should mean.
unibind_beds <- function(tf_id) {
  u <- sprintf("https://unibind.uio.no/api/v1/datasets/%s/?format=json", tf_id)
  j <- tryCatch(jsonlite::fromJSON(u, simplifyVector = FALSE), error = function(e) NULL)
  if (is.null(j)) {
    log_warn("UniBind API failed for ", tf_id)
    return(character(0))
  }
  damo <- unlist(lapply(j$tfbs, function(x) x$DAMO), recursive = FALSE)
  vapply(damo, function(m) m$bed_url, character(1))
}

read_peaks <- function(urls) {
  grs <- lapply(urls, function(u) {
    dest <- tempfile(fileext = if (grepl("\\.gz$", u)) ".bed.gz" else ".bed")
    ok <- tryCatch(
      {
        utils::download.file(u, dest, mode = "wb", quiet = TRUE)
        TRUE
      },
      error = function(e) {
        log_warn("download failed: ", u, " (", conditionMessage(e), ")")
        FALSE
      }
    )
    if (!ok) return(NULL)
    dt <- tryCatch(data.table::fread(dest, header = FALSE), error = function(e) NULL)
    unlink(dest)
    if (is.null(dt) || nrow(dt) == 0) return(NULL)
    GRanges(
      seqnames = dt[[1]],
      ranges = IRanges(start = dt[[2]] + 1L, end = dt[[3]])
    )
  })
  grs <- Filter(Negate(is.null), grs)
  if (length(grs) == 0) return(NULL)
  grs
}

normalise_peaks <- function(gr) {
  gr <- keepSeqlevels(gr, intersect(seqlevels(gr), keep.seqlevels), pruning.mode = "coarse")
  if (!is.na(peak.halfwidth)) {
    mid <- start(gr) + floor(width(gr) / 2)
    gr <- GRanges(seqnames(gr), IRanges(mid - peak.halfwidth, mid + peak.halfwidth))
  }
  sort(reduce(gr))
}

# Peaks present in at least k of the supplied per-dataset GRanges. With k = 1
# this is the union, which is what a single-dataset entry falls back to.
consensus_peaks <- function(grs, k) {
  k <- min(k, length(grs))
  if (length(grs) == 1L) return(normalise_peaks(grs[[1]]))
  # each dataset is reduced first, so one dataset contributes at most 1 to the
  # depth at any base and the depth is literally "how many datasets cover this"
  grs <- lapply(grs, normalise_peaks)
  pooled <- unlist(GRangesList(grs), use.names = FALSE)
  cov <- GenomicRanges::coverage(pooled)
  keep <- IRanges::slice(cov, lower = k, rangesOnly = TRUE)
  gr <- GRanges(keep)
  sort(reduce(gr))
}

#####################################################################
# Manifest
#
# One row per peak set that analyses A, B and C are allowed to use. Cell
# lines are in here too but they carry match = "proxy" and every figure and
# table must keep that label. urls is a semicolon separated list; several
# urls mean several datasets that are collapsed by consensus_peaks.
#####################################################################

manifest <- data.table(
  set_id = character(), factor = character(), bp_celltype = character(),
  chip_label = character(), match = character(), source = character(),
  urls = character()
)

add <- function(set_id, factor, bp_celltype, chip_label, match, source, urls) {
  manifest <<- rbind(manifest, data.table(
    set_id = set_id, factor = factor, bp_celltype = bp_celltype,
    chip_label = chip_label, match = match, source = source,
    urls = paste(urls, collapse = ";")
  ))
}

## --- primary myeloid, the part of the manifest the argument rests on -----
add("SPI1_monocyte_ReMap", "SPI1", "mono", "primary CD14+ monocyte", "primary", "ReMap2022",
    c(remap_url("GSE31621", "SPI1", "monocyte"),
      remap_url("GSE98367", "SPI1", "monocyte")))
add("CEBPB_monocyte_ReMap", "CEBPB", "mono", "primary CD14+ monocyte", "primary", "ReMap2022",
    c(remap_url("GSE31621", "CEBPB", "monocyte"),
      remap_url("GSE98367", "CEBPB", "monocyte")))
add("SPI1_macrophage_ReMap", "SPI1", "Mf", "monocyte-derived macrophage", "primary", "ReMap2022",
    c(remap_url("GSE31621", "SPI1", "monocyte_MACROPHAGE"),
      remap_url("GSE47188", "SPI1", "macrophage_IFNG"),
      remap_url("GSE47188", "SPI1", "macrophage_IL4")))
add("SPI1_DC_ReMap", "SPI1", "DC", "dendritic cell", "primary", "ReMap2022",
    c(remap_url("GSE123347", "SPI1", "DC"),
      remap_url("GSE58864", "SPI1", "dendrite")))

## --- primary myeloid, GSE128834. Read the caveat in the docs before use --
## This study transfects PU.1 mRNA. Some arms are endogenous PU.1 in the
## stated cell type and some are ectopic. Verify per GSM. It is kept in the
## manifest because it is the ONLY TF ChIP in human neutrophils that exists.
add("SPI1_neutrophil_ReMap", "SPI1", "gran", "primary neutrophil", "primary (verify)", "ReMap2022",
    c(remap_url("GSE128834", "SPI1", "primary-neutrophil_donorE"),
      remap_url("GSE128834", "SPI1", "primary-neutrophil_donorF")))
add("SPI1_primaryB_ReMap", "SPI1", "Bcell", "primary B cell", "primary (verify)", "ReMap2022",
    c(remap_url("GSE128834", "SPI1", "primary-B-cell_donorA"),
      remap_url("GSE128834", "SPI1", "primary-B-cell_donorC")))

## --- primary T cells, the reciprocal arm of analysis C ------------------
add("ETS1_CD4T_UniBind", "ETS1", "Tcell", "primary CD4+ T cell", "primary", "UniBind",
    unlist(lapply(c("EXP032565.CD4_CD25_CD45RA_T-cells.ETS1",
                    "EXP032566.CD4_CD25_CD45RA_T-cells.ETS1",
                    "EXP032567.CD4_CD25_CD45RA_T-cells.ETS1",
                    "EXP032568.CD4_CD25-_T-cells.ETS1",
                    "EXP032569.CD4_CD25-_T-cells.ETS1"), unibind_beds)))
add("RELA_CD4Th1_UniBind", "RELA", "Tcell", "primary CD4+ Th1", "primary", "UniBind",
    unibind_beds("GSE62482.CD4_TH1_DMSO.RELA"))

## --- primary macrophage, UniBind ----------------------------------------
add("CEBPB_MDM_UniBind", "CEBPB", "Mf", "monocyte-derived macrophage", "primary", "UniBind",
    unlist(lapply(c("EXP038588.monocyte-derived_macrophages.CEBPB",
                    "EXP038589.monocyte-derived_macrophages.CEBPB",
                    "EXP040652.monocyte-derived_macrophages.CEBPB"), unibind_beds)))
add("RELA_MDM_UniBind", "RELA", "Mf", "monocyte-derived macrophage", "primary", "UniBind",
    unlist(lapply(c("EXP036831.monocyte-derived_macrophages.RELA",
                    "EXP036832.monocyte-derived_macrophages.RELA"), unibind_beds)))

## --- GM12878, the B cell PROXY. Labelled as such everywhere.
## K562 is not here: it is a chronic myeloid leukaemia line with an erythroid
## phenotype, and cancer lines carry globally distorted methylomes ----------
gm <- list(SPI1 = "ENCFF492ZRZ", EBF1 = "ENCFF895MHN", POU2F2 = "ENCFF934JFA",
           ELF1 = "ENCFF146SYU", ETS1 = "ENCFF568AZT", GABPA = "ENCFF093KLR",
           BATF = "ENCFF728KFD", CEBPB = "ENCFF955YFB", JUNB = "ENCFF912OPT",
           JUND = "ENCFF134BQO", FOS = "ENCFF571DGT", RELA = "ENCFF141SMI")
for (tf in names(gm)) {
  add(paste0(tf, "_GM12878_ENCODE"), tf, "Bcell", "GM12878", "proxy", "ENCODE",
      encode_url(gm[[tf]]))
}

fwrite(manifest, manifest.file)
log_info("Manifest: ", nrow(manifest), " peak sets across ",
         uniqueN(manifest$factor), " factors")

#####################################################################
# Fetch
#####################################################################

summ <- rbindlist(lapply(seq_len(nrow(manifest)), function(i) {
  m <- manifest[i]
  out.file <- file.path(cache.dir, paste0(m$set_id, ".rds"))

  if (file.exists(out.file)) {
    gr <- readRDS(out.file)
    log_info("cached: ", m$set_id, " (", length(gr), " peaks)")
    return(data.table(m, n_datasets = NA_integer_, n_peaks = length(gr), status = "cached"))
  }
  if (offline) {
    log_warn("missing and offline: ", m$set_id)
    return(data.table(m, n_datasets = NA_integer_, n_peaks = NA_integer_, status = "missing"))
  }

  urls <- strsplit(m$urls, ";", fixed = TRUE)[[1]]
  urls <- urls[nzchar(urls)]
  if (length(urls) == 0) {
    log_warn("no urls resolved for ", m$set_id)
    return(data.table(m, n_datasets = 0L, n_peaks = NA_integer_, status = "no_urls"))
  }

  grs <- read_peaks(urls)
  if (is.null(grs)) {
    return(data.table(m, n_datasets = 0L, n_peaks = NA_integer_, status = "download_failed"))
  }

  # the minimum peak rule is applied per dataset, before any merging
  sizes <- vapply(grs, length, integer(1))
  drop <- sizes < min.peaks.per.dataset
  if (any(drop)) {
    log_warn(m$set_id, ": dropping ", sum(drop), " dataset(s) below ",
             min.peaks.per.dataset, " peaks")
    grs <- grs[!drop]
  }
  if (length(grs) == 0) {
    return(data.table(m, n_datasets = 0L, n_peaks = NA_integer_, status = "all_below_min_peaks"))
  }

  gr <- consensus_peaks(grs, min.datasets.for.consensus)
  saveRDS(gr, out.file)
  log_success(m$set_id, ": ", length(gr), " peaks from ", length(grs), " dataset(s)")
  data.table(m, n_datasets = length(grs), n_peaks = length(gr),
             status = if (length(grs) >= min.datasets.for.consensus) "consensus" else "single_dataset")
}))

summ[, urls := NULL]
fwrite(summ, summary.out)
log_success("Wrote ", summary.out)
print(summ[, .(set_id, factor, bp_celltype, match, n_datasets, n_peaks, status)])

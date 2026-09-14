#!/usr/bin/env Rscript

#####################################################################
# 13_bp_tf_chip_validation.R
# created on 10-09-2026 by Irem B Gunduz
#
# Does the methylTFR deviation score track transcription factor occupancy?
#
# For each motif, the distal binding sites are split into those with ChIP-seq
# support in a given cell type and those without, the unsupported ones are
# matched on CpG count and GC content, and the deviation is computed for each
# group. Every peak set is scored against a WGBS methylome of its own cell
# type, so the comparison is within one cell type rather than across two.
#
# The control is the mismatched methylome: the identical site sets rescored
# against another lineage. If the split is driven by open chromatin rather
# than occupancy the two give the same answer. Set n.permutations above zero
# for a within-pair label permutation as well.
#####################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(logger)
  library(GenomicRanges)
  library(Biostrings)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(methylTFR)
  library(methylTFRAnnotationHg38)
  library(methylTFRAnnotationBuilder)
})
set.seed(42)

#####################################################################
# Settings
#####################################################################

motifSet <- "jaspar2020_distal"
tfSet <- "jaspar2020"

overlap.slop <- 0L
min.sites.per.motif <- 500L
min.sites.per.group <- 100L

match.window <- 250L
match.gc.tolerance <- 0.02
match.ratio <- 1L

cov.threshold <- 5L
n.permutations <- 0L

# Analysis B, restricted to motifs with ChIP in the same cell type
min.factors.for.B <- 6L
b.peak.budget <- 20000L

analysis.dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/blueprint/"
peaks.dir <- file.path(analysis.dir, "tf_chip_peaks")
meth.dir <- file.path(analysis.dir, "tf_chip_methylomes")
cache.dir <- file.path(analysis.dir, "tf_chip_validation")
for (d in c(meth.dir, cache.dir)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

github.dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR_manuscript/"
out.dir <- file.path(github.dir, "tables", "chip_validation")
if (!dir.exists(out.dir)) dir.create(out.dir, recursive = TRUE)
manifest.file <- file.path(out.dir, "chip_manifest.csv")

distal.file <- "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFRAnnotationHg38_old/inst/extdata/distal_regions.RDS"

options(timeout = 3600)

#####################################################################
# Paired methylomes
#
# ENCODE WGBS, GRCh38, bedMethyl, the preferred_default file of each
# experiment. There is no human neutrophil, macrophage or dendritic cell WGBS
# in ENCODE, so those peak sets have no partner and are skipped.
#
# The three primary methylomes are all donor ENCDO661BYS, so any contrast
# between them is donor controlled.
#####################################################################

methylomes <- data.table(
  mkey = c("GM12878", "monocyte", "Bcell", "Tcell"),
  label = c("GM12878", "CD14+ monocyte", "B cell", "T cell"),
  experiment = c("ENCSR890UQO", "ENCSR017BUL", "ENCSR284TCU", "ENCSR663MXB"),
  file = c("ENCFF570TIL", "ENCFF689TNG", "ENCFF703XLD", "ENCFF953DKC")
)
methylomes[, url := sprintf(
  "https://www.encodeproject.org/files/%s/@@download/%s.bed.gz", file, file)]

# Which methylome each peak set belongs with, and which one is its negative
# control. The control is always a different lineage and always a normal
# methylome: the lymphoid sets are answered by monocyte, which is the largest
# identity gap available among the four.
pairing <- data.table(
  set_id = c(
    "SPI1_monocyte_ReMap", "CEBPB_monocyte_ReMap", "SPI1_primaryB_ReMap",
    "ETS1_CD4T_UniBind", "RELA_CD4Th1_UniBind",
    paste0(c("SPI1", "EBF1", "POU2F2", "ELF1", "ETS1", "GABPA", "BATF",
             "CEBPB", "JUNB", "JUND", "FOS", "RELA"), "_GM12878_ENCODE")
  )
)
pairing[, methylome := fifelse(grepl("_GM12878_", set_id), "GM12878",
                fifelse(grepl("monocyte", set_id), "monocyte",
                fifelse(grepl("primaryB", set_id), "Bcell", "Tcell")))]
pairing[, control := fifelse(methylome == "monocyte", "Tcell", "monocyte")]

#####################################################################
# Helpers
#####################################################################

motif_to_tfs <- function(motif) {
  parts <- unlist(strsplit(motif, "::", fixed = TRUE))
  parts <- sub("\\s*\\(var\\.[0-9]+\\)$", "", parts)
  toupper(trimws(parts))
}

tf_to_motifs <- function(tf, motifs) {
  motifs[vapply(motifs, function(m) tf %in% motif_to_tfs(m), logical(1))]
}

subset_gcfreq <- function(motif, sites) {
  setNames(list(processMotifs2Matrix(
    motif = motif, gc_bin = gc.bin, genome = BSgenome.Hsapiens.UCSC.hg38,
    tf_bindsites = setNames(list(sites), motif), enhancer = NULL
  )), motif)
}

score_sites <- function(motif, sites, key, gcf) {
  if (length(sites) < 20L) return(NULL)
  r <- tryCatch(
    computeDeviation(
      motif = motif, msites = msites[[key]],
      tf_bindsites = setNames(list(sites), motif), gcfreqs = gcf,
      enhancer = enhancer, ignoreStrand = TRUE, binMsites = bin_msites[[key]]
    ),
    error = function(e) {
      log_warn(motif, " / ", key, ": ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(r) || !is.finite(r$dev)) return(NULL)
  data.table(current = r$dev, expected_only = r$exp_dev,
             observed_only = r$dev + r$exp_dev, n_sites = length(sites))
}

site_composition <- function(sites, halfwidth = match.window) {
  mid <- start(sites) + floor(width(sites) / 2)
  sq <- getSeq(BSgenome.Hsapiens.UCSC.hg38,
               trim(GRanges(seqnames(sites), IRanges(mid - halfwidth, mid + halfwidth))))
  data.table(idx = seq_along(sites),
             cpg = as.integer(vcountPattern("CG", sq)),
             gc = as.numeric(letterFrequency(sq, "GC", as.prob = TRUE)))
}

match_background <- function(comp_sup, comp_unsup, ratio = 1L) {
  sup_keep <- integer(0); unsup_pick <- integer(0)
  for (k in intersect(unique(comp_sup$cpg), unique(comp_unsup$cpg))) {
    s <- comp_sup[cpg == k][order(gc)]
    u <- comp_unsup[cpg == k][order(gc)]
    if (nrow(u) == 0L) next
    used <- rep(FALSE, nrow(u))
    for (i in seq_len(nrow(s))) {
      cand <- which(!used & abs(u$gc - s$gc[i]) <= match.gc.tolerance)
      if (length(cand) < ratio) next
      take <- cand[order(abs(u$gc[cand] - s$gc[i]))][seq_len(ratio)]
      used[take] <- TRUE
      sup_keep <- c(sup_keep, s$idx[i])
      unsup_pick <- c(unsup_pick, u$idx[take])
    }
  }
  list(supported = sup_keep, matched = unsup_pick)
}

#####################################################################
# Inputs
#####################################################################

tf_bindsites <- getTFbindsites(motifSet = tfSet)
gcfreqs <- getGCfreq(motifSet = motifSet)
gc_dist <- getGenomeGC()
enhancer <- readRDS(distal.file)

gc.bin <- S4Vectors::metadata(gc_dist)$gc_breaks
if (is.list(gc.bin)) gc.bin <- gc.bin[[1]]
if (is.null(gc.bin)) gc.bin <- gcBreaks(gc_dist$GC_bias)

manifest <- fread(manifest.file)
manifest <- merge(manifest, pairing, by = "set_id")
manifest[, peak_file := file.path(peaks.dir, paste0(set_id, ".rds"))]
manifest <- manifest[file.exists(peak_file)]
log_info(nrow(manifest), " peak sets with a paired methylome")

#####################################################################
# Fetch and read the methylomes
#####################################################################

read_encode_bedmethyl <- function(filename, cov_threshold = 1) {
  dt <- fread(filename, header = FALSE, skip = 1, showProgress = FALSE)
  if (ncol(dt) < 11) stop(filename, " has ", ncol(dt), " columns, expected 11")
  gr <- GRanges(
    seqnames = dt$V1,
    ranges = IRanges(start = dt$V2 + 1L, end = dt$V3),
    strand = dt$V6,
    score = dt$V11 / 100,
    coverage = dt$V10
  )
  gr <- gr[!is.na(gr$score) & gr$coverage >= cov_threshold]
  sort(gr)
}

msites <- list()
for (i in seq_len(nrow(methylomes))) {
  m <- methylomes[i]
  if (!m$mkey %in% manifest$methylome && !m$mkey %in% manifest$control) next
  rds <- file.path(meth.dir, paste0(m$mkey, "_cov", cov.threshold, ".rds"))
  if (file.exists(rds)) {
    msites[[m$mkey]] <- readRDS(rds)
    log_info("cached methylome: ", m$mkey, " (", length(msites[[m$mkey]]), " sites)")
    next
  }
  bed <- file.path(meth.dir, paste0(m$file, ".bed.gz"))
  if (!file.exists(bed)) {
    log_info("downloading ", m$label, " (", m$file, ") ...")
    utils::download.file(m$url, bed, mode = "wb", quiet = TRUE)
  }
  gr <- read_encode_bedmethyl(bed, cov_threshold = cov.threshold)
  rng <- range(gr$score, na.rm = TRUE)
  if (rng[2] > 1.001 || rng[1] < 0) {
    stop("Methylation scores for ", m$mkey, " span ", signif(rng[1], 3), " to ",
         signif(rng[2], 3), ", which is not a fraction.")
  }
  saveRDS(gr, rds)
  msites[[m$mkey]] <- gr
  log_success(m$label, ": ", length(gr), " CpGs at coverage >= ", cov.threshold)
}

bin_msites <- lapply(msites, addGCBintoMethylome, gcdist = gc_dist)

#####################################################################
# Preflight
#####################################################################

pf.key <- manifest$methylome[1]
pf.motif <- names(which.max(vapply(tf_bindsites, length, integer(1))))
pf.sites <- subsetByOverlaps(tf_bindsites[[pf.motif]], enhancer)
pf.comp <- site_composition(pf.sites)
pf.n <- floor(length(pf.sites) / 4)
pf.hi <- pf.sites[pf.comp[order(-gc), idx][seq_len(pf.n)]]
pf.lo <- pf.sites[pf.comp[order(gc), idx][seq_len(pf.n)]]

pf <- rbindlist(list(
  cbind(split = "gc_high_motif_table", score_sites(pf.motif, pf.hi, pf.key, gcfreqs[pf.motif])),
  cbind(split = "gc_low_motif_table", score_sites(pf.motif, pf.lo, pf.key, gcfreqs[pf.motif])),
  cbind(split = "gc_high_subset_table", score_sites(pf.motif, pf.hi, pf.key, subset_gcfreq(pf.motif, pf.hi))),
  cbind(split = "gc_low_subset_table", score_sites(pf.motif, pf.lo, pf.key, subset_gcfreq(pf.motif, pf.lo)))
), fill = TRUE)
pf[, `:=`(motif = pf.motif, methylome = pf.key)]
fwrite(pf, file.path(out.dir, "preflight_expected_model.csv"))
print(pf)

gap <- function(a, b, col) abs(pf[split == a][[col]] - pf[split == b][[col]])
exp.gap <- gap("gc_high_subset_table", "gc_low_subset_table", "expected_only")
obs.gap <- gap("gc_high_subset_table", "gc_low_subset_table", "observed_only")
log_info("Expected gap across GC quartiles: ", signif(exp.gap, 3),
         " with subset tables, ",
         signif(gap("gc_high_motif_table", "gc_low_motif_table", "expected_only"), 3),
         " with the per motif table; observed gap ", signif(obs.gap, 3))

if (!is.finite(exp.gap) || exp.gap < 0.1 * obs.gap) {
  stop("Preflight failed: the expected term moves ", signif(exp.gap, 3),
       " while the observed term moves ", signif(obs.gap, 3),
       ", so the correction is inert at site level.")
}
log_success("Preflight passed: expected gap is ",
            signif(100 * exp.gap / obs.gap, 3), "% of the observed gap")

#####################################################################
# Site level split, matched methylome against mismatched
#####################################################################

run_split <- function(row) {
  cache <- file.path(cache.dir, paste0(row$set_id, ".rds"))
  if (file.exists(cache)) return(readRDS(cache))
  if (!row$methylome %in% names(msites)) return(NULL)

  peaks <- readRDS(row$peak_file)
  motifs <- tf_to_motifs(row$factor, names(tf_bindsites))

  res <- rbindlist(lapply(motifs, function(motif) {
    sites <- subsetByOverlaps(tf_bindsites[[motif]], enhancer)
    if (length(sites) < min.sites.per.motif) return(NULL)
    mid <- start(sites) + floor(width(sites) / 2)
    is_sup <- overlapsAny(
      GRanges(seqnames(sites), IRanges(mid - overlap.slop, mid + overlap.slop)), peaks)
    if (sum(is_sup) < min.sites.per.group || sum(!is_sup) < min.sites.per.group) return(NULL)

    comp <- site_composition(sites)
    mm <- match_background(comp[is_sup], comp[!is_sup], match.ratio)
    if (length(mm$supported) < min.sites.per.group) return(NULL)

    sup <- sites[mm$supported]; mat <- sites[mm$matched]
    gcf_sup <- subset_gcfreq(motif, sup)
    gcf_mat <- subset_gcfreq(motif, mat)

    s_sup <- score_sites(motif, sup, row$methylome, gcf_sup)
    s_mat <- score_sites(motif, mat, row$methylome, gcf_mat)
    if (is.null(s_sup) || is.null(s_mat)) return(NULL)

    has_ctrl <- row$control %in% names(msites)
    c_sup <- if (has_ctrl) score_sites(motif, sup, row$control, gcf_sup) else NULL
    c_mat <- if (has_ctrl) score_sites(motif, mat, row$control, gcf_mat) else NULL

    delta <- s_sup$current - s_mat$current
    delta_ctrl <- if (is.null(c_sup) || is.null(c_mat)) NA_real_ else c_sup$current - c_mat$current

    data.table(
      set_id = row$set_id, factor = row$factor, motif = motif,
      chip_label = row$chip_label, match = row$match,
      methylome = row$methylome, control_methylome = row$control,
      n_sites_total = length(sites),
      n_supported = sum(is_sup), n_matched_pairs = length(mm$supported),
      frac_supported = sum(is_sup) / length(sites),
      delta = delta,
      delta_uncorrected = s_sup$observed_only - s_mat$observed_only,
      delta_control = delta_ctrl,
      specificity = delta_ctrl - delta,
      cpg_sup = mean(comp[is_sup][idx %in% mm$supported, cpg]),
      cpg_mat = mean(comp[idx %in% mm$matched, cpg]),
      gc_sup = mean(comp[is_sup][idx %in% mm$supported, gc]),
      gc_mat = mean(comp[idx %in% mm$matched, gc])
    )
  }))
  saveRDS(res, cache)
  res
}

A <- rbindlist(lapply(seq_len(nrow(manifest)), function(i) run_split(manifest[i])), fill = TRUE)
fwrite(A, file.path(out.dir, "A_site_level_split.csv"))
log_success("A: ", nrow(A), " comparisons across ", uniqueN(A$methylome), " methylomes")

#####################################################################
# Motif level, deviation taken from the same methylome as the ChIP
#####################################################################

B <- rbindlist(lapply(unique(manifest$methylome), function(key) {
  sub <- manifest[methylome == key]
  rbindlist(lapply(seq_len(nrow(sub)), function(i) {
    row <- sub[i]
    peaks <- readRDS(row$peak_file)
    peaks_b <- peaks[seq_len(min(length(peaks), b.peak.budget))]
    rbindlist(lapply(tf_to_motifs(row$factor, names(tf_bindsites)), function(motif) {
      sites <- subsetByOverlaps(tf_bindsites[[motif]], enhancer)
      if (length(sites) < min.sites.per.motif) return(NULL)
      d <- score_sites(motif, sites, key, gcfreqs[motif])
      if (is.null(d)) return(NULL)
      mid <- start(sites) + floor(width(sites) / 2)
      data.table(
        methylome = key, set_id = row$set_id, factor = row$factor, motif = motif,
        match = row$match, n_peaks = length(peaks), n_sites = length(sites),
        frac_supported = mean(overlapsAny(
          GRanges(seqnames(sites), IRanges(mid, mid)), peaks_b)),
        deviation = d$current
      )
    }))
  }))
}), fill = TRUE)

B_cor <- B[, if (.N >= min.factors.for.B) {
  ct <- suppressWarnings(cor.test(frac_supported, deviation, method = "spearman"))
  .(n_factors = .N, rho = unname(ct$estimate), p = ct$p.value)
} else .(n_factors = .N, rho = NA_real_, p = NA_real_), by = methylome]

fwrite(B, file.path(out.dir, "B_motif_level.csv"))
fwrite(B_cor, file.path(out.dir, "B_motif_level_correlation.csv"))
print(B_cor)

#####################################################################

if (nrow(A)) {
  log_info("Median delta, matched methylome: ", signif(median(A$delta, na.rm = TRUE), 3),
           " ; mismatched methylome: ", signif(median(A$delta_control, na.rm = TRUE), 3))
  log_info(sum(A$specificity > 0, na.rm = TRUE), " of ",
           sum(is.finite(A$specificity)),
           " comparisons are stronger in the matched methylome")
}
log_success("Done.")
#!/usr/bin/env Rscript

#####################################################################
# 11_bp_chip_validation.R
# created on 06-09-2026 by Irem B Gunduz
# Histone ChIP signal at the binding sites of the motifs the methylTFR
# and VIPER integration calls active. This script only computes and
# stores, the panels are drawn in 10_bp_viper_supplementary.R
#
# The bigwigs are large and read over a shared mount, so the per sample
# signal is cached and the script can be stopped and resumed
#####################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(logger)
  library(GenomicRanges)
  library(rtracklayer)
  library(parallel)
  library(methylTFR)
  library(methylTFRAnnotationHg38)
})
set.seed(42)

#####################################################################
# Settings
#####################################################################

motifSet <- "jaspar2020_distal"
tfSet <- "jaspar2020"

marks <- c("H3K27ac", "H3K4me1", "H3K4me3", "H3K36me3", "H3K27me3", "H3K9me3")
mark.primary <- "H3K27ac"

# Sites are subsampled and trimmed to a window around their centre, otherwise
# a single sample means reading a few million ranges off the shared mount
max.sites.per.motif <- 2000L
site.window <- 200L

# fc.signal is scaled per experiment, so every motif in a sample shares that
# sample's global level. Standardising each sample across the motifs removes it
# and leaves the motif to motif structure the correlations are meant to test
standardise.per.sample <- TRUE

# Every motif that carries a deviation and a VIPER regulon, rather than the
# integration hits against a sample of the rest. The set still records which
# is which, it is simply no longer the thing that defines the pool
use.all.motifs <- TRUE
background.n <- 50L

# The aggregate signal profile: the conventional way to show a ChIP result,
# and the first thing a reviewer asks for. Computed only for the example
# motifs plus a set of random distal regions that normalises every sample
profile.window <- 2000L
profile.bins <- 40L
profile.sites <- 500L
profile.random.n <- 500L
# Ordered preference. The first available ones are taken, anything missing is
# reported and the remaining slots fall back to the automatic pick
profile.motifs <- c(
  "CEBPB", "SPI1", "JUNB", "RELA",
  "FOSL2::JUNB", "SPIB", "MAFB", "MAF", "NR1H3::RXRA"
)
profile.motifs.n <- 4L

min.samples.per.celltype <- 5L
drop.cell.types <- c("other", "thymocyte")

n.cores <- 4L
# Set to a small number to try the pipeline on a few samples first
max.samples <- NULL

analysis.dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/blueprint/"
dev.tag <- "mTFR_devs_230826"
mofa.dir <- file.path(analysis.dir, "mofa_230826")

dev.file <- file.path(analysis.dir, dev.tag, paste0(motifSet, "_deviations.RDS"))
viper.file <- file.path(mofa.dir, "bp_viper_activity.RDS")
map.file <- file.path(mofa.dir, "bp_rna_wgbs_map.tsv")
distal.file <- "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFRAnnotationHg38_old/inst/extdata/distal_regions.RDS"

cache.dir <- file.path(analysis.dir, "chip_signal")
if (!dir.exists(cache.dir)) dir.create(cache.dir, recursive = TRUE)

github.dir <- "/icbb_triton/scratch/igunduz/methylTFR_manuscript/github/methylTFR_manuscript/"
annot.file <- file.path(github.dir, "tables", "annotated_wgbs_with_rna_chip_matches.csv")
integration.file <- file.path(github.dir, "tables", "mixed", "bp_mtfr_viper_integration.csv")
mixed.dir <- file.path(github.dir, "tables", "mixed")
if (!dir.exists(mixed.dir)) dir.create(mixed.dir, recursive = TRUE)

signal.out <- file.path(mofa.dir, "bp_chip_signal.RDS")

#####################################################################
# Helpers
#####################################################################

zscore_vec <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) {
    return(rep(0, length(x)))
  }
  (x - mean(x, na.rm = TRUE)) / s
}

# JASPAR writes heterodimers as FOS::JUNB and variants as JUN(var.2), so a
# motif can map to several TFs
motif_to_tfs <- function(motif) {
  parts <- unlist(strsplit(motif, "::", fixed = TRUE))
  parts <- sub("\\s*\\(var\\.[0-9]+\\)$", "", parts)
  toupper(trimws(parts))
}

cor_across <- function(a, b) {
  ok <- is.finite(a) & is.finite(b)
  if (sum(ok) < 4) {
    return(NA_real_)
  }
  suppressWarnings(cor(a[ok], b[ok], method = "pearson"))
}

# Cell type means removed from both sides, so what is left is the agreement
# within a cell type rather than the shared lineage signal
centre_within <- function(v, g) {
  m <- tapply(v, g, mean, na.rm = TRUE)
  v - as.numeric(m[as.character(g)])
}

#####################################################################
# The ChIP inventory
#####################################################################

for (f in c(annot.file, integration.file, dev.file, viper.file, map.file)) {
  if (!file.exists(f)) stop("Missing input: ", f)
}

annot <- fread(annot.file)
chip <- annot[has_matching_chip == TRUE & has_matching_rna == TRUE]
if (nrow(chip) == 0) stop("No sample carries both ChIP and RNA")

# matching_chip_files_with_mark reads "MARK: path;MARK: path"
inventory <- rbindlist(lapply(seq_len(nrow(chip)), function(i) {
  parts <- trimws(unlist(strsplit(chip$matching_chip_files_with_mark[i], ";", fixed = TRUE)))
  parts <- parts[nzchar(parts)]
  data.table(
    bedFile = chip$bedFile[i],
    cellTypeGroup = chip$cellTypeGroup[i],
    mark = trimws(sub(":.*$", "", parts)),
    path = trimws(sub("^[^:]*:\\s*", "", parts))
  )
}))
inventory <- inventory[mark %in% marks]
log_info(
  "ChIP inventory: ", uniqueN(inventory$bedFile), " samples, ",
  paste(names(table(inventory$mark)), table(inventory$mark), sep = " = ", collapse = ", ")
)

#####################################################################
# The two other modalities, on the samples the ChIP covers
#####################################################################

dev_obj <- readRDS(dev.file)
mtfr_z <- deviationZScores(dev_obj)
viper_act <- readRDS(viper.file)
map <- as.data.table(read.delim(map.file, stringsAsFactors = FALSE))

map <- map[bedFile %in% colnames(mtfr_z) & rna_id %in% colnames(viper_act) &
  bedFile %in% inventory$bedFile & !cellTypeGroup %in% drop.cell.types]
keep_types <- names(which(table(map$cellTypeGroup) >= min.samples.per.celltype))
map <- map[cellTypeGroup %in% keep_types]
if (nrow(map) < 6) {
  stop(
    "Only ", nrow(map), " samples carry ChIP, RNA and methylation in a cell type ",
    "with at least ", min.samples.per.celltype, " of them"
  )
}
if (!is.null(max.samples)) map <- head(map, max.samples)

inventory <- inventory[bedFile %in% map$bedFile]
log_info(
  nrow(map), " samples with all three modalities across ",
  length(unique(map$cellTypeGroup)), " cell types: ",
  paste(names(table(map$cellTypeGroup)), table(map$cellTypeGroup),
    sep = " = ", collapse = ", "
  )
)

#####################################################################
# The motifs to test
#####################################################################

integration <- fread(integration.file)
hits <- unique(integration[in_heatmap == TRUE]$motif)
pool <- setdiff(unique(integration$motif), hits)
background <- if (use.all.motifs) pool else sample(pool, min(background.n, length(pool)))
log_info(
  length(hits), " integration motifs, ", length(background), " other motifs",
  if (use.all.motifs) " (every paired motif)" else ""
)

motif_set <- data.table(
  motif = c(hits, background),
  set = c(
    rep("Integration motifs", length(hits)),
    rep("Background motifs", length(background))
  )
)

#####################################################################
# Binding sites, restricted the same way the deviations were
#####################################################################

tf_bindsites <- getTFbindsites(motifSet = tfSet)
missing <- setdiff(motif_set$motif, names(tf_bindsites))
if (length(missing) > 0) {
  log_warn(
    length(missing), " motif(s) have no binding sites and are dropped: ",
    paste(head(missing, 5), collapse = ", ")
  )
  motif_set <- motif_set[!motif %in% missing]
}

enhancer <- readRDS(distal.file)
log_info("Restricting to ", length(enhancer), " distal regions")

site_list <- lapply(motif_set$motif, function(m) {
  gr <- subsetByOverlaps(tf_bindsites[[m]], enhancer)
  if (length(gr) == 0) {
    return(NULL)
  }
  if (length(gr) > max.sites.per.motif) gr <- gr[sample(length(gr), max.sites.per.motif)]
  gr <- resize(gr, width = 2L * site.window, fix = "center")
  mcols(gr) <- DataFrame(motif = rep(m, length(gr)))
  gr
})
site_list <- Filter(Negate(is.null), site_list)
if (length(site_list) == 0) stop("No motif has a binding site inside the distal regions")
sites <- do.call(c, unname(site_list))
seqlevelsStyle(sites) <- "UCSC"
sites <- sort(trim(sites))
site_motif <- factor(mcols(sites)$motif, levels = motif_set$motif)
log_info(
  length(sites), " binding site windows over ", uniqueN(site_motif), " motifs, ",
  "median ", median(as.integer(table(site_motif))), " per motif"
)

#####################################################################
# Mean signal per motif per sample, cached
#####################################################################

# The cached signal is only valid for the motif list it was computed on, so
# the list goes into the file name. Without it, changing the motifs silently
# returns the previous run's numbers
motif_tag <- function(motifs) {
  substr(
    paste0(
      length(motifs), "m",
      sum(utf8ToInt(paste(sort(motifs), collapse = "")))
    ),
    1, 24
  )
}
sites_tag <- motif_tag(motif_set$motif)
log_info("Signal cache tag: ", sites_tag)

signal_for <- function(path, mk, sample_id) {
  cache <- file.path(
    cache.dir, paste0(sample_id, "_", mk, "_", motifSet, "_", sites_tag, ".rds")
  )
  if (file.exists(cache)) {
    return(readRDS(cache))
  }
  if (!file.exists(path)) {
    log_warn("Missing bigwig, skipped: ", path)
    return(NULL)
  }
  bw <- BigWigFile(path)
  keep <- seqlevels(sites) %in% seqlevels(seqinfo(bw))
  query <- keepSeqlevels(sites, seqlevels(sites)[keep], pruning.mode = "coarse")
  vals <- import(bw, which = query, as = "NumericList")
  per_site <- vapply(vals, function(v) {
    if (length(v) == 0) NA_real_ else mean(v, na.rm = TRUE)
  }, numeric(1))
  out <- tapply(per_site, factor(mcols(query)$motif, levels = levels(site_motif)),
    mean,
    na.rm = TRUE
  )
  saveRDS(out, cache)
  out
}

signal_matrix <- function(mk) {
  inv <- inventory[mark == mk][!duplicated(bedFile)]
  inv <- inv[match(map$bedFile, bedFile)]
  inv <- inv[!is.na(bedFile)]
  if (nrow(inv) == 0) {
    return(NULL)
  }
  log_info("Reading ", nrow(inv), " ", mk, " bigwigs")
  cols <- mclapply(seq_len(nrow(inv)), function(i) {
    signal_for(inv$path[i], mk, inv$bedFile[i])
  }, mc.cores = n.cores)
  ok <- !vapply(cols, is.null, logical(1))
  if (!any(ok)) {
    return(NULL)
  }
  mat <- do.call(cbind, cols[ok])
  colnames(mat) <- inv$bedFile[ok]
  mat
}

signals <- setNames(lapply(marks, signal_matrix), marks)
signals <- signals[!vapply(signals, is.null, logical(1))]
if (length(signals) == 0) stop("No ChIP signal could be read")
if (!mark.primary %in% names(signals)) {
  stop("The primary mark ", mark.primary, " could not be read")
}
if (standardise.per.sample) {
  signals <- lapply(signals, function(mat) {
    out <- apply(mat, 2, zscore_vec)
    dimnames(out) <- dimnames(mat)
    out
  })
}

#####################################################################
# Correlations of the signal with each modality
#####################################################################

viper_tf <- setNames(integration$viper_tf, integration$motif)

cor_tab <- rbindlist(lapply(names(signals), function(mk) {
  mat <- signals[[mk]]
  samples <- colnames(mat)
  rna_ids <- map$rna_id[match(samples, map$bedFile)]
  groups <- map$cellTypeGroup[match(samples, map$bedFile)]
  rbindlist(lapply(intersect(rownames(mat), motif_set$motif), function(m) {
    chip_v <- as.numeric(mat[m, ])
    mtfr_v <- as.numeric(mtfr_z[m, samples])
    tf <- viper_tf[[m]]
    has_viper <- !is.null(tf) && tf %in% rownames(viper_act)
    viper_v <- if (has_viper) as.numeric(viper_act[tf, rna_ids]) else rep(NA_real_, length(samples))
    data.table(
      motif = m,
      viper_tf = if (has_viper) tf else NA_character_,
      mark = mk,
      r_mtfr = cor_across(chip_v, mtfr_v),
      r_viper = cor_across(chip_v, viper_v),
      # The motif's own relationship between the two modalities, on the same
      # samples. It sets which way the chromatin is expected to point
      r_mtfr_viper = cor_across(mtfr_v, viper_v),
      r_mtfr_within = cor_across(centre_within(chip_v, groups), centre_within(mtfr_v, groups)),
      r_viper_within = cor_across(centre_within(chip_v, groups), centre_within(viper_v, groups))
    )
  }))
}))
cor_tab <- merge(cor_tab, motif_set, by = "motif")
write.csv(cor_tab, file.path(mixed.dir, "bp_chip_validation_correlations.csv"), row.names = FALSE)

primary <- cor_tab[mark == mark.primary]
log_info(
  mark.primary, " median |r| with methylTFR, integration ",
  round(median(abs(primary[set == "Integration motifs"]$r_mtfr), na.rm = TRUE), 3),
  " against background ",
  round(median(abs(primary[set == "Background motifs"]$r_mtfr), na.rm = TRUE), 3)
)
log_info(
  mark.primary, " median |r| within cell type, methylTFR ",
  round(median(abs(primary$r_mtfr_within), na.rm = TRUE), 3),
  ", VIPER ", round(median(abs(primary$r_viper_within), na.rm = TRUE), 3)
)

#####################################################################
# The example motifs
#
# JASPAR is full of AP-1 dimers and ranking on agreement alone returns
# nothing else, so a TF family is allowed one example
#####################################################################

tf_family <- function(tf) {
  fam <- sub("[0-9]+$", "", tf)
  fam <- sub("^FOSL$", "FOS", fam)
  fam <- sub("^JUN[BD]$", "JUN", fam)
  fam
}

families_of <- function(motifs) {
  unique(unlist(lapply(motifs, function(m) {
    vapply(motif_to_tfs(m), tf_family, character(1), USE.NAMES = FALSE)
  })))
}

auto_examples <- function(n, exclude = character(0)) {
  if (n <= 0) {
    return(primary[0])
  }
  ranked <- primary[r_mtfr < 0 & r_viper > 0][order(-(abs(r_mtfr) + abs(r_viper)))]
  taken <- families_of(exclude)
  keep <- integer(0)
  for (i in seq_len(nrow(ranked))) {
    if (ranked$motif[i] %in% exclude) next
    fams <- families_of(ranked$motif[i])
    if (any(fams %in% taken)) next
    taken <- c(taken, fams)
    keep <- c(keep, i)
    if (length(keep) == n) break
  }
  ranked[keep]
}

# The footprint needs only binding sites, so a requested motif is shown even
# without a paired ChIP/mRNA/VIPER correlation; its sites must overlap the distal set
has_distal_sites <- function(m) {
  m %in% names(tf_bindsites) && length(subsetByOverlaps(tf_bindsites[[m]], enhancer)) > 0
}
example_row <- function(m) if (m %in% primary$motif) primary[match(m, primary$motif)] else data.table(motif = m)

examples <- primary[0]
if (!is.null(profile.motifs)) {
  ok <- profile.motifs[vapply(profile.motifs, has_distal_sites, logical(1))]
  miss <- setdiff(profile.motifs, ok)
  if (length(miss) > 0) log_warn("No distal binding sites, skipped as a profile motif: ", paste(miss, collapse = ", "))
  examples <- rbindlist(lapply(head(ok, profile.motifs.n), example_row), fill = TRUE)
} else {
  examples <- auto_examples(profile.motifs.n)
}
if (nrow(examples) == 0) stop("No example motif could be chosen")
log_info("Example motifs: ", paste(examples$motif, collapse = ", "))

#####################################################################
# Aggregate signal profiles around the motif centres
#
# Every window is the same width so the bins can be averaged as a matrix.
# Each sample is divided by its own mean over the random distal regions,
# which puts every sample on an enrichment scale and removes the per
# experiment scaling of fc.signal
#####################################################################

profile_windows <- function(gr, label) {
  gr <- resize(gr, width = 1L, fix = "center")
  gr <- resize(gr, width = 2L * profile.window, fix = "center")
  seqlevelsStyle(gr) <- "UCSC"
  gr <- gr[width(gr) == 2L * profile.window]
  gr <- gr[start(gr) > 0]
  mcols(gr) <- DataFrame(group = rep(label, length(gr)))
  gr
}

set.seed(42)
profile_list <- lapply(seq_len(nrow(examples)), function(i) {
  m <- examples$motif[i]
  gr <- subsetByOverlaps(tf_bindsites[[m]], enhancer)
  if (length(gr) > profile.sites) gr <- gr[sample(length(gr), profile.sites)]
  profile_windows(gr, m)
})

# Distal regions are regulatory elements, so centring the control windows on
# their midpoints centres them on their own peaks and the null comes out with
# a peak in it. Random positions drawn uniformly inside the regions have no
# such centre, which is what a null needs
rand_regions <- enhancer[sample(length(enhancer), min(profile.random.n, length(enhancer)))]
rand <- GRanges(
  seqnames(rand_regions),
  IRanges(
    start = start(rand_regions) + floor(runif(length(rand_regions)) * width(rand_regions)),
    width = 1
  )
)
profile_list <- c(profile_list, list(profile_windows(rand, "Random distal positions")))
profile_ranges <- sort(do.call(c, unname(profile_list)))
profile_group <- factor(mcols(profile_ranges)$group,
  levels = c(examples$motif, "Random distal positions")
)
profile_tag <- motif_tag(c(examples$motif, "rand"))
log_info(
  length(profile_ranges), " profile windows of ", 2L * profile.window, " bp: ",
  paste(names(table(profile_group)), table(profile_group), sep = " = ", collapse = ", ")
)

bin_width <- (2L * profile.window) %/% profile.bins
bin_index <- rep(seq_len(profile.bins), each = bin_width)
bin_pos <- (seq_len(profile.bins) - 0.5) * bin_width - profile.window

profile_for <- function(path, sample_id) {
  cache <- file.path(
    cache.dir, paste0(sample_id, "_", mark.primary, "_profile_", profile_tag, ".rds")
  )
  if (file.exists(cache)) {
    return(readRDS(cache))
  }
  if (!file.exists(path)) {
    return(NULL)
  }
  bw <- BigWigFile(path)
  keep <- seqlevels(profile_ranges) %in% seqlevels(seqinfo(bw))
  query <- keepSeqlevels(profile_ranges, seqlevels(profile_ranges)[keep],
    pruning.mode = "coarse"
  )
  grp <- factor(mcols(query)$group, levels = levels(profile_group))
  out <- rbindlist(lapply(levels(grp), function(g) {
    sel <- query[grp == g]
    if (length(sel) == 0) {
      return(NULL)
    }
    vals <- import(bw, which = sel, as = "NumericList")
    ok <- lengths(vals) == 2L * profile.window
    if (!any(ok)) {
      return(NULL)
    }
    mat <- matrix(as.numeric(unlist(vals[ok])), nrow = 2L * profile.window)
    mat[!is.finite(mat)] <- 0
    binned <- rowsum(mat, bin_index) / bin_width
    data.table(group = g, pos = bin_pos, value = rowMeans(binned), n = sum(ok))
  }))
  if (is.null(out) || nrow(out) == 0) {
    return(NULL)
  }
  # Every profile of this sample is divided by the same number: the signal in
  # the flanks of the random windows. That puts the sample on an enrichment
  # scale without flattening the differences between cell types, which are
  # the point, and leaves the random facet sitting at one
  base <- mean(
    out[group == "Random distal positions" & abs(pos) >= profile.window * 0.75]$value,
    na.rm = TRUE
  )
  out[, value := if (is.finite(base) && base > 0) value / base else NA_real_]
  out[, sample := sample_id]
  saveRDS(out, cache)
  out
}

inv_primary <- inventory[mark == mark.primary][!duplicated(bedFile)]
inv_primary <- inv_primary[match(map$bedFile, bedFile)][!is.na(bedFile)]
log_info("Reading ", nrow(inv_primary), " ", mark.primary, " bigwigs for the profiles")
profiles <- rbindlist(Filter(Negate(is.null), mclapply(seq_len(nrow(inv_primary)), function(i) {
  profile_for(inv_primary$path[i], inv_primary$bedFile[i])
}, mc.cores = n.cores)))

if (nrow(profiles) == 0) {
  log_warn("No profile could be computed")
} else {
  profiles[, cellTypeGroup := map$cellTypeGroup[match(sample, map$bedFile)]]
  log_info(
    "Profiles over ", uniqueN(profiles$sample), " samples, peak enrichment ",
    round(max(profiles[group != "Random distal regions"]$value, na.rm = TRUE), 2)
  )
  write.csv(profiles, file.path(mixed.dir, "bp_chip_profiles.csv"), row.names = FALSE)
}

#####################################################################
# Everything the figure script needs
#####################################################################

saveRDS(
  list(
    signals = signals,
    map = map,
    motif_set = motif_set,
    correlations = cor_tab,
    examples = examples,
    profiles = profiles,
    profile.window = profile.window,
    mark.primary = mark.primary,
    marks = names(signals),
    standardised = standardise.per.sample
  ),
  signal.out
)
log_success("Wrote ", signal.out)
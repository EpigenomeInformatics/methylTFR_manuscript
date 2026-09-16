#!/usr/bin/env Rscript
# 07_scECHO_mTFR_ind.R
# Created on 15-09-26 by Irem Begum Gunduz
# Standalone methylTFR UMAPs, one file per TF, plus diagnostics.
# Purpose: isolate what makes the methylTFR panel look "weird" in the merged
# figure. No ATAC, no patchwork, no MAGIC -- just the deviations on the LSI UMAP.

suppressPackageStartupMessages({
  library(ggplot2)
  library(methylTFR)
  library(logger)
})
set.seed(42)

#####################################################################
# Config  (edit paths if needed)
#####################################################################
features   <- c("POU2F2", "SPIB", "EOMES")
motifSet   <- "jaspar2020_distal"
dev.dir    <- file.path("/scratch/icbb/igunduz/methylTFR_manuscript/echo",
                        paste0("mTFR_sc_230826_", motifSet))
lsi.file   <- "/icbb/projects/igunduz/DARPA_analysis/artemis_031023/itLSI_res/itLSI_res_sub11kpc30.rds"
cell.id.col <- "Cell_UID"

out.dir    <- "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR_manuscript/figures/echo/methylTFR_individual"
if (!dir.exists(out.dir)) dir.create(out.dir, recursive = TRUE)

# Display choices -- flip these to see which one matches your good standalone plot
zscore        <- TRUE   # scale each motif to unit variance (comparable across TFs)
zscore.center <- FALSE  # FALSE = x / sd  -> keeps grey at biological zero (recommended;
                        #   matches the good raw look, just rescaled).
                        # TRUE  = (x - mean) / sd -> classic z, but recenters grey on the
                        #   per-motif MEAN, which flips a broadly-negative motif to orange.
clip.quantile <- 0.98   # symmetric colour clip at quantile(abs(value), clip.quantile)
point.size   <- 0.35    # a bit larger than the merged panel; sparse cells read faint
plot.na.grey <- TRUE    # draw unmatched / NA cells as light grey underneath (context)
pal.low <- "#5E3C99"; pal.mid <- "grey92"; pal.high <- "#E66101"

#####################################################################
# 1. Load UMAP coordinates + cell IDs
#####################################################################
log_info("Loading LSI result: ", lsi.file)
res <- readRDS(lsi.file)
umap <- as.data.frame(res$umapCoord)[, 1:2]
colnames(umap) <- c("UMAP1", "UMAP2")
umap[[cell.id.col]] <- rownames(umap)
log_info("UMAP has ", nrow(umap), " cells")

#####################################################################
# 2. Load methylTFR deviations
#####################################################################
dev.files <- list.files(dev.dir, pattern = "_deviations\\.RDS$", full.names = TRUE)
log_info("Found ", length(dev.files), " deviation file(s)")
mats <- lapply(dev.files, function(f) as.matrix(methylTFR::deviations(readRDS(f))))
motif.rows <- Reduce(intersect, lapply(mats, rownames))
dev.mat <- do.call(cbind, lapply(mats, function(m) m[motif.rows, , drop = FALSE]))
dev.mat <- dev.mat[, !duplicated(colnames(dev.mat)), drop = FALSE]
dev.mat[!is.finite(dev.mat)] <- NA
log_info("dev.mat: ", nrow(dev.mat), " motifs x ", ncol(dev.mat), " cells")
log_info("motifs available: ", paste(rownames(dev.mat), collapse = ", "))

#####################################################################
# 3. DIAGNOSTICS  -- this is what we are here to see
#####################################################################
cat("\n=========== DIAGNOSTICS ===========\n")
cat("UMAP cell IDs (head):    ", paste(head(umap[[cell.id.col]], 3), collapse = " | "), "\n")
cat("dev.mat cell IDs (head): ", paste(head(colnames(dev.mat), 3), collapse = " | "), "\n")

n.match <- sum(!is.na(match(umap[[cell.id.col]], colnames(dev.mat))))
cat(sprintf("MATCHED cells: %d / %d  (%.1f%% of the UMAP)\n",
            n.match, nrow(umap), 100 * n.match / nrow(umap)))
if (n.match < 0.5 * nrow(umap)) {
  cat(">>> LOW OVERLAP: the ID strings almost certainly differ in format.\n",
      ">>> Compare the two 'head' lines above (prefixes, '#', trailing '-1').\n")
}
cat("dev.mat cells NOT in UMAP: ", sum(is.na(match(colnames(dev.mat), umap[[cell.id.col]]))), "\n")

for (g in intersect(features, rownames(dev.mat))) {
  v <- dev.mat[g, ]
  cat(sprintf("%-8s  n=%d  NA=%d  range=[%.3g, %.3g]  sd=%.3g\n",
              g, length(v), sum(is.na(v)),
              suppressWarnings(min(v, na.rm = TRUE)),
              suppressWarnings(max(v, na.rm = TRUE)),
              stats::sd(v, na.rm = TRUE)))
}
cat("===================================\n\n")

#####################################################################
# 4. Optional per-motif z-score
#####################################################################
if (zscore) {
  sdv <- apply(dev.mat, 1, stats::sd, na.rm = TRUE); sdv[!is.finite(sdv) | sdv == 0] <- 1
  mu  <- if (zscore.center) rowMeans(dev.mat, na.rm = TRUE) else 0
  dev.mat <- (dev.mat - mu) / sdv
  log_info("Applied per-motif scaling (center = ", zscore.center, ")")
}

#####################################################################
# 5. Theme + plot one UMAP per TF
#####################################################################
theme_umap <- theme_classic(base_size = 9) + theme(
  plot.title  = element_text(hjust = 0.5, size = 12, face = "bold"),
  axis.text   = element_blank(), axis.ticks = element_blank(),
  legend.key.width = unit(0.3, "cm"), legend.key.height = unit(0.6, "cm")
)

plot_features <- intersect(features, rownames(dev.mat))
if (length(plot_features) == 0) stop("None of `features` are present in dev.mat")

for (gene in plot_features) {
  score <- rep(NA_real_, nrow(umap))
  hit <- match(umap[[cell.id.col]], colnames(dev.mat))
  score[!is.na(hit)] <- dev.mat[gene, hit[!is.na(hit)]]

  df <- data.frame(UMAP1 = umap$UMAP1, UMAP2 = umap$UMAP2, val = score)

  lim <- stats::quantile(abs(df$val), clip.quantile, na.rm = TRUE)
  if (!is.finite(lim) || lim == 0) lim <- max(abs(df$val), na.rm = TRUE)

  p <- ggplot()
  # NA / unmatched cells as context so the full UMAP shape is always visible
  if (plot.na.grey) {
    p <- p + geom_point(data = df[is.na(df$val), ],
                        aes(UMAP1, UMAP2), colour = "grey88", size = point.size)
  }
  # coloured cells on top, high-|value| last so signal is not buried
  dcol <- df[!is.na(df$val), ]
  dcol <- dcol[order(abs(dcol$val)), ]
  p <- p +
    geom_point(data = dcol, aes(UMAP1, UMAP2, colour = val), size = point.size) +
    scale_colour_gradient2(low = pal.low, mid = pal.mid, high = pal.high,
                           midpoint = 0, limits = c(-lim, lim),
                           oob = scales::squish, name = "Score") +
    labs(title = paste(gene, "methylTFR"), x = "UMAP1", y = "UMAP2") +
    coord_fixed() + theme_umap

  f.pdf <- file.path(out.dir, paste0("methylTFR_", gene, ".pdf"))
  f.png <- file.path(out.dir, paste0("methylTFR_", gene, ".png"))
  ggsave(f.pdf, p, width = 5.5, height = 5, bg = "white")
  ggsave(f.png, p, width = 5.5, height = 5, dpi = 200, bg = "white")
  log_success("Saved ", f.png, "  (", sum(!is.na(df$val)), " coloured cells)")
}

log_success("Done. Individual UMAPs in: ", out.dir)
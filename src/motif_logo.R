#!/usr/bin/env Rscript

#####################################################################
# motif_logos.R
# Fetch PWMs from JASPAR2020 for a list of TFs and draw motif logos.
#####################################################################

suppressPackageStartupMessages({
  library(JASPAR2020); library(TFBSTools); library(ggplot2); library(ggseqlogo)
})

fig_dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR_manuscript/figures/"   # change if you want the PDF elsewhere
out_pdf <- file.path(fig_dir, "motif_logos.pdf")

tfs <- c(
  "TFAP2B", "PAX5", "PAX9", "PAX1", "NHLH2", "ASCL1", "NHLH1", "BHLHE22", "FERD3L", "PAX6",
  "EMX1", "PAX4", "EN1", "LHX1", "ELF4", "ELF2", "ETV5", "ETV6", "ELF5", "SPIC",
  "SPIB", "SPI1", "EHF", "ELF3", "IKZF1", "TFAP2C", "TFAP2B", "TFAP2A", "TFAP2E", "TFAP2C",
  "TGIF1", "CREB3L4", "PBX3", "TAL1::TCF3", "MYOG", "ATOH1", "MYF5", "BHLHA15", "ZBTB18", "EBF3",
  "EBF1", "TFAP4", "NEUROD1", "VSX2", "EGR4", "DPRX", "SOX8", "CUX1", "CUX2", "TEAD3", "CEBPB",
  "FOSL1::JUND", "CEBPB", "BATF", "TBX21", "SPIB", "JUN", "FOS2", "ETV5"
)
tfs <- unique(tfs)                       # drop duplicates

# --- pull the full human CORE set once, index by name -----------------------
pfms <- getMatrixSet(JASPAR2020, list(species = 9606, collection = "CORE"))
nm   <- vapply(pfms, name, character(1))

# case-insensitive match; keep first hit per requested TF
idx  <- match(toupper(tfs), toupper(nm))
found_tf  <- tfs[!is.na(idx)]
missing   <- tfs[is.na(idx)]
sel_pfms  <- pfms[idx[!is.na(idx)]]

cat(sprintf("Requested (unique): %d\n", length(tfs)))
cat(sprintf("Matched in JASPAR2020 CORE: %d\n", length(found_tf)))
if (length(missing) > 0)
  cat("NOT found (check spelling / not in JASPAR2020):\n  ",
      paste(missing, collapse = ", "), "\n")
if (length(sel_pfms) == 0) stop("No motifs matched; nothing to draw.")

# --- convert PFMs to position-probability matrices --------------------------
mats <- lapply(sel_pfms, function(p) { m <- as.matrix(p); sweep(m, 2, colSums(m), "/") })
names(mats) <- make.unique(found_tf)     # unique panel titles

# --- draw logos (a few columns, height scales with row count) ---------------
ncol_grid <- 5
nrow_grid <- ceiling(length(mats) / ncol_grid)
logo <- ggseqlogo(mats, method = "bits", ncol = ncol_grid)

ggsave(out_pdf, logo, width = 2.4 * ncol_grid, height = 1.8 * nrow_grid, limitsize = FALSE)
cat("Saved", out_pdf, "with", length(mats), "logos.\n")

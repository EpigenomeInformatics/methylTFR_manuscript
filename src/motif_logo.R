#!/usr/bin/env Rscript

#####################################################################
# motif_logo.R
# Fetch PWMs from JASPAR2020 for a list of TFs and draw motif logos.
#####################################################################

suppressPackageStartupMessages({
  library(JASPAR2020)
  library(TFBSTools)
  library(ggplot2)
  library(ggseqlogo)
})

#####################################################################
# Settings
#####################################################################

# CORE vertebrates rather than species 9606. Several vertebrate CORE
# matrices are annotated to the species they were derived in, so a human
# only query drops motifs that the analysis does use, among them the
# mouse derived bHLH matrices.
collection <- "CORE"
tax.group <- "vertebrates"

# Panels per row of the logo grid
ncol.grid <- 5

# Directories
github.dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR_manuscript"
fig.dir <- file.path(github.dir, "figures")
if (!dir.exists(fig.dir)) dir.create(fig.dir, recursive = TRUE)

out.pdf <- file.path(fig.dir, "motif_logos.pdf")
out.table <- file.path(github.dir, "tables", "motif_logo_matrices.csv")
if (!dir.exists(dirname(out.table))) dir.create(dirname(out.table), recursive = TRUE)

tfs <- c(
  "TFAP2B", "PAX5", "PAX9", "PAX1", "NHLH2", "ASCL1", "NHLH1", "BHLHE22", "FERD3L", "PAX6",
  "EMX1", "PAX4", "EN1", "LHX1", "ELF4", "ELF2", "ETV5", "ETV6", "ELF5", "SPIC",
  "SPIB", "SPI1", "EHF", "ELF3", "IKZF1", "TFAP2C", "TFAP2B", "TFAP2A", "TFAP2E", "TFAP2C",
  "TGIF1", "CREB3L4", "PBX3", "POU2F3", "MYOG", "ATOH1", "MYF5", "BHLHA15", "ZBTB18", "EBF3",
  "EBF1", "TFAP4", "NEUROD1", "VSX2", "EGR4", "DPRX", "SOX8", "CUX1", "CUX2", "TEAD3", "CEBPB",
  "FOSL1::JUND", "CEBPB", "BATF", "TBX21", "SPIB", "JUN", "FOS2", "ETV5", "FOSL2"
)
tfs <- unique(tfs)

#####################################################################
# Match the requested names against JASPAR
#####################################################################

pfms <- getMatrixSet(JASPAR2020, list(collection = collection, tax_group = tax.group))
nm <- vapply(pfms, name, character(1))
ids <- vapply(pfms, ID, character(1))

# Case insensitive, first matrix per name
idx <- match(toupper(tfs), toupper(nm))
found <- tfs[!is.na(idx)]
missing <- tfs[is.na(idx)]
sel_pfms <- pfms[idx[!is.na(idx)]]

cat(sprintf("Requested (unique): %d\n", length(tfs)))
cat(sprintf(
  "Matched in JASPAR2020 %s %s: %d\n",
  collection, tax.group, length(found)
))

#####################################################################
# Report what is still missing
#####################################################################

# A name absent on its own may exist only as half of a heterodimer, so
# the candidates are reported rather than substituted. A dimer motif is
# not the same motif as its monomer and swapping one in silently would
# mislabel the panel.
if (length(missing) > 0) {
  cat("NOT found:\n")
  for (tf in missing) {
    partners <- nm[grepl("::", nm, fixed = TRUE) &
      grepl(paste0("\\b", tf, "\\b"), nm, ignore.case = TRUE)]
    if (length(partners) > 0) {
      cat(sprintf(
        "  %-14s absent alone, present as: %s\n",
        tf, paste(unique(partners), collapse = ", ")
      ))
    } else {
      cat(sprintf("  %-14s no matrix under this name\n", tf))
    }
  }
}

if (length(sel_pfms) == 0) stop("No motifs matched, nothing to draw")

#####################################################################
# Draw
#####################################################################

# Position probability matrices, one column summing to one
mats <- lapply(sel_pfms, function(p) {
  m <- as.matrix(p)
  sweep(m, 2, colSums(m), "/")
})
names(mats) <- make.unique(found)

nrow.grid <- ceiling(length(mats) / ncol.grid)
logo <- ggseqlogo(mats, method = "bits", ncol = ncol.grid)

ggsave(out.pdf, logo,
  width = 2.4 * ncol.grid, height = 1.8 * nrow.grid,
  limitsize = FALSE
)
cat("Saved", out.pdf, "with", length(mats), "logos\n")


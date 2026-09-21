#!/usr/bin/env Rscript

#####################################################################
# 00_supplementary_tables.R
# created on 21-09-2026 by Irem B Gunduz
# Assemble the manuscript supplementary tables into one workbook
#   S1  methylTFR deviation Z-scores of human immune cells
#   S2  MOFA feature weights per factor and modality (VIPER, methylTFR)
#   S3  Differentially methylated distal regions, TN vs TEM (RnBeads)
#   S4  Differentially methylated distal regions, TN vs TCM (RnBeads)
#   S5  Differentially active TFs, TN vs TEM and TCM (methylTFR)
#   S6  LOLA enrichment between TN, TCM and TEM
#   S7  JASPAR2020 distal TFs: methylTFR activity, differential, VIPER
#   S8  MOFA feature weights per factor and modality (chromVAR, methylTFR)
#####################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(logger)
  library(openxlsx)
  library(RnBeads)
  library(LOLA)
  library(methylTFR)
  library(MOFA2)
})
set.seed(42)

#####################################################################
# Settings and paths
#####################################################################

motifSet <- "jaspar2020_distal"
region.type <- "distal"

# Repository and shared helper functions
github.dir <- "/icbb_triton/scratch/igunduz/methylTFR_manuscript/github/methylTFR_manuscript"
table.dir <- file.path(github.dir, "tables")
src.dir <- file.path(github.dir, "src")
source(file.path(src.dir, "lola_utils.R"))

# Blueprint immune cells: deviations and the VIPER MOFA model
bp.dir <- "/icbb_triton/scratch/igunduz/methylTFR_manuscript/blueprint"
bp.dev.file <- file.path(bp.dir, "mTFR_devs_230826", paste0(motifSet, "_deviations.RDS"))
bp.viper.model <- file.path(bp.dir, "mofa_230826", "mtfr_viper_model_bp.rds")

# Memory T cells: RnBeads differential methylation and LOLA enrichment
mt.dir <- "/icbb_triton/scratch/igunduz/methylTFR_manuscript/memoryTcells"
mt.rnb.set.path <- file.path(mt.dir, "reports", "data_import_data", "rnb.set_preprocessed")
mt.diffmeth.dir <- file.path(mt.dir, "reports", "differential_methylation_data", "differential_rnbDiffMeth")
mt.lola.file <- file.path(mt.diffmeth.dir, "TF_motifs_lola.rds")

# Distal region set the RnBeads run was built against, re-registered below so
# annotation() can return the coordinates of this custom region type
distal.file <- "/icbb_triton/scratch/igunduz/methylTFR_manuscript/github/methylTFRAnnotationHg38_old/inst/extdata/distal_regions.RDS"

# ECHO single-cell MOFA weights, already exported by src/echo/01_mofa_integration.R
echo.weights.file <- "/icbb_triton/scratch/igunduz/methylTFR_manuscript/echo/mofa_230826/MOFA_weights_long.csv"

# Exported tables read as-is
diff.em.file <- file.path(table.dir, paste0("diff_em_", motifSet, ".csv"))
diff.cm.file <- file.path(table.dir, paste0("diff_cm_", motifSet, ".csv"))
integ.file <- file.path(table.dir, "mixed", "bp_mtfr_viper_integration.csv")

# Output workbook
out.file <- file.path(table.dir, "supplementary_tables.xlsx")

#####################################################################
# S1  methylTFR deviation Z-scores of human immune cells
#####################################################################

log_info("S1: methylTFR deviation Z-scores (Blueprint immune cells)")
dev.bp <- readRDS(bp.dev.file)
z.bp <- deviationZScores(dev.bp)
S1 <- data.frame(motif = rownames(z.bp), z.bp, check.names = FALSE, row.names = NULL)

#####################################################################
# S2  MOFA feature weights per factor and modality (VIPER, methylTFR)
#####################################################################

log_info("S2: MOFA weights (VIPER + methylTFR)")
mofa.viper <- readRDS(bp.viper.model)
S2 <- get_weights(mofa.viper, as.data.frame = TRUE)
S2$view <- as.character(S2$view)

#####################################################################
# S3 & S4  Differentially methylated distal regions (RnBeads)
#####################################################################

log_info("S3/S4: RnBeads differentially methylated distal regions")

# Custom region types are not restored to the annotation registry on load, so
# re-register the distal set before annotation() can resolve its coordinates
distal <- as.data.frame(readRDS(distal.file))[, c("seqnames", "start", "end")]
colnames(distal) <- c("Chromosome", "Start", "End")
rnb.set.annotation(type = region.type, regions = distal, assembly = "hg38")

diffMeth <- load.rnb.diffmeth(mt.diffmeth.dir)
rnb.set <- load.rnb.set(mt.rnb.set.path)
distal.annot <- annotation(rnb.set, type = region.type)[, c("Chromosome", "Start", "End")]

cmps <- get.comparisons(diffMeth)
labels <- comparison_labels(cmps)

# One differential distal table per naive-vs-memory contrast, coordinates prepended
diffmeth_distal <- function(test, ref) {
  hit <- resolve_comparison(labels, test = test, ref = ref, what = "RnBeads comparison")
  tbl <- get.table(diffMeth, cmps[[hit$index]], region.type, return.data.frame = TRUE)
  stopifnot(nrow(tbl) == nrow(distal.annot))
  cbind(distal.annot, tbl)
}

S3 <- diffmeth_distal("TEM", "TN")
S4 <- diffmeth_distal("TCM", "TN")

#####################################################################
# S5  Differentially active TFs, TN vs TEM and TCM (methylTFR)
#####################################################################

log_info("S5: methylTFR differential TFs (TN vs TEM and TCM)")
em <- fread(diff.em.file)
em$comparison <- "TEM vs TN"
cm <- fread(diff.cm.file)
cm$comparison <- "TCM vs TN"
S5 <- rbindlist(list(cm, em))
setcolorder(S5, c("comparison", setdiff(names(S5), "comparison")))

#####################################################################
# S6  LOLA enrichment between TN, TCM and TEM
#####################################################################

log_info("S6: LOLA enrichment (TN, TCM, TEM)")
lola.res <- readRDS(mt.lola.file)
lola.comps <- names(lola.res$region)

# Flatten the distal-region enrichment of every pairwise contrast into one table
S6 <- rbindlist(lapply(lola.comps, function(cmp) {
  reg <- resolve_lola_region(lola.res, cmp, c(region.type, "tiling1kb"))
  d <- as.data.frame(lola.res$region[[cmp]][[reg]])
  d$comparison <- cmp
  d$region <- reg
  d
}), fill = TRUE)

#####################################################################
# S7  JASPAR2020 distal TFs: methylTFR activity, differential, VIPER
#####################################################################

log_info("S7: methylTFR / VIPER integration table")
S7 <- fread(integ.file)

#####################################################################
# S8  MOFA feature weights per factor and modality (chromVAR, methylTFR)
#####################################################################

log_info("S8: MOFA weights (chromVAR + methylTFR)")
S8 <- fread(echo.weights.file)

#####################################################################
# Write the workbook
#####################################################################

log_info("Writing workbook to ", out.file)
overview <- data.frame(
  Sheet = paste0("S", 1:8),
  Table = paste0("Supplementary Table ", 1:8),
  Title = c(
    "methylTFR deviation Z-scores of human immune cells",
    "Feature weights per factor and modality (VIPER, methylTFR deviations); MOFA results",
    "Pairwise differentially methylated distal regions between naive T cells and TEM (RnBeads)",
    "Pairwise differentially methylated distal regions between naive T cells and TCM (RnBeads)",
    "Pairwise differentially active TFs between naive T cells and TEM and TCM (methylTFR)",
    "LOLA enrichment statistics between TN, TCM and TEM",
    "JASPAR2020 distal TFs with methylTFR activity, differential status and VIPER scores",
    "Feature weights per factor and modality (chromVAR deviations, methylTFR deviations); MOFA results"
  )
)

sheets <- list(
  Overview = overview,
  S1 = S1, S2 = S2, S3 = S3, S4 = S4,
  S5 = S5, S6 = S6, S7 = S7, S8 = S8
)

write.xlsx(sheets, file = out.file, overwrite = TRUE)
log_success("Wrote ", out.file)

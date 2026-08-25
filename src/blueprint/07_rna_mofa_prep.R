#!/usr/bin/env Rscript

#####################################################################
# 07_rna_mofa_prep.R
# created on 07-04-2026 by Regina Nitsch
# Updated by IBG on 25-08-2026
# Build the RNA counts matrix for the Blueprint samples that have
# matching WGBS data, and write the WGBS to RNA sample map that 08
# reads. Deriving that map twice, with two different regular
# expressions, was what silently emptied the MOFA input.
#####################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(logger)
  library(stringr)
  library(SummarizedExperiment)
  library(tibble)
})
set.seed(12)

#####################################################################
# Settings
#####################################################################

# Disease samples are excluded here for the same reason as in 01
drop.disease <- TRUE

# Directories
github.dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR_manuscript"
annot.file <- file.path(github.dir, "tables", "annotated_wgbs_with_rna_chip_matches.csv")
rna.se.file <- "/icbb/projects/share/datasets/blueprint/rna_se.rds"

analysis.dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/blueprint"
out.dir <- file.path(analysis.dir, "mofa_230826")
if (!dir.exists(out.dir)) dir.create(out.dir, recursive = TRUE)

counts.file <- file.path(out.dir, "bp_rawRNAcounts.RDS")
map.file <- file.path(out.dir, "bp_rna_wgbs_map.tsv")

#####################################################################
# Sample annotation
#####################################################################

if (!file.exists(annot.file)) stop("Annotation not found: ", annot.file)
annot <- read.csv(annot.file, stringsAsFactors = FALSE)
log_info(nrow(annot), " samples in ", basename(annot.file))

# has_matching_rna reads back as logical or as the string "TRUE" depending
# on the writer, so it is normalised rather than compared to one of them
has_rna <- toupper(as.character(annot$has_matching_rna)) %in% c("TRUE", "T", "YES", "1")
annot <- annot[has_rna, ]
log_info(nrow(annot), " samples with matching RNA")

if (drop.disease) {
  keep <- !is.na(annot$DISEASE) & annot$DISEASE == "None"
  if (any(!keep)) {
    log_info("Dropping ", sum(!keep), " disease samples: ",
      paste(unique(annot$DISEASE[!keep]), collapse = ", "))
  }
  annot <- annot[keep, ]
}

# A sample can match several RNA files, separated by ";". The first is
# taken, as before, but the split is explicit rather than folded into the
# regular expression.
annot$rna_file <- sub(";.*$", "", annot$matching_rna_files)
annot <- annot[!is.na(annot$rna_file) & nzchar(annot$rna_file), ]

# Two identifier forms live in these file names:
#   ...ihec-grapenf-containerv1.1.0.IHECRE00000279.2.<uuid>.genes.results
# the IHEC reference id with its version, and the trailing uuid. Which one
# the SummarizedExperiment is keyed on is decided below by counting hits,
# because the two prep scripts previously assumed different ones.
annot$id_uuid <- sub(".*\\.([a-f0-9-]+)\\.genes\\.results$", "\\1", annot$rna_file)
annot$id_ihec <- sub(".*\\.(IHECRE[0-9]+\\.[0-9]+)\\..*$", "\\1", annot$rna_file)

# Doublets: several WGBS samples can point at the same RNA file
annot <- annot[!duplicated(annot$rna_file), ]
log_info(nrow(annot), " unique RNA files")

#####################################################################
# RNA SummarizedExperiment
#####################################################################

if (!file.exists(rna.se.file)) stop("RNA object not found: ", rna.se.file)
log_info("Loading ", rna.se.file)
rna_se <- readRDS(rna.se.file)

# Candidate keys on the object side, candidate ids on the annotation side.
# The pair with the most hits wins and is reported, so a naming change in
# either shows up as a number instead of an empty result.
se_keys <- list(colnames = colnames(rna_se))
for (nm in colnames(colData(rna_se))) {
  se_keys[[nm]] <- as.character(colData(rna_se)[[nm]])
}
annot_ids <- list(uuid = annot$id_uuid, ihec = annot$id_ihec)

hits <- expand.grid(
  se_key = names(se_keys), annot_id = names(annot_ids),
  stringsAsFactors = FALSE
)
hits$n <- mapply(function(sk, ai) {
  sum(se_keys[[sk]] %in% annot_ids[[ai]])
}, hits$se_key, hits$annot_id)
hits <- hits[order(-hits$n), ]
log_info("Best identifier matches:")
print(utils::head(hits, 5))

if (hits$n[1] == 0) {
  stop(
    "No RNA sample identifier matched the annotation. Annotation ids look ",
    "like '", annot$id_uuid[1], "' or '", annot$id_ihec[1],
    "'; the object carries ", paste(names(se_keys), collapse = ", "), "."
  )
}
se_key <- hits$se_key[1]
annot_id <- hits$annot_id[1]
log_info("Matching colData column '", se_key, "' against the ", annot_id, " ids")

key_vec <- se_keys[[se_key]]
annot$rna_id <- annot[[paste0("id_", annot_id)]]

keep_cols <- which(key_vec %in% annot$rna_id)
if (length(keep_cols) == 0) stop("No RNA columns left after matching")
rna_se <- rna_se[, keep_cols]
key_vec <- key_vec[keep_cols]
log_info(ncol(rna_se), " RNA samples retained")

# Column names of the counts matrix, and the map that 08 reads
colnames(rna_se) <- key_vec
annot <- annot[annot$rna_id %in% key_vec, ]
annot <- annot[match(key_vec, annot$rna_id), ]

#####################################################################
# Counts matrix, one row per gene symbol
#####################################################################

rna <- assay(rna_se, "count")

ids2gene <- as.data.frame(rowData(rna_se)[, c("gene_name", "gene_id")])
ids2gene$gene_id_short <- rownames(ids2gene)

rna <- rna %>%
  as.data.frame() %>%
  rownames_to_column("gene_id_short") %>%
  left_join(ids2gene, by = "gene_id_short")

# Genes with no symbol cannot be matched to a motif and would collapse
# into a single NA row
n_before <- nrow(rna)
rna <- rna[!is.na(rna$gene_name) & nzchar(rna$gene_name), ]
if (nrow(rna) < n_before) {
  log_info("Dropped ", n_before - nrow(rna), " genes without a symbol")
}

# Duplicate symbols are averaged. across() takes a function rather than
# extra dots, which dplyr 1.1 no longer forwards.
rna <- rna %>%
  select(-gene_id, -gene_id_short) %>%
  group_by(gene_name) %>%
  summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)), .groups = "drop") %>%
  as.data.frame()

rownames(rna) <- rna$gene_name
rna$gene_name <- NULL
rna <- as.matrix(rna)
log_info("Counts matrix: ", nrow(rna), " genes x ", ncol(rna), " samples")

saveRDS(rna, counts.file)
log_success("Wrote ", counts.file)

map <- data.frame(
  bedFile = annot$bedFile,
  rna_id = annot$rna_id,
  rna_file = annot$rna_file,
  cellTypeGroup = annot$cellTypeGroup,
  cellTypeShort = annot$cellTypeShort,
  stringsAsFactors = FALSE
)
write.table(map, map.file, sep = "\t", quote = FALSE, row.names = FALSE)
log_success("Wrote ", map.file, " with ", nrow(map), " paired samples")

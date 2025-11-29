## Viper analysis using RNAseq expression data
library(dorothea)
library(dplyr)
library(stringr)
library(SummarizedExperiment)
set.seed(12)

# Read in annotation file
annot <- read.csv("methylTFR/sample_annotation/annotated_wgbs_with_rna_chip_matches.csv", stringsAsFactors=FALSE)

# Subset to samples that have matching RNA
annot <- subset(annot, has_matching_rna=="TRUE")

# Remove doublets
annot_unique <- annot[!duplicated(annot$matching_rna_files),]

# Extract IDs from filenames
#(some samples have more files because donor_id matched several rna_paths)
rna_ids_list <- strsplit(annot_unique$matching_rna_files, ";") |>
  lapply(function(paths) {
    sub(".*\\.[0-9]+\\.([a-z0-9\\-]+)\\.genes\\.results$", "\\1", paths)
  })
rna_ids <- sub(".*\\.[0-9]+\\.([a-z0-9\\-]+)\\.genes\\.results$", "\\1", annot_unique$matching_rna_files)

# Load summarized experiments
rna_se <- readRDS("/icbb/projects/share/datasets/blueprint/rna_se.rds")

# Filter for samples in annot (i just took the last one if there are several samples per id but idk)
keep <- colData(rna_se)$id2 %in% rna_ids
rna_se <- rna_se[,keep]

# Expression matrix (fpkm or counts?)
rna_fpkm <- assay(rna_se, "fpkm")

# Link gene_id to gene names
ids2gene <-as.data.frame(rowData(rna_se)[,c("gene_name", "gene_id")])
ids2gene$gene_id_short <- rownames(lookup)


## Dorothea
# Get the dorothea regulons
dorothea_hs <- decoupleR::get_dorothea(levels = c('A', 'B', 'C', 'D'))

# Keep only high-confidence interactions
dorothea_regulon <- subset(dorothea_hs, confidence %in% c("A", "B", 'C', 'D'))
dorothea_regulon$target <- toupper(dorothea_regulon$target)
dorothea_regulon$tf <- toupper(dorothea_regulon$source)

# Merge to get tracking_ids for DoRothEA targets
dorothea_with_tracking <- left_join(dorothea_regulon, ids2gene, 
                                    by = c("target" = "gene_name"))

# Remove entries with no match
dorothea_with_tracking <- dorothea_with_tracking[!is.na(dorothea_with_tracking$gene_id_short), ]

# TODO: Viper on rna
# TODO: Load blueprint deviations scores and subsett that it contains same samples as rna
# TODO: Heatmap rna + mtfr
# TODO: row correlation


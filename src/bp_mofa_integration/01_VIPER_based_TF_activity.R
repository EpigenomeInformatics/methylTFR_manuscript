#!/usr/bin/env Rscript

#####################################################################
# 01_VIPER_based_TF_activity.R
# created on 29.08.2026 
# Creating TF activity based on gene expression data using VIPER
####################################################################### 
suppressPackageStartupMessages({
  library(dorothea)
  library(dplyr)
  library(stringr)
  library(SummarizedExperiment)
  library(viper)
})

set.seed(13)

## Read in RNA files
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

# Expression matrix
rna_fpkm <- assay(rna_se, "fpkm")

# Link gene_id to gene names
ids2gene <-as.data.frame(rowData(rna_se)[,c("gene_name", "gene_id")])
ids2gene$gene_id_short <- sub("\\..*$", "", ids2gene$gene_id)


## Dorothea
# Get the dorothea regulons
data(dorothea_hs, package = "dorothea")

# Keep only high-confidence interactions
dorothea_regulon <- subset(dorothea_hs, confidence %in% c("A", "B", 'C', 'D'))
dorothea_regulon$target <- toupper(dorothea_regulon$target)
dorothea_regulon$tf <- toupper(dorothea_regulon$tf)
viper_regulon <- df2regulon(dorothea_regulon)

# Merge to get tracking_ids for DoRothEA targets
dorothea_with_tracking <- left_join(dorothea_regulon, ids2gene, 
                                    by = c("target" = "gene_name"))

# Remove entries with no match
dorothea_with_tracking <- dorothea_with_tracking[!is.na(dorothea_with_tracking$gene_id_short), ]

## Prepare expression matrix for VIPER
rna_fpkm <- as.data.frame(rna_fpkm)
# Gene id as own column
rna_fpkm$gene_id_short <- rownames(rna_fpkm)
# Match gene ids with corresponding gene names
rna_fpkm_with_genes <- left_join(rna_fpkm, ids2gene[, c("gene_id_short", "gene_name")], by = "gene_id_short")

# Take median of IDs with same gene_name
rna_fpkm_clean <- rna_fpkm_with_genes %>%
  group_by(gene_name) %>%
  summarise(across(where(is.numeric), function(x) median(x, na.rm = TRUE))) %>%
  ungroup()

# Set rownames to gene_name column and remove it
rna_mat <- as.matrix(rna_fpkm_clean[,-1])
rownames(rna_mat) <- rna_fpkm_clean$gene_name

# Log transformation
rna_mat_log <- log2(rna_mat +1)

## Viper analysis
tf_activity <- viper(rna_mat_log, viper_regulon, method="scale", minsize=10)

saveRDS(tf_activity,"/icbb/projects/nitschre/methylTFR/r_objects/viper_activity_bp.RDS")
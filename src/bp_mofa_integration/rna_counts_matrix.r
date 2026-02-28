## Heatmap direct TF expression of bp
library(dplyr)
library(stringr)
library(SummarizedExperiment)
library(tibble)

set.seed(12)
plot_dir <- "/icbb/projects/nitschre/methylTFR/figures/blueprint/Heatmaps_differentials/"

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
rna_ids <- sapply(rna_ids_list, `[`, 1)

# Load summarized experiments
rna_se <- readRDS("/icbb/projects/share/datasets/blueprint/rna_se.rds")

# Filter for samples in annot (i just took the last one if there are several samples per id but idk)
keep <- colData(rna_se)$id2 %in% rna_ids
rna_se <- rna_se[,keep]

# Expression matrix 
rna <- assay(rna_se, "count")

# Link gene_id to gene names
ids2gene <-as.data.frame(rowData(rna_se)[,c("gene_name", "gene_id")])
ids2gene$gene_id_short <- rownames(ids2gene)

# Annotate gene names
rna <- rna %>%
  as.data.frame() %>%
  rownames_to_column("gene_id_short") %>%
  left_join(ids2gene, by="gene_id_short")

# Duplicate gene names replaced with mean gene expression
rna <- rna %>%
  group_by(gene_name) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
  as.data.frame()

# Gene ids as rownames
rownames(rna) <- rna$gene_name
rna$gene_name <- NULL

## Change sample names
annotation_with_ids <- read.csv("methylTFR/sample_annotation/full_annotation.csv", stringsAsFactors=FALSE)

#Indices of ids in table
indices <- match(colData(rna_se)$epirrId, annotation_with_ids$epirr_id_without_version)

# Extract corresponding names
colnames(rna) <- annotation_with_ids$cellTypeGroup[indices]

saveRDS(rna, "/icbb/projects/nitschre/methylTFR/r_objects/bp_rawRNAcounts.RDS")
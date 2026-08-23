#!/usr/bin/env Rscript

#####################################################################
# mfoot_echo_cell_type_jaspar.R
# created on 25-10-25 by Irem Gunduz
# Plot footprints for the JASPAR motifs for mapped cell-types
#####################################################################

# load libraries
suppressPackageStartupMessages({
  library(ArchR)
  library(dplyr)
  library(muLogR)
  library(GenomicRanges)
  library(BSgenome.Hsapiens.UCSC.hg38)
})
cell_type_colors <- c(
  "B-cell" = "#AE017E",
  "Monocyte" = "#CC4C02",
  "NK-cell" = "#A65628",
  "Th-Mem" = "#41B6C4",
  "Tc-Mem" = "#4292C6",
  "Tc-Naive" = "#888FB5",
  "Th-Naive" = "#C7E9B4"
)
set.seed(12) 
outputDir <- "/icbb/projects/igunduz/archr_projects/icbb/projects/igunduz/archr_project_011023"
project <- ArchR::loadArchRProject(outputDir, showLogo = FALSE)
addArchRThreads(threads = 30) # set the cores

# Get the overlapped samples between scATAC-seq and scWGBS
annot_acc <- readRDS("/icbb/projects/igunduz/DARPA_analysis/artemis_031023/cell_matching/annot_acc_cellType.rds")
cell_annot_atac <- read.delim("/icbb/projects/igunduz/DARPA_analysis/artemis_031023/sample_annot_atac.tsv")
rownames(cell_annot_atac) <- cell_annot_atac$cellId_archr
annot_acc <- readRDS("/icbb/projects/igunduz/DARPA_analysis/artemis_031023/cell_matching/annot_acc_cellType.rds")
annot_acc$archr_id <- rownames(annot_acc)
new_annot <- merge(cell_annot_atac, annot_acc, by = "row.names")
new_annot$atac_pb_id <- paste0(new_annot$class, "_", new_annot$sample_sampleId_cminid)

# Subset the cells
idxSample <- BiocGenerics::which(project$cellNames %in% new_annot$archr_id)
cellsSample <- project$cellNames[idxSample]
project <- project[cellsSample, ]


# Add mapped cell-type information
project <- addCellColData(
  ArchRProj = project, data = new_annot$atac_pb_id,
  name = "atac_pb_id", cells = new_annot$archr_id
)
project <- addCellColData(
  ArchRProj = project, data = new_annot$mapped_cell_class,
  name = "cell_type_meth", cells = new_annot$archr_id, force = TRUE
)

# Remove the Other-cell
idxSample <- BiocGenerics::which(project$cell_type_meth %in% "Other-cell")
cellsSample <- project$cellNames[-idxSample]
project <- project[cellsSample, ]

# get motif positions
motifPositions <- getPositions(project, name = "jaspar2020")

# Motifs from MOFA2 analysis
markerMotifs <- c(
  "OLIG2", "NEUROG1", "HES1", "HES6", "HEY2", "HEY1", "KLF17", "KLF16",
  "TFAP2A", "ZNF652", "TCF21(var.2)", "TFAP4", "MSC", "TBXT", "OSR2",
  "ETV5", "ELF4", "ELF5", "SPIC", "SPIB", "SPI1", "EHF", "PRDM4", "GFI1",
  "MXI1", "ZBTB14", "E2F2", "FOXN3", "ZFP57", "ESR1", "TCF4", "ASCL1(var.2)",
  "EBF3", "POU2F3", "POU2F2", "POU5F1", "POU3F4", "LHX6", "EMX2", "TBX21",
  "TBX2", "TBR1", "EOMES", "TBX20", "TBX1", "MGA", "HIC2"
)

# add group coverages
project <- addGroupCoverages(ArchRProj = project, groupBy = "cell_type_meth")

# get footprints
seFoot <- getFootprints(
  ArchRProj = project,
  positions = motifPositions[markerMotifs],
  groupBy = "cell_type_meth"
)

plotFootprints(
  seFoot = seFoot,
  ArchRProj = project,
  normMethod = "divide",
  pal = cell_type_colors,
  plotName = "Integrative_Footprints_JASPAR_by_cellType",
  addDOC = FALSE
)

#####################################################################

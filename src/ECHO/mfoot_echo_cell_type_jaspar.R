
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

set.seed(12) # set seed
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
  ArchRProj = project, data = new_annot$class,
  name = "cell_type_meth", cells = new_annot$archr_id
)

# get motif positions
motifPositions <- getPositions(project,name="jaspar2020")

# Get the marker motifs
markerMotifs <- names(motifPositions)

# add group coverages
project <- addGroupCoverages(ArchRProj = project, groupBy = "cell_type_meth")

# get footprints
seFoot <- getFootprints(
  ArchRProj = project, 
  positions = motifPositions, 
  groupBy = "cell_type_meth"
)

plotFootprints(
  seFoot = seFoot,
  ArchRProj = project, 
  normMethod = "divide",
  plotName = "Integrative_Footprints_JASPAR_by_cellType",
  addDOC = FALSE
)

#####################################################################
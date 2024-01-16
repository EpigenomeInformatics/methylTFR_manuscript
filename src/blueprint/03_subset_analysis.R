set.seed(42)
muLogR::logger.info("Loading libraries...")
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(logger)
  library(muLogR)
  library(RnBeads)
})
sample_dir <- "/icbb/projects/igunduz/methylTFR_manuscript/results/BLUEPRINT/"
logger.info("Loading RnBeads object...")
rnb_set <- RnBeads::load.rnb.set(paste0(sample_dir,"reports/data_import_data/rnb.set_preprocessed"))
report.dir <- "/icbb/projects/igunduz/methylTFR_manuscript/results/BLUEPRINT/"

# Multiprocess
parallel.setup(30)

# Set options
rnb.options(
  analysis.name = "Blueprint Bcell VS Tcell",
  assembly = "hg38",
  import.table.separator = "\t",
  region.aggregation = "sum",
  region.types = c("cpgislands"),
  import.default.data.type = "data.dir",
  import.bed.style = "EPP",
  analyze.sites = FALSE,
  disk.dump.big.matrices = TRUE,
  strand.specific = FALSE,
  filtering.sex.chromosomes.removal = TRUE,
  differential.enrichment.lola = TRUE,
  differential.enrichment.lola.dbs = "/icbb/projects/share/annotations/lolaDB",
  identifiers.column = "bedFile",
  differential.comparison.columns = "cellTypeShort" #  exclusive cell-types
)


#Subset to only include B cells
idx <- rnb_set@pheno[rnb_set@pheno$cellTypeGroup != "Bcell", ]$bedFile
idxt <- rnb_set@pheno[rnb_set@pheno$cellTypeGroup == "Bcell", ]$bedFile
rnb_set_bcells <- remove.samples(rnb_set,idx)
rnb_set_tcells <- remove.samples(rnb_set,idxt)

# Set report directories
breport.dir <- paste0(report.dir,"/bcell")
treport.dir <- paste0(report.dir,"/tcell")
if(!dir.exists(breport.dir)){dir.create(breport.dir)}
if(!dir.exists(treport.dir)){dir.create(treport.dir)}

# Differential methylation
rnb.run.differential(rnb_set_bcells, breport.dir)
rnb.run.differential(rnb_set_tcells, treport.dir)
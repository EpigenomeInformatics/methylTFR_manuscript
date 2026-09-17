#!/usr/bin/env Rscript
# 06_scECHO_plot.R
# Created on 15-09-26 by Irem Begum Gunduz
# Merged scECHO (06) and ATAC (07) Plotting Script
# Includes optional MAGIC imputation for methylTFR, caching for ATAC/ArchR,
# and a unified side-by-side patchwork layout.

suppressPackageStartupMessages({
  library(ArchR)
  library(SummarizedExperiment)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(logger)
  library(methylTFR)
  library(Rmagic) # For optional MAGIC imputation
})
set.seed(42)

# Point reticulate at the Python that has magic-impute, BEFORE any Python is
# initialised. Set python.bin to the interpreter you pip-installed magic-impute
# into (e.g. `which python` inside the conda env). Only needed if run.magic = TRUE.
python.bin <- "/icbb/projects/share/software/packages/miniconda3/envs/methyltfr/bin/python"
if (file.exists(python.bin)) {
  Sys.setenv(RETICULATE_PYTHON = python.bin)
  try(suppressWarnings(reticulate::use_python(python.bin, required = FALSE)), silent = TRUE)
}
addArchRThreads(threads = 30)

# Rasterise the point clouds so the figure stays small. Prefer per-layer raster
# (ggrastr, or the lighter scattermore: axes/text stay vector). If neither is
# installed, fall back to writing the whole figure as a high-res PNG (fully
# raster, no extra packages).
raster.method <- if (requireNamespace("ggrastr", quietly = TRUE)) {
  "ggrastr"
} else if (requireNamespace("scattermore", quietly = TRUE)) {
  "scattermore"
} else {
  "png"
}
log_info("Rasterisation method: ", raster.method)

#####################################################################
# Configuration & Settings
#####################################################################

# -- General --
features <- c("POU2F3", "SPIB", "EOMES")
point.size <- 0.25
raster.dpi <- 300
github.dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR_manuscript"
out.dir <- file.path(github.dir, "figures", "echo", "scECHO_230826_merged")
if (!dir.exists(out.dir)) dir.create(out.dir, recursive = TRUE)

# -- Methylation (scECHO) Settings --
motifSet <- "jaspar2020_distal"
dev.dir <- file.path("/scratch/icbb/igunduz/methylTFR_manuscript/echo", paste0("mTFR_sc_150926_", motifSet))
lsi.file <- "/icbb/projects/igunduz/DARPA_analysis/artemis_031023/itLSI_res/itLSI_res_sub11kpc30.rds"
cell.id.col <- "Cell_UID"
pal.low <- "#5E3C99"; pal.mid <- "grey92"; pal.high <- "#E66101"
clip.quantile <- 0.98

# -- MAGIC imputation (methylTFR) --
# OFF by default. MAGIC smooths values over a cell-cell graph it BUILDS FROM the
# feature matrix it is given. These deviation files carry only a couple of motifs
# (the run log shows "... on N cells and 2 genes" + a graphtools zero-distance
# warning), which is far too few to build a meaningful graph -- MAGIC then just
# collapses everything toward the mean and the panel looks worse, not better.
# Enable this ONLY if dev.mat below holds the FULL motif set (i.e. you recomputed
# methylTFR deviations over all JASPAR motifs) AND Python magic-impute is on the
# interpreter above. With the current thin deviation files, leave it FALSE.
run.magic  <- FALSE
magic.knn  <- 5      # neighbours for the affinity graph
magic.t    <- 3      # diffusion time; explicit value beats "auto" for sparse data
magic.npca <- 50     # PCA dims before the graph (clamped to n_motifs - 1)

# -- ATAC (ArchR) Settings --
outputDir <- "/icbb/projects/igunduz/archr_projects/icbb/projects/igunduz/archr_project_011023/"
embedding <- "UMAPHarmony"
gene.matrix <- "GeneScoreMatrix"
motif.matrix <- "jaspar2020Matrix"
gene.pal <- ArchRPalettes[["horizonExtra"]]
motif.pal <- ArchRPalettes[["solarExtra"]]
gene.quant <- c(0.01, 0.99)
motif.quant <- c(0.02, 0.98)

# Cache file to prevent reloading ArchR repeatedly
atac.cache.file <- file.path(out.dir, "atac_cache.rds")

# -- Color Palettes --
ClusterCellTypes_colors = c(
  B_mem = "#AE017E", B_naive = "#F768A1", DC = "#67000D", Mono_CD14 = "#FE9929",
  Mono_CD16 = "#CC4C02", NK_CD16 = "#A65628", Plasma = "#A106BD", T_mait = "#41B6C4",
  T_mem_CD4 = "#4292c6", T_mem_CD8 = "#0074cc", T_mix = "#888FB5", T_naive = "#C7E9B4"
)
cell_type_colors = c(
  "B-cell" = "#AE017E", "Monocyte" = "#CC4C02", "NK-cell" = "#A65628",
  "Th-Mem" = "#41B6C4", "Tc-Mem" = "#4292C6", "Tc-Naive" = "#888FB5",
  "Th-Naive" = "#C7E9B4", "Other-cell" = "#CCCCCC"
)

# Shared Raster point layer helper
point_layer <- function(data, aes_mapping, size = point.size, ...) {
  if (raster.method == "ggrastr") {
    return(ggrastr::rasterise(
      ggplot2::geom_point(data = data, mapping = aes_mapping, size = size, ...),
      dpi = raster.dpi))
  }
  if (raster.method == "scattermore") {
    return(scattermore::geom_scattermore(
      data = data, mapping = aes_mapping, pointsize = max(1, size * 4), ...))
  }
  # PNG fallback: points stay vector in the object; the whole figure is saved raster
  ggplot2::geom_point(data = data, mapping = aes_mapping, size = size, ...)
}

#####################################################################
# 1. ATAC Data Processing (ArchR with Caching)
#####################################################################

atac_data <- if (file.exists(atac.cache.file)) readRDS(atac.cache.file) else NULL
if (!is.null(atac_data) && !is.null(atac_data$motif)) {
  log_info("Loading cached ATAC data from ", atac.cache.file)
} else {
  if (!is.null(atac_data)) log_warn("Cached ATAC has no chromVAR; rebuilding the cache")
  log_info("Loading ArchR project and extracting matrices (will be cached)...")

  load_archr_lenient <- function(path) {
    p <- normalizePath(path, mustWork = TRUE)
    proj <- ArchR::recoverArchRProject(readRDS(file.path(p, "Save-ArchR-Project.rds")))
    arrows <- file.path(p, "ArrowFiles", basename(proj@sampleColData$ArrowFiles))
    proj@sampleColData$ArrowFiles <- arrows
    proj@projectMetadata$outputDirectory <- p
    proj
  }

  project <- tryCatch(
    ArchR::loadArchRProject(outputDir, showLogo = FALSE),
    error = function(e) { load_archr_lenient(outputDir) }
  )

  subset.cell.types <- c("B_mem", "B_naive", "T_naive", "T_mem_CD4", "T_mem_CD8", "NK_CD16", "Mono_CD14", "Mono_CD16")
  cd <- getCellColData(project, select = "ClusterCellTypes")
  keep.cells <- rownames(cd)[cd$ClusterCellTypes %in% subset.cell.types]
  project <- ArchR::subsetCells(project, cellNames = keep.cells)

  # Get UMAP & Merge Cell Types
  umap_atac <- as.data.frame(getEmbedding(project, embedding = embedding, returnDF = TRUE))[, 1:2]
  colnames(umap_atac) <- c("UMAP1", "UMAP2")
  umap_atac$cell <- rownames(umap_atac)
  umap_atac$ClusterCellTypes <- as.character(cd[rownames(umap_atac), "ClusterCellTypes"])

  # Extract matrices (rows subset to `features`); z_scores keeps the chromVAR z
  # rows (assay "z", else rowData$seqnames == "z").
  get_values <- function(useMatrix, wanted, z_scores = FALSE) {
    se <- tryCatch(getMatrixFromProject(project, useMatrix = useMatrix),
      error = function(e) { log_warn(useMatrix, ": ", conditionMessage(e)); NULL })
    if (is.null(se)) return(NULL)

    if (z_scores) {
      if ("z" %in% assayNames(se)) {
        m <- assay(se, "z"); rownames(m) <- rowData(se)$name
      } else if ("seqnames" %in% colnames(rowData(se))) {
        se <- se[rowData(se)$seqnames == "z", ]; m <- assay(se); rownames(m) <- rowData(se)$name
      } else {
        m <- assay(se); rownames(m) <- rowData(se)$name
      }
    } else {
      m <- assay(se); rownames(m) <- rowData(se)$name
    }
    colnames(m) <- colnames(se)


    hit <- intersect(wanted, rownames(m))
    if (length(hit) == 0) { log_warn(useMatrix, ": none of the features found"); return(NULL) }
    m[hit, , drop = FALSE]  # few features x all cells

  }

  gene.vals <- get_values(gene.matrix, features, z_scores = FALSE)
  motif.vals <- get_values(motif.matrix, features, z_scores = TRUE)

  # Impute on the cells shared by both matrices and the embedding. imputeMatrix
  # needs every impute-weight cell present in the matrix, and getMatrixFromProject
  # can return a slightly different cell set, so the weights are rebuilt on the
  # shared cells. Wrapped in a fallback so a residual mismatch degrades to the
  # raw (un-smoothed) values instead of aborting the run.
  cell.sets <- list(umap_atac$cell)
  if (!is.null(gene.vals)) cell.sets <- c(cell.sets, list(colnames(gene.vals)))
  if (!is.null(motif.vals)) cell.sets <- c(cell.sets, list(colnames(motif.vals)))
  common.cells <- Reduce(intersect, cell.sets)
  log_info(length(common.cells), " cells shared by the embedding and the matrices")

  proj.i <- ArchR::subsetCells(project, cellNames = common.cells)
  proj.i <- addImputeWeights(proj.i)
  impute.weights <- getImputeWeights(proj.i)

  smooth_vals <- function(m) {
    if (is.null(m)) return(NULL)
    m <- m[, common.cells, drop = FALSE]
    tryCatch(imputeMatrix(mat = m, imputeWeights = impute.weights),
      error = function(e) { log_warn("imputeMatrix failed (", conditionMessage(e), "); using raw values"); m })
  }
  gene.vals <- smooth_vals(gene.vals)
  motif.vals <- smooth_vals(motif.vals)

  atac_data <- list(umap = umap_atac, gene = gene.vals, motif = motif.vals)
  saveRDS(atac_data, atac.cache.file)
  log_success("Saved ATAC cache to ", atac.cache.file)
}

#####################################################################
# 2. Methylation Data Processing (scECHO)
#####################################################################

log_info("Loading Methylation (LSI) data...")
res <- readRDS(lsi.file)
umap_meth <- as.data.frame(res$umapCoord)[, 1:2]
colnames(umap_meth) <- c("UMAP1", "UMAP2")
umap_meth[[cell.id.col]] <- rownames(umap_meth)

# Methylation cell types: the LSI result has none, so map Cell_UID -> cell_type
# from the sample annotation used by 05_scECHO.R.
sannot.file <- "/icbb/projects/igunduz/DARPA_analysis/artemis_031023/sample_annot.tsv"
umap_meth$cell_type <- "Other-cell"
if (file.exists(sannot.file)) {
  sa <- utils::read.delim(sannot.file, stringsAsFactors = FALSE)
  if (all(c("Cell_UID", "cell_type") %in% colnames(sa))) {
    umap_meth$cell_type <- sa$cell_type[match(umap_meth[[cell.id.col]], sa$Cell_UID)]
    n.lab <- sum(!is.na(umap_meth$cell_type))
    log_info("Mapped cell_type for ", n.lab, " of ", nrow(umap_meth), " methylation cells")
    umap_meth$cell_type[is.na(umap_meth$cell_type)] <- "Other-cell"
  } else {
    log_warn("Cell_UID/cell_type columns not in ", sannot.file, "; using dummy labels")
  }
} else {
  log_warn("Sample annotation not found (", sannot.file, "); using dummy labels")
}

# Load Deviations
dev.files <- list.files(dev.dir, pattern = "_deviations\\.RDS$", full.names = TRUE)
mats <- lapply(dev.files, function(f) as.matrix(methylTFR::deviations(readRDS(f))))
motif.rows <- Reduce(intersect, lapply(mats, rownames))
dev.mat <- do.call(cbind, lapply(mats, function(m) m[motif.rows, , drop = FALSE]))
dev.mat <- dev.mat[, !duplicated(colnames(dev.mat)), drop = FALSE]

# methylTFR emits Inf/NaN for motifs uncovered in a cell. Keep these as NA so they
# are dropped from the SD and simply not plotted; forcing them to 0 would pull
# real signal toward grey. (The MAGIC branch below 0-fills its own input copy.)
dev.mat[!is.finite(dev.mat)] <- NA

#####################################################################
# 2b. MAGIC imputation (optional; OFF by default -- see config note)
#
# Runs MAGIC on the full dev.mat (all `motif.rows` x cells) so the affinity graph
# is built from the full methylation feature space, then returns only `features`.
# With the current thin deviation files this graph is degenerate (see config), so
# run.magic is FALSE and this block is skipped -- dev.mat passes through unchanged.
#####################################################################

magic.done <- FALSE
if (run.magic) {
  py.ok <- requireNamespace("reticulate", quietly = TRUE) &&
    tryCatch(reticulate::py_module_available("magic"), error = function(e) FALSE)

  if (!py.ok) {
    log_warn("Python 'magic' module not importable for reticulate; plotting RAW ",
             "deviations. Install into the interpreter above with:")
    log_warn("    ", python.bin, " -m pip install magic-impute")
  } else {
    input <- t(dev.mat)                                   # cells x motifs (full space)
    input[!is.finite(input)] <- 0                         # MAGIC needs a finite matrix
    if (is.null(rownames(input))) rownames(input) <- colnames(dev.mat)
    wanted <- intersect(features, colnames(input))
    npca_cfg <- if (exists("magic.npca")) magic.npca else 50L   # tolerate missing config
    npca <- max(2L, min(npca_cfg, ncol(input) - 1L))
    log_info("MAGIC input: ", nrow(input), " cells x ", ncol(input),
             " motifs; imputing: ", paste(wanted, collapse = ", "),
             " (knn=", magic.knn, ", t=", magic.t, ", npca=", npca, ")")
    if (ncol(input) < 10L)
      log_warn("Only ", ncol(input), " motifs in dev.mat -- MAGIC's graph will be ",
               "degenerate. See the run.magic note in the config.")

    imp <- tryCatch({
      mg <- Rmagic::magic(input, genes = wanted,
                          knn = magic.knn, t = magic.t, npca = npca)
      m <- t(as.matrix(mg$result))                        # wanted x cells
      m[, rownames(input), drop = FALSE]                  # align cells BY NAME
    }, error = function(e) {
      log_warn("MAGIC failed (", conditionMessage(e), "); plotting RAW deviations")
      NULL
    })

    if (!is.null(imp)) {
      for (g in rownames(imp)) {
        raw.sd <- stats::sd(dev.mat[g, ], na.rm = TRUE)
        imp.sd <- stats::sd(imp[g, ], na.rm = TRUE)
        log_info(sprintf("  %-8s raw sd=%.4g  ->  MAGIC sd=%.4g", g, raw.sd, imp.sd))
      }
      dev.mat <- imp        # now holds ONLY the plotted motifs, MAGIC-smoothed
      magic.done <- TRUE
      log_success("MAGIC imputation applied to methylTFR deviations")
    }
  }
}

# Per-motif scaling (NOT a full z-score): divide each motif by its SD but DO NOT
# subtract the mean. methylTFR deviations have a meaningful zero (= no methylation
# preference), so grey must sit at zero. Subtracting the per-motif mean would
# recenter grey on the AVERAGE cell, which paints a broadly-negative motif like
# SPIB orange and flips the whole map. Scaling only keeps the good, zero-centred
# look while putting every motif on a comparable unit-SD scale.
row_zscore <- function(m) {
  sdv <- apply(m, 1, stats::sd, na.rm = TRUE)
  sdv[!is.finite(sdv) | sdv == 0] <- 1
  z <- m / sdv
  dimnames(z) <- dimnames(m)
  z
}
dev.mat <- row_zscore(dev.mat)
plot_features <- intersect(features, rownames(dev.mat))
if (length(plot_features) == 0)
  log_warn("None of `features` present in dev.mat after processing")

#####################################################################
# 3. Plotting Setup & Functions
#####################################################################

theme_umap <- theme_classic(base_size = 8) + theme(
  plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
  axis.text = element_blank(), axis.ticks = element_blank(),
  legend.key.width = unit(0.25, "cm"), legend.key.height = unit(0.5, "cm"),
  legend.text = element_text(size = 7),
  plot.background = element_rect(fill = "transparent", colour = NA),
  panel.background = element_rect(fill = "transparent", colour = NA),
  legend.background = element_rect(fill = "transparent", colour = NA)
)

# --- Cell Type Plots ---
p_meth_ct <- ggplot() +
  point_layer(umap_meth, aes(UMAP1, UMAP2, colour = cell_type), size = point.size) +
  scale_colour_manual(values = cell_type_colors, na.value = "#CCCCCC", name = "Cell Type") +
  labs(title = "Methylation Cell Types") + coord_fixed() + theme_umap

p_atac_ct <- ggplot() +
  point_layer(atac_data$umap, aes(UMAP1, UMAP2, colour = ClusterCellTypes), size = point.size) +
  scale_colour_manual(values = ClusterCellTypes_colors, na.value = "#CCCCCC", name = "ATAC Cell Types") +
  labs(title = "ATAC Cell Types") + coord_fixed() + theme_umap

# --- Value UMAP Plotter (Used for all 3 metrics) ---
plot_val_umap <- function(df, val_col, title, palette, limits, symmetric = FALSE) {
  if (is.null(df)) return(plot_spacer())
  df <- df[!is.na(df[[val_col]]), ]            # never plot NA cells (grey) on top
  if (nrow(df) == 0) return(plot_spacer())
  df <- df[order(abs(df[[val_col]])), ]
  p <- ggplot() +
    point_layer(df, aes(UMAP1, UMAP2, colour = .data[[val_col]]), size = point.size) +
    labs(title = title, x = "UMAP1", y = "UMAP2") +
    coord_fixed() + theme_umap

  if (symmetric) {
    p <- p + scale_colour_gradient2(
      low = palette[1], mid = palette[2], high = palette[3], midpoint = 0,
      limits = limits, oob = scales::squish, name = "Score"
    )
  } else {
    p <- p + scale_colour_gradientn(
      colours = palette, limits = limits, oob = scales::squish, name = "Score"
    )
  }
  p
}

#####################################################################
# 4. Constructing the Unified Layout
#####################################################################

# Row 1: The two Cell Type Plots + Spacer (to align with 3-panel Motif rows)
row1 <- wrap_plots(p_meth_ct, p_atac_ct, plot_spacer(), nrow = 1)
all_rows <- list(row1)

for (gene in plot_features) {

  # A. ATAC chromVAR
  cvar_df <- if (!is.null(atac_data$motif) && gene %in% rownames(atac_data$motif)) {
    data.frame(UMAP1 = atac_data$umap$UMAP1, UMAP2 = atac_data$umap$UMAP2,
               val = atac_data$motif[gene, match(atac_data$umap$cell, colnames(atac_data$motif))])
  } else NULL
  if (!is.null(cvar_df)) {
    lim <- max(abs(stats::quantile(cvar_df$val, motif.quant, na.rm = TRUE)))
    p_cvar <- plot_val_umap(cvar_df, "val", paste(gene, "chromVAR"), motif.pal, c(-lim, lim), FALSE)
  } else p_cvar <- plot_spacer()

  # B. MethylTFR
  score <- rep(NA_real_, nrow(umap_meth))
  hit <- match(umap_meth[[cell.id.col]], colnames(dev.mat))
  score[!is.na(hit)] <- dev.mat[gene, hit[!is.na(hit)]]
  meth_df <- data.frame(UMAP1 = umap_meth$UMAP1, UMAP2 = umap_meth$UMAP2, val = score)
  lim_meth <- stats::quantile(abs(meth_df$val), clip.quantile, na.rm = TRUE)
  if (!is.finite(lim_meth) || lim_meth == 0) lim_meth <- max(abs(meth_df$val), na.rm = TRUE)
  p_mtfr <- plot_val_umap(meth_df[!is.na(meth_df$val), ], "val", paste(gene, "methylTFR"),
                          c(pal.low, pal.mid, pal.high), c(-lim_meth, lim_meth), TRUE)

  # C. ATAC Gene Score
  gene_df <- if (!is.null(atac_data$gene) && gene %in% rownames(atac_data$gene)) {
    data.frame(UMAP1 = atac_data$umap$UMAP1, UMAP2 = atac_data$umap$UMAP2,
               val = atac_data$gene[gene, match(atac_data$umap$cell, colnames(atac_data$gene))])
  } else NULL
  if (!is.null(gene_df)) {
    g_lim <- stats::quantile(gene_df$val, gene.quant, na.rm = TRUE)
    p_gene <- plot_val_umap(gene_df, "val", paste(gene, "Gene Score"), gene.pal, g_lim, FALSE)
  } else p_gene <- plot_spacer()

  # Combine row: chromVAR | methylTFR | Gene Score
  row_m <- wrap_plots(p_cvar, p_mtfr, p_gene, nrow = 1)
  all_rows[[gene]] <- row_m
}

# Wrap all rows into one master plot vertically
final_grid <- wrap_plots(all_rows, ncol = 1)

grid.file <- file.path(out.dir, "merged_scECHO_ATAC_grid.pdf")
if (raster.method == "png") {
  grid.file <- sub("\\.pdf$", ".png", grid.file)
  ggsave(grid.file, final_grid, width = 14, height = 4 * length(all_rows),
    limitsize = FALSE, bg = "transparent", dpi = raster.dpi)
} else {
  ggsave(grid.file, final_grid, width = 14, height = 4 * length(all_rows),
    limitsize = FALSE, bg = "transparent")
}
log_success("Successfully plotted and saved merged visualization: ", grid.file)
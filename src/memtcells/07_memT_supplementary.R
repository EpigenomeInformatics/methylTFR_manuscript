#!/usr/bin/env Rscript

#####################################################################
# 07_memT_supplementary.R
# created on 06-09-2026 by Irem B Gunduz
# The memory T cell supplementary figure
#   A  Mixed footprints (Obs vs Exp for JUN/FOSL2, Diff for BATF/SPI1)
#   B  Scatter plots of methylTFR vs TF expression (TCM vs TN, TEM vs TN)
#   C  Scatter plots of methylTFR vs VIPER activity (TCM vs TN, TEM vs TN)
#   D  Boxplots of cross-modality concordance across cell types with 
#      comparison brackets and exact padj values.
#
# LOLA and MA plots removed to conserve memory.
#####################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(GenomicRanges)
  library(logger)
  library(RnBeads)
  library(methylTFR)
  library(methylTFRAnnotationHg38)
  library(org.Hs.eg.db)
})
set.seed(13)

src.dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR_manuscript/src"
source(file.path(src.dir, "utils.R"))

#####################################################################
# Settings
#####################################################################

motifSet <- "jaspar2020_distal"
tfSet <- "jaspar2020"

comparisons <- list(
  CM = list(test = "TCM", ref = "TN"),
  EM = list(test = "TEM", ref = "TN")
)
drop.cell.types <- "TEMRA"
group_order <- c("TN", "TCM", "TEM")

# Color matching the provided screenshot
cell_type_colors <- c(
  "TN" = "#C5E1A5",  # Light Green
  "TCM" = "#2CA25F", # Medium Green 
  "TEM" = "#005A32"  # Dark Green
)

# Colors for Expected vs Observed footprints
tint <- function(hex, amount = 0.55) {
  rgb_vals <- grDevices::col2rgb(hex)[, 1]
  mixed <- rgb_vals + (255 - rgb_vals) * amount
  grDevices::rgb(mixed[1], mixed[2], mixed[3], maxColorValue = 255)
}
base_colors <- unlist(lapply(names(cell_type_colors), function(ct) {
  setNames(
    c(tint(cell_type_colors[[ct]]), unname(cell_type_colors[[ct]])),
    paste0(c("Expected_", "Observed_"), ct)
  )
}))

# Canvas
fig.width <- 16
fig.height <- 22
base.size <- 12
flank.norm <- 50

distal.file <- "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFRAnnotationHg38_old/inst/extdata/distal_regions.RDS"

# Directories
analysis.dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/memoryTcells"
dev.tag <- "mTFR_devs_230826"
dev.file <- file.path(analysis.dir, dev.tag, paste0(motifSet, "_deviations.RDS"))
rnb.set.path <- file.path(analysis.dir, "reports", "data_import_data", "rnb.set_preprocessed")

# RNA & VIPER data paths
rna.dir <- "/icbb/projects/share/datasets/memoryTcells/rna_raw_081025"
meth.dir <- "/icbb/projects/share/datasets/memoryTcells"
viper.out <- file.path(analysis.dir, "memT_viper_activity.RDS")

cache.dir <- file.path(analysis.dir, "debug")
if (!dir.exists(cache.dir)) dir.create(cache.dir, recursive = TRUE)
cache.file <- file.path(cache.dir, "msites_merged_cellType_230826.rds")

github.dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR_manuscript"
plot.dir <- file.path(github.dir, "figures", "memtcells", "diff_TFs_230826")
if (!dir.exists(plot.dir)) dir.create(plot.dir, recursive = TRUE)

# Enforce perfect transparency across patchwork and overall device
no_bg <- theme(
  plot.background = element_rect(fill = "transparent", colour = NA),
  panel.background = element_rect(fill = "transparent", colour = NA),
  legend.background = element_rect(fill = "transparent", colour = NA),
  legend.box.background = element_rect(fill = "transparent", colour = NA),
  legend.key = element_rect(fill = "transparent", colour = NA),
  strip.background = element_blank()
)
tag_theme <- theme(plot.tag = element_text(face = "bold", size = 16))

#####################################################################
# Helper functions
#####################################################################

ensure_msites_cols <- function(gr, label) {
  if (!"score" %in% names(mcols(gr))) {
    alt <- intersect(c("beta", "meth", "methylation", "avg_methyl"), names(mcols(gr)))
    mcols(gr)$score <- as.numeric(mcols(gr)[[alt[1]]])
  }
  if (!"coverage" %in% names(mcols(gr))) {
    alt <- intersect(c("covg", "cov", "reads"), names(mcols(gr)))
    mcols(gr)$coverage <- if (length(alt) > 0) as.numeric(mcols(gr)[[alt[1]]]) else NA_real_
  }
  gr
}

load_cell_type_deviations <- function(file) {
  if (!file.exists(file)) return(NULL)
  obj <- readRDS(file)
  cd <- as.data.frame(colData(obj), stringsAsFactors = FALSE)
  mat <- deviations(obj)
  split(seq_len(ncol(mat)), as.character(cd$cellType)) |>
    lapply(function(idx) mat[, idx, drop = FALSE])
}

motif_mean_deviations <- function(motif, dev_list, groups) {
  vapply(groups, function(g) {
    mat <- dev_list[[g]]
    idx <- which(rownames(mat) == motif)
    if (length(idx) == 0) return(NA_real_)
    vals <- mat[idx[1], ]
    if (all(is.na(vals))) NA_real_ else mean(vals, na.rm = TRUE)
  }, numeric(1))
}

make_label <- function(group, dev_val) {
  if (is.na(dev_val)) return(group)
  paste0(group, ": ", format(round(dev_val, 2), nsmall = 2))
}

clean_gex_names <- function(x) {
  x <- sub("\\..*$", "", x)
  x <- sub("^[0-9]+_", "", x)
  x <- sub("_Ct$", "", x)
  x <- sub("_BlCM$", "_TCM", x)
  x <- sub("_BlEM$", "_TEM", x)
  x <- sub("_BlTN$", "_TN", x)
  x <- sub("_BlTR$", "_TEMRA", x)
  x
}

sample_donors <- function(ids) {
  donors <- rep(NA_character_, length(ids))
  hit <- grepl("Hf[0-9]+", ids)
  donors[hit] <- sub(".*?(Hf[0-9]+).*", "\\1", ids[hit])
  donors
}

#####################################################################
# Panel A: Footprinting
#####################################################################

if (!file.exists(cache.file)) {
  log_info("Loading RnBeads object from ", rnb.set.path)
  rnb_set <- RnBeads::load.rnb.set(rnb.set.path)
  msites <- rnb.RnBSet.to.GRangesList(mergeSamples(rnb_set, "cellType"))
  saveRDS(msites, cache.file)
} else {
  msites <- readRDS(cache.file)
}

msites <- as.list(msites)
msites <- mapply(ensure_msites_cols, msites, names(msites), SIMPLIFY = FALSE)
msites <- msites[intersect(group_order, setdiff(names(msites), drop.cell.types))]

dev_list <- load_cell_type_deviations(dev.file)
tf_bindsites <- getTFbindsites(motifSet = tfSet)
gcfreqs <- getGCfreq(motifSet = motifSet)
gc_dist <- getGenomeGC()

enhancer <- NULL
if (motifSet == "jaspar2020_distal" && file.exists(distal.file)) {
  enhancer <- readRDS(distal.file)
}

footprint_panel <- function(motif, type = "diff") {
  if (!motif %in% names(tf_bindsites)) return(NULL)
  log_info("Generating footprint panel (", type, ") for motif: ", motif)
  
  groups_here <- names(msites)
  
  per_group <- lapply(groups_here, function(g) {
    df <- plotExpectedFootprint(
        motif = motif, tf_bindsites = tf_bindsites, msites = msites[[g]],
        sample_name = g, gc_dist = gc_dist, gcfreqs = gcfreqs,
        enhancer = enhancer, returnPlotData = TRUE
      )$plotDF
    if (is.null(df)) return(NULL)
    
    if (type == "diff") {
        diff_df <- df[, .(avg_methyl = avg_methyl[type == "Observed"] - avg_methyl[type == "Expected"]), by = x]
        flank <- max(abs(diff_df$x), na.rm = TRUE)
        diff_df[, avg_methyl := avg_methyl - mean(avg_methyl[abs(x) >= flank - flank.norm], na.rm = TRUE)]
        diff_df[, group := g]
        return(diff_df)
    } else {
        eo_df <- copy(df)
        eo_df[, type := paste(type, g, sep = "_")]
        return(eo_df)
    }
  })
  
  df <- rbindlist(Filter(Negate(is.null), per_group))
  if(nrow(df) == 0) return(NULL)
  
  mean_devs <- motif_mean_deviations(motif, dev_list, groups_here)
  labels <- vapply(groups_here, function(g) make_label(g, mean_devs[[g]]), character(1))
  
  if (type == "diff") {
      present <- intersect(group_order, unique(df$group))
      colours_here <- setNames(unname(cell_type_colors[present]), labels[present])
      df[, group := factor(labels[match(group, present)], levels = labels[present])]

      p <- ggplot(df, aes(x = x, y = avg_methyl, colour = group)) +
        geom_line(linewidth = 0.8) +
        scale_colour_manual(values = colours_here, name = NULL) +
        coord_cartesian(xlim = c(-200, 200)) +
        labs(title = motif, x = "Distance from motif center", y = "Methylation difference (Obs. - Exp.)") +
        theme_classic(base_size = base.size) +
        theme(
          legend.position = "bottom", 
          plot.title = element_text(face = "bold"),
          plot.background = element_rect(fill = "transparent", colour = NA),
          panel.background = element_rect(fill = "transparent", colour = NA)
        )
      return(p)
  } else {
      mapped_colors <- c()
      df$type_label <- ""
      
      for (g in groups_here) {
          exp_lbl <- paste("Expected", g)
          obs_lbl <- paste("Observed", labels[[g]])
          
          df[type == paste0("Expected_", g), type_label := exp_lbl]
          df[type == paste0("Observed_", g), type_label := obs_lbl]
          
          mapped_colors[exp_lbl] <- base_colors[[paste0("Expected_", g)]]
          mapped_colors[obs_lbl] <- base_colors[[paste0("Observed_", g)]]
      }
      
      ordered_labels <- c()
      for(g in group_order) {
          if(g %in% groups_here) {
              ordered_labels <- c(ordered_labels, paste("Expected", g), paste("Observed", labels[[g]]))
          }
      }
      
      present_labels <- intersect(ordered_labels, unique(df$type_label))
      df[, type_label := factor(type_label, levels = present_labels)]
      
      p <- ggplot(df, aes(x = x, y = avg_methyl, color = type_label)) +
        geom_line(linewidth = 0.8) +
        scale_color_manual(values = mapped_colors, name = NULL) +
        coord_cartesian(xlim = c(-200, 200)) +
        labs(title = motif, x = "Distance from motif center", y = "Average methylation") +
        theme_classic(base_size = base.size) +
        theme(
          legend.position = "bottom", 
          plot.title = element_text(face = "bold"),
          plot.background = element_rect(fill = "transparent", colour = NA),
          panel.background = element_rect(fill = "transparent", colour = NA)
        )
      return(p)
  }
}

p_jun <- footprint_panel("JUN", "eo")
p_fosl2 <- footprint_panel("FOSL2", "eo")
p_batf <- footprint_panel("BATF", "diff")
p_spi1 <- footprint_panel("SPI1", "diff")
p_foot <- Filter(Negate(is.null), list(p_jun, p_fosl2, p_batf, p_spi1))

#####################################################################
# Load Data for Scatters & Boxplots: methylTFR, RNA, VIPER
#####################################################################

log_info("Loading deviations for scatter plots and boxplots...")
dev_obj <- readRDS(dev.file)
cd <- as.data.frame(colData(dev_obj))
cell_types <- as.character(cd$cellType)
keep <- which(!is.na(cell_types) & !cell_types %in% drop.cell.types)
dev_obj <- dev_obj[, keep]
cell_types <- cell_types[keep]
sample_labels <- paste(sample_donors(colnames(dev_obj)), cell_types, sep = "_")
deviations_mat <- deviations(dev_obj)
colnames(deviations_mat) <- sample_labels

log_info("Loading matched RNA-seq data from ", rna.dir)
rna_files <- list.files(rna.dir, pattern = "IHECrsem.20190606.hs38.genes.results$", full.names = TRUE)
meth_files <- list.files(meth.dir, pattern = "bed$", full.names = TRUE)
meth_sample_ids <- gsub(".*?/(\\d+_Hf\\d+_Bl\\w+_Ct)_.*", "\\1", meth_files)
rna_sample_ids <- sub("\\..*", "", basename(rna_files))
matched_rna_files <- rna_files[rna_sample_ids %in% meth_sample_ids]

rna_list <- lapply(matched_rna_files, function(file) {
  df <- read.delim(file, stringsAsFactors = FALSE)[, c("gene_id", "expected_count")]
  colnames(df)[2] <- tools::file_path_sans_ext(basename(file))
  df
})
expr_matrix <- Reduce(function(x, y) full_join(x, y, by = "gene_id"), rna_list)
expr_matrix$gene_id_clean <- sub("\\..*$", "", expr_matrix$gene_id)

mapping <- AnnotationDbi::select(org.Hs.eg.db, keys = expr_matrix$gene_id_clean, keytype = "ENSEMBL", columns = "SYMBOL")
mapping <- mapping[!is.na(mapping$SYMBOL), ]
mapping <- mapping[!duplicated(mapping$ENSEMBL), ]
colnames(mapping) <- c("gene_id_clean", "gene_short_name")

expr <- merge(expr_matrix, mapping, by = "gene_id_clean")
rownames(expr) <- make.unique(expr$gene_short_name)
expr <- expr[, !colnames(expr) %in% c("gene_id_clean", "gene_id", "gene_short_name")]
colnames(expr) <- clean_gex_names(colnames(expr))
gex_cell_types <- sub("^.*_", "", colnames(expr))
expr <- expr[, !gex_cell_types %in% drop.cell.types, drop = FALSE]

rna_mat <- log2(t(t(as.matrix(expr)) / colSums(as.matrix(expr), na.rm = TRUE)) * 1e6 + 1)
viper_mat <- readRDS(viper.out)

#####################################################################
# Panel B & C: Scatter plots with top 10 labeling per category
#####################################################################

create_scatter <- function(nm, comp_info, x_mat, y_mat, x_label, y_label, x_name, y_name, label_motifs=TRUE) {
  test <- comp_info$test
  ref <- comp_info$ref
  
  log_info("Building scatter plot for ", test, " vs ", ref, " (", nm, ")")
  
  shared_motifs <- intersect(rownames(x_mat), rownames(y_mat))
  shared_samples <- intersect(colnames(x_mat), colnames(y_mat))
  
  x_sub <- x_mat[shared_motifs, shared_samples, drop=FALSE]
  y_sub <- y_mat[shared_motifs, shared_samples, drop=FALSE]
  ct_sub <- sub("^.*_", "", shared_samples)
  
  diff_x <- rowMeans(x_sub[, ct_sub == test, drop=FALSE]) - rowMeans(x_sub[, ct_sub == ref, drop=FALSE])
  diff_y <- rowMeans(y_sub[, ct_sub == test, drop=FALSE]) - rowMeans(y_sub[, ct_sub == ref, drop=FALSE])
  
  # Significance X
  sig_x <- c()
  if(x_name == "methylTFR") {
      test_motifs <- intersect(shared_motifs, rownames(deviations_mat))
      idx <- which(cell_types %in% c(test, ref))
      diff_test <- differential_deviation_test(
          deviations_mat[test_motifs, idx, drop=FALSE], 
          groups = cell_types[idx], 
          alternative = "two.sided", 
          parametric = TRUE
      )
      sig_x <- diff_test$motifs[diff_test$p_value_adjusted < 0.05]
  } else {
      pvals_x <- apply(x_sub, 1, function(row) {
         tryCatch(t.test(row[ct_sub == test], row[ct_sub == ref])$p.value, error=function(e) 1)
      })
      sig_x <- shared_motifs[p.adjust(pvals_x, method="BH") < 0.05]
  }
  
  # Significance Y (TF Expr / VIPER)
  pvals_y <- apply(y_sub, 1, function(row) {
     tryCatch(t.test(row[ct_sub == test], row[ct_sub == ref], var.equal = TRUE)$p.value, error=function(e) 1)
  })
  sig_y <- shared_motifs[pvals_y < 0.05]
  
  df <- data.frame(Motif = shared_motifs, DiffX = diff_x, DiffY = diff_y)
  
  # 4-color status assignment
  df$Status <- "Not differential"
  df$Status[df$Motif %in% sig_x] <- paste("Differential in", x_name)
  df$Status[df$Motif %in% sig_y] <- paste("Differential in", y_name)
  df$Status[df$Motif %in% sig_x & df$Motif %in% sig_y] <- "Differential in both"
  
  df$Status <- factor(df$Status, levels = c("Not differential", paste("Differential in", x_name), paste("Differential in", y_name), "Differential in both"))
  
  # Label top 10 motifs per category
  df$Label <- ""
  for (cat in levels(df$Status)) {
    if (cat == "Not differential") next
    rows <- which(df$Status == cat)
    if (length(rows) == 0) next
    
    if (cat == "Differential in both") {
      rank_by <- sqrt(df$DiffX[rows]^2 + df$DiffY[rows]^2)
    } else if (cat == paste("Differential in", x_name)) {
      rank_by <- abs(df$DiffX[rows])
    } else if (cat == paste("Differential in", y_name)) {
      rank_by <- abs(df$DiffY[rows])
    }
    
    keep <- rows[order(rank_by, decreasing = TRUE)[seq_len(min(10, length(rows)))]]
    df$Label[keep] <- df$Motif[keep]
  }
  
  df <- df[order(df$Status), ]
  
  scatter_colors <- c(
    "Not differential" = "#BDBDBD",
    "Differential in both" = "#E41A1C",            # Red
    "Differential in methylTFR" = "#4DAF4A",       # Green
    "Differential in TF Expression" = "#377EB8",   # Blue
    "Differential in VIPER" = "#377EB8"            # Blue
  )
  
  cor_val <- cor(df$DiffX, df$DiffY, use="complete.obs", method="pearson")
  
  max_x <- max(abs(df$DiffX), na.rm = TRUE) * 1.1
  max_y <- max(abs(df$DiffY), na.rm = TRUE) * 1.1
  
  p <- ggplot(df, aes(x = DiffX, y = DiffY)) +
    geom_hline(yintercept = 0, linetype="dashed", color="gray") +
    geom_vline(xintercept = 0, linetype="dashed", color="gray") +
    geom_point(aes(color = Status), size = 2, alpha=0.7) +
    scale_color_manual(values = scatter_colors) +
    annotate("text", x = -max_x * 0.95, y = max_y * 0.95, 
             label = sprintf("r = %.2f\nn = %d", cor_val, nrow(df)), hjust=0, vjust=1, fontface="bold") +
    labs(title = paste(test, "vs", ref), x = x_label, y = y_label) +
    coord_cartesian(xlim = c(-max_x, max_x), ylim = c(-max_y, max_y)) + 
    theme_classic(base_size = base.size) +
    theme(
      legend.title = element_blank(),
      plot.background = element_rect(fill = "transparent", colour = NA),
      panel.background = element_rect(fill = "transparent", colour = NA)
    )
    
  if(label_motifs) {
    p <- p + geom_text_repel(
      data = subset(df, Label != ""), 
      aes(label = Label, color = Status), 
      size = 5, 
      max.overlaps = Inf,
      min.segment.length = Inf, 
      show.legend = FALSE
    )
  }
  return(p)
}

# Generate scatter plots
log_info("Generating methylTFR vs TF Expression scatters...")
p_scat_mtfr_expr_cm <- create_scatter("CM", comparisons$CM, deviations_mat, rna_mat, "methylTFR Diff (Z-score)", "TF Expr Diff (Log2FC)", "methylTFR", "TF Expression")
p_scat_mtfr_expr_em <- create_scatter("EM", comparisons$EM, deviations_mat, rna_mat, "methylTFR Diff (Z-score)", "TF Expr Diff (Log2FC)", "methylTFR", "TF Expression")

log_info("Generating methylTFR vs VIPER scatters...")
p_scat_mtfr_viper_cm <- create_scatter("CM", comparisons$CM, deviations_mat, viper_mat, "methylTFR Diff (Z-score)", "VIPER Activity Diff", "methylTFR", "VIPER")
p_scat_mtfr_viper_em <- create_scatter("EM", comparisons$EM, deviations_mat, viper_mat, "methylTFR Diff (Z-score)", "VIPER Activity Diff", "methylTFR", "VIPER")

#####################################################################
# Panel D: Concordance Boxplots across cell types
#####################################################################

create_concordance_boxplot <- function(mat1, mat2, title_text, y_label) {
  shared_motifs <- intersect(rownames(mat1), rownames(mat2))
  shared_samples <- intersect(colnames(mat1), colnames(mat2))
  
  m1 <- mat1[shared_motifs, shared_samples, drop=FALSE]
  m2 <- mat2[shared_motifs, shared_samples, drop=FALSE]
  ct <- sub("^.*_", "", shared_samples)
  
  # Z-score standardization across all 6 samples
  m1_z <- t(apply(m1, 1, scale))
  m2_z <- t(apply(m2, 1, scale))
  colnames(m1_z) <- shared_samples
  colnames(m2_z) <- shared_samples
  
  df_list <- lapply(unique(ct), function(g) {
    samps <- which(ct == g)
    val1 <- rowMeans(m1_z[, samps, drop=FALSE], na.rm=TRUE)
    val2 <- rowMeans(m2_z[, samps, drop=FALSE], na.rm=TRUE)
    concordance <- val1 * val2
    data.frame(Motif = shared_motifs, CellType = g, Concordance = concordance)
  })
  
  df <- rbindlist(df_list)
  df$CellType <- factor(df$CellType, levels = group_order)
  n_motifs <- length(shared_motifs)
  
  # Calculate pairwise statistics
  pairs <- list(c(1, 2), c(2, 3), c(1, 3))
  pvals <- vapply(pairs, function(p) {
    val1 <- df$Concordance[as.numeric(df$CellType) == p[1]]
    val2 <- df$Concordance[as.numeric(df$CellType) == p[2]]
    if(length(val1) < 2 || length(val2) < 2) return(1)
    wilcox.test(val1, val2, exact = FALSE)$p.value
  }, numeric(1))
  padj <- p.adjust(pvals, method = "BH")
  
  # Clean padj formatting (4 decimals, or scientific for very small values)
  fmt_p <- function(p) {
    if (p < 0.0001) sprintf("p.adj = %.2e", p)
    else sprintf("p.adj = %.4f", p)
  }
  padj_labels <- vapply(padj, fmt_p, character(1))
  
  # Dynamically map out coordinates for the comparison brackets
  y_min <- quantile(df$Concordance, 0.01, na.rm = TRUE)
  y_max <- quantile(df$Concordance, 0.99, na.rm = TRUE)
  step <- (y_max - y_min) * 0.1
  
  y_bracket1 <- y_max + step * 0.6  # TN vs TCM
  y_bracket2 <- y_max + step * 1.4  # TCM vs TEM
  y_bracket3 <- y_max + step * 2.2  # TN vs TEM
  y_top_text <- y_max + step * 3.0  # Top motif label
  y_lim_top  <- y_top_text + step
  
  tick <- step * 0.15 # little vertical downtick for the edges of the bracket
  
  # Data frame guiding the drawing of brackets
  brackets <- data.frame(
    x = c(1, 2, 1),
    xend = c(2, 3, 3),
    y = c(y_bracket1, y_bracket2, y_bracket3),
    label = padj_labels
  )
  
  ggplot(df, aes(x = CellType, y = Concordance)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_jitter(aes(color = CellType), width = 0.2, height = 0, size = 0.8, alpha = 0.4) +
    geom_boxplot(aes(fill = CellType), alpha = 0.8, outlier.shape = NA, width = 0.4, color = "black", linewidth = 0.6) +
    
    # 1. Main horizontal lines of the brackets
    geom_segment(data = brackets, aes(x = x, xend = xend, y = y, yend = y), color = "black", linewidth = 0.5) +
    # 2. Left side vertical ticks
    geom_segment(data = brackets, aes(x = x, xend = x, y = y, yend = y - tick), color = "black", linewidth = 0.5) +
    # 3. Right side vertical ticks
    geom_segment(data = brackets, aes(x = xend, xend = xend, y = y, yend = y - tick), color = "black", linewidth = 0.5) +
    # 4. Actual padj Text on top of brackets
    geom_text(data = brackets, aes(x = (x + xend) / 2, y = y + tick * 1.5, label = label), size = 3, vjust = 0) +
    
    scale_fill_manual(values = cell_type_colors) +
    scale_color_manual(values = cell_type_colors) +
    coord_cartesian(ylim = c(y_min, y_lim_top)) + 
    labs(title = title_text, x = "Cell Type", y = y_label) +
    
    # Number of motifs on top center, standard font (not bold)
    annotate("text", x = 2, y = y_top_text, label = sprintf("n = %s motifs", format(n_motifs, big.mark=",")), size = 3.5, fontface = "plain") +
    
    theme_classic(base_size = base.size) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      # Strict transparency enforcement to prevent clipping mask artifacts in Adobe Illustrator
      plot.background = element_rect(fill = "transparent", colour = NA),
      panel.background = element_rect(fill = "transparent", colour = NA)
    )
}

log_info("Generating concordance boxplot panels...")
p_box_mtfr_expr <- create_concordance_boxplot(deviations_mat, rna_mat, "methylTFR vs. TF Expression", "Concordance (Z-score Product)")
p_box_mtfr_viper <- create_concordance_boxplot(deviations_mat, viper_mat, "methylTFR vs. VIPER Activity", "Concordance (Z-score Product)")
p_box_expr_viper <- create_concordance_boxplot(rna_mat, viper_mat, "TF Expression vs. VIPER Activity", "Concordance (Z-score Product)")

#####################################################################
# Final Assembly
#####################################################################
log_info("Assembling final supplementary figure...")

p_foot[[1]] <- p_foot[[1]] + labs(tag = "A")
p_scat_mtfr_expr_cm <- p_scat_mtfr_expr_cm + labs(tag = "B")
p_scat_mtfr_viper_cm <- p_scat_mtfr_viper_cm + labs(tag = "C")
p_box_mtfr_expr <- p_box_mtfr_expr + labs(tag = "D")

scatter_grid <- wrap_plots(
  p_scat_mtfr_expr_cm, p_scat_mtfr_expr_em,
  p_scat_mtfr_viper_cm, p_scat_mtfr_viper_em,
  ncol = 2
) + plot_layout(guides = "collect") & theme(legend.position = "bottom")

box_grid <- wrap_plots(p_box_mtfr_expr, p_box_mtfr_viper, p_box_expr_viper, ncol = 3)

supplementary <- wrap_plots(
  wrap_plots(p_foot, nrow = 1),
  scatter_grid,
  box_grid,
  ncol = 1, heights = c(1, 3, 1.5)
) & tag_theme

file <- file.path(plot.dir, "memT_supplementary.pdf")
ggsave(file, supplementary & no_bg, width = fig.width, height = fig.height, bg = "transparent", limitsize = FALSE)
log_success("Wrote supplementary figure to ", file)
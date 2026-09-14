#!/usr/bin/env Rscript

#####################################################################
# test_altius.R
# created on 08-09-2026 by Irem B Gunduz
# Footprints of the AP-1 archetypes of the Altius set. Top row is the
# observed profile against the expected one, bottom row is observed
# minus expected with the flanks set to zero, the way the main figure
# draws it. Nothing is read from the deviations.
#####################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(logger)
  library(GenomicRanges)
  library(methylTFR)
  library(methylTFRAnnotationHg38)
})

motifSet <- "altius"
motifs <- c("ap1_1", "ap1_2","ccaat_cebp")
plot.window <- 200L

# The flanks are the baseline of the difference curve, so it starts from
# zero away from the motif, as in the main figure
flank.norm <- 30L
# The number printed next to each cell type is the mean of the difference
# curve over this centre window
centre.window <- 25L

group_colors <- c(
  "Bcell_naive" = "#C46B8C",
  "Bcell_mem" = "#7B1E3D",
  "Tcell" = "#4FC3D9",
  "Monocytes" = "#C2703D",
  "Granulocytes" = "#E8A33D"
)

analysis.dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/blueprint/"
rnb.tag <- "RnBeads_230826"
cache.file5 <- file.path(
  analysis.dir, "debug", paste0("msites_merged_celltype5G_", rnb.tag, ".rds")
)
fig.dir <- file.path(
  "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR_manuscript",
  "figures", "blueprint"
)
if (!dir.exists(fig.dir)) dir.create(fig.dir, recursive = TRUE)

#####################################################################
# Annotation and merged methylation
#####################################################################

tf_bindsites <- getTFbindsites(motifSet = motifSet)
gcfreqs <- getGCfreq(motifSet = motifSet)
gc_dist <- getGenomeGC()

absent <- setdiff(motifs, names(tf_bindsites))
if (length(absent) > 0) {
  stop(
    paste(absent, collapse = ", "), " not in ", motifSet, ". Names look like: ",
    paste(head(grep("ap", names(tf_bindsites), value = TRUE), 10), collapse = ", ")
  )
}

if (!file.exists(cache.file5)) stop("Merged methylation not found: ", cache.file5)
msites <- as.list(readRDS(cache.file5))
msites <- msites[intersect(names(group_colors), names(msites))]
if (length(msites) == 0) stop("No known cell type in the cache")
log_info("Cell types: ", paste(names(msites), collapse = ", "))

#####################################################################
# One footprint per motif, drawn two ways
#####################################################################

footprint <- function(motif) {
  df <- rbindlist(Filter(Negate(is.null), lapply(names(msites), function(g) {
    d <- tryCatch(
      plotExpectedFootprint(
        motif = motif, tf_bindsites = tf_bindsites, msites = msites[[g]],
        sample_name = g, gc_dist = gc_dist, gcfreqs = gcfreqs,
        returnPlotData = TRUE
      )$plotDF,
      error = function(e) {
        log_warn(motif, " / ", g, ": ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(d)) {
      return(NULL)
    }
    d <- as.data.table(d)
    d[, group := g]
    d
  })))
  if (nrow(df) == 0) {
    return(NULL)
  }
  df[, group := factor(group, levels = intersect(names(group_colors), unique(group)))]
  df
}

base_theme <- theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(hjust = 0, face = "plain", size = 10),
    legend.position = "bottom",
    legend.key.size = unit(3.5, "mm")
  )

raw_panel <- function(df, motif) {
  ggplot(df, aes(x = x, y = avg_methyl, colour = group, linetype = type)) +
    geom_line(linewidth = 0.4) +
    scale_colour_manual(values = group_colors[levels(df$group)], name = NULL) +
    scale_linetype_manual(
      values = c("Observed" = "solid", "Expected" = "dashed"), name = NULL
    ) +
    coord_cartesian(xlim = c(-plot.window, plot.window)) +
    labs(
      title = paste0(motif, " (", length(tf_bindsites[[motif]]), " sites)"),
      x = "Distance from motif centre", y = "Average methylation"
    ) +
    base_theme
}

diff_panel <- function(df, motif) {
  d <- df[, .(
    avg_methyl = avg_methyl[type == "Observed"] - avg_methyl[type == "Expected"]
  ), by = .(x, group)]

  flank <- max(abs(d$x), na.rm = TRUE)
  d[, avg_methyl := avg_methyl -
    mean(avg_methyl[abs(x) >= flank - flank.norm], na.rm = TRUE), by = group]

  centre <- d[abs(x) <= centre.window, .(dip = mean(avg_methyl, na.rm = TRUE)), by = group]
  log_info(
    motif, " centre of the difference curve: ",
    paste(centre$group, round(centre$dip, 3), sep = " = ", collapse = ", ")
  )

  present <- levels(droplevels(d$group))
  labels <- paste0(present, " (", format(round(centre$dip[match(present, centre$group)], 2),
    nsmall = 2
  ), ")")
  d[, group := factor(labels[match(as.character(group), present)], levels = labels)]
  colours <- setNames(unname(group_colors[present]), labels)

  ggplot(d, aes(x = x, y = avg_methyl, colour = group)) +
    geom_hline(yintercept = 0, linetype = "dotted", colour = "grey55") +
    geom_line(linewidth = 0.4) +
    scale_colour_manual(values = colours, name = NULL) +
    coord_cartesian(xlim = c(-plot.window, plot.window)) +
    labs(
      title = paste0(motif, ", bias corrected"),
      x = "Distance from motif centre",
      y = "Methylation difference (Observed - Expected)"
    ) +
    base_theme
}

profiles <- lapply(motifs, footprint)
names(profiles) <- motifs
profiles <- Filter(Negate(is.null), profiles)
if (length(profiles) == 0) stop("No footprint could be drawn")

panels <- c(
  lapply(names(profiles), function(m) raw_panel(profiles[[m]], m)),
  lapply(names(profiles), function(m) diff_panel(profiles[[m]], m))
)

file <- file.path(fig.dir, paste0("test_", motifSet, "_ap1_footprints.pdf"))
ggsave(file, wrap_plots(panels, ncol = length(profiles)),
  width = 6 * length(profiles), height = 9, bg = "transparent"
)
log_success("Wrote ", file)
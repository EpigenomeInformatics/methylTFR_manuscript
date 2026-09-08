#!/usr/bin/env Rscript

#####################################################################
# test_altius.R
# created on 08-09-2026 by Irem B Gunduz
# Footprints of the AP-1 archetypes of the Altius set, one panel per
# motif, one line per cell type. Nothing is read from the deviations.
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
motifs <- c("ap1_1", "ap1_2")
plot.window <- 200L

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
# One panel per motif
#####################################################################

panel <- function(motif) {
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
    log_warn(motif, ": no footprint")
    return(NULL)
  }
  present <- intersect(names(group_colors), unique(df$group))
  df[, group := factor(group, levels = present)]

  ggplot(df, aes(x = x, y = avg_methyl, colour = group, linetype = type)) +
    geom_line(linewidth = 0.4) +
    scale_colour_manual(values = group_colors[present], name = NULL) +
    scale_linetype_manual(
      values = c("Observed" = "solid", "Expected" = "dashed"), name = NULL
    ) +
    coord_cartesian(xlim = c(-plot.window, plot.window)) +
    labs(
      title = paste0(motif, " (", length(tf_bindsites[[motif]]), " sites)"),
      x = "Distance from motif centre", y = "Average methylation"
    ) +
    theme_classic(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0, face = "plain", size = 10),
      legend.position = "bottom",
      legend.key.size = unit(3.5, "mm")
    )
}

panels <- Filter(Negate(is.null), lapply(motifs, panel))
if (length(panels) == 0) stop("No footprint could be drawn")

file <- file.path(fig.dir, paste0("test_", motifSet, "_ap1_footprints.pdf"))
ggsave(file,
  wrap_plots(panels, ncol = length(panels), guides = "collect") &
    theme(legend.position = "bottom"),
  width = 6 * length(panels), height = 4.5, bg = "transparent"
)
log_success("Wrote ", file)
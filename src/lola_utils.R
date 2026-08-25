#####################################################################
# lola_utils.R
# LOLA enrichment plotting helpers for the methylTFR manuscript
#####################################################################

# Motifs enriched on the naive side are the lighter green, the memory side
# the darker one, matching CELL_TYPE_COLORS in utils.R
LOLA_COLORS <- c("loss" = "#A1D99B", "gain" = "#2E8B57")

#' Direction bar drawn above a volcano panel.
#'
#' A standalone strip rather than an annotation, because the volcano is
#' facetted into a loss and a gain half and an annotation would be repeated
#' in each facet instead of spanning both.
#'
#' @param grp1 group on the right, the one "gain" refers to
#' @param grp2 group on the left
#' @param cell_colors named colour vector covering both groups
directionHeader <- function(grp1, grp2, cell_colors = CELL_TYPE_COLORS) {
  missing <- setdiff(c(grp1, grp2), names(cell_colors))
  if (length(missing) > 0) {
    stop("No colour for ", paste(missing, collapse = ", "))
  }
  ggplot(data.frame(x = 0, y = 0), aes(x, y)) +
    annotate("point", x = 0.03, y = 0, colour = cell_colors[[grp2]], size = 3) +
    annotate("text", x = 0.07, y = 0, label = grp2, hjust = 0, size = 4.2) +
    annotate("segment",
      x = 0.18, xend = 0.82, y = 0, yend = 0,
      arrow = arrow(ends = "both", length = unit(1.8, "mm")), linewidth = 0.4
    ) +
    annotate("text", x = 0.93, y = 0, label = grp1, hjust = 1, size = 4.2) +
    annotate("point", x = 0.97, y = 0, colour = cell_colors[[grp1]], size = 3) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-0.5, 0.5), expand = c(0, 0)) +
    theme_void()
}

#' Volcano plot of a LOLA enrichment result.
#'
#' The previous version read the table from a global object called `res`
#' rather than from its own `lolaRes` argument, so it silently plotted
#' whatever happened to be in the calling environment, or failed with
#' "object 'res' not found". It now uses the argument.
#'
#' @param lolaRes the LOLA result object, res$region is indexed below
#' @param outputDir directory the pdf is written to
#' @param comparison name or index of the element of lolaRes$region
#' @param region region type inside that element, e.g. "tiling1kb"
#' @param label a label used in the file name, defaults to the comparison
#' @param grp1 group that "hyper" refers to, drawn on the right of the header
#' @param grp2 the reference group, drawn on the left
#' @param signifCol significance column, "qValue" is converted to -log10
#' @param n.label number of motifs labelled per direction
#' @param userSets the two userSet names holding the hyper and hypo sets
#' @param cell_colors colours of the two groups in the header
#' @return the ggplot object, invisibly
lolaVolcanoPlot <- function(lolaRes, outputDir, comparison, region,
                            label = NULL, grp1 = NULL, grp2 = NULL,
                            signifCol = "qValue", n.label = 10,
                            userSets = c("rankCut_1000_hyper", "rankCut_1000_hypo"),
                            cell_colors = CELL_TYPE_COLORS) {
  if (is.null(label)) label <- as.character(comparison)

  if (is.null(lolaRes$region[[comparison]])) {
    stop(
      "lolaRes$region[[", comparison, "]] is missing. Comparisons available: ",
      paste(names(lolaRes$region), collapse = " | ")
    )
  }
  if (is.null(lolaRes$region[[comparison]][[region]])) {
    stop(
      "Region '", region, "' is missing from '", comparison,
      "'. Regions available: ",
      paste(names(lolaRes$region[[comparison]]), collapse = ", ")
    )
  }
  # LOLA results are data.tables. Their [[ ]] errors with "subscript out of
  # bounds" where a data.frame returns NULL, and $<- on them copies with a
  # warning, so the table is converted once here.
  df <- as.data.frame(lolaRes$region[[comparison]][[region]])
  needed <- c("userSet", "description", "oddsRatio", "qValue")
  missing_cols <- setdiff(needed, colnames(df))
  if (length(missing_cols) > 0) {
    stop(
      "Column(s) missing from the LOLA table: ",
      paste(missing_cols, collapse = ", "), ". Present: ",
      paste(colnames(df), collapse = ", ")
    )
  }
  df <- df[df$userSet %in% userSets, ]
  if (nrow(df) == 0) {
    stop("No rows for userSets ", paste(userSets, collapse = ", "))
  }

  # Rows with no description cannot be named
  n_before <- nrow(df)
  df <- df[!is.na(df$description), ]
  if (nrow(df) < n_before) {
    message(n_before - nrow(df), " LOLA rows without a description dropped")
  }

  df$condition <- ifelse(grepl("hyper$", df$userSet), "gain", "loss")

  if (signifCol == "qValue") {
    signifCol <- "qValueLog"
    df$qValueLog <- -log10(df$qValue)
    df$differential <- ifelse(df$qValueLog > 1.3, "Differential", "Non-differential")
  }
  oddsRatioCol <- "log2OR"
  df$log2OR <- log2(df$oddsRatio)

  pick_top <- function(cond) {
    d <- df[df$condition == cond & df$differential == "Differential", ]
    d <- d[order(-d[[signifCol]]), ]
    d <- head(d, n = n.label)
    dplyr::select(d, description, condition, log2OR, qValueLog)
  }
  top_motifs_gain <- pick_top("gain")
  top_motifs_loss <- pick_top("loss")
  top_motifs_loss$log2OR <- -top_motifs_loss$log2OR
  top_motifs <- rbind(top_motifs_gain, top_motifs_loss)

  pp <- plotVolcano(df, top_motifs, oddsRatioCol, signifCol)

  # Direction bar on top, when the two groups are known
  if (!is.null(grp1) && !is.null(grp2)) {
    if (requireNamespace("patchwork", quietly = TRUE)) {
      pp <- patchwork::wrap_plots(
        directionHeader(grp1, grp2, cell_colors), pp,
        ncol = 1, heights = c(1, 14)
      )
    } else {
      message("patchwork is not installed, the direction header is skipped")
    }
  }

  if (!dir.exists(outputDir)) dir.create(outputDir, recursive = TRUE)
  file <- file.path(outputDir, paste0("lolaVolcanoPlot_", region, "_", label, ".pdf"))
  pdf(file, width = 7, height = 6)
  print(pp)
  dev.off()

  invisible(pp)
}

#' Volcano panel used by lolaVolcanoPlot.
#'
#' Points are black and only the labelled motifs carry colour, which is what
#' keeps a few hundred region sets readable. Loss sits on the negative side
#' of the axis, gain on the positive one, each in its own facet.
plotVolcano <- function(df, top_motifs, oddsRatioCol, signifCol,
                        colors = LOLA_COLORS) {
  df$condition <- factor(df$condition, levels = c("loss", "gain"))
  top_motifs$condition <- factor(top_motifs$condition, levels = c("loss", "gain"))

  df_combined <- dplyr::bind_rows(
    df %>% dplyr::filter(condition == "loss") %>% dplyr::mutate(log2OR = pmin(log2OR, 0)),
    df %>% dplyr::filter(condition == "gain") %>% dplyr::mutate(log2OR = pmax(log2OR, 0))
  )
  top_motifs_combined <- dplyr::bind_rows(
    top_motifs %>% dplyr::filter(condition == "loss") %>% dplyr::mutate(log2OR = pmin(log2OR, 0)),
    top_motifs %>% dplyr::filter(condition == "gain") %>% dplyr::mutate(log2OR = pmax(log2OR, 0))
  )

  ggplot(df_combined) +
    aes_string(oddsRatioCol, signifCol) +
    geom_point(colour = "black", size = 1.1) +
    ggrepel::geom_text_repel(aes(label = description),
      data = top_motifs_combined,
      size = 3.2,
      colour = colors[as.character(top_motifs_combined$condition)],
      min.segment.length = 0,
      segment.size = 0.25,
      segment.alpha = 0.6,
      seed = 42,
      box.padding = 0.5,
      force = 4,
      max.overlaps = Inf
    ) +
    ylab(expression(-log[10] * "(Q-Value)")) +
    xlab(expression(log[2] * "(Odds-Ratio)")) +
    theme_classic(base_size = 13) +
    theme(
      axis.text = element_text(size = 11),
      legend.position = "none",
      strip.text = element_blank(),
      strip.background = element_blank(),
      panel.spacing = unit(0, "lines")
    ) +
    facet_wrap(~condition, scales = "free_x")
}

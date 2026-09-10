#####################################################################
# lola_utils.R
# LOLA enrichment plotting helpers for the methylTFR manuscript
#####################################################################

# Motifs enriched on the naive side are the lighter green, the memory side
# the darker one, matching CELL_TYPE_COLORS in utils.R
LOLA_COLORS <- c("loss" = "#A1D99B", "gain" = "#2E8B57")

#' Direction bar drawn above a volcano panel.
#'
#' A standalone strip rather than an annotation: the volcano is facetted into
#' a loss and a gain half, and an annotation is repeated in each facet instead
#' of spanning both.
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
#' @param lolaRes the LOLA result object, res$region is indexed below
#' @param outputDir directory the pdf is written to
#' @param comparison name or index of the element of lolaRes$region
#' @param region region type inside that element, e.g. "tiling1kb"
#' @param label a label used in the file name, defaults to the comparison
#' @param grp1 group that "hyper" refers to, drawn on the right of the header
#' @param grp2 the reference group, drawn on the left
#' @param signifCol significance column, "qValue" is converted to -log10
#' @param n.label number of motifs labelled per direction, used only when
#'   motifs is NULL
#' @param motifs character vector of motif names to label. When given, these
#'   are labelled on whichever side they fall and the ranking is ignored, so
#'   the same motifs appear in every contrast. Names are matched after the
#'   JASPAR accession is stripped, so "FOSL1::JUN" rather than
#'   "MA1128.1_FOSL1::JUN".
#' @param userSets the two userSet names holding the hyper and hypo sets
#' @param cell_colors colours of the two groups in the header
#' @return the ggplot object, invisibly
lolaVolcanoPlot <- function(lolaRes, outputDir, comparison, region,
                            label = NULL, grp1 = NULL, grp2 = NULL,
                            signifCol = "qValue", n.label = 10, motifs = NULL,
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
  # LOLA results are data.tables, whose [[ ]] and $<- semantics differ from
  # data.frames, so the table is converted once here.
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

  df$motif <- .strip_motif_id(df$description)

  if (is.null(motifs)) {
    # Rank by significance within each direction
    pick_top <- function(cond) {
      d <- df[df$condition == cond & df$differential == "Differential", ]
      d <- d[order(-d[[signifCol]]), ]
      head(d, n = n.label)
    }
    top_gain <- pick_top("gain")
    top_loss <- pick_top("loss")
  } else {
    # A fixed set, so the same motifs are named in every contrast. A motif
    # can appear in both userSets, in which case it is labelled on the side
    # where it is enriched.
    absent <- setdiff(motifs, df$motif)
    if (length(absent) > 0) {
      message(
        length(absent), " requested motif(s) are not in this LOLA table: ",
        paste(utils::head(absent, 15), collapse = ", ")
      )
    }
    wanted <- df[df$motif %in% motifs, ]
    if (nrow(wanted) == 0) {
      stop(
        "None of the requested motifs are in the LOLA table. Descriptions ",
        "look like '", df$description[1], "'."
      )
    }
    # One row per motif and direction, the most significant if repeated
    wanted <- wanted[order(-wanted[[signifCol]]), ]
    wanted <- wanted[!duplicated(paste(wanted$motif, wanted$condition)), ]
    top_gain <- wanted[wanted$condition == "gain", ]
    top_loss <- wanted[wanted$condition == "loss", ]
  }

  keep_cols <- c("description", "condition", "log2OR", "qValueLog")
  top_gain <- top_gain[, keep_cols, drop = FALSE]
  top_loss <- top_loss[, keep_cols, drop = FALSE]
  top_loss$log2OR <- -top_loss$log2OR
  top_motifs <- rbind(top_gain, top_loss)
  if (nrow(top_motifs) == 0) {
    message("No motif could be labelled on this volcano")
  }

  pp <- plotVolcano(df, top_motifs, oddsRatioCol, signifCol)

  # Direction bar on top, when the two groups are known
  if (!is.null(grp1) && !is.null(grp2)) {
    if (requireNamespace("patchwork", quietly = TRUE)) {
      pp <- patchwork::wrap_plots(
        directionHeader(grp1, grp2, cell_colors),
        patchwork::wrap_elements(full = pp),
        ncol = 1, heights = c(1, 14)
      )
    } else {
      message("patchwork is not installed, the direction header is skipped")
    }
  }

  if (!dir.exists(outputDir)) dir.create(outputDir, recursive = TRUE)
  file <- file.path(outputDir, paste0("lolaVolcanoPlot_", region, "_", label, ".pdf"))
  pdf(file, width = 9, height = 5.5)
  print(pp)
  dev.off()

  invisible(pp)
}

#' Strip the JASPAR accession from a LOLA description.
#' "MA0637.1_CENPB" becomes "CENPB", "MA1128.1_FOSL1::JUN" becomes
#' "FOSL1::JUN". Same rule the LOLA against methylTFR panel uses.
.strip_motif_id <- function(x) gsub(".*_", "", x)

#' Volcano panel used by lolaVolcanoPlot.
#'
#' Two panels that mirror each other around zero: the loss side runs
#' outwards to the left on a reversed axis, the gain side outwards to the
#' right, so the distance from the centre is the enrichment in both. Built
#' with patchwork rather than facet_wrap because ggplot cannot reverse the
#' scale of one facet only.
#'
#' Points are black and only the labelled motifs carry colour, which is what
#' keeps a few hundred region sets readable.
plotVolcano <- function(df, top_motifs, oddsRatioCol, signifCol,
                        colors = LOLA_COLORS) {
  df$condition <- as.character(df$condition)
  top_motifs$condition <- as.character(top_motifs$condition)

  # Distance from the centre, so both sides run outwards from zero. The
  # clamp keeps a loss set with a positive odds ratio (or the reverse) from
  # crossing into the other panel.
  df$xplot <- ifelse(df$condition == "loss",
    abs(pmin(df[[oddsRatioCol]], 0)),
    pmax(df[[oddsRatioCol]], 0)
  )
  top_motifs$xplot <- abs(top_motifs[[oddsRatioCol]])
  df$yplot <- df[[signifCol]]
  top_motifs$yplot <- top_motifs[[signifCol]]
  top_motifs$lab <- .strip_motif_id(top_motifs$description)

  y_max <- max(df$yplot, na.rm = TRUE) * 1.05
  x_max <- max(df$xplot, na.rm = TRUE) * 1.05

  panel <- function(cond, reverse) {
    d <- df[df$condition == cond, ]
    tm <- top_motifs[top_motifs$condition == cond, ]
    p <- ggplot(d, aes(x = xplot, y = yplot)) +
      geom_point(colour = "black", size = 1.1) +
      ggrepel::geom_text_repel(
        data = tm, aes(label = lab),
        size = 3.2, colour = colors[[cond]],
        min.segment.length = 0, segment.size = 0.25, segment.alpha = 0.6,
        box.padding = 0.5, force = 4, max.overlaps = Inf, seed = 42
      ) +
      ylab(expression(-log[10] * "(Q-Value)")) +
      xlab(NULL) +
      theme_classic(base_size = 13) +
      theme(axis.text = element_text(size = 11), legend.position = "none")
    p <- p + if (reverse) {
      scale_x_reverse(limits = c(x_max, 0))
    } else {
      scale_x_continuous(limits = c(0, x_max))
    }
    p + scale_y_continuous(limits = c(0, y_max))
  }

  p_loss <- panel("loss", reverse = TRUE)
  p_gain <- panel("gain", reverse = FALSE) +
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank()
    )

  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("patchwork is required for the mirrored volcano layout")
  }
  patchwork::wrap_plots(p_loss, p_gain, nrow = 1) +
    patchwork::plot_annotation(
      caption = expression(log[2] * "(Odds-Ratio)"),
      theme = ggplot2::theme(
        plot.caption = element_text(hjust = 0.5, size = 12, vjust = 1)
      )
    )
}

#####################################################################
# Resolving a contrast to an index
#
# RnBeads labels comparisons "<group1> vs. <group2> (based on <column>)" but
# returns them in a vector whose names are cmp1, cmp2 and so on, so the
# description is the value and not the name. Matching the wrong one of the
# two silently resolves every contrast to the same element, which is how two
# contrasts end up sharing a figure.
#####################################################################

# Index of the single label carrying both group names as whole words.
# Errors rather than guessing, and never returns more than one.
resolve_comparison <- function(labels, test, ref, override = NA, what = "comparison") {
  if (length(labels) == 0) stop("No ", what, " labels to match against")
  if (!is.na(override)) {
    idx <- if (is.character(override)) match(override, labels) else as.integer(override)
    if (is.na(idx) || idx < 1 || idx > length(labels)) {
      stop("Override '", override, "' does not match any ", what)
    }
  } else {
    hit <- grepl(paste0("\\b", test, "\\b"), labels) &
      grepl(paste0("\\b", ref, "\\b"), labels)
    if (sum(hit) == 0) {
      stop(
        "No ", what, " matches ", test, " and ", ref, ". Available: ",
        paste(labels, collapse = " | ")
      )
    }
    if (sum(hit) > 1) {
      stop(
        "Several ", what, "s match ", test, " and ", ref, ": ",
        paste(labels[hit], collapse = " | "), ". Set an override."
      )
    }
    idx <- which(hit)
  }
  parts <- strsplit(sub(" \\(based on.*", "", labels[idx]), "\\s*vs\\.\\s*")[[1]]
  list(
    index = idx, name = labels[idx],
    grp1 = trimws(parts[1]), grp2 = trimws(parts[2])
  )
}

# The descriptive labels of an RnBeads comparison vector, which live in the
# values when the names are the generic cmp identifiers
comparison_labels <- function(cmps) {
  nms <- names(cmps)
  if (is.null(nms) || all(grepl("^cmp[0-9]+$", nms))) as.character(cmps) else nms
}

# The same for a LOLA result, whose region list is named by comparison
resolve_lola_comparison <- function(res_lola, test, ref, override = NA) {
  labels <- names(res_lola$region)
  if (is.null(labels)) stop("res_lola$region has no names")
  resolve_comparison(labels, test, ref, override, what = "LOLA comparison")
}

# Two contrasts resolving to one index means the figures would be identical.
# Called once per script, after every contrast has been resolved.
assert_distinct_comparisons <- function(hits) {
  idx <- vapply(hits, function(h) h$index, numeric(1))
  if (anyDuplicated(idx)) {
    stop(
      "Contrasts resolved to the same comparison: ",
      paste(names(hits), idx, sep = " -> ", collapse = ", "),
      ". Their figures would be identical."
    )
  }
  invisible(TRUE)
}

# First of the preferred region names that this comparison actually carries
resolve_lola_region <- function(res_lola, index, preferred) {
  available <- names(res_lola$region[[index]])
  hit <- intersect(preferred, available)
  if (length(hit) == 0) {
    stop(
      "None of ", paste(preferred, collapse = ", "), " found. Available: ",
      paste(available, collapse = ", ")
    )
  }
  hit[1]
}

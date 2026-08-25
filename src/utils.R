#####################################################################
# utils.R
# Shared plotting helpers for the methylTFR manuscript
#####################################################################

# Green scheme shared by every figure that splits by T cell subtype
CELL_TYPE_COLORS <- c(
  "TN" = "#C7E9C0",
  "TCM" = "#41AB5D",
  "TEM" = "#00441B",
  "TEMRA" = "#7FCDBB"
)

#' Density scatter of the group mean methylation for one comparison.
#'
#' The previous version looped over every comparison and kept only the one
#' the loop happened to end on, so the plot never matched the comparison the
#' caller wanted. The comparison is now selected explicitly, either by name
#' or by index, and the loop is gone.
#'
#' @param diffMeth an RnBeads differential methylation object
#' @param region region type, e.g. "tiling" or "sites"
#' @param comparison name or index of the comparison to plot
#' @param p.cut adjusted p-value cut used to highlight DMPs
rnbeadsDensityScatter <- function(diffMeth, region, comparison = 1, p.cut = 0.05) {
  cmps <- get.comparisons(diffMeth)
  if (is.character(comparison)) {
    if (!comparison %in% names(cmps)) {
      stop(
        "Unknown comparison '", comparison, "'. Available: ",
        paste(names(cmps), collapse = ", ")
      )
    }
    comparison <- match(comparison, names(cmps))
  }
  if (comparison < 1 || comparison > length(cmps)) {
    stop("comparison index out of range, ", length(cmps), " comparisons available")
  }

  cc <- names(cmps)[comparison]
  ccc <- cmps[cc]
  df2p <- get.table(diffMeth, ccc, region, return.data.frame = TRUE)
  grp.names <- get.comparison.grouplabels(diffMeth)[ccc, ]

  if (!is.element("comb.p.adj.fdr", colnames(df2p))) {
    stop("No comb.p.adj.fdr column in the ", region, " table of ", cc)
  }
  df2p$isDMP <- df2p[, "comb.p.adj.fdr"] < p.cut

  pp <- create.densityScatter(df2p[, c("mean.mean.g2", "mean.mean.g1")],
    is.special = df2p$isDMP,
    dens.subsample = TRUE, sparse.points = 0.001, add.text.cor = TRUE
  ) +
    labs(
      x = paste0("Mean-Methylation ( ", grp.names[2], " )"),
      y = paste0("Mean-Methylation ( ", grp.names[1], " )")
    ) +
    coord_fixed() +
    theme_classic() +
    theme(legend.position = "none")

  return(pp)
}

#' Resolve the group mean columns of an RnBeads differential table.
#'
#' Region tables carry mean.mean.g1 / mean.mean.g2, site tables carry
#' mean.g1 / mean.g2. Hardcoding the region names made the site tables fail.
#' @return a character vector, group 1 first then group 2
.diffmeth_mean_cols <- function(tbl) {
  if (all(c("mean.mean.g1", "mean.mean.g2") %in% colnames(tbl))) {
    return(c("mean.mean.g1", "mean.mean.g2"))
  }
  if (all(c("mean.g1", "mean.g2") %in% colnames(tbl))) {
    return(c("mean.g1", "mean.g2"))
  }
  stop(
    "No group mean columns found. Present: ",
    paste(colnames(tbl), collapse = ", ")
  )
}

#' Resolve a comparison by name or index.
.diffmeth_comparison <- function(diffMeth, comparison) {
  cmps <- get.comparisons(diffMeth)
  if (is.character(comparison)) {
    if (!comparison %in% names(cmps)) {
      stop(
        "Unknown comparison '", comparison, "'. Available: ",
        paste(names(cmps), collapse = " | ")
      )
    }
    comparison <- match(comparison, names(cmps))
  }
  if (comparison < 1 || comparison > length(cmps)) {
    stop("comparison index out of range, ", length(cmps), " available")
  }
  list(index = comparison, name = names(cmps)[comparison], cmp = cmps[names(cmps)[comparison]])
}

#' Density scatter with the differential regions taken from the RANK CUT.
#'
#' The companion of rnbeadsDensityScatter, which highlights regions by
#' adjusted p-value. Here a region counts as differential when its RnBeads
#' combinedRank falls within the top `rank.cut`, which is the same criterion
#' the LOLA enrichment uses through its rankCut_<n>_hyper / _hypo user sets.
#' Using it here keeps the scatter and the enrichment panels describing the
#' same set of regions.
#'
#' @param diffMeth an RnBeads differential methylation object
#' @param region region type, e.g. "tiling1kb" or "sites"
#' @param comparison name or index of the comparison to plot
#' @param rank.cut how many top ranking regions count as differential
#' @param auto.rank.cut if TRUE, RnBeads picks the cut itself and rank.cut
#' is only the fallback when that fails
rnbeadsDensityScatterRankCut <- function(diffMeth, region, comparison = 1,
                                         rank.cut = 1000, auto.rank.cut = FALSE) {
  cc <- .diffmeth_comparison(diffMeth, comparison)
  df2p <- get.table(diffMeth, cc$cmp, region, return.data.frame = TRUE)
  grp.names <- get.comparison.grouplabels(diffMeth)[cc$cmp, ]

  if (!"combinedRank" %in% colnames(df2p)) {
    stop(
      "No combinedRank column in the ", region, " table of ", cc$name,
      ". Present: ", paste(colnames(df2p), collapse = ", ")
    )
  }

  cut_used <- rank.cut
  if (auto.rank.cut) {
    p_col <- intersect(c("comb.p.val", "comb.p.adj.fdr"), colnames(df2p))[1]
    cut_used <- tryCatch(
      RnBeads:::auto.select.rank.cut(df2p[[p_col]], df2p$combinedRank, alpha = 0.1),
      error = function(e) {
        message("auto rank cut failed (", conditionMessage(e), "), using ", rank.cut)
        rank.cut
      }
    )
  }
  n_special <- sum(!is.na(df2p$combinedRank) & df2p$combinedRank <= cut_used)
  message(
    cc$name, " / ", region, ": rank cut ", cut_used, ", ",
    n_special, " of ", nrow(df2p), " regions highlighted"
  )
  df2p$isDMP <- !is.na(df2p$combinedRank) & df2p$combinedRank <= cut_used

  mcols <- .diffmeth_mean_cols(df2p)
  create.densityScatter(df2p[, c(mcols[2], mcols[1])],
    is.special = df2p$isDMP,
    dens.subsample = TRUE, sparse.points = 0.001, add.text.cor = TRUE
  ) +
    labs(
      x = paste0("Mean-Methylation ( ", grp.names[2], " )"),
      y = paste0("Mean-Methylation ( ", grp.names[1], " )")
    ) +
    coord_fixed() +
    theme_classic() +
    theme(legend.position = "none")
}

#' MA plot of the group mean methylation for one comparison.
#'
#' A is the average methylation of the two groups, M their difference, so a
#' point far from the horizontal zero line is a region that differs between
#' the groups and its x position says whether that happens in a lowly or
#' highly methylated part of the genome.
#'
#' @param diffMeth an RnBeads differential methylation object
#' @param region region type, e.g. "tiling" or "sites"
#' @param comparison name or index of the comparison to plot
#' @param p.cut adjusted p-value cut used to colour the significant points
#' @param point.colors named vector for the significant and background points
maPlot <- function(diffMeth, region, comparison = 1,
                   rank.cut = 1000, auto.rank.cut = FALSE, p.cut = NULL,
                   point.colors = c(
                     "Differential" = "#00441B",
                     "Not differential" = "grey75"
                   )) {
  cc <- .diffmeth_comparison(diffMeth, comparison)
  tbl <- get.table(diffMeth, cc$cmp, region, return.data.frame = TRUE)
  grp.names <- get.comparison.grouplabels(diffMeth)[cc$cmp, ]
  mcols <- .diffmeth_mean_cols(tbl)

  df <- data.frame(
    A = (tbl[[mcols[1]]] + tbl[[mcols[2]]]) / 2,
    M = tbl[[mcols[1]]] - tbl[[mcols[2]]]
  )

  # Differential by rank cut, matching rnbeadsDensityScatterRankCut and the
  # rankCut_<n> user sets of the LOLA enrichment. Pass p.cut instead to fall
  # back to the adjusted p-value.
  if (!is.null(p.cut)) {
    if (!"comb.p.adj.fdr" %in% colnames(tbl)) {
      stop("No comb.p.adj.fdr column in the ", region, " table of ", cc$name)
    }
    df$sig <- ifelse(tbl$comb.p.adj.fdr < p.cut, "Differential", "Not differential")
    sub <- paste0("adjusted p < ", p.cut)
  } else {
    if (!"combinedRank" %in% colnames(tbl)) {
      stop("No combinedRank column in the ", region, " table of ", cc$name)
    }
    cut_used <- rank.cut
    if (auto.rank.cut) {
      p_col <- intersect(c("comb.p.val", "comb.p.adj.fdr"), colnames(tbl))[1]
      cut_used <- tryCatch(
        RnBeads:::auto.select.rank.cut(tbl[[p_col]], tbl$combinedRank, alpha = 0.1),
        error = function(e) rank.cut
      )
    }
    df$sig <- ifelse(!is.na(tbl$combinedRank) & tbl$combinedRank <= cut_used,
      "Differential", "Not differential"
    )
    sub <- paste0("top ", cut_used, " by combined rank")
  }
  df$sig <- factor(df$sig, levels = names(point.colors))
  df <- df[is.finite(df$A) & is.finite(df$M), ]

  # Differential points last so they are not buried under the background
  df <- df[order(df$sig, decreasing = TRUE), ]

  ggplot(df, aes(x = A, y = M, colour = sig)) +
    geom_point(size = 0.6, alpha = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "black", linewidth = 0.3) +
    scale_colour_manual(values = point.colors, name = NULL) +
    labs(
      x = paste0("Mean methylation ( ", grp.names[1], " + ", grp.names[2], " ) / 2"),
      y = paste0("Methylation difference ( ", grp.names[1], " - ", grp.names[2], " )"),
      title = paste0(grp.names[1], " vs ", grp.names[2], ", ", region),
      subtitle = sub
    ) +
    theme_classic(base_size = 13) +
    theme(legend.position = "bottom") +
    guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))
}

lolaBarPlot.hyp <- function(lolaDb, lolaRes, signifCol = "qValue", scoreCol = "pValueLog",
                            orderCol = scoreCol, includedCollections = c(), pvalCut = 0.05,
                            maxTerms = 50, colorpanel = c(brewer.pal(2, "Set2")),
                            groupByCollection = TRUE, orderDecreasing = NULL) {
  isHypo <- grepl("hypo", lolaRes[["userSet"]], ignore.case = TRUE)
  isHyper <- grepl("hyper", lolaRes[["userSet"]], ignore.case = TRUE)

  if (sum(isHypo) < 1) logger.error("LOLA Result does not contain hypOmethylated userSet")
  if (sum(isHyper) < 1) logger.error("LOLA Result does not contain hypERmethylated userSet")

  lolaRes$userSet[isHypo] <- "hypomethylated"
  lolaRes$userSet[isHyper] <- "hypermethylated"

  if (sum(isHypo) + sum(isHyper) != nrow(lolaRes)) {
    lolaRes <- lolaRes[isHypo | isHyper, ]
    logger.warning("LOLA Result contains userSets not annotated as hyper or hypo")
  }

  if (scoreCol == "logOddsRatio") {
    logger.warning("In newer versions of LOLA the odds ratio column is called 'oddsRatio'")
    scoreCol <- "oddsRatio"
  }

  df2p <- muRtools:::lolaPrepareDataFrameForPlot(lolaDb, lolaRes,
    signifCol = signifCol, scoreCol = scoreCol, orderCol = orderCol,
    includedCollections = includedCollections, pvalCut = pvalCut,
    maxTerms = maxTerms, perUserSet = TRUE,
    groupByCollection = groupByCollection, orderDecreasing = orderDecreasing
  )
  if (is.null(df2p)) {
    return(rnb.message.plot("No significant association found"))
  }
  df2p$name <- sub("\\s\\[.*", "", df2p$name)

  # The hypo and hyper bars now use the shared green scheme instead of
  # viridis, so this panel matches the rest of the figure
  cpanel <- c("hypomethylated" = "#C7E9C0", "hypermethylated" = "#00441B")

  # The trailing + before the if() used to append the value of the if block
  # to the plot, which only worked by accident. The facet is added normally.
  pp <- ggplot(df2p) +
    aes_string("name", scoreCol, fill = "userSet") +
    geom_col(position = "dodge") +
    scale_x_discrete(name = "") +
    scale_y_continuous(name = "Odds-Ratio") +
    scale_fill_manual(values = cpanel) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      legend.position = "bottom"
    )

  if (groupByCollection) {
    pp <- pp + facet_grid(. ~ collection, scales = "free", space = "free")
  }

  return(pp)
}
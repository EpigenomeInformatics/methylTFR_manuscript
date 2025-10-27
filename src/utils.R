rnbeadsDensityScatter <- function(diffMeth, region) {
  # the rank cuts
  rank.cuts.auto <- 0

  for (i in 1:length(get.comparisons(diffMeth))) {
    cc <- names(get.comparisons(diffMeth))[i]
    ccc <- get.comparisons(diffMeth)[cc]
    dmt <- get.table(diffMeth, ccc, region, return.data.frame = TRUE)
    ccn <- ifelse(RnBeads:::is.valid.fname(cc), cc, paste("cmp", i, sep = ""))
    grp.names <- get.comparison.grouplabels(diffMeth)[ccc, ]
    auto.rank.cut <- rank.cuts.auto[[i]]
    df2p <- dmt # data frame to plot
    ChrAccR:::cleanMem()
  }
  # scatterplot based on adjusted p-value significance
  if (is.element("comb.p.adj.fdr", colnames(df2p))) {
    df2p$isDMP <- df2p[, "comb.p.adj.fdr"] < 0.05
    sparse.points <- 0.001
    dens.subsample <- TRUE
    pp <- create.densityScatter(df2p[, c("mean.mean.g2", "mean.mean.g1")],
      is.special = df2p$isDMP,
      dens.subsample = dens.subsample, sparse.points = sparse.points, add.text.cor = TRUE
    ) +
      labs(
        x = paste("Mean-Methylation (", grp.names[2], ")", sep = " "),
        y = paste("Mean-Methylation (", grp.names[1], ")", sep = "  ")
      ) +
      coord_fixed() +
      theme_classic() +
      theme(legend.position = "none")
  }
  return(pp)
}

lolaBarPlot.hyp <- function(lolaDb, lolaRes, signifCol = "qValue", scoreCol = "pValueLog", orderCol = scoreCol, includedCollections = c(), pvalCut = 0.05, maxTerms = 50, colorpanel = c(brewer.pal(2, "Set2")), groupByCollection = TRUE, orderDecreasing = NULL) {
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
    logger.warning("In newer versions of LOLA the odds ratio column is called 'oddsRatio' (no longer 'logOddsRatio')")
    scoreCol <- "oddsRatio"
  }

  df2p <- muRtools:::lolaPrepareDataFrameForPlot(lolaDb, lolaRes, signifCol = signifCol, scoreCol = scoreCol, orderCol = orderCol, includedCollections = includedCollections, pvalCut = pvalCut, maxTerms = maxTerms, perUserSet = TRUE, groupByCollection = groupByCollection, orderDecreasing = orderDecreasing)
  df2p$name <- sub("\\s\\[.*", "", df2p$name)

  if (is.null(df2p)) {
    return(rnb.message.plot("No significant association found"))
  }

  cpanel <- viridis::viridis(2)
  aesObj <- aes_string("name", scoreCol, fill = "userSet")

  pp <- ggplot(df2p) +
    aesObj +
    geom_col(position = "dodge") +
    scale_x_discrete(name = "") +
    scale_fill_manual(values = cpanel) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      legend.position = "bottom"
    ) +
    scale_y_continuous(name = "Odds-Ratio") + # Add this line


    if (groupByCollection) {
      pp <- pp + facet_grid(. ~ collection, scales = "free", space = "free")
    }

  return(pp)
}

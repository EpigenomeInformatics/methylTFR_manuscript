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
        labs(x = paste("Mean-Methylation (", grp.names[2], ")", sep = " "),
        y = paste("Mean-Methylation (", grp.names[1], ")", sep = "  ")) +
        coord_fixed() +
        theme_classic() +
        theme(legend.position = "none")
    }
  return(pp)
}

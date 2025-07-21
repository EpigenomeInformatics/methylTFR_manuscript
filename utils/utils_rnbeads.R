
rnbeadsDensityScatter <- function(diffMeth, region) {
    # Initialize a list to store plots
    plot_list <- list()
    rank.cuts.auto <- 0 
    # Iterate over each comparison
    for (i in 1:length(get.comparisons(diffMeth))) {
        cc <- names(get.comparisons(diffMeth))[i]
        ccc <- get.comparisons(diffMeth)[cc]
        dmt <- get.table(diffMeth, ccc, region, return.data.frame = TRUE)
        ccn <- ifelse(RnBeads:::is.valid.fname(cc), cc, paste("cmp", i, sep = ""))
        grp.names <- get.comparison.grouplabels(diffMeth)[ccc, ]
        auto.rank.cut <- rank.cuts.auto
        df2p <- dmt # data frame to plot
        ChrAccR:::cleanMem()

        # Generate scatterplot if the required column exists
        if (is.element("comb.p.adj.fdr", colnames(df2p))) {
            df2p$isDMP <- df2p[, "comb.p.adj.fdr"] <= 0.05
            sparse.points <- 0.01
            dens.subsample <- FALSE
            pp <- create.densityScatter(df2p[, c("mean.mean.g2", "mean.mean.g1")],
                is.special = df2p$isDMP,
                dens.subsample = dens.subsample, sparse.points = sparse.points, add.text.cor = TRUE
            ) +
                labs(x = paste("Mean-Methylation (", grp.names[2], ")", sep = " "),
                     y = paste("Mean-Methylation (", grp.names[1], ")", sep = "  ")) +
                coord_fixed() +
                theme_classic() +
                theme(legend.position = "none")

            # Add the plot to the list
            plot_list[[ccn]] <- pp
        }
    }

    # Return the list of plots
    return(plot_list)
}

rnbeadsDensityScatterRankCut <- function(diffMeth, region = "distal", alpha = 0.1, sparse.points = 0.01,
                                         dens.subsample = TRUE) {
    # Initialize a list to store plots
    plot_list <- list()
    comps <- get.comparisons(diffmeth)
    rank.cuts.auto <- lapply(1:length(comps),FUN=function(i){
			dmt <- get.table(diffmeth,comps[i],region,return.data.frame=TRUE)
			res <- RnBeads:::auto.select.rank.cut(dmt$comb.p.adj.fdr,dmt$combinedRank,alpha=0.1)
			return(as.integer(res))
	})
    tt <- data.frame(unlist(rank.cuts.auto))
colnames(tt) <- c("Rank Cutoff")
rownames(tt) <- get.comparisons(diffmeth)

    # Iterate over each comparison
    for (i in 1:length(get.comparisons(diffMeth))) {
        cc <- names(get.comparisons(diffMeth))[i]
        ccc <- get.comparisons(diffMeth)[cc]
        dmt <- get.table(diffMeth, ccc, region, return.data.frame = TRUE)
        ccn <- ifelse(RnBeads:::is.valid.fname(cc), cc, paste("cmp", i, sep = ""))
        grp.names <- get.comparison.grouplabels(diffMeth)[ccc, ]
        df2p <- dmt # data frame to plot
        ChrAccR:::cleanMem()
        autoRankCut <- tt[i,]

        # Generate scatterplot if the required column exists
        if (is.element("combinedRank", colnames(df2p))) {
            df2p$isDMP <- df2p[,"combinedRank"] <= autoRankCut
            sparse.points <- 0.01
            dens.subsample <- FALSE
            pp <- create.densityScatter(df2p[, c("mean.mean.g2", "mean.mean.g1")],
                is.special = df2p$isDMP,
                dens.subsample = dens.subsample, sparse.points = sparse.points, add.text.cor = TRUE
            ) +
                labs(x = paste("Mean-Methylation (", grp.names[2], ")", sep = " "),
                     y = paste("Mean-Methylation (", grp.names[1], ")", sep = "  ")) +
                coord_fixed() +
                theme_classic() +
                theme(legend.position = "none")

            # Add the plot to the list
            plot_list[[ccn]] <- pp
        }
    }

    # Return the list of plots
    return(plot_list)
}

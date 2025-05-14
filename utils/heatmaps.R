collapseMotifMatrix <- function(X, motifClust=NULL, assembly="hg38", motifs="jaspar", aggrFun=mean){
	if (motifs != "jaspar") logger.error(c("Currently motif clustering is only supported for JASPAR motifs"))

	if (is.null(motifClust)){
        logger.error("motifClust must be provided")
	}
	names(motifClust$clustAssign) <- sub(".*_", "", names(motifClust$clustAssign))
	# Filter the motifs present in motifClust
	X <- X[rownames(X) %in% names(motifClust$clustAssign),]
	if (is.null(rownames(X)) || !all(rownames(X) %in% names(motifClust$clustAssign))){
		logger.error("X must have valid rownames that are all contained in the clustering")
	}

	df <- reshape2::melt(X)
	colnames(df)[1:2] <- c("motif", "sample")
	df[,"cluster"] <- motifClust$clustAssign[as.character(df[,"motif"])]

	Xc <- reshape2::acast(df, formula=cluster~sample, fun.aggregate=aggrFun, value.var="value")

	if (!is.null(motifClust$clustNames)){
		rownames(Xc) <- motifClust$clustNames[rownames(Xc)]
	}

	res <- list(
		collapsed = Xc,
		clustering = motifClust
	)
	class(res) <- "CollapsedMotifMatrix"
	return(res)
}

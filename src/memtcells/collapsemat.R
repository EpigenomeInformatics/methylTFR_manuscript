library(methylTFR)
library(ComplexHeatmap)
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

# Example usage
mtfr <- readRDS("/icbb/projects/nitschre/methylTFR/r_objects/jaspar2020_distal_deviations.RDS")
motifClust <- readRDS("/icbb/projects/igunduz/irem_github/methylTFR_manuscript/utils/motifclust.RDS")	
X <- deviations(mtfr)
Xc <- collapseMotifMatrix(X, motifClust=motifClust, aggrFun=mean)
head(Xc$collapsed)

# Rowwise zscoring
Xc_zscores <- methylTFR:::computeRowZScore(Xc$collapsed) #all 110
Xc_zscores <- Xc_zscores[1:80,]

# Heatmap
colnames(Xc_zscores) <- c("Hf03_CM", "Hf03_EM", "Hf03_TN", "Hf04_CM", "Hf04_EM",
                                "Hf04_TN", "Hf04_TEMRA")

ha <- HeatmapAnnotation(
  celltypes=c("CM","EM", "TN", "CM", "EM","TN","TEMRA"),
  col = list(celltypes = c("TN" = "#ff0000", "TEMRA" = "#ff6600", "CM"="#008cff", "EM"="#0037ff"))
)

rownames(Xc_zscores) <- sub(":.*", "", rownames(Xc_zscores))

plot_dir <- "/icbb/projects/nitschre/methylTFR/figures"
path <- file.path(plot_dir, "heatmap_collapsed.pdf")

pdf(path)
Heatmap(
  Xc_zscores,
  row_names_gp = gpar(fontsize = 5),
  top_annotation = ha,
  column_title = "Z-Scores of collapsed matrix (80 of 110)",
  show_row_names = TRUE,
  column_names_side = "top")
dev.off()
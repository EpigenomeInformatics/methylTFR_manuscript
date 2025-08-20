library(MOFA2)
library(ggplot2)

set.seed(12)

plot_dir <- "/icbb/projects/nitschre/methylTFR/figures/"

# Load pseudobulk TF activity matrices
chromVar <- readRDS(
  gzfile("/icbb/projects/nitschre/methylTFR/scripts/other/cvar_zscores_adj_psuedobulk.R")
)

mtfr <- readRDS(
    gzfile("/icbb/projects/nitschre/methylTFR/scripts/other/mtfr_zscores_adj_psuedobulk.R")
)

# Create Mofa object
data <- list(chromVar = chromVar, mtfr = mtfr)
MOFAobject <- create_mofa(data)

# Add meta data
metadata <- data.frame(
    celltype = sub("_.*", "", colnames(mtfr)),
    sample = colnames(mtfr)
    )

samples_metadata(MOFAobject) <- metadata

# Train model
MOFAobject <- prepare_mofa(MOFAobject)
outfile = file.path("/icbb/projects/nitschre/methylTFR/r_objects","model.hdf5")
MOFAobject.trained <- run_mofa(MOFAobject, outfile, use_basilisk = TRUE)

# Plot factors
p <- plot_factor(MOFAobject.trained, 
  factor = c(1,2,3,4,7),
  color_by = "celltype",
)
ggsave(paste0(plot_dir, "factors.pdf"), p)

# ANOVA
factors <- get_factors(MOFAobject.trained, as.data.frame=TRUE)
factors$celltype <- metadata$celltype 

res <- data.frame(
  Factor = character(),
  R2 = numeric()  
)

for(fa in unique(factors$factor)){
  df <- subset(factors,factor == fa)
  anova <- aov(value ~ celltype, data=df)

  # Caluclate R2 for each factor
  sum_sq_factor <- summary(anova)[[1]]["celltype", "Sum Sq"]
  sum_sq_res <- summary(anova)[[1]]["Residuals", "Sum Sq"]
  r2 <- sum_sq_factor/(sum_sq_res+sum_sq_factor)

  # Save in df
  res <- rbind(res, data.frame(Factor = fa, R2 = r2))
}
res <- res[order(res$R2, decreasing = TRUE),]
write.csv(res, "/icbb/projects/nitschre/methylTFR/r_objects/r2Values.csv")


# Variance explained by each modality
var_expl <- get_variance_explained(MOFAobject.trained)$r2_per_factor[[1]]
var_top <- subset(var_expl, rownames(var_expl) %in% res$Factor[1:5])
write.csv(var_expl, "/icbb/projects/nitschre/methylTFR/r_objects/var_expl_all.csv")

# Top TFs for each modality
weights <- get_weights(MOFAobject.trained, as.data.frame = TRUE)
weights$value_abs <- abs(weights$value)

topTfs <- data.frame(
  feature = character(),
  value = numeric(),
  factor = character()
)

for(f in res$Factor[1:5]){
  for(modality in unique(weights$view)){
    df <- subset(weights, view == modality & factor == f)
    df <- df[order(df$value_abs, decreasing = TRUE),]

    new_df <- df[1:5, c("feature", "value", "factor")]

    topTfs <- rbind(topTfs, new_df)
  }
}
write.csv(topTfs, "/icbb/projects/nitschre/methylTFR/r_objects/topTFs.csv")
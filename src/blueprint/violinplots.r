suppressPackageStartupMessages({
    library(RnBeads)
    library(dplyr)
    library(methylTFR)})


set.seed(12)
plotdir <- "/icbb/projects/nitschre/methylTFR/figures/blueprint/Violinplots/"
# Load Selex data
selex <- read.csv("/icbb/projects/igunduz/exposure_atlas_manuscript/sample_annots/Selex_data.csv", skip = 20, header = 21, sep = ";")[, c(1, 2, 3, 4, 6)]
selex <- selex[, c("TF.name", "methyl.SELEX.call")]

mtfr <- readRDS("/icbb/projects/nitschre/methylTFR/r_objects/blueprintjaspar2020_distal_deviations.RDS")
mtfr <- deviations(mtfr)
mtfr <- as.data.frame(mtfr)
mtfr$TF.name <- rownames(mtfr)


# Match Selex labels with Tfs
df <- merge(mtfr, selex, by="TF.name")

# To long format
df_long <- df %>%
  pivot_longer(
    cols = -c(TF.name, methyl.SELEX.call),   
    names_to = "bedFile",
    values_to = "Activity"
  ) %>%
  rename(SelexLabel = methyl.SELEX.call)

# Celltypes added
sannot <- read.csv("/icbb/projects/igunduz/methylTFR_manuscript/results/BLUEPRINT_080725/reports/data_import_data/annotation.csv", stringsAsFactors = FALSE)
df_long <- merge(df_long, sannot[, c("bedFile", "cellTypeGroup")], by="bedFile")

# Violin Plot
for(celltype in unique(df_long$cellTypeGroup)){
    df <- subset(df_long, cellTypeGroup==celltype & SelexLabel %in% c("MethylMinus", "MethylPlus"))
    
    p <- ggplot(df, aes(x=SelexLabel, y=Activity, fill=SelexLabel))+
            geom_violin()+
            theme_classic()+
            ggtitle(paste("TF activity by SELEX label across ", celltype))

    ggsave(filename = paste0(plotdir,celltype, "_violin_plot.pdf"), plot = p)

}

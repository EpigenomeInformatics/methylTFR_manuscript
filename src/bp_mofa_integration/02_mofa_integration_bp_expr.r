#!/usr/bin/env Rscript

#####################################################################
# mofa_integration.R
# Created by IBG on 25-10-2025
# Script to perform TF activity integration using MOFA2
#####################################################################
suppressPackageStartupMessages({
library(MOFA2)
library(dplyr)
library(tidyr)
library(viridis)
library(tibble)
library(patchwork)
library(ggrepel)
library(scales)
library(methylTFR)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
})

set.seed(12)

source("/icbb/projects/nitschre/methylTFR/scripts/other/helpers.R")

# Define custom colors for cell types
cell_type_colors <- c(
  "Bcell" = "#980043",
  "plasma" = "#CD2990",
  "DC" = "#EED5B7",
  "other" = "#8B8682",
  "gran" = "#ff7f50",
  "eryt" = "#67000d",
  "Mf" = "#864a38",
  "mono" = "#CD7054",
  "megK" = "#4a0221",
  "NK" = "#bf812d",
  "osteoclast" = "#DEB887",
  "Tcell" = "#40E0D0",
  "thymocyte" = "#74c476",
  "Progenitors" = "#df65b0"
)

# Paths
plot_dir <- "/icbb/projects/nitschre/methylTFR/figures/figure3_with_expr/"
r_objects_dir <- "/icbb/projects/nitschre/methylTFR/r_objects/"
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(r_objects_dir, recursive = TRUE, showWarnings = FALSE)

# Load two modalities: mtfr and rna 
mtfr <- readRDS("/scratch/icbb/igunduz/methylTFR_manuscript/blueprint/mTFR_devs_311025/jaspar2020_distal_deviations.RDS")
rna <- readRDS("/icbb/projects/nitschre/methylTFR/r_objects/bp_rawRNAcounts.RDS")

# Subset that rna + mtfr have same samples
annot <- read.csv("/icbb/projects/nitschre/methylTFR/sample_annotation/annotated_wgbs_with_rna_chip_matches.csv", stringsAsFactors=FALSE)
annot <- annot %>%
  select(bedFile, matching_rna_files, cellTypeGroup) %>%                                 
  mutate(ids = sub(".*\\.([A-Z0-9]+\\.[0-9]+\\.[a-z0-9-]+)\\.genes\\.results$", "\\1", matching_rna_files), # Get IDs
  pos_in_rna = match(ids, colnames(rna))) %>%  # Get position of ids in rna                      
  filter(!is.na(pos_in_rna)) %>% # Only keep samples that are in RNA
  distinct(ids, .keep_all = TRUE) %>% # Remove duplicate ids
  arrange(pos_in_rna) %>% # Same order as in rna
  select(-pos_in_rna) # remove pos column again


# MTFR preprocessing
mtfr <- computeRowZScore(deviations(mtfr))
mtfr <- mtfr[, colnames(mtfr) %in% annot$bedFile] # Only keep samples that are in subsetted annotation file
mtfr <- mtfr[,match(annot$bedFile, colnames(mtfr))] # Same order as in RNA samples
colnames(mtfr) <- colnames(rna) # Same sample names as in RNA

# Subset 20 000 TFs in RNA to have only TFs also in mtfr
rna <- rna[intersect(rownames(rna), rownames(mtfr)),]

# RNA zscores
rna <- computeRowZScore(as.matrix(rna))

# Quick sanity checks
stopifnot(is.matrix(mtfr) | is.data.frame(mtfr))
stopifnot(is.matrix(rna) | is.data.frame(rna))
if(!all(colnames(mtfr) == colnames(rna))){
  warning("Column names (samples) do not match exactly between rna and mtfr.")
}

# Create MOFA object and metadata
data_list <- list(mtfr = as.matrix(mtfr), expr = as.matrix(rna))
MOFAobject <- create_mofa(data_list)

metadata <- data.frame(
  sample = colnames(mtfr),
  celltype = annot$cellTypeGroup,
  stringsAsFactors = FALSE
)
samples_metadata(MOFAobject) <- metadata

# Reduce number of factors (because not enough samples)
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 14

# Prepare and train MOFA (save HDF5)
MOFAobject <- prepare_mofa(MOFAobject, model_options = model_opts)
outfile <- file.path(r_objects_dir, "mtfr_expr_model_bp.hdf5")

# If model file exists, load trained object; otherwise run training
if(file.exists(outfile)){
  MOFAobject.trained <- load_model(outfile)
} else {
  MOFAobject.trained <- run_mofa(MOFAobject, outfile, use_basilisk = TRUE)
  saveRDS(MOFAobject.trained, file.path(r_objects_dir, "mtfr_expr_model_bp.hdf5"))

}

# Remove sample "others"
metadata <- samples_metadata(MOFAobject.trained)
metadata <- metadata[metadata$celltype != "other",]
samples_to_keep <- metadata$sample
MOFAobject.trained <- subset_samples(MOFAobject.trained, samples_to_keep)


# Extract factors and compute ANOVA R2 for celltype
factors_long <- get_factors(MOFAobject.trained, as.data.frame = TRUE)
factors_long$celltype <- metadata$celltype

r2_df <- data.frame(Factor = character(), R2 = numeric(), stringsAsFactors = FALSE)

for(fa in unique(factors_long$factor)){
  df <- subset(factors_long, factor == fa)
  mod <- aov(value ~ celltype, data = df)
  ss_tab <- summary(mod)[[1]]
  ss_factor <- ss_tab["celltype", "Sum Sq"]
  ss_resid  <- ss_tab["Residuals", "Sum Sq"]
  r2_val <- ss_factor / (ss_factor + ss_resid)
  r2_df <- base::rbind(r2_df, data.frame(Factor = fa, R2 = r2_val, stringsAsFactors = FALSE))
}
r2_df <- r2_df %>% arrange(desc(R2))
write.csv(r2_df, file.path(plot_dir, "r2Values_bp.csv"), row.names = FALSE)


# Choose top factors (by R2) for plotting
r2_df <- r2_df[r2_df$Factor!="Factor1",] # Remove Factor 1 because it is highly correlated with total numbers of expressed features
top_n <- 7
top_factors <- head(r2_df$Factor, top_n)

# Plot 1: Factor scatter (example F1 vs F2) colored by celltype
factors_wide <- factors_long %>%
  pivot_wider(names_from = factor, values_from = value) %>%
  subset(celltype != "other") # Remove "other" cell type for plotting

# default pair: top two factors by R2 (if <2 fallback to 1,2)
f_x <- ifelse(length(top_factors) >= 1, top_factors[1], "1")
f_y <- ifelse(length(top_factors) >= 2, top_factors[2], "2")

if(!(f_x %in% colnames(factors_wide)) | !(f_y %in% colnames(factors_wide))){
  # fallback to numeric factor names if needed
  available_factors <- setdiff(colnames(factors_wide), "celltype")
  f_x <- available_factors[1]; f_y <- available_factors[2]
}

p_factors_scatter <- ggplot(factors_wide, aes_string(x = f_x, y = f_y, color = "celltype")) +
  geom_point(alpha = 0.9, size = 2) +
  stat_ellipse(aes(group = celltype), linetype = 2, alpha = 0.5, size = 0.4) +
  scale_color_manual(values = cell_type_colors) +  # Use custom colors
  labs(x = paste0("Factor ", f_x), y = paste0("Factor ", f_y),
       title = paste0("MOFA: Factor ", f_x, " vs ", f_y, " (colored by celltype)")) +
  theme_classic(base_size = 13) +
  theme(legend.position = "right")

ggsave(filename = file.path(plot_dir, "factor_scatter_Fx_Fy.pdf"), plot = p_factors_scatter,
       width = 7, height = 5)

# Variance explained per view × factor (from MOFA)
ve <- get_variance_explained(MOFAobject.trained)

# Convert to long table (robust for MOFA2 output)
r2_mat <- as.data.frame(ve$r2_per_factor[[1]])
r2_mat$view <- rownames(r2_mat)
r2_mat <- r2_mat[r2_mat$view != "Factor1",] # Because factor1 highly correlated with num of features
r2_long <- r2_mat %>% pivot_longer(cols = -view, names_to = "factor", values_to = "r2")

# Keep only top factors (from R2 ranking)
r2_df <- r2_df[r2_df$Factor != "Factor1",]
top_factors <- head(r2_df$Factor, 7)
r2_long <- r2_long %>% filter(view %in% top_factors)
r2_long$factor <- factor(r2_long$view, levels = top_factors)

# Get weights and prepare modality contribution plots
weights_df <- get_weights(MOFAobject.trained, as.data.frame = TRUE)
weights_df <- weights_df %>% mutate(value_abs = abs(value), value_signed = value)

# Aggregate absolute weights per view × factor and compute fraction per factor
agg <- weights_df %>%
  filter(factor %in% top_factors) %>%
  group_by(factor, view) %>%
  summarise(sum_abs = sum(value_abs, na.rm = TRUE), .groups = "drop") %>%
  group_by(factor) %>%
  mutate(frac = sum_abs / sum(sum_abs)) %>%
  ungroup()

agg$factor <- factor(agg$factor, levels = c("Factor3", "Factor8", "Factor2", "Factor6", "Factor12", "Factor7", "Factor5"))

# Plot stacked contribution (fraction)
p_modality_frac <- ggplot(agg, aes(x = factor, y = frac, fill = view)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = ifelse(frac > 0.03, paste0(round(frac*100,1), "%"), "")),
            position = position_stack(vjust = 0.5), size = 3) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_fill_manual(values = c("#ED4B4A", "#6BC75A")) +
  labs(x = "Factor", y = "Relative contribution (|weights| fraction)", fill = "View",
       title = "Modality contribution per factor (weights-based)") +
  theme_classic(base_size = 12)

ggsave(filename = file.path(plot_dir, "mtfr_modality_contribution_fraction_topFactors.pdf"), plot = p_modality_frac,
       width = 8, height = 4)

# Filter to top factors (optional)
#top_factors <- c("Factor1","Factor2","Factor3","Factor4","Factor7") 
factors_long <- factors_long %>% 
  filter(factor %in% top_factors) %>%
  subset(celltype != "other") # Remove "other" cell type for plotting

# Plot: strip / dot plot per factor
factors_long$factor <- factor(factors_long$factor, levels = c("Factor3", "Factor8", "Factor2", "Factor6", "Factor12", "Factor7", "Factor5"))

p <- ggplot(factors_long, aes(x = factor, y = value, color = celltype)) +
  geom_jitter(width = 0.2, height = 0, size = 2, alpha = 0.8) +
  scale_color_manual(values = cell_type_colors) +  # Use custom colors
  labs(x = "Factor", y = "Factor value", color = "Cell type",
       title = "MOFA factor values per sample by cell type") +
  theme_classic(base_size = 13) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())
ggsave(filename = file.path(plot_dir, "factors_stripplot_by_celltype.pdf"), plot = p,
       width = 8, height = 5)


# Top TFs for each modality
weights <- get_weights(MOFAobject.trained, as.data.frame = TRUE)
weights$value_abs <- abs(weights$value)

topTfs <- data.frame(
  feature = character(),
  value = numeric(),
  factor = character(),
  view = character()
)

# Create DF with top  5 TFs for each factor and view
for(f in r2_df$Factor[1:5]){
  for(modality in unique(weights$view)){
    df <- subset(weights, view == modality & factor == f)
    df <- df[order(df$value_abs, decreasing = TRUE),]

    new_df <- df[1:5, c("feature", "value", "factor", "view")]
    new_df$feature <- gsub("_.*$", "", new_df$feature)

    topTfs <- base::rbind(topTfs, new_df)
  }
}
write.csv(topTfs, file.path(plot_dir, "topTFs.csv"), row.names=FALSE)

## Heatmaps
# Filter for unique TFs
mofaTFs <- unique(topTfs$feature)

# Filter for only Mofa TFs that are present in mtfr and rna
common_TFs <- Reduce(intersect, list(mofaTFs, rownames(rna), rownames(mtfr)))
rna_filtered <- rna[common_TFs,]
mtfr_filtered <- mtfr[common_TFs,]

# Extract cell types 
groups <- annot$cellTypeGroup

# Change colnames
colnames(mtfr_filtered) <- groups
colnames(rna_filtered) <- groups

# Remove "other" sample + plasma
mtfr_filtered <- mtfr_filtered[, !colnames(mtfr_filtered) %in% c("other", "plasma")]
rna_filtered <- rna_filtered[, !colnames(rna_filtered) %in% c("other", "plasma")]
groups <- groups[groups != "other" & groups != "plasma"]

# Annotation
ha <- HeatmapAnnotation(
  celltypes=colnames(mtfr_filtered),
  col = list(celltypes = cell_type_colors
))
# Column order
column_order <- c("megK", "eryt", "gran", "mono", "Mf", "DC", "osteoclast", "NK", "Tcell", "thymocyte", "Bcell")

# Column groups
column_split_factor <- factor(groups, levels = column_order)

# MTFR
p<-Heatmap(
  mtfr_filtered,
  row_names_gp = gpar(fontsize = 4),
  top_annotation = ha,
  column_title = "mtfr Z-Scores of TFs from Mofa analysis",
  show_row_names = TRUE,
  show_column_names = FALSE,
  column_split = column_split_factor, # Split columns by the desired order
  cluster_column_slices = FALSE # Prevent clustering of the groups themselves)
)
pdf(file.path(plot_dir, "heatmap_mofa_mtfr_bp.pdf"))
draw(p)
dev.off()

# expr
# Get row order
p <- draw(p)
row_order_indices <- row_order(p)  

# Set color scheme
col_fun <- colorRamp2(seq(-2,2, length.out = 100), viridis(100))

pdf(file.path(plot_dir, "heatmap_mofa_expr_bp.pdf"))
Heatmap(
  rna_filtered,
  row_names_gp = gpar(fontsize = 4),
  top_annotation = ha,
  column_title = "rna Z-Scores of TFs from Mofa analysis",
  show_row_names = TRUE,
  show_column_names = FALSE,
  column_split = column_split_factor, # Split columns by the desired order
  cluster_column_slices = FALSE, # Prevent clustering of the groups themselves),
  row_order = row_order_indices,
  col = col_fun)
dev.off()

# Calculate row-wise correlation
row_correlation <- sapply(seq_len(nrow(mtfr_filtered)), function(i) {
  cor(mtfr_filtered[i, ], rna_filtered[i, ]) # ), use = "complete.obs")  # Only use rows with complete observations
})

# Convert row-wise correlation vector to a matrix for heatmap plotting
correlation_matrix <- as.matrix(row_correlation)
rownames(correlation_matrix) <- rownames(mtfr_filtered) # Add row names for better interpretation

# Create the heatmap
cols <- colorRampPalette(brewer.pal(11,"PiYG"))(200)
col_fun <- colorRamp2(
  seq(-1, 1, length.out = length(cols)),
  cols
)

cm <- Heatmap(
  correlation_matrix,
  name = "Correlation", # Name for the heatmap legend
  cluster_rows = FALSE, # Cluster rows
  show_row_names = TRUE, # Show row names
  cluster_columns = FALSE, # No clustering for columns since it's a single column
  col = col_fun,
  row_order = row_order_indices,
  heatmap_legend_param = list(title = "Row Correlation") # Legend settings
)

pdf(paste0(plot_dir, "correlation_heatmap.pdf"))
draw(cm)
dev.off()

## Factor column
mofa_filt <- topTfs[!duplicated(topTfs$feature),]

# Set levels for each factor 
levels <- levels(factor(mofa_filt$factor))
level_col <- c("#E41A1C", "#4DAF4A", "#377EB8", "#FF7F00", "#984EA3")
names(level_col) <- levels

# Add col column
mofa_filt$col <- level_col[mofa_filt$factor]

# Order mofa df same as zscores_filtered/chromVar
mofa_filt <- mofa_filt[match(row.names(mtfr_filtered), mofa_filt$feature),]

path <- file.path(plot_dir, "factor_column.pdf")
pdf(path, width=1.5)
Heatmap(
  mofa_filt$factor, column_names_gp = gpar(fontsize = 9),
  cluster_rows = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  row_order = row_order_indices,
  heatmap_legend_param = list(title = "Factors"),
  col = level_col)
dev.off()

# SELEX annotation column
selex <- read.csv("/icbb/projects/igunduz/exposure_atlas_manuscript/sample_annots/Selex_data.csv", skip = 20, header = 21, sep = ";")[, c(1, 2, 3, 4, 6)]
motifs <- rownames(mtfr_filtered)
table(motifs %in% selex$TF.name)

# FALSE  TRUE
#   24    19

# Add a new annotation column based on selex$Call
# Ensure the order of `motifs` matches the heatmap rows
selex_annotation <- data.frame(
  Call = ifelse(motifs %in% selex$TF.name,
    selex$Call[match(motifs, selex$TF.name)],
    "Inconclusive"
  ) # Replace NA with "Inconclusive"
)

# Replace NA values with "Inconclusive"
rownames(selex_annotation) <- motifs

# Set levels for each call 
levels <- levels(factor(selex_annotation$Call))
level_col <- c("#CCCCCC", "#BDB76B", "#8B0000", "#008080")
names(level_col) <- levels

# Add col column
selex_annotation$col <- level_col[selex_annotation$Call]

# Order annotation df same as zscores_filtered/chromVar
selex_annotation <- selex_annotation[match(row.names(mtfr_filtered), rownames(selex_annotation)),]

path <- file.path(plot_dir, "selex_annotation_heatmap.pdf")
pdf(path)
Heatmap(
  selex_annotation$Call, column_names_gp = gpar(fontsize = 9),
  cluster_rows = FALSE,
  show_column_names = FALSE,
  row_order = row_order_indices,
  heatmap_legend_param = list(title = "SELEX Call"), 
  col = level_col,
  width = unit(1, "cm"))
dev.off()
#######################################################################
# Correlation plots for SELEX annotation
#######################################################################

# Sort to match the order of rna
mtfr <- mtfr[rownames(rna), ]

# Calculate row-wise correlation
row_correlation <- sapply(seq_len(nrow(mtfr)), function(i) {
  cor(mtfr[i, ], rna[i, ]) # ), use = "complete.obs")  # Only use rows with complete observations
})

# Convert to a data frame for better readability
row_correlation_df <- data.frame(
  TF.name = rownames(rna),
  Correlation = row_correlation
)

# Add SELEX annotation and filter methylPlus and methylMinus
row_correlation_df <- merge(row_correlation_df, selex, by = "TF.name")
row_correlation_df <- row_correlation_df[row_correlation_df$methyl.SELEX.call %in% c("MethylPlus", "MethylMinus"), ]

# Ensure row_correlation_df is a proper data frame
row_correlation_df <- data.frame(row_correlation_df)

# Count the number of motifs in each class using base R
motif_counts <- as.data.frame(table(row_correlation_df$methyl.SELEX.call))
colnames(motif_counts) <- c("methyl.SELEX.call", "count")


# Create the boxplot with jittered points and annotations
p <- ggplot(row_correlation_df, aes(x = methyl.SELEX.call, y = Correlation, fill = methyl.SELEX.call)) +
  geom_boxplot(color = "black", outlier.shape = NA, width = 0.6) + # Boxplot without outliers
  geom_jitter(aes(color = methyl.SELEX.call), width = 0.2, size = 1.5, alpha = 0.7) + # Add jittered points
  scale_fill_manual(values = c("MethylPlus" = "darkgreen", "MethylMinus" = "darkred")) + # Custom colors for boxplot
  scale_color_manual(values = c("MethylPlus" = "darkgreen", "MethylMinus" = "darkred")) + # Custom colors for points
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.2)) + # Set y-axis limits and breaks
  # Add whisker caps
  stat_boxplot(geom = "errorbar", width = 0.3, size = 0.8) + # Horizontal caps at whisker ends
  labs(
    x = "SELEX Group",
    y = "Correlation"
  ) +
  theme_classic(base_size = 14) + # Use a larger base font size for readability
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"), # Center and bold the title
    axis.title = element_text(face = "bold"), # Bold axis titles
    legend.position = "none" # Remove legend if unnecessary
  ) +
  # Annotate the number of motifs in each class
  geom_text(
    data = motif_counts,
    aes(x = methyl.SELEX.call, y = 1, label = count), # Position annotations above the boxplots
    vjust = -0.5,
    size = 4
  )

# Save the plot
ggsave(
  filename = paste0(plot_dir, "correlation_boxplot_with_caps.pdf"),
  plot = p,
  width = 8,
  height = 6
)

#######################################################################


# Barplot R2 values
p <- ggplot(r2_df,aes(x=reorder(Factor, -R2), y=R2))+
        geom_bar(stat="identity")+
        theme_classic()+
        ylim(0,1)+
        labs(
            title="R2 values of MOFA analysis",
            x="Factors")+
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(paste0(plot_dir, "barplots_R2_2.pdf"), p)



# Save intermediate tables for reproducibility
write.csv(agg, file.path(plot_dir, "modality_contribution_by_factor.csv"), row.names = FALSE)
write.csv(weights_df, file.path(plot_dir, "MOFA_weights_long.csv"), row.names = FALSE)
write.csv(r2_df, file.path(plot_dir, "factor_celltype_R2.csv"), row.names = FALSE)

message("All done. Figures saved to: ", plot_dir, "  |  R objects / tables saved to: ", r_objects_dir)


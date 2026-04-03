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

# Load pseudobulk TF activity matrices 
mtfr     <- readRDS("/scratch/icbb/igunduz/methylTFR_manuscript/blueprint/mTFR_devs_311025/jaspar2020_distal_deviations.RDS")
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

# zscores
rna <- computeRowZScore(as.matrix(rna))

# Quick sanity checks
stopifnot(is.matrix(rnaViper) | is.data.frame(rnaViper))
stopifnot(is.matrix(mtfr) | is.data.frame(mtfr))
stopifnot(is.matrix(rna) | is.data.frame(rna))
if(!all(colnames(rnaViper) == colnames(rna))){
  warning("Column names (samples) do not match exactly between rnaViper and mtfr.")
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
outfile <- file.path(r_objects_dir, "mtfr_expr_model.hdf5")

# If model file exists, load trained object; otherwise run training
if(file.exists(outfile)){
  MOFAobject.trained <- load_model(outfile)
} else {
  MOFAobject.trained <- run_mofa(MOFAobject, outfile, use_basilisk = TRUE)
}

# Save trained object in rds 
saveRDS(MOFAobject.trained, file.path(r_objects_dir, "mtfr_expr_model.hdf5_bp.rds"))

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
  pivot_wider(names_from = factor, values_from = value)

# default pair: top two factors by R2 (if <2 fallback to 1,2)
f_x <- ifelse(length(top_factors) >= 1, top_factors[5], "1")
f_y <- ifelse(length(top_factors) >= 2, top_factors[1], "2")

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

ggsave(filename = file.path(plot_dir, "factor_scatter_Fx_Fy_Factor3_12.pdf"), plot = p_factors_scatter,
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

# Plot
pvar <- ggplot(r2_long, aes(x = factor, y = r2, fill = view)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_fill_viridis_d(option = "C") +
  labs(x = "Factor", y = "Variance explained (R²)", fill = "View",
       title = "Variance explained per factor by view (top factors)") +
  theme_classic(base_size = 12)

ggsave(filename = file.path(plot_dir, "variance_explained_topFactors.pdf"), plot = pvar,
       width = 8, height = 4)

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

# Also plot absolute sums (dodged), with optional log-scale if large dynamic range
p_modality_abs <- ggplot(agg, aes(x = factor, y = sum_abs, fill = view)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_text(aes(label = round(sum_abs, 1)), position = position_dodge(width = 0.8), vjust = -0.5, size = 3) +
  scale_y_continuous(trans = "log10", labels = scales::comma_format()) +
  scale_fill_manual(values =  c("#ED4B4A", "#6BC75A")) +
  labs(x = "Factor", y = "Sum of |weights| (log10 scale)", fill = "View",
       title = "Absolute modality contribution per factor (sum |weights|)") +
  theme_classic(base_size = 12)

ggsave(filename = file.path(plot_dir, "modality_contribution_expr_abs_topFactors_log.pdf"), plot = p_modality_abs,
       width = 9, height = 4)

# Top TF dotplot (signed + magnitude)
dot_topN <- 6
top_by_factor <- weights_df %>%
  group_by(factor) %>%
  arrange(desc(value_abs)) %>%
  slice_head(n = dot_topN) %>%
  ungroup() %>%
  mutate(feature = factor(feature, levels = unique(feature)))

p_dot <- ggplot(top_by_factor, aes(x = factor, y = feature, size = value_abs, color = value_signed)) +
  geom_point(alpha = 0.9) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(min(top_by_factor$value_signed), max(top_by_factor$value_signed))) +
  scale_size_continuous(range = c(2, 6)) +
  labs(x = "Factor", y = "TF (feature)", color = "Signed weight", size = "|weight|",
       title = paste0("Top ", dot_topN, " TFs per factor (signed weights)")) +
  theme_classic(base_size = 11) +
  theme(axis.text.y = element_text(size = 8))

ggsave(filename = file.path(plot_dir, "topTFs_dotplot_signed.pdf"), plot = p_dot, width = 10, height = 6)

# Combine panels into a single figure for main text (example)
# Use p_factors_scatter, p_modality_frac, p_variance_expl, p_dot
combined_fig <- (p_factors_scatter + p_modality_frac) / (p_variance_expl + p_dot) + plot_annotation(tag_levels = "A")
ggsave(filename = file.path(plot_dir, "figure_combined_panels.pdf"), plot = combined_fig, width = 12, height = 10)

# Save intermediate tables for reproducibility
write.csv(agg, file.path(plot_dir, "modality_contribution_by_factor.csv"), row.names = FALSE)
write.csv(weights_df, file.path(plot_dir, "MOFA_weights_long.csv"), row.names = FALSE)
write.csv(r2_df, file.path(plot_dir, "factor_celltype_R2.csv"), row.names = FALSE)

message("All done. Figures saved to: ", plot_dir, "  |  R objects / tables saved to: ", r_objects_dir)

# Extract factors
factors_long <- get_factors(MOFAobject.trained, as.data.frame = TRUE)
factors_long$celltype <- metadata$celltype  # ensure metadata aligned

# Filter to top factors (optional)
#top_factors <- c("Factor1","Factor2","Factor3","Factor4","Factor7") 
factors_long <- factors_long %>% filter(factor %in% top_factors)

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

# Total variance explained per view
r2_per_factor_df <- get_variance_explained(MOFAobject.trained)$r2_per_factor$group1
r2_per_factor_df <- as.data.frame(r2_per_factor_df)

# Add Factor as a column
r2_per_factor_df <- r2_per_factor_df %>%
  rownames_to_column(var = "Factor")

# Convert to long format and use fractions
r2_long <- r2_per_factor_df %>%
  pivot_longer(cols = -Factor, names_to = "View", values_to = "R2_percent") %>%
  mutate(R2_fraction = R2_percent / 100)

r2_long <- r2_long[r2_long$Factor != "Factor1",]
# Filter for the top factors (e.g., Factor 1-7, where contribution drops off)
top_n_factors <- 7
top_factors_levels <- paste0("Factor", 1:top_n_factors)

r2_long_top <- r2_long %>%
  filter(Factor %in% top_factors_levels) %>%
  mutate(Factor = factor(Factor, levels = top_factors_levels))

# Define custom colors for the views
view_colors <- c("green", "blue")

# Plot 2: Variance Explained per Factor by View

pvar_per_factor <- ggplot(r2_long_top, aes(x = Factor, y = R2_fraction, fill = View)) +
  geom_bar(stat = "identity", position = "dodge") +  # Use dodge for side-by-side bars
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_fill_manual(values = view_colors) +  # Use custom green and blue colors
  labs(x = "Factor", y = "Variance Explained (R²)", fill = "View",
       title = paste0("Variance Explained per Factor by View (Top ", top_n_factors, " Factors)")) +
  theme_classic(base_size = 14)

ggsave(filename = file.path(plot_dir, "variance_explained_per_factor_topFactors.pdf"), plot = pvar_per_factor,
       width = 9, height = 5)

# Plot a correlation for single top TF between chromVar and mtfr (e.g., BATF)
viper_BATF <- rnaViper["BATF",]
mtfr_BATF <- mtfr["BATF",]
data_scatter <- data.frame(
  Viper = viper_BATF,
  mtfr = mtfr_BATF,
  celltype = metadata$celltype
)

p_data_scatter <- ggplot(data_scatter, aes(x = Viper, y = mtfr, color = celltype)) +
  geom_point() +
  geom_smooth(method = "lm", color = "black", se = FALSE) +
  scale_color_manual(values = cell_type_colors) +
  theme_classic(base_size = 13) +
  labs(title = "TF Activity Correlation (chromVar vs mtfr) for BATF")
ggsave(filename = file.path(plot_dir,"BATF_chromVar_vs_mtfr_scatter.pdf"),
 plot = p_data_scatter, width = 6, height = 5)

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

# Annotation
ha <- HeatmapAnnotation(
  celltypes=colnames(mtfr_filtered),
  col = list(celltypes = cell_type_colors
))
# Column order
column_order <- c("megK", "eryt", "gran", "mono", "Mf", "DC", "osteoclast", "NK", "Tcell", "thymocyte", "Bcell", "plasma", "other")

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
pdf(file.path(plot_dir, "mofa_mtfr_bp.pdf"))
draw(p)
dev.off()

# expr
# Get row order
p <- draw(p)
row_order_indices <- row_order(p)  

# Set color scheme
col_fun <- colorRamp2(seq(-2,2, length.out = 100), viridis(100))

pdf(file.path(plot_dir, "mofa_expr_bp.pdf"))
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
level_col <- c("#331E36", "#41337A", "#6EA4BF", "#C2EFEB", "#ECFEE8")
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

# Correlation plots for all tfs
for (motif in common_TFs){
  rna_subs <- rna[motif,]
  mtfr_subs <- mtfr[motif,]
  data_scatter <- data.frame(
    rna = rna_subs,
    mtfr = mtfr_subs,
    celltype = metadata$celltype
)
cor <- cor(data_scatter$rna, data_scatter$mtfr)

p_data_scatter <- ggplot(data_scatter, aes(x = rna, y = mtfr, color = celltype)) +
  geom_point() +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  scale_color_manual(values = cell_type_colors) +
  theme_classic(base_size = 13) +
  labs(caption = paste0("Cor. Coef.:",round(cor, 2))) +
  theme(plot.caption = element_text(hjust = 0.5, size = 10))+
  labs(title = paste0("TF Activity Correlation (rna vs mtfr) for ", motif))
ggsave(filename = file.path(paste0(plot_dir,"scatterplots/",motif,"_rna_vs_mtfr_scatter.pdf")),
 plot = p_data_scatter, width = 6, height = 5)
}

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
ggsave(paste0(plot_dir, "barplots_R2.pdf"), p)


##### Heatmap with less celltypes 
annot_subs <- subset(annot, cellTypeGroup %in% c("Tcell", "gran", "mono", "Mf"))

# Extract new cell types 
groups <- annot_subs$cellTypeGroup

# Change colnames
colnames(mtfr_filtered) <- groups
colnames(rna_filtered) <- groups

# Annotation
ha <- HeatmapAnnotation(
  celltypes=colnames(mtfr_filtered),
  col = list(celltypes = cell_type_colors
))
# Column groups
column_split_factor <- factor(groups, levels = c("Tcell", "gran", "mono", "Mf"))


# MTFR
p<-Heatmap(
  mtfr_filtered,
  row_names_gp = gpar(fontsize = 4),
  top_annotation = ha,
  column_title = "mtfr Z-Scores of TFs from Mofa analysis",
  show_row_names = TRUE,
  show_column_names = FALSE,
  column_split = column_split_factor, # Split columns by the desired order
  cluster_column_slices = FALSE, # Prevent clustering of the groups themselves)
)
pdf(file.path(plot_dir, "mofa_mtfr_bp_v2.pdf"))
draw(p)
dev.off()

# expr
# Color scheme
col_fun <- colorRamp2(seq(-2,2, length.out = 100), viridis(100))

# Row order
p <- draw(p)
row_order_indices <- row_order(p)  

pdf(file.path(plot_dir, "mofa_expr_bp_v2.pdf"))
Heatmap(
  rna_filtered,
  row_names_gp = gpar(fontsize = 4),
  top_annotation = ha,
  col = col_fun,
  column_title = "expr  Z-Scores of TFs from Mofa analysis",
  show_row_names = TRUE,
  show_column_names = FALSE,
  row_order = row_order_indices,
  column_split = column_split_factor, # Split columns by the desired order
  cluster_column_slices = FALSE, # Prevent clustering of the groups themselves)
)
dev.off()

# Calculate row-wise correlation
row_correlation <- sapply(seq_len(nrow(mtfr_filtered)), function(i) {
  cor(mtfr_filtered[i, ], viper_filtered[i, ]) # ), use = "complete.obs")  # Only use rows with complete observations
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

pdf(paste0(plot_dir, "correlation_heatmap_v2.pdf"))
draw(cm)
dev.off()

## Factor column
# Set levels for each factor 
levels <- levels(factor(mofa_filt$factor))
level_col <- c("#331E36", "#41337A", "#6EA4BF", "#C2EFEB", "#ECFEE8")
names(level_col) <- levels

# Add col column
mofa_filt$col <- level_col[mofa_filt$factor]

# Order mofa df same as zscores_filtered/chromVar
mofa_filt <- mofa_filt[match(row.names(mtfr_filtered), mofa_filt$feature),]

path <- file.path(plot_dir, "factor_column_v2.pdf")
pdf(path)
Heatmap(
  mofa_filt$factor, column_names_gp = gpar(fontsize = 9),
  cluster_rows = FALSE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  row_order = row_order_indices,
  heatmap_legend_param = list(title = "Factors"),
  col = level_col)
dev.off()

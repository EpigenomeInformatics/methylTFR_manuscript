#!/usr/bin/env Rscript

#####################################################################
# 02_mofa_integration.R
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
  library(ComplexHeatmap)
})
set.seed(12)

# Define custom colors for cell types
cell_type_colors <- c(
  "B-cell" = "#AE017E",
  "Monocyte" = "#CC4C02",
  "NK-cell" = "#A65628",
  "Th-Mem" = "#41B6C4",
  "Tc-Mem" = "#4292C6",
  "Tc-Naive" = "#888FB5",
  "Th-Naive" = "#C7E9B4"
)

# Paths
plot_dir <- "/icbb/projects/nitschre/methylTFR/figures/figure5"
r_objects_dir <- "/icbb/projects/nitschre/methylTFR/r_objects/"
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(r_objects_dir, recursive = TRUE, showWarnings = FALSE)

# Load pseudobulk TF activity matrices
chromVar <- readRDS(gzfile("/icbb/projects/nitschre/methylTFR/scripts/other/cvar_zscores_adj_psuedobulk.R"))
mtfr <- readRDS(gzfile("/icbb/projects/nitschre/methylTFR/scripts/other/mtfr_zscores_adj_psuedobulk.R"))

# Quick sanity checks
stopifnot(is.matrix(chromVar) | is.data.frame(chromVar))
stopifnot(is.matrix(mtfr) | is.data.frame(mtfr))
if (!all(colnames(chromVar) == colnames(mtfr))) {
  warning("Column names (samples) do not match exactly between chromVar and mtfr.")
}

# Create MOFA object and metadata
data_list <- list(chromVar = as.matrix(chromVar), mtfr = as.matrix(mtfr))
MOFAobject <- create_mofa(data_list)

metadata <- data.frame(
  sample = colnames(mtfr),
  celltype = sub("_.*", "", colnames(mtfr)),
  stringsAsFactors = FALSE
)
samples_metadata(MOFAobject) <- metadata


# Prepare and train MOFA (save HDF5)
MOFAobject <- prepare_mofa(MOFAobject)
outfile <- file.path(r_objects_dir, "MOFAobject_trained_Echo.rds")

# If model file exists, load trained object; otherwise run training
if (file.exists(outfile)) {
  MOFAobject.trained <- load_model(outfile)
} else {
  MOFAobject.trained <- run_mofa(MOFAobject, outfile, use_basilisk = TRUE)
}

# Extract factors and compute ANOVA R2 for celltype
factors_long <- get_factors(MOFAobject.trained, as.data.frame = TRUE)
factors_long$celltype <- metadata$celltype

r2_df <- data.frame(Factor = character(), R2 = numeric(), stringsAsFactors = FALSE)

for (fa in unique(factors_long$factor)) {
  df <- subset(factors_long, factor == fa)
  mod <- aov(value ~ celltype, data = df)
  ss_tab <- summary(mod)[[1]]
  ss_factor <- ss_tab["celltype", "Sum Sq"]
  ss_resid <- ss_tab["Residuals", "Sum Sq"]
  r2_val <- ss_factor / (ss_factor + ss_resid)
  r2_df <- rbind(r2_df, data.frame(Factor = fa, R2 = r2_val, stringsAsFactors = FALSE))
}
r2_df <- r2_df %>% arrange(desc(R2))
write.csv(r2_df, file.path(plot_dir, "r2Values.csv"), row.names = FALSE)

# Keep only top factors (from R2 ranking)
top_factors <- head(r2_df$Factor, 5)#

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

agg$factor <- factor(agg$factor, levels = top_factors)

# Plot stacked contribution (fraction)
p_modality_frac <- ggplot(agg, aes(x = factor, y = frac, fill = view)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = ifelse(frac > 0.03, paste0(round(frac * 100, 1), "%"), "")),
    position = position_stack(vjust = 0.5), size = 3
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_fill_manual(values = c("#ED4B4A", "#6BC75A")) +
  labs(
    x = "Factor", y = "Relative contribution (|weights| fraction)", fill = "View",
    title = "Modality contribution per factor (weights-based)"
  ) +
  theme_classic(base_size = 12)

ggsave(
  filename = file.path(plot_dir, "modality_contribution_fraction_topFactors.pdf"), plot = p_modality_frac,
  width = 8, height = 4
)

# Also plot absolute sums (dodged), with optional log-scale if large dynamic range
p_modality_abs <- ggplot(agg, aes(x = factor, y = sum_abs, fill = view)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_text(aes(label = round(sum_abs, 1)), position = position_dodge(width = 0.8), vjust = -0.5, size = 3) +
  scale_y_continuous(trans = "log10", labels = scales::comma_format()) +
  scale_fill_manual(values = c("#ED4B4A", "#6BC75A")) +
  labs(
    x = "Factor", y = "Sum of |weights| (log10 scale)", fill = "View",
    title = "Absolute modality contribution per factor (sum |weights|)"
  ) +
  theme_classic(base_size = 12)

ggsave(
  filename = file.path(plot_dir, "modality_contribution_abs_topFactors_log.pdf"), plot = p_modality_abs,
  width = 9, height = 4
)

# Filter to top factors (optional)
top_factors <- c("Factor1", "Factor2", "Factor3", "Factor4", "Factor7")
factors_long <- factors_long %>% filter(factor %in% top_factors)

# Plot: strip / dot plot per factor
p <- ggplot(factors_long, aes(x = factor, y = value, color = celltype)) +
  geom_jitter(width = 0.2, height = 0, size = 2, alpha = 0.8) +
  scale_color_manual(values = cell_type_colors) + # Use custom colors
  labs(
    x = "Factor", y = "Factor value", color = "Cell type",
    title = "MOFA factor values per sample by cell type"
  ) +
  theme_classic(base_size = 13) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
ggsave(
  filename = file.path(plot_dir, "factors_stripplot_by_celltype.pdf"), plot = p,
  width = 8, height = 5
)

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

# Rename Tcells
samples_metadata(MOFAobject.trained)$celltype <- case_when(
    grepl("Naive", samples_metadata(MOFAobject.trained)$celltype) ~"T-naive",
    grepl("Mem", samples_metadata(MOFAobject.trained)$celltype) ~ "T-mem",
    TRUE ~ samples_metadata(MOFAobject.trained)$celltype
    )

# Top TFs for each modality
topTfs <- data.frame(
  feature = character(),
  value = numeric(),
  factor = character(),
  view = character()
)

# Create DF with top  5 TFs for each factor and view
for(f in r2_df$Factor[1:5]){
  for(modality in unique(weights_df$view)){
    df <- subset(weights_df, view == modality & factor == f)
    df <- df[order(df$value_abs, decreasing = TRUE),]

    new_df <- df[1:5, c("feature", "value", "factor", "view")]
    new_df$feature <- gsub("_.*$", "", new_df$feature)

    topTfs <- base::rbind(topTfs, new_df)
  }
}
write.csv(topTfs, file.path(plot_dir, "topTFs.csv"), row.names=FALSE)


#### Heatmap of TFs from integration
# Filter for unique TFs
mofaTFs <- unique(topTfs$feature)

# Filter for only Mofa TFs that are present in mtfr and chromVar
common_TFs <- Reduce(intersect, list(mofaTFs, rownames(mtfr), rownames(chromVar)))
chromVar_filtered <- chromVar[common_TFs,]
mtfr_filtered <- mtfr[common_TFs,]

# Extract cell types 
groups <- sapply(strsplit(colnames(mtfr_filtered), "_"), `[`, 1)

# Change colnames
colnames(mtfr_filtered) <- groups
colnames(chromVar_filtered) <- groups

# Plot
path <- file.path(plot_dir, "heatmap_MofaTFs_mtfr.pdf")

# Annotation
ha <- HeatmapAnnotation(
  celltypes=colnames(mtfr_filtered),
  col = list(celltypes = cell_type_colors
))

# Column groups
column_split_factor <- factor(groups, levels = unique(groups))

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
pdf(path)
draw(p)
dev.off()

# Chromvar
# Annotation
ha <- HeatmapAnnotation(
  celltypes=colnames(chromVar_filtered),
  col = list(celltypes = cell_type_colors)
)

path <- file.path(plot_dir, "heatmap_MofaTFs_chromVar.pdf")

# Get row order
p <- draw(p)
row_order_indices <- row_order(p)    

# Set color scheme
colors.cv <- ChrAccR::getConfigElement("colorSchemesCont")
colors.cv <- colors.cv[[".default.div"]]
c <- grDevices::colorRampPalette(colors.cv)(nrow(chromVar_filtered))

col_fun <- colorRamp2(
  seq(-5, 5, length.out = length(c)),
  c
)
# Column groups
column_split_factor <- factor(groups, levels = unique(groups))

pdf(path)
Heatmap(
  chromVar_filtered,
  row_names_gp = gpar(fontsize = 4),
  top_annotation = ha,
  column_title = "chromVar Z-Scores of TFs from Mofa analysis",
  show_row_names = TRUE,
  show_column_names = FALSE,
  column_split = column_split_factor, # Split columns by the desired order
  cluster_column_slices = FALSE, # Prevent clustering of the groups themselves),
  row_order = row_order_indices,
  col = col_fun)
dev.off()

# Calculate row-wise correlation
row_correlation <- sapply(seq_len(nrow(mtfr_filtered)), function(i) {
  cor(mtfr_filtered[i, ], chromVar_filtered[i, ]) # ), use = "complete.obs")  # Only use rows with complete observations
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

# Column with Factors
#  Removing duplicates
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

# Plot a correlation plots for all TFs
for (motif in common_TFs){
  chromVar_subs <- chromVar[motif,]
  mtfr_subs <- mtfr[motif,]
  data_scatter <- data.frame(
    chromVar = chromVar_subs,
    mtfr = mtfr_subs,
    celltype = metadata$celltype
)
cor <- cor(data_scatter$chromVar, data_scatter$mtfr)

p_data_scatter <- ggplot(data_scatter, aes(x = chromVar, y = mtfr, color = celltype)) +
  geom_point() +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  scale_color_manual(values = cell_type_colors) +
  theme_classic(base_size = 13) +
  labs(caption = paste0("Cor. Coef.:",round(cor, 1))) +
  theme(plot.caption = element_text(hjust = 0.5, size = 10))+
  labs(title = paste0("TF Activity Correlation (chromVar vs mtfr) for ", motif))
ggsave(filename = file.path(paste0(plot_dir,motif,"_chromVar_vs_mtfr_scatter.pdf")),
 plot = p_data_scatter, width = 6, height = 5)
}

# Save intermediate tables for reproducibility
write.csv(agg, file.path(r_objects_dir, "modality_contribution_by_factor.csv"), row.names = FALSE)
write.csv(weights_df, file.path(r_objects_dir, "MOFA_weights_long.csv"), row.names = FALSE)
write.csv(r2_df, file.path(r_objects_dir, "factor_celltype_R2.csv"), row.names = FALSE)

message("All done. Figures saved to: ", plot_dir, "  |  R objects / tables saved to: ", r_objects_dir)

#####################################################################
